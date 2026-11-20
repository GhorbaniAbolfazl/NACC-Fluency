"""Inventory the course files and audit the variables used by the project."""

import argparse
import csv
import json
import math
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path

UDS_FILE = "udslong_imputed.csv"
MRI_FILE = "mriscanclean.csv"
MATRIX_FILE = "combined_adjacency_matrix.csv"
PARTICIPANT_ID = "NACCID"

# These groups document the delivered UDS columns by content.
UDS_GROUPS = {
    "identifier": "NACCID",
    "demographics": "SEX HISPANIC HISPOR RACE PRIMLANG EDUC MARISTAT NACCLIVS INDEPEND RESIDENC NACCNIHR NACCAGE",
    "family_history": "NACCFAM NACCMOM NACCDAD NACCAM NACCAMS NACCFM NACCFMS NACCOM NACCOMS NACCFADM NACCFFTD",
    "medications": (
        "ANYMEDS NACCAAAS NACCAANX NACCAC NACCACEI NACCADEP NACCAHTN NACCAMD NACCANGI "
        "NACCAPSY NACCBETA NACCCCBS NACCDBMD NACCDIUR NACCEMD NACCEPMD NACCHTNC NACCLIPL "
        "NACCNSD NACCPDMD NACCVASD"
    ),
    "tobacco": "TOBAC30 TOBAC100 SMOKYRS PACKSPER QUITSMOK",
    "cardiovascular_history": "CVHATT CVAFIB CVANGIO CVBYPASS CVPACDEF CVCHF CVOTHR",
    "neurological_history": "CBSTROKE CBTIA SEIZURES NACCTBI TBI TBIBRIEF TBIEXTEN TBIWOLOS",
    "other_medical_history": (
        "DIABETES DIABTYPE HYPERTEN HYPERCHO B12DEF THYROID ARTHRIT ARTHTYPE ARTHUPEX "
        "ARTHLOEX ARTHSPIN ARTHUNK INCONTU INCONTF"
    ),
    "sleep": "APNEA RBD INSOMN OTHSLEEP",
    "psychiatry_and_substance_use": "ALCOHOL ABUSOTHR PTSD BIPOLAR SCHIZ DEP2YRS DEPOTHR ANXIETY OCD NPSYDEV PSYCDIS",
    "physical_exam_and_senses": "HEIGHT WEIGHT BPSYS BPDIAS HRATE VISION VISCORR VISWCORR HEARING HEARAID HEARWAID NACCBMI",
    "neurological_exam": "NACCNREX NORMEXAM",
    "behavioral_symptoms": (
        "NPIQINF DEL DELSEV HALL HALLSEV AGIT AGITSEV DEPD DEPDSEV ANX ANXSEV ELAT ELATSEV "
        "APA APASEV DISN DISNSEV IRR IRRSEV MOT MOTSEV NITE NITESEV APP APPSEV"
    ),
    "depression_scale": (
        "NOGDS SATIS DROPACT EMPTY BORED SPIRITS AFRAID HAPPY HELPLESS STAYHOME MEMPROB "
        "WONDRFUL WRTHLESS ENERGY HOPELESS BETTER NACCGDS"
    ),
    "everyday_function": "BILLS TAXES SHOPPING GAMES STOVE MEALPREP EVENTS PAYATTN REMDATES TRAVEL",
    "subjective_memory_decline": "DECSUB",
    "cognitive_tests": "ANIMALS VEG TRAILA TRAILARR TRAILALI TRAILB TRAILBRR TRAILBLI MOCATOTS",
    "genetics": "NACCAPOE",
    "clinical_diagnosis": "NACCUDSD",
    "clinical_dementia_rating": "MEMORY ORIENT JUDGMENT COMMUN HOMEHOBB PERSCARE CDRSUM CDRGLOB",
}

CATEGORY_LABELS = {
    "SEX": {1: "male", 2: "female"},
    "NACCUDSD": {1: "normal cognition", 2: "impaired without MCI", 3: "MCI", 4: "dementia"},
}


@dataclass(frozen=True)
class FieldRule:
    minimum: int
    maximum: int | None = None
    integer: bool = False
    strictly_positive: bool = False
    missing_codes: tuple[int, ...] = ()


# The variable-specific references are linked in PROJECT_CARD.md.
FLUENCY_RULE = FieldRule(0, 77, integer=True, missing_codes=(-4, 95, 96, 97, 98))
THICKNESS_RULE = FieldRule(0, strictly_positive=True)
UDS_RULES = {
    "ANIMALS": FLUENCY_RULE,
    "VEG": FLUENCY_RULE,
    "NACCAGE": FieldRule(18, 120),
    "SEX": FieldRule(1, 2, integer=True),
    "EDUC": FieldRule(0, 36, integer=True, missing_codes=(99,)),
    "NACCUDSD": FieldRule(1, 4, integer=True),
}
MRI_RULES = {
    "LH_INFERIORPARIETAL_AVGTH": THICKNESS_RULE,
    "RH_INFERIORPARIETAL_AVGTH": THICKNESS_RULE,
    "NACCUDSD": UDS_RULES["NACCUDSD"],
}


def assess_value(raw: str, rule: FieldRule) -> tuple[int | float | None, str | None]:
    """Return the usable value or a single reason the raw value is unavailable."""
    raw = raw.strip()
    if not raw:
        return None, "blank"
    try:
        number = Decimal(raw)
    except InvalidOperation:
        return None, "non_numeric"
    if not number.is_finite():
        return None, "non_finite"
    if number in rule.missing_codes:
        return None, f"special_code_{int(number)}"
    if number < rule.minimum or (rule.maximum is not None and number > rule.maximum):
        return None, "out_of_range"
    if rule.integer:
        if number != number.to_integral_value():
            return None, "non_integer"
        return int(number), None
    value = float(number)
    if not math.isfinite(value):
        return None, "non_finite"
    if rule.strictly_positive and value <= 0:
        return None, "non_positive"
    return value, None


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a rectangular CSV and keep every input value as text."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, [])
        if not header or len(header) != len(set(header)):
            raise ValueError(f"{path.name}: missing or duplicate column names")
        rows = []
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{path.name}, line {reader.line_num}: unexpected number of fields")
            rows.append(row)
    if not rows:
        raise ValueError(f"{path.name}: no data rows")
    return header, rows


def read_participant_file(path: Path, rules: dict[str, FieldRule]) -> tuple[list[str], list[dict[str, str]]]:
    header, rows = read_csv(path)
    missing = sorted(({PARTICIPANT_ID} | rules.keys()) - set(header))
    if missing or "" in header:
        raise ValueError(f"{path.name}: missing columns {missing} or blank column name")
    return header, [dict(zip(header, row)) for row in rows]


def group_uds(columns: list[str]) -> dict[str, list[str]]:
    group_by_column = {
        column: group for group, members in UDS_GROUPS.items() for column in members.split()
    }
    groups: dict[str, list[str]] = {group: [] for group in UDS_GROUPS}
    for column in columns:
        group = group_by_column.get(column, "uncategorized")
        groups.setdefault(group, []).append(column)
    return groups


def group_mri(columns: list[str]) -> dict:
    metadata = {PARTICIPANT_ID, "NACCUDSD"}
    measures = {"GVOL": "volume_mm3", "AVGTH": "thickness_mm"}
    regions: dict[str, dict[str, str]] = {}
    for column in columns:
        if column in metadata:
            continue
        parts = column.split("_")
        if len(parts) != 3 or parts[0] not in {"LH", "RH"} or parts[2] not in measures:
            raise ValueError(f"Unmapped MRI column: {column}")
        hemisphere, region, measure = parts
        regions.setdefault(f"{hemisphere}_{region}", {})[measures[measure]] = column
    return {"participant_fields": [column for column in columns if column in metadata], "regions": regions}


def inspect_matrix(path: Path, mri_regions: dict) -> dict:
    header, rows = read_csv(path)
    labels = header[1:]
    if header[0] or [row[0] for row in rows] != labels:
        raise ValueError(f"{path.name}: row and column region labels do not match")
    matrix = [row[1:] for row in rows]
    sides = {"left": "LH", "right": "RH"}
    normalized_labels = set()
    for label in labels:
        hemisphere, separator, region = label.partition("_")
        if hemisphere not in sides or not separator or not region:
            raise ValueError(f"{path.name}: unrecognized region label {label}")
        normalized_labels.add(f"{sides[hemisphere]}_{region.replace('_', '').upper()}")
    return {
        "rows": len(rows),
        "data_columns": len(labels),
        "binary": all(value in {"0", "1"} for row in matrix for value in row),
        "symmetric": matrix == [list(column) for column in zip(*matrix)],
        "matches_mri_region_labels": len(normalized_labels) == len(labels) and normalized_labels == set(mri_regions),
    }


def summarize_field(rows: list[dict[str, str]], column: str, rule: FieldRule) -> dict:
    values = []
    issues: Counter[str] = Counter()
    for row in rows:
        value, issue = assess_value(row[column], rule)
        if issue is None:
            values.append(value)
        else:
            issues[issue] += 1
    missing_count = issues["blank"] + sum(count for key, count in issues.items() if key.startswith("special_code_"))
    result = {
        "total": len(rows),
        "usable": len(values),
        "unusable": len(rows) - len(values),
        "unusable_percent": 100 * (len(rows) - len(values)) / len(rows),
        "missing_or_special_code": missing_count,
        "missing_or_special_code_percent": 100 * missing_count / len(rows),
        "issues": dict(sorted(issues.items())),
    }
    if column in CATEGORY_LABELS:
        result["usable_category_counts"] = {
            label: values.count(code) for code, label in CATEGORY_LABELS[column].items()
        }
    else:
        result["usable_minimum"] = min(values, default=None)
        result["usable_maximum"] = max(values, default=None)
    return result


def summarize_file(columns: list[str], rows: list[dict[str, str]], rules: dict[str, FieldRule]) -> dict:
    ids = Counter(row[PARTICIPANT_ID].strip() for row in rows)
    blank_ids = ids.pop("", 0)
    return {
        "rows": len(rows),
        "columns": len(columns),
        "unique_nonblank_ids": len(ids),
        "blank_ids": blank_ids,
        "extra_rows_with_repeated_id": sum(count - 1 for count in ids.values()),
        "selected_variables": {column: summarize_field(rows, column, rule) for column, rule in rules.items()},
    }


def audit(data_dir: Path) -> dict:
    uds_columns, uds_rows = read_participant_file(data_dir / UDS_FILE, UDS_RULES)
    mri_columns, mri_rows = read_participant_file(data_dir / MRI_FILE, MRI_RULES)
    uds_map = group_uds(uds_columns)
    mri_map = group_mri(mri_columns)
    return {
        "uds": {"groups": uds_map, **summarize_file(uds_columns, uds_rows, UDS_RULES)},
        "mri": {"map": mri_map, **summarize_file(mri_columns, mri_rows, MRI_RULES)},
        "adjacency_matrix": inspect_matrix(data_dir / MATRIX_FILE, mri_map["regions"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args()
    print(json.dumps(audit(args.data_dir), ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

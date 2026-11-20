"""Join UDS and MRI by participant and produce the phase-one sample flow."""

import argparse
import csv
import json
from collections import Counter
from hashlib import sha256
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from .data_audit import (
    CATEGORY_LABELS,
    MRI_FILE,
    MRI_RULES,
    PARTICIPANT_ID,
    UDS_FILE,
    UDS_RULES,
    assess_value,
    read_participant_file,
)

STATUS = "NACCUDSD"
MRI_STATUS = "MRI_NACCUDSD"
NORMAL_COGNITION = 1
LEFT_THICKNESS = "LH_INFERIORPARIETAL_AVGTH"
PAIR_FIELDS = ("ANIMALS", LEFT_THICKNESS)
MODEL_FIELDS = (*PAIR_FIELDS, "NACCAGE", "SEX", "EDUC")

# Required fields and whether analysis is restricted to normal cognition.
SAMPLES = {
    "status_and_mri": ((STATUS, LEFT_THICKNESS), False),
    "paired_all": (PAIR_FIELDS, False),
    "paired_normal": (PAIR_FIELDS, True),
    "model_all": (MODEL_FIELDS, False),
    "model_normal": (MODEL_FIELDS, True),
}

MRI_COLUMNS = {column: MRI_STATUS if column == STATUS else column for column in MRI_RULES}
ANALYSIS_COLUMNS = (*UDS_RULES, *MRI_COLUMNS.values())
OUTPUT_COLUMNS = (
    PARTICIPANT_ID,
    *ANALYSIS_COLUMNS,
    *(f"eligible_{sample}" for sample in SAMPLES),
    "status_check",
    *(f"{column}_RAW" for column in ANALYSIS_COLUMNS),
    *(f"{column}_ISSUE" for column in ANALYSIS_COLUMNS),
)


def index_participants(rows: list[dict[str, str]], filename: str) -> dict[str, dict[str, str]]:
    indexed = {}
    for row_number, row in enumerate(rows, start=1):
        participant_id = row[PARTICIPANT_ID].strip()
        if not participant_id:
            raise ValueError(f"{filename}: blank NACCID in data row {row_number}")
        if participant_id in indexed:
            raise ValueError(f"{filename}: repeated NACCID in data row {row_number}")
        indexed[participant_id] = row
    return indexed


def assess_participant(participant_id: str, uds_row: dict[str, str], mri_row: dict[str, str]) -> dict:
    row = {PARTICIPANT_ID: participant_id}
    for column, rule in UDS_RULES.items():
        raw = uds_row[column]
        row[f"{column}_RAW"] = raw
        row[column], row[f"{column}_ISSUE"] = assess_value(raw, rule)
    for column, rule in MRI_RULES.items():
        target = MRI_COLUMNS[column]
        raw = mri_row[column]
        row[f"{target}_RAW"] = raw
        row[target], row[f"{target}_ISSUE"] = assess_value(raw, rule)
    if row[STATUS] is None or row[MRI_STATUS] is None:
        row["status_check"] = "unavailable"
    elif row[STATUS] == row[MRI_STATUS]:
        row["status_check"] = "equal"
    else:
        row["status_check"] = "different"
    return row


def summarize_sample(rows: list[dict], name: str, fields: tuple[str, ...], normal_only: bool) -> dict:
    categories: Counter[int | None] = Counter()
    exclusions: Counter[str] = Counter()
    for row in rows:
        if normal_only and row[STATUS] is None:
            reason = "cognitive_status_unavailable"
        elif normal_only and row[STATUS] != NORMAL_COGNITION:
            reason = "outside_normal_cognition"
        else:
            reason = next((column for column in fields if row[column] is None), None)

        row[f"eligible_{name}"] = int(reason is None)
        if reason is None:
            categories[row[STATUS]] += 1
        else:
            exclusions[reason] += 1

    return {
        "required_values": list(fields),
        "normal_cognition_only": normal_only,
        "included": sum(categories.values()),
        "excluded": sum(exclusions.values()),
        "first_exclusion_reason": dict(exclusions),
        "included_by_cognitive_status": {
            **{label: categories[code] for code, label in CATEGORY_LABELS[STATUS].items()},
            "unavailable": categories[None],
        },
    }


def build_dataset(data_dir: Path) -> tuple[list[dict], dict]:
    _, uds_rows = read_participant_file(data_dir / UDS_FILE, UDS_RULES)
    _, mri_rows = read_participant_file(data_dir / MRI_FILE, MRI_RULES)
    uds = index_participants(uds_rows, UDS_FILE)
    mri = index_participants(mri_rows, MRI_FILE)
    shared_ids = sorted(uds.keys() & mri.keys())
    if not shared_ids:
        raise ValueError("The UDS and MRI files have no matching participant identifiers")

    rows = [assess_participant(participant_id, uds[participant_id], mri[participant_id]) for participant_id in shared_ids]
    status_checks = dict.fromkeys(("equal", "unavailable", "different"), 0)
    status_checks.update(Counter(row["status_check"] for row in rows))
    if status_checks["different"]:
        raise ValueError(
            f"UDS and MRI cognitive status disagree for {status_checks['different']} matched participants"
        )

    flow = {
        "linkage": {
            "uds_participants": len(uds),
            "mri_participants": len(mri),
            "matched_participants": len(rows),
            "uds_only": len(uds) - len(rows),
            "mri_only": len(mri) - len(rows),
            "matched_with_usable_left_thickness": sum(row[LEFT_THICKNESS] is not None for row in rows),
            "matched_with_unusable_left_thickness": sum(row[LEFT_THICKNESS] is None for row in rows),
            "cognitive_status_agreement": status_checks,
        },
        "matched_normal_cognition_before_eligibility": sum(row[STATUS] == NORMAL_COGNITION for row in rows),
        "unusable_values_in_matched_table": {
            column: sum(row[column] is None for row in rows) for column in ANALYSIS_COLUMNS
        },
        "samples": {
            name: summarize_sample(rows, name, fields, normal_only)
            for name, (fields, normal_only) in SAMPLES.items()
        },
    }
    return rows, flow


def write_output(path: Path, data: bytes) -> None:
    """Replace each complete output file after writing it on the same filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent) as directory:
        staged = Path(directory) / path.name
        staged.write_bytes(data)
        staged.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args()

    rows, flow = build_dataset(args.data_dir)
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    csv_bytes = buffer.getvalue().encode("utf-8")
    flow["participants_csv_sha256"] = sha256(csv_bytes).hexdigest()
    report = json.dumps(flow, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(Path("data/processed/participants.csv"), csv_bytes)
    write_output(Path("reports/phase1/sample_flow.json"), report.encode("utf-8"))


if __name__ == "__main__":
    main()

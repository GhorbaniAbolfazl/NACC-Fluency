"""Summarize participants, model samples, and cognitive groups for phase one."""

import csv
import json
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from io import StringIO
from pathlib import Path
from statistics import fmean, quantiles, stdev

from .build_dataset import LEFT_THICKNESS, NORMAL_COGNITION, OUTPUT_COLUMNS, SAMPLES, STATUS, write_output
from .data_audit import CATEGORY_LABELS

NUMERIC_FIELDS = {
    "ANIMALS": int,
    LEFT_THICKNESS: float,
    "NACCAGE": float,
    "EDUC": int,
}
SAMPLE_FLAGS = tuple(f"eligible_{sample}" for sample in SAMPLES)
PARTICIPANTS_FILE = Path("data/processed/participants.csv")
FLOW_FILE = Path("reports/phase1/sample_flow.json")
REPORT_FILE = Path("reports/phase1/descriptive_statistics.json")


def load_participants(
    table_path: Path = PARTICIPANTS_FILE, flow_path: Path = FLOW_FILE
) -> tuple[list[dict], dict]:
    """Read analysis values and flags from a table that matches its sample flow."""
    table_bytes = table_path.read_bytes()
    flow = json.loads(flow_path.read_text(encoding="utf-8"))
    if sha256(table_bytes).hexdigest() != flow["participants_csv_sha256"]:
        raise ValueError("Participant table and sample flow differ. Rebuild with python3 -m src.build_dataset.")

    types = {**NUMERIC_FIELDS, **dict.fromkeys(CATEGORY_LABELS, int)}
    with StringIO(table_bytes.decode("utf-8"), newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        if reader.fieldnames != list(OUTPUT_COLUMNS):
            raise ValueError("Unexpected columns in the participant table")
        rows = []
        for raw in reader:
            row = {column: convert(raw[column]) if raw[column] else None for column, convert in types.items()}
            for flag in SAMPLE_FLAGS:
                if raw[flag] not in {"0", "1"}:
                    raise ValueError(f"Invalid sample flag at line {reader.line_num}: {flag}")
                row[flag] = raw[flag] == "1"
            rows.append(row)

    if len(rows) != flow["linkage"]["matched_participants"]:
        raise ValueError("Row count differs from the matched sample")
    return rows, flow


def availability(usable: int, total: int) -> dict:
    unavailable = total - usable
    return {
        "total": total,
        "usable": usable,
        "unavailable": unavailable,
        "unavailable_percent": 100 * unavailable / total if total else None,
    }


def numeric_summary(rows: list[dict], column: str) -> dict:
    values = [row[column] for row in rows if row[column] is not None]
    size = len(values)
    # Keep values exactly on a fence inside it, including decimal MRI measurements.
    decimal_values = [Decimal(str(value)) for value in values]
    if size > 1:
        q1, median, q3 = quantiles(decimal_values, n=4, method="inclusive")
    else:
        q1 = median = q3 = next(iter(decimal_values), None)
    iqr = q3 - q1 if size else None
    review = None
    if size >= 4 and iqr > 0:
        margin = Decimal("1.5") * iqr
        lower = q1 - margin
        upper = q3 + margin
        review = {
            "lower_fence": float(lower),
            "upper_fence": float(upper),
            "below": sum(value < lower for value in decimal_values),
            "above": sum(value > upper for value in decimal_values),
        }
    return {
        **availability(size, len(rows)),
        "mean": fmean(values) if values else None,
        "sample_sd": stdev(values) if size > 1 else None,
        "minimum": min(values, default=None),
        "q1": float(q1) if size else None,
        "median": float(median) if size else None,
        "q3": float(q3) if size else None,
        "maximum": max(values, default=None),
        "iqr": float(iqr) if size else None,
        "iqr_review": review,
    }


def category_summary(rows: list[dict], column: str) -> dict:
    counts = Counter(row[column] for row in rows)
    usable = len(rows) - counts[None]
    return {
        **availability(usable, len(rows)),
        "levels": [
            {
                "code": code,
                "label": label,
                "count": counts[code],
                "percent_of_all": 100 * counts[code] / len(rows) if rows else None,
                "percent_of_usable": 100 * counts[code] / usable if usable else None,
            }
            for code, label in CATEGORY_LABELS[column].items()
        ],
    }


def empirical_frequencies(rows: list[dict], column: str) -> dict:
    counts = Counter(row[column] for row in rows if row[column] is not None)
    usable = sum(counts.values())
    cumulative = 0
    levels = []
    for value, count in sorted(counts.items()):
        cumulative += count
        levels.append(
            {
                "value": value,
                "count": count,
                "percent_of_usable": 100 * count / usable,
                "cumulative_percent_of_usable": 100 * cumulative / usable,
            }
        )
    return {**availability(usable, len(rows)), "levels": levels}


def summarize_group(rows: list[dict]) -> dict:
    return {
        "participants": len(rows),
        "numeric": {column: numeric_summary(rows, column) for column in NUMERIC_FIELDS},
        "categorical": {column: category_summary(rows, column) for column in CATEGORY_LABELS},
    }


def describe(rows: list[dict], flow: dict) -> dict:
    by_status = {
        label: summarize_group([row for row in rows if row[STATUS] == code])
        for code, label in CATEGORY_LABELS[STATUS].items()
    }
    by_status["unavailable"] = summarize_group([row for row in rows if row[STATUS] is None])

    normal = [row for row in rows if row[STATUS] == NORMAL_COGNITION]
    if len(normal) != flow["matched_normal_cognition_before_eligibility"]:
        raise ValueError("Normal-cognition count differs from the sample flow")
    samples = {"matched_all": rows, "matched_normal": normal}
    for sample in SAMPLES:
        eligible = [row for row in rows if row[f"eligible_{sample}"]]
        if len(eligible) != flow["samples"][sample]["included"]:
            raise ValueError(f"Count for {sample} differs from the sample flow")
        samples[sample] = eligible

    return {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "methods": {
            "quantitative_denominator": "Available values per variable within each stated sample or group",
            "sample_sd_ddof": 1,
            "quartiles": "Linear interpolation at ranks (n - 1) * p, zero-indexed",
            "iqr_review": "Flag values beyond 1.5 IQR fences when n >= 4 and IQR > 0; only for review",
            "group_denominator": "All matched participants in each UDS cognitive status before sample eligibility",
            "units": {"ANIMALS": "correct responses in 60 seconds", LEFT_THICKNESS: "mm", "NACCAGE": "years", "EDUC": "years"},
        },
        "samples": {name: summarize_group(members) for name, members in samples.items()},
        "by_cognitive_status": by_status,
        "discrete_frequencies": {
            name: {column: empirical_frequencies(samples[name], column) for column in ("ANIMALS", "EDUC")}
            for name in ("matched_all", "matched_normal")
        },
    }


def main() -> None:
    rows, flow = load_participants()
    report = json.dumps(describe(rows, flow), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    write_output(REPORT_FILE, report.encode("utf-8"))


if __name__ == "__main__":
    main()

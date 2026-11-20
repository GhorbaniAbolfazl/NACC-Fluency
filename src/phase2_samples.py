"""Validate the participant samples for Phase 2 and write their aggregate definitions."""

import json
from pathlib import Path

from .build_dataset import (
    LEFT_THICKNESS,
    NORMAL_COGNITION,
    SAMPLES,
    STATUS,
    summarize_sample,
    write_output,
)
from .data_audit import MRI_FILE, PARTICIPANT_ID, UDS_FILE
from .descriptive_statistics import FLOW_FILE, PARTICIPANTS_FILE, load_participants

REPORT_FILE = Path("reports/phase2/sample_definitions.json")
PHASE2_SAMPLES = ("status_and_mri", "paired_normal")


def load_phase2_samples() -> tuple[dict[str, list[dict]], dict]:
    """Use the phase-one table after checking both phase-two membership rules."""
    rows, flow = load_participants()
    checked_rows = [row.copy() for row in rows]
    samples = {}

    for name in PHASE2_SAMPLES:
        flag = f"eligible_{name}"
        summary = summarize_sample(checked_rows, name, *SAMPLES[name])
        if summary != flow["samples"][name] or any(
            recorded[flag] != checked[flag] for recorded, checked in zip(rows, checked_rows)
        ):
            raise ValueError(
                f"Membership for {name} differs from the phase-one rules. "
                "Rebuild with python3 -m src.build_dataset."
            )

        members = [row for row in rows if row[flag]]
        if not members:
            raise ValueError(f"The phase-two sample {name} is empty.")
        samples[name] = members

    return samples, flow


def main() -> None:
    _, flow = load_phase2_samples()
    report = {
        "inputs": {
            "participant_table": str(PARTICIPANTS_FILE),
            "participant_table_sha256": flow["participants_csv_sha256"],
            "phase_one_sample_flow": str(FLOW_FILE),
            "source_files": {"uds": UDS_FILE, "mri": MRI_FILE},
        },
        "observation": {
            "unit": "One participant with matched UDS and MRI records",
            "join_key": PARTICIPANT_ID,
            "index_record_rule": (
                "Use the single record supplied per participant in each file. "
                "The phase-one build rejects blank identifiers and repeated participant rows before matching."
            ),
            "visit_and_scan_timing": (
                "The supplied fields do not identify the clinical visit date or number, "
                "the MRI scan date, or the interval between the two measurements."
            ),
        },
        "linkage": flow["linkage"],
        "normal_cognition": {"uds_column": STATUS, "code": NORMAL_COGNITION},
        "samples": {
            "status_and_mri": {
                "use": "UDS cognitive-status and MRI thickness event probabilities in all cognitive groups",
                **flow["samples"]["status_and_mri"],
            },
            "paired_normal": {
                "use": "Score and thickness distributions and covariance among participants with normal cognition",
                "normal_cognition_before_eligibility": flow["matched_normal_cognition_before_eligibility"],
                **flow["samples"]["paired_normal"],
            },
        },
        "planned_mri_grouping": {
            "column": LEFT_THICKNESS,
            "unit": "mm",
            "reference_sample": "status_and_mri",
            "cutoff_rule": "Pooled median of usable thickness values in the reference sample",
            "groups": ["below median", "at or above median"],
            "purpose": "Descriptive grouping for the probability exercise",
        },
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

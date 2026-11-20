"""Summarize empirical probabilities for the Phase 2 UDS and MRI events."""

import json
from collections import Counter
from decimal import Decimal
from pathlib import Path
from statistics import median

from .build_dataset import LEFT_THICKNESS, NORMAL_COGNITION, STATUS, write_output
from .phase2_samples import load_phase2_samples

REPORT_FILE = Path("reports/phase2/event_probabilities.json")
SAMPLE_NAME = "status_and_mri"


def probability(numerator: int, denominator: int) -> dict:
    """Keep the original group denominator, including for an undefined conditional."""
    return {
        "numerator": numerator,
        "denominator": denominator,
        "probability": numerator / denominator if denominator else None,
    }


def summarize_events(rows: list[dict]) -> dict:
    """Use the pooled thickness median within the validated status and MRI sample."""
    thicknesses = [Decimal(str(row[LEFT_THICKNESS])) for row in rows]
    cutoff = median(thicknesses)

    counts = Counter(
        (row[STATUS] == NORMAL_COGNITION, thickness >= cutoff)
        for row, thickness in zip(rows, thicknesses)
    )
    a_and_b = counts[True, True]
    a_and_not_b = counts[True, False]
    not_a_and_b = counts[False, True]
    not_a_and_not_b = counts[False, False]
    a = a_and_b + a_and_not_b
    not_a = not_a_and_b + not_a_and_not_b
    b = a_and_b + not_a_and_b
    not_b = a_and_not_b + not_a_and_not_b
    total = len(rows)

    a_given_b = probability(a_and_b, b)
    a_given_not_b = probability(a_and_not_b, not_b)
    rate_at_or_above = a_given_b["probability"]
    rate_below = a_given_not_b["probability"]
    gap = None
    if rate_at_or_above is not None and rate_below is not None:
        gap = 100 * (rate_at_or_above - rate_below)

    return {
        "sample_size": total,
        "events": {
            "A": {
                "definition": "Cognitively normal",
                "source": "UDS",
                "column": STATUS,
                "code": NORMAL_COGNITION,
            },
            "B": {
                "definition": "Left inferior parietal cortical thickness at or above the pooled sample median",
                "source": "MRI",
                "column": LEFT_THICKNESS,
                "median_mm": float(cutoff),
                "equal_to_median": thicknesses.count(cutoff),
            },
        },
        "contingency_table": {
            "A": {"B": a_and_b, "not_B": a_and_not_b, "total": a},
            "not_A": {"B": not_a_and_b, "not_B": not_a_and_not_b, "total": not_a},
            "total": {"B": b, "not_B": not_b, "total": total},
        },
        "marginal_probabilities": {
            "A": probability(a, total),
            "not_A": probability(not_a, total),
            "B": probability(b, total),
            "not_B": probability(not_b, total),
        },
        "joint_probabilities": {
            "A_and_B": probability(a_and_b, total),
            "A_and_not_B": probability(a_and_not_b, total),
            "not_A_and_B": probability(not_a_and_b, total),
            "not_A_and_not_B": probability(not_a_and_not_b, total),
        },
        "conditional_probabilities": {
            "A_given_B": a_given_b,
            "A_given_not_B": a_given_not_b,
        },
        "conditional_difference_percentage_points": gap,
    }


def main() -> None:
    samples, flow = load_phase2_samples()
    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "sample": SAMPLE_NAME,
        **summarize_events(samples[SAMPLE_NAME]),
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

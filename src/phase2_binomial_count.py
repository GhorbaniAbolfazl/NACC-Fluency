"""Model the number of cognitively normal participants in 20 hypothetical draws."""

import json
from fractions import Fraction
from math import comb
from pathlib import Path

from .build_dataset import NORMAL_COGNITION, STATUS, write_output
from .phase2_event_probabilities import SAMPLE_NAME, probability
from .phase2_samples import load_phase2_samples

REPORT_FILE = Path("reports/phase2/binomial_count.json")
DRAWS = 20
TAIL_MINIMUM = 15  # Three quarters of the fixed number of draws.


def summarize_count(population_size: int, normal_count: int) -> dict:
    """Compare uniform draws with and without replacement from the same finite set."""
    if population_size < DRAWS or not 0 <= normal_count <= population_size:
        raise ValueError("The reference sample must contain at least 20 participants with a valid success count.")

    success_rate = Fraction(normal_count, population_size)
    mean = DRAWS * success_rate
    variance_with_replacement = DRAWS * success_rate * (1 - success_rate)
    variance_multiplier = Fraction(population_size - DRAWS, population_size - 1)
    variance_without_replacement = variance_with_replacement * variance_multiplier

    binomial = [
        comb(DRAWS, k) * success_rate**k * (1 - success_rate)**(DRAWS - k)
        for k in range(DRAWS + 1)
    ]
    other_count = population_size - normal_count
    combinations = comb(population_size, DRAWS)
    hypergeometric = [
        Fraction(comb(normal_count, k) * comb(other_count, DRAWS - k), combinations)
        for k in range(DRAWS + 1)
    ]

    binomial_tail = sum(binomial[TAIL_MINIMUM:])
    hypergeometric_tail = sum(hypergeometric[TAIL_MINIMUM:])
    return {
        "reference_sample_size": population_size,
        "scope": "Hypothetical draws from the fixed course participant set, with uniform selection.",
        "success_event": {"label": "Cognitively normal", "uds_column": STATUS, "code": NORMAL_COGNITION},
        "single_draw_success_probability": probability(normal_count, population_size),
        "random_variable": "Number of draws yielding a participant with normal cognition",
        "number_of_draws": DRAWS,
        "tail_minimum": TAIL_MINIMUM,
        "with_replacement": {
            "distribution": "binomial",
            "draw_rule": (
                "Return each selected participant before the next draw. "
                "Draws are independent with the same success probability."
            ),
            "expected_count": float(mean),
            "variance": float(variance_with_replacement),
            "probability_at_least_tail_minimum": float(binomial_tail),
        },
        "without_replacement": {
            "distribution": "hypergeometric",
            "draw_rule": (
                "Select distinct participants. Each conditional success probability "
                "uses the counts remaining after earlier draws."
            ),
            "expected_count": float(mean),
            "variance": float(variance_without_replacement),
            "variance_multiplier_vs_with_replacement": float(variance_multiplier),
            "probability_at_least_tail_minimum": float(hypergeometric_tail),
        },
        "tail_difference_percentage_points_with_minus_without_replacement": float(
            100 * (binomial_tail - hypergeometric_tail)
        ),
        "probability_mass": [
            {
                "count": k,
                "with_replacement": float(with_replacement),
                "without_replacement": float(without_replacement),
            }
            for k, (with_replacement, without_replacement) in enumerate(zip(binomial, hypergeometric))
        ],
    }


def main() -> None:
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    normal_count = sum(row[STATUS] == NORMAL_COGNITION for row in rows)
    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "sample": SAMPLE_NAME,
        **summarize_count(len(rows), normal_count),
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

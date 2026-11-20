"""Compare empirical event probabilities with independence and apply Bayes' rule."""

import json
from fractions import Fraction
from pathlib import Path

from .build_dataset import write_output
from .phase2_event_probabilities import SAMPLE_NAME, probability, summarize_events
from .phase2_samples import load_phase2_samples

REPORT_FILE = Path("reports/phase2/independence_bayes.json")


def compare_probabilities(event_summary: dict) -> dict:
    """Use counts from the Stage 2 table; evaluate both identities before rounding."""
    table = event_summary["contingency_table"]
    n = table["total"]["total"]
    a = table["A"]["total"]
    b = table["total"]["B"]
    both = table["A"]["B"]

    prior = Fraction(a, n)
    evidence = Fraction(b, n)
    observed_joint = Fraction(both, n)
    independent_joint = prior * evidence
    expected_count = n * independent_joint

    posterior = None
    matches_direct = None
    update = None
    if a and b:
        likelihood = Fraction(both, a)
        posterior = likelihood * prior / evidence
        matches_direct = posterior == Fraction(both, b)
        update = float(100 * (posterior - prior))

    return {
        "sample_size": n,
        "events": event_summary["events"],
        "independence": {
            "observed_joint_A_and_B": probability(both, n),
            "joint_if_independent": probability(a * b, n * n),
            "observed_minus_expected_percentage_points": float(100 * (observed_joint - independent_joint)),
            "observed_count": both,
            "expected_count_if_independent": float(expected_count),
            "observed_minus_expected_count": float(both - expected_count),
            "empirical_joint_equals_product": observed_joint == independent_joint,
        },
        "bayes": {
            "prior_A": probability(a, n),
            "likelihood_B_given_A": probability(both, a),
            "evidence_B": probability(b, n),
            "posterior_A_given_B_via_bayes": float(posterior) if posterior is not None else None,
            "posterior_A_given_B_direct": probability(both, b),
            "matches_direct_exactly": matches_direct,
            "posterior_minus_prior_percentage_points": update,
        },
    }


def main() -> None:
    samples, flow = load_phase2_samples()
    event_summary = summarize_events(samples[SAMPLE_NAME])
    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "sample": SAMPLE_NAME,
        **compare_probabilities(event_summary),
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

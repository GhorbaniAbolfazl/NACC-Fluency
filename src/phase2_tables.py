"""Generate the shared numeric values for both Phase 2 LaTeX reports."""

import json
from pathlib import Path

from .build_dataset import LEFT_THICKNESS, PAIR_FIELDS, write_output
from .phase2_samples import load_phase2_samples

REPORT_DIR = Path("reports/phase2")
TABLE_FILE = REPORT_DIR / "tables.tex"


def load_results() -> dict:
    """Read aggregate results from one verified participant table."""
    samples, flow = load_phase2_samples()
    fingerprint = flow["participants_csv_sha256"]
    results = {}
    for name in ("sample_definitions", "event_probabilities", "independence_bayes",
                 "binomial_count", "distributions", "paired_covariance"):
        report = json.loads((REPORT_DIR / f"{name}.json").read_text(encoding="utf-8"))
        if name == "sample_definitions":
            source_fingerprint = report["inputs"]["participant_table_sha256"]
            for sample, members in samples.items():
                if report["samples"][sample]["included"] != len(members):
                    raise ValueError(f"Inconsistent sample count in {name}. Regenerate the Phase 2 results.")
        else:
            source_fingerprint = report["participants_csv_sha256"]
            size_key = "reference_sample_size" if name == "binomial_count" else "sample_size"
            if report[size_key] != len(samples[report["sample"]]):
                raise ValueError(f"Inconsistent sample count in {name}. Regenerate the Phase 2 results.")
        if source_fingerprint != fingerprint:
            raise ValueError(f"Stale Phase 2 result: {name}. Regenerate the analysis before building the reports.")
        results[name] = report
    return results


def percent(value: float) -> str:
    return f"{100 * value:.2f}"


def table_values(results: dict) -> dict:
    event = results["event_probabilities"]
    cell = event["contingency_table"]
    dependence = results["independence_bayes"]
    independence, bayes = dependence["independence"], dependence["bayes"]
    count = results["binomial_count"]
    binomial, hypergeometric = count["with_replacement"], count["without_replacement"]
    distributions = results["distributions"]
    score = distributions["summaries"][PAIR_FIELDS[0]]
    thickness = distributions["summaries"][LEFT_THICKNESS]
    covariance = results["paired_covariance"]
    score_flags, thickness_flags = score["iqr_review"], thickness["iqr_review"]
    score_mode = max(distributions["animal_score_frequency"]["levels"], key=lambda level: level["count"])
    low_tail = distributions["thickness_normal_reference"]["quantile_checks"][0]

    values = {
        "EventN": event["sample_size"], "PairedN": distributions["sample_size"],
        "NormalN": cell["A"]["total"], "OtherN": cell["not_A"]["total"],
        "UpperN": cell["total"]["B"], "LowerN": cell["total"]["not_B"],
        "MissingScoreN": results["sample_definitions"]["samples"]["paired_normal"]["first_exclusion_reason"].get("ANIMALS", 0),
        "ThicknessCut": f"{event['events']['B']['median_mm']:.3f}",
        "TiedN": event["events"]["B"]["equal_to_median"],
        "CellAB": cell["A"]["B"], "CellAnB": cell["A"]["not_B"],
        "CellnAB": cell["not_A"]["B"], "CellnAnB": cell["not_A"]["not_B"],
        "ConditionalGap": f"{event['conditional_difference_percentage_points']:.2f}",
        "ExpectedJoint": percent(independence["joint_if_independent"]["probability"]),
        "ExpectedCount": f"{independence['expected_count_if_independent']:.1f}",
        "JointGap": f"{independence['observed_minus_expected_percentage_points']:.2f}",
        "BayesLikelihood": percent(bayes["likelihood_B_given_A"]["probability"]),
        "BayesPosterior": percent(bayes["posterior_A_given_B_via_bayes"]),
        "DrawN": count["number_of_draws"], "TailMin": count["tail_minimum"],
        "DrawFraction": percent(count["number_of_draws"] / count["reference_sample_size"]),
        "CountMean": f"{binomial['expected_count']:.2f}",
        "BinomialVar": f"{binomial['variance']:.2f}",
        "HypergeomVar": f"{hypergeometric['variance']:.2f}",
        "BinomialTail": percent(binomial["probability_at_least_tail_minimum"]),
        "HypergeomTail": percent(hypergeometric["probability_at_least_tail_minimum"]),
        "TailGap": f"{count['tail_difference_percentage_points_with_minus_without_replacement']:.2f}",
        "ScoreMin": score["minimum"], "ScoreMax": score["maximum"],
        "ScoreMean": f"{score['mean']:.2f}", "ScoreMedian": f"{score['median']:.0f}",
        "ScoreFirstQuartile": f"{score['q1']:.0f}", "ScoreThirdQuartile": f"{score['q3']:.0f}",
        "ScoreEmpVar": f"{score['empirical_variance']:.3f}", "ScoreSampleVar": f"{score['sample_variance']:.3f}",
        "ScoreSkew": f"{score['empirical_skewness']:.2f}",
        "ScoreLowFlags": score_flags["below"], "ScoreHighFlags": score_flags["above"],
        "ScoreMode": score_mode["value"], "ScoreModeCount": score_mode["count"],
        "ScoreModePercent": percent(score_mode["probability"]),
        "ThicknessMin": f"{thickness['minimum']:.3f}", "ThicknessMax": f"{thickness['maximum']:.3f}",
        "ThicknessMean": f"{thickness['mean']:.3f}", "ThicknessMedian": f"{thickness['median']:.3f}",
        "ThicknessFirstQuartile": f"{thickness['q1']:.4f}", "ThicknessThirdQuartile": f"{thickness['q3']:.4f}",
        "ThicknessEmpVar": f"{thickness['empirical_variance']:.6f}",
        "ThicknessSampleVar": f"{thickness['sample_variance']:.6f}",
        "ThicknessSkew": f"{thickness['empirical_skewness']:.2f}",
        "ThicknessLowFlags": thickness_flags["below"], "ThicknessHighFlags": thickness_flags["above"],
        "FirstPercentileObserved": f"{low_tail['empirical_mm']:.3f}",
        "FirstPercentileNormal": f"{low_tail['fitted_normal_mm']:.3f}",
        "CovSample": f"{covariance['covariance']['sample']['value']:.6f}",
        "CovEmpirical": f"{covariance['covariance']['empirical']['value']:.6f}",
        "DistinctPoints": covariance["scatter"]["distinct_coordinates"],
    }
    for section, keys in (
        ("marginal_probabilities", {"A": "MargA", "not_A": "MargnA", "B": "MargB", "not_B": "MargnB"}),
        ("joint_probabilities", {"A_and_B": "JointAB", "A_and_not_B": "JointAnB",
                                 "not_A_and_B": "JointnAB", "not_A_and_not_B": "JointnAnB"}),
        ("conditional_probabilities", {"A_given_B": "ConditionAB", "A_given_not_B": "ConditionAnB"}),
    ):
        values.update({macro: percent(event[section][key]["probability"]) for key, macro in keys.items()})
    return values


def main() -> None:
    values = table_values(load_results())
    lines = ["% Generated by python3 -m src.phase2_tables."]
    lines.extend(rf"\newcommand{{\{name}}}{{{value}}}" for name, value in values.items())
    write_output(TABLE_FILE, ("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()

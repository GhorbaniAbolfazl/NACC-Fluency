"""Compare animal-naming means across the two registered MRI thickness groups."""

import json
from decimal import Decimal
from hashlib import sha256
from importlib.metadata import version
from math import isfinite, sqrt
from pathlib import Path
from statistics import fmean, variance

from scipy.stats import ttest_ind

from .build_dataset import write_output
from .descriptive_statistics import numeric_summary
from .phase2_samples import load_phase2_samples
from .phase3_mean_intervals import PLAN_FILE, SAMPLE_NAME, load_registered_plan

REPORT_FILE = Path("reports/phase3/score_difference_by_mri.json")
ESTIMAND = "score_difference_by_mri"
GROUP_ORDER = ("at_or_above_cutoff", "below_cutoff")


def load_mri_groups() -> tuple[dict, dict[str, list[dict]], dict]:
    """Preserve complete rows and validate the fixed cutoff and group counts."""
    plan = load_registered_plan()
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    specification = plan["estimands"][ESTIMAND]
    if (
        specification["sample"] != SAMPLE_NAME
        or len(rows) != plan["samples"][SAMPLE_NAME]["included"]
        or flow["participants_csv_sha256"] != plan["inputs"]["participants_csv_sha256"]
    ):
        raise ValueError("The MRI contrast sample differs from the registered Phase 3 plan.")

    grouping = plan["mri_grouping"]
    cutoff = Decimal(str(grouping["cutoff_mm"]))
    groups = {name: [] for name in GROUP_ORDER}
    for row in rows:
        thickness = Decimal(str(row[grouping["column"]]))
        group = "at_or_above_cutoff" if thickness >= cutoff else "below_cutoff"
        groups[group].append(row)
    if {name: len(members) for name, members in groups.items()} != plan["samples"][SAMPLE_NAME]["mri_groups"]:
        raise ValueError("The MRI contrast group counts differ from the registered Phase 3 plan.")
    return plan, groups, flow


def welch_contrast(higher: list[float], lower: list[float], confidence_level: float) -> dict:
    """Estimate the higher-minus-lower mean difference with a two-sided Welch test."""
    if min(len(higher), len(lower)) < 2 or not all(isfinite(v) for v in (*higher, *lower)):
        raise ValueError("A Welch contrast requires at least two finite scores per group.")
    if not 0 < confidence_level < 1:
        raise ValueError("The confidence level must be strictly between zero and one.")
    standard_error = sqrt(variance(higher) / len(higher) + variance(lower) / len(lower))
    if standard_error == 0:
        raise ValueError("A Welch contrast requires positive score variance in at least one group.")
    result = ttest_ind(higher, lower, equal_var=False, alternative="two-sided", nan_policy="raise")
    interval = result.confidence_interval(confidence_level=confidence_level)
    return {
        "mean_difference": fmean(higher) - fmean(lower),
        "standard_error": standard_error,
        "t_statistic": float(result.statistic),
        "degrees_of_freedom": float(result.df),
        "p_value_raw": float(result.pvalue),
        "confidence_interval": {
            "confidence_level": confidence_level,
            "lower": float(interval.low),
            "upper": float(interval.high),
            "interval_width": float(interval.high - interval.low),
        },
    }


def build_report() -> dict:
    plan, groups, flow = load_mri_groups()
    score_column = plan["estimands"]["mean_score"]["column"]
    scores = {name: [row[score_column] for row in groups[name]] for name in GROUP_ORDER}
    estimate = welch_contrast(scores[GROUP_ORDER[0]], scores[GROUP_ORDER[1]], plan["confidence_level"])
    grouping = plan["mri_grouping"]
    cutoff = Decimal(str(grouping["cutoff_mm"]))
    equal_to_cutoff = sum(
        Decimal(str(row[grouping["column"]])) == cutoff for row in groups["at_or_above_cutoff"]
    )
    alpha = plan["testing_policy"]["familywise_alpha"]
    return {
        "stage": "Phase 3, Step 4",
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "analysis_plan_sha256": sha256(PLAN_FILE.read_bytes()).hexdigest(),
        "sample": SAMPLE_NAME,
        "sample_size": sum(map(len, groups.values())),
        "estimand": plan["estimands"][ESTIMAND],
        "score_column": score_column,
        "mri_grouping": {**grouping, "participants_equal_to_cutoff": equal_to_cutoff},
        "group_order": list(GROUP_ORDER),
        "groups": {
            name: {"participants": len(groups[name]), "score_summary": numeric_summary(groups[name], score_column)}
            for name in GROUP_ORDER
        },
        "method": {
            "name": plan["hypothesis_tests"][ESTIMAND]["method"],
            "scipy_call": "scipy.stats.ttest_ind(higher, lower, equal_var=False, alternative='two-sided', nan_policy='raise')",
            "scipy_version": version("scipy"),
            "confidence_interval": plan["estimands"][ESTIMAND]["standard_interval"],
            "standard_error": "sqrt(s_higher^2 / n_higher + s_lower^2 / n_lower), using sample variances with n - 1 denominators",
            "degrees_of_freedom": "Welch-Satterthwaite approximation",
            "interval_scope": plan["testing_policy"]["interval_scope"],
            "assumptions": {
                "independence": "Each participant contributes to one MRI group. The sampling model treats different participants as independent.",
                "large_sample": "Group means use the large-sample Welch approximation. Individual score distributions need not be exactly normal.",
                "variances": "The method allows unequal population variances and unequal group sizes.",
                "value_rules": plan["sampling_model"]["value_rules"],
                "cutoff": "Inference is conditional on the fixed Phase 2 cutoff; uncertainty in selecting that cutoff is not propagated.",
            },
            "scope": "Unadjusted exploratory comparison among eligible participants with normal cognition",
        },
        "hypothesis_test": {
            **plan["hypothesis_tests"][ESTIMAND],
            "null_mean_difference": 0,
            "unadjusted_reference_alpha": alpha,
            "reject_unadjusted": estimate["p_value_raw"] < alpha,
        },
        "estimate": estimate,
        "multiplicity": {
            **plan["testing_policy"],
            "holm_adjusted_p_value": None,
            "reject_holm": None,
            "status": "Pending the supplementary thickness-by-status test in Step 6. Calculate Holm adjustment once both planned raw p-values are available.",
        },
        "diagnostics": {
            "group_distributions": "Group summaries include the score range, quartiles, sample SD and IQR review counts; all flagged values are retained.",
            "group_sizes": "Both group counts are checked against the saved plan before testing.",
            "pairing": "Linked UDS and MRI records belong to the same person. The two MRI groups contain different people, so the contrast uses independent groups.",
        },
        "interpretation": {
            "confidence_interval": plan["interpretation"]["confidence_interval"],
            "p_value": plan["interpretation"]["p_value"],
            "effect_units": "The mean difference and its interval are in correct animal names in 60 seconds, higher MRI group minus lower MRI group.",
            "generalization": plan["sampling_model"]["scope"],
            "scope": "The contrast describes an unadjusted association. Age, sex and education are addressed in Phase 4.",
            "practical_importance": "Interpret the effect and its interval in score units. A practical equivalence margin has not been specified.",
        },
    }


def main() -> None:
    output = json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

"""Measure Cook-distance sensitivity and bootstrap uncertainty of the reference slope."""

import json
from pathlib import Path

import numpy as np

from .build_dataset import write_output
from .phase4_design import ModelDesign, canonical_hash, prepare_designs
from .phase4_diagnostics import diagnostic_arrays, subset
from .phase4_plan import fingerprint, require
from .phase4_regression import OLSFit, coefficient_table, fit_ols, fit_summary
from .phase5_plan import REPORT_FILE as PHASE5_PLAN, load_registered_plan

REPORT_FILE = Path("reports/final/influence_bootstrap.json")
THICKNESS_TERM = "thickness_0_1mm"
METHODS = ("model_based", "hc3")


def coefficient_intervals(fit: OLSFit, confidence: float) -> dict:
    """Return estimates and usual-formula intervals without diagnostic p-values."""
    return {
        term: {"estimate": item["estimate"], **{
            method: {"standard_error": item[method]["standard_error"],
                     "confidence_interval": item[method]["confidence_interval"]}
            for method in METHODS
        }}
        for term, item in coefficient_table(fit, confidence).items()
    }


def influence_comparison(design: ModelDesign, full_fit: OLSFit, full_coefs: dict,
                         phase4: dict, plan: dict) -> dict:
    arrays = diagnostic_arrays(full_fit, phase4["diagnostics"]["review_thresholds"])
    cooks = arrays["cooks_distance"]
    threshold = arrays["limits"]["cooks_distance"]
    flagged = np.flatnonzero(arrays["masks"]["cooks_distance"])
    require(len(design.columns) == 5 and threshold == 4 / len(design.response),
            "The Cook-distance screen differs from the registered reference.")
    ranked = sorted(flagged, key=lambda i: (-cooks[i], design.participant_ids[i]))
    omitted = set(ranked[:plan["maximum_joint_omissions"]])
    keep = tuple(i for i in range(len(design.response)) if i not in omitted)
    reduced = fit_ols(subset(design, keep)) if omitted else full_fit
    reduced_coefs = coefficient_intervals(reduced, phase4["estimation"]["confidence_level"])
    selected_values = [float(cooks[i]) for i in omitted]
    return {
        "screen": plan["screen"], "ranking": plan["order"],
        "threshold_strict_gt": threshold,
        "flagged_n": len(flagged), "maximum_joint_omissions": plan["maximum_joint_omissions"],
        "omitted_n": len(omitted), "retained_n": len(keep),
        "omitted_cook_distance_range": (
            {"minimum": min(selected_values), "maximum": max(selected_values)} if omitted else None
        ),
        "temporary_fit": {"coefficients": reduced_coefs, "fit": fit_summary(reduced)},
        "coefficient_changes_temporary_minus_full": {
            term: reduced_coefs[term]["estimate"] - item["estimate"]
            for term, item in full_coefs.items()
        },
        "interpretation": "The full model supplies the ranking once. Diagnostics flag records for review. The usual model-based and HC3 intervals describe the reduced fit; nominal coverage is not guaranteed after outcome-dependent omission. The full fit remains the primary result.",
    }


def bootstrap_comparison(design: ModelDesign, full_thickness: dict, plan: dict) -> dict:
    """Draw whole participant rows with replacement from the full reference."""
    matrix = np.asarray(design.matrix)
    response = np.asarray(design.response)
    n, p = matrix.shape
    draws = plan["draws"]
    require(plan["draw_size"] == n and plan["replacement"]
            and design.participant_ids == tuple(sorted(design.participant_ids)),
            "Bootstrap sample size, replacement or row order differs from the registered plan.")
    rng = np.random.Generator(np.random.PCG64(plan["seed"]))
    estimates = np.empty(draws)
    invalid = {"rank_deficient": 0, "nonfinite": 0, "solver_failure": 0}
    index = design.columns.index(THICKNESS_TERM)

    for draw in range(draws):
        indices = rng.integers(0, n, size=n, dtype=np.int64)
        try:
            coefficients, _, rank, _ = np.linalg.lstsq(matrix[indices], response[indices], rcond=None)
        except np.linalg.LinAlgError:
            invalid["solver_failure"] += 1
            continue
        if not np.all(np.isfinite(coefficients)):
            invalid["nonfinite"] += 1
        elif rank != p:
            invalid["rank_deficient"] += 1
        else:
            estimates[draw] = coefficients[index]

    invalid_n = sum(invalid.values())
    require(invalid_n == 0,
            f"Bootstrap withheld after {draws} scheduled draws: {invalid_n} invalid fits ({invalid}).")
    settings = plan["interval"]
    require(settings["method"] == "percentile", "The registered bootstrap interval must use percentiles.")
    low, high = (float(v) for v in np.quantile(estimates, settings["quantiles"], method=settings["quantile_method"]))
    width = high - low
    se = float(np.std(estimates, ddof=1))
    comparisons = {}
    for method in METHODS:
        reference = full_thickness[method]["confidence_interval"]
        comparisons[method] = {
            "lower_endpoint_change": low - reference["lower"],
            "upper_endpoint_change": high - reference["upper"],
            "width_change": width - (reference["upper"] - reference["lower"]),
            "standard_error_change": se - full_thickness[method]["standard_error"],
        }
    return {
        "protocol": plan,
        "draws_attempted_n": draws, "valid_fits_n": draws - invalid_n,
        "invalid_fits_n": invalid_n, "invalid_fits_by_reason": invalid,
        "thickness_mean_across_draws": float(np.mean(estimates)),
        "empirical_standard_error": se,
        "confidence_interval": {"confidence_level": settings["confidence_level"],
                                "lower": low, "upper": high, "width": width},
        "bootstrap_minus_full_fit": comparisons,
        "interpretation": plan["limits"],
    }


def build_report() -> dict:
    sources = (PHASE5_PLAN, Path("src/phase5_influence.py"),
               Path("src/phase4_diagnostics.py"))
    hashes = {str(path): fingerprint(path) for path in sources}
    plan = load_registered_plan()
    phase4, designs = prepare_designs()
    design = designs["integrated"]
    full = fit_ols(design)
    full_coefs = coefficient_intervals(full, phase4["estimation"]["confidence_level"])
    require(full.participants == plan["reference"]["n"]
            and list(design.columns) == plan["reference"]["column_order"],
            "The full fit differs from the registered reference.")
    report = {
        "stage": "Phase 5, Step 5", "course_requirements": ["general_4", "general_5"],
        "source_files_sha256": hashes,
        "reference": {
            "participants": full.participants, "formula": plan["reference"]["model"],
            "membership_sha256": plan["reference"]["membership_sha256"],
            "matrix_sha256": canonical_hash(design.matrix),
            "response_sha256": canonical_hash(design.response),
            "coefficients": full_coefs, "fit": fit_summary(full),
        },
        "influence": influence_comparison(design, full, full_coefs, phase4, plan["comparisons"]["influence"]),
        "participant_bootstrap": bootstrap_comparison(design, full_coefs[THICKNESS_TERM], plan["comparisons"]["bootstrap"]),
        "privacy": "Only aggregate summaries and hashes are persisted. Individual diagnostics, participant IDs, resample indices and bootstrap replicates stay in memory.",
    }
    require(hashes == {str(path): fingerprint(path) for path in sources},
            "An input changed during the influence and bootstrap checks.")
    return report


def main() -> None:
    data = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != data:
        write_output(REPORT_FILE, data)


if __name__ == "__main__":
    main()

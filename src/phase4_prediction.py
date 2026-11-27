"""Record the Phase 4 model decision and calculate intervals at the registered profile."""

import json
from math import fsum, isclose, isfinite, sqrt
from pathlib import Path

from scipy.stats import t

from .build_dataset import write_output
from .phase4_design import FEATURE_ORDER, encode_predictors
from .phase4_diagnostics import REPORT_FILE as DIAGNOSTIC_FILE, analyze
from .phase4_plan import REPORT_FILE as PLAN_FILE, fingerprint, require
from .phase4_regression import OLSFit, coefficient_table

REPORT_FILE = Path("reports/phase4/profile_prediction.json")


def model_decision(plan: dict, diagnostics: dict) -> dict:
    """Record the judgment from the reviewed Step 6 evidence before computing intervals."""
    model = diagnostics["models"]["integrated"]
    review = diagnostics["graphical_review"]
    followup = diagnostics["integrated_influence_followup"]
    bins = model["fitted_value_bins"]
    return {
        "decision": "Retain the registered additive linear OLS specifications.",
        "extension_fitted": False,
        "registered_trigger": plan["conditional_model_extension"]["trigger"],
        "graphical_review_date": review["review_date"],
        "mean_form": {
            "evidence": review["mean_form"],
            "integrated_fitted_bin_mean_residual_range": [min(b["residual_mean"] for b in bins),
                                                          max(b["residual_mean"] for b in bins)],
            "judgment": "The central residual pattern provides insufficient evidence to add a quadratic thickness term. Sparse tails limit what can be concluded about curvature.",
        },
        "residual_spread_and_tails": {
            "evidence": [review["variance"], review["normal_reference"]],
            "integrated_fitted_bin_residual_SD_range": [min(b["residual_SD"] for b in bins),
                                                        max(b["residual_SD"] for b in bins)],
            "integrated_residual_skewness": model["residual_skewness_bias_corrected"],
            "integrated_residual_excess_kurtosis": model["residual_excess_kurtosis_bias_corrected"],
            "judgment": "Keep the score in correct-response units and report the prespecified HC3 coefficient companion. The spread and tail patterns qualify the classical interval calibration and give insufficient direction for a specific response transformation.",
        },
        "influence": {
            "flagged_union": model["review_flags"]["union_count"],
            "omission_slope_range": [followup["individual_leave_one_out"]["estimate"][key]
                                     for key in ("minimum", "maximum")],
            "joint_omission_slope": followup["joint_union_omission"]["thickness"]["estimate"],
            "judgment": "Use the full eligible sample for prediction. Interpret the omission fits as sensitivity checks because their membership depends on the observed fit.",
        },
        "scope": "A judgment about model form for this sample and course analysis. The graphical evidence leaves uncertainty about linearity, constant conditional variance and normal errors, especially in sparse regions. The primary model-based test and its HC3 companion come from the full registered fit.",
    }


def profile_intervals(fit: OLSFit, row: tuple[float, ...], confidence: float,
                      usable_range: tuple[float, float]) -> dict:
    """Use the classical mean and single-new-observation variance formulas."""
    require(len(row) == len(fit.columns) and all(map(isfinite, row)), "Invalid prediction row.")
    require(0 < confidence < 1 and fit.residual_df > 0, "Invalid interval level or residual degrees of freedom.")
    fitted_mean = fsum(value * beta for value, beta in zip(row, fit.coefficients))
    row_variance_factor = fsum(x * fit.unscaled_covariance[i][j] * z
                               for i, x in enumerate(row) for j, z in enumerate(row))
    residual_variance = fit.sse / fit.residual_df
    critical = float(t.isf((1 - confidence) / 2, fit.residual_df))
    require(all(isfinite(v) and v > 0 for v in (row_variance_factor, residual_variance, critical)),
            "Prediction requires finite positive variance factors and t critical value.")
    score_min, score_max = usable_range

    def interval(variance: float) -> dict:
        se = sqrt(variance)
        half_width = critical * se
        lower, upper = fitted_mean - half_width, fitted_mean + half_width
        return {"variance": variance, "standard_error": se, "lower": lower, "upper": upper,
                "width": upper - lower, "outside_usable_score_range": lower < score_min or upper > score_max}

    mean_variance = residual_variance * row_variance_factor
    mean_ci = interval(mean_variance)
    prediction_interval = interval(mean_variance + residual_variance)
    require(prediction_interval["lower"] < mean_ci["lower"] < fitted_mean < mean_ci["upper"] < prediction_interval["upper"],
            "The mean interval must be inside the individual prediction interval.")
    return {"point_estimate": fitted_mean, "confidence_level": confidence,
            "residual_df": fit.residual_df, "t_critical": critical,
            "residual_variance": residual_variance, "residual_SD": sqrt(residual_variance),
            "profile_variance_factor": row_variance_factor,
            "mean_confidence_interval": mean_ci, "individual_prediction_interval": prediction_interval,
            "prediction_to_mean_width_ratio": prediction_interval["width"] / mean_ci["width"],
            "point_estimate_outside_usable_score_range": not score_min <= fitted_mean <= score_max}


def build_report() -> dict:
    plan = json.loads(PLAN_FILE.read_text(encoding="utf-8"))
    saved = json.loads(DIAGNOSTIC_FILE.read_text(encoding="utf-8"))
    evidence_paths = [item["path"] for item in saved["figures"]]
    tracked = {**saved["source_files_sha256"], str(DIAGNOSTIC_FILE): fingerprint(DIAGNOSTIC_FILE),
               "src/phase4_prediction.py": fingerprint(Path("src/phase4_prediction.py")),
               **{name: fingerprint(Path(name)) for name in evidence_paths}}
    diagnostics, context = analyze()
    require(saved == diagnostics, "Saved Step 6 diagnostics differ from the current inputs or code.")
    require(plan == json.loads(PLAN_FILE.read_text(encoding="utf-8")), "The registered plan changed during prediction.")
    decision = model_decision(plan, diagnostics)
    prediction = plan["prediction"]
    profile = prediction["profile"]
    require(prediction["model"] == "integrated" and profile["jointly_observed_predictor_profile"]
            and not profile["response_used_for_selection"], "Expected the registered predictor-only integrated profile.")
    design, fit = context["designs"]["integrated"], context["fits"]["integrated"]
    encoded = dict(zip(FEATURE_ORDER, encode_predictors(profile["values"], plan)))
    row = tuple(encoded[name] for name in design.columns)
    registered = tuple(float(prediction["design_row"][name]) for name in design.columns)
    require(row == registered == design.reference_row == (1.0, 0.0, 0.0, 0.0, 0.0)
            and row in design.matrix, "The reference profile must be jointly observed and use the registered coding.")
    require(fit.columns == design.columns and fit.participants == len(design.response) == diagnostics["participants"],
            "The profile must use the full retained integrated fit in registered column order.")
    confidence = plan["estimation"]["confidence_level"]
    usable_range = tuple(plan["variables"]["response"]["usable_range"])
    intervals = profile_intervals(fit, row, confidence, usable_range)
    intercept = coefficient_table(fit, confidence)["intercept"]
    mean_ci, individual = intervals["mean_confidence_interval"], intervals["individual_prediction_interval"]
    require(isclose(intervals["point_estimate"], intercept["estimate"], rel_tol=0, abs_tol=1e-12)
            and all(isclose(mean_ci[key], intercept["model_based"]["confidence_interval"][key], rel_tol=0, abs_tol=1e-12)
                    for key in ("lower", "upper")), "The centered profile mean interval must match the classical intercept interval.")
    require(isclose(individual["variance"] - mean_ci["variance"], intervals["residual_variance"], rel_tol=0, abs_tol=1e-12),
            "The individual interval must add the residual variance.")
    require(tracked == {name: fingerprint(Path(name)) for name in tracked}, "Inputs or sources changed during prediction.")
    return {
        "stage": "Phase 4, Step 7", "status": "Model retention documented and both registered profile intervals calculated.",
        "sample": "model_normal", "participants": fit.participants, "model": fit.name,
        "source_files_sha256": tracked, "environment": diagnostics["environment"],
        "model_decision": decision,
        "profile": {"predictor_values": profile["values"], "selection_rule": profile["selection_rule"],
                    "response_used_for_selection": False, "joint_support_verified": True,
                    "column_order": list(design.columns), "design_row": list(row),
                    "target": prediction["target"], "privacy": profile["privacy"]},
        "intervals": intervals,
        "methods": {"mean_interval": prediction["mean_interval"],
                    "individual_prediction_interval": prediction["individual_prediction_interval"],
                    "usable_score_range": list(usable_range), "scale_rule": prediction["scale_rule"],
                    "fitted_scores_outside_usable_range": diagnostics["models"]["integrated"]["fitted_outside_usable_score_range"],
                    "response_unit": plan["variables"]["response"]["unit"]},
        "interpretation": {
            "mean": "The confidence interval quantifies uncertainty about the conditional mean score for this predictor profile.",
            "individual": "The prediction interval targets one new independent eligible participant with the same predictors. Its width includes the residual variation of an individual score.",
            "calibration": "Both intervals are nominal 95% classical OLS intervals. Exact calibration relies on an adequate linear mean and independent, homoskedastic normal errors, including an independent new observation from the same conditional distribution for prediction. The observed asymmetry and variation in residual spread warrant treating the limits as model-based approximations, particularly for an individual score.",
            "hc3": "HC3 remains a companion for coefficient uncertainty. The classical individual prediction interval uses the pooled residual variance as specified in the plan.",
            "scope": "An educational application to a fixed profile on observed support. Actual predictive coverage in new data remains unverified.",
        },
        "primary_test": diagnostics["primary_test_unchanged"],
        "validation": {"step_6_reproduced": True, "registered_observed_profile_verified": True,
                       "mean_interval_matches_model_based_intercept": True,
                       "individual_variance_adds_residual_variance": True, "individual_interval_contains_mean_interval": True},
        "requirement_evidence": {"general_8": "Retained model form and rationale recorded against Step 6 evidence.",
                                 "general_9": "Conditional mean and single-observation prediction intervals reported separately."},
        "progress": {"completed_through_step": 7, "remaining": "Step 8 final verification, LaTeX report and executed notebook."},
        "sources": {"prediction_intervals": "https://www.itl.nist.gov/div898/handbook/pmd/section5/pmd512.htm",
                    "prediction_targets": "https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.get_prediction.html"},
    }


def load_prediction_results() -> dict:
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_report(), "Saved profile intervals differ from current inputs or code. Review and rebuild Step 7.")
    return saved


def main() -> None:
    output = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != output:
        write_output(REPORT_FILE, output)


if __name__ == "__main__":
    main()

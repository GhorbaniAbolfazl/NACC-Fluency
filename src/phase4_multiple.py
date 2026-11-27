"""Fit the registered UDS-only and integrated OLS models on the common sample."""

import json
from importlib.metadata import version
from math import fsum, isclose
from pathlib import Path
from platform import python_version

from scipy.linalg import lstsq

from .build_dataset import write_output
from .phase4_design import REPORT_FILE as DESIGN_FILE, canonical_hash, load_prepared_designs
from .phase4_plan import REPORT_FILE as PLAN_FILE, fingerprint, require
from .phase4_regression import OLSFit, coefficient_table, fit_ols, fit_summary
from .phase4_simple import REPORT_FILE as SIMPLE_FILE, load_simple_results

REPORT_FILE = Path("reports/phase4/multiple_regression.json")
MODEL_NAMES = ("uds_only", "integrated")


def coefficient_interpretations(name: str) -> dict[str, str]:
    """Describe between-person conditional contrasts in the registered units."""
    integrated = name == "integrated"
    conditioning = "left thickness, " if integrated else ""
    result = {
        "intercept": ("Fitted mean animal-naming score for recorded males aged 68 with 16 years of education"
                      + (" and left thickness 2.344 mm." if integrated else "; thickness is not included in this model.")),
        "age_c": f"Mean-score difference per one year older age at the same {conditioning}recorded sex and years of education; not within-person aging.",
        "female": f"Mean-score difference for recorded female (SEX=2) minus recorded male (SEX=1), at the same {conditioning}age and years of education.",
        "education_c": f"Mean-score difference per one additional year of education at the same {conditioning}age and recorded sex.",
    }
    if integrated:
        result["thickness_0_1mm"] = "Mean-score difference per 0.1 mm greater left inferior parietal thickness at the same age, recorded sex and years of education; the primary adjusted association."
    return result


def analyze() -> tuple[dict, dict[str, OLSFit]]:
    """Validate earlier stages, fit both multiple models, and retain private arrays in memory."""
    simple = load_simple_results()
    plan, designs = load_prepared_designs()
    integrated = designs["integrated"]
    require(tuple(plan["models"]) == ("uds_only", "mri_only", "integrated"), "Unexpected registered model set or order.")
    primary = plan["testing_policy"]["primary"]
    require(primary["model"] == "integrated" and primary["term"] == "thickness_0_1mm",
            "The registered primary estimand has changed.")
    tracked = {**simple["source_files_sha256"], str(PLAN_FILE): fingerprint(PLAN_FILE),
               str(DESIGN_FILE): fingerprint(DESIGN_FILE), str(SIMPLE_FILE): fingerprint(SIMPLE_FILE),
               "src/phase4_multiple.py": fingerprint(Path("src/phase4_multiple.py"))}
    require(simple["participants"] == plan["sample"]["counts"]["model_normal"] == len(integrated.response)
            and simple["design"]["membership_sha256"] == plan["sample"]["membership_sha256"]
            and simple["design"]["response_sha256"] == canonical_hash(integrated.response),
            "Simple and multiple regression must use identical registered participants and responses.")
    require(simple["design"]["matrix_sha256"] == canonical_hash(designs["mri_only"].matrix),
            "The previous MRI-only matrix differs from the prepared design.")
    fits, models = {}, {}
    for name in MODEL_NAMES:
        design = designs[name]
        specification = plan["models"][name]
        require(design.participant_ids == integrated.participant_ids and design.response == integrated.response,
                "Multiple models must have exactly the same participant/response order.")
        require(design.columns == tuple(specification["columns"]), "Registered columns or their order changed.")
        fit = fit_ols(design)
        require(fit.residual_df == specification["residual_df_if_full_rank"], "Unexpected residual degrees of freedom.")
        # An SVD-based LAPACK solve independently checks the pivoted-QR coefficient order.
        oracle, _, rank, _ = lstsq(design.matrix, design.response, lapack_driver="gelsd", check_finite=True)
        require(rank == len(design.columns), "The independent SVD solution is not full rank.")
        require(all(isclose(value, float(expected), rel_tol=1e-10, abs_tol=1e-11)
                    for value, expected in zip(fit.coefficients, oracle)), "QR and SVD coefficient estimates disagree.")
        oracle_sse = fsum((y - fsum(x * float(b) for x, b in zip(row, oracle))) ** 2
                          for row, y in zip(design.matrix, design.response))
        require(isclose(fit.sse, oracle_sse, rel_tol=1e-12, abs_tol=1e-9), "QR and SVD residual sums of squares disagree.")
        coefficients = coefficient_table(fit, plan["estimation"]["confidence_level"])
        fits[name] = fit
        models[name] = {
            "role": specification["role"], "formula": specification["formula"],
            "column_order": list(design.columns), "matrix_sha256": canonical_hash(design.matrix),
            "response_sha256": canonical_hash(design.response), "membership_sha256": plan["sample"]["membership_sha256"],
            "reference_row": list(design.reference_row), "coefficients": coefficients,
            "coefficient_interpretations": coefficient_interpretations(name),
            "coefficient_test_roles": {term: ("Primary registered model-based t test; HC3 is a nominal companion."
                                               if name == primary["model"] and term == primary["term"]
                                               else "Supporting coefficient description; p-values are nominal.")
                                       for term in design.columns},
            "covariance": {"column_order": list(fit.columns),
                           "unscaled": [list(row) for row in fit.unscaled_covariance],
                           "model_based": [list(row) for row in fit.model_based_covariance],
                           "hc3": [list(row) for row in fit.hc3_covariance]},
            "fit": fit_summary(fit),
            "numerical_verification": {"svd_rank": int(rank), "qr_svd_coefficient_agreement": True,
                                       "maximum_absolute_coefficient_difference": max(abs(value - float(expected)) for value, expected in zip(fit.coefficients, oracle)),
                                       "absolute_sse_difference": abs(fit.sse - oracle_sse),
                                       "solver": "scipy.linalg.lstsq, LAPACK gelsd; separate from the estimation QR solver."},
        }
    slope = models[primary["model"]]["coefficients"][primary["term"]]
    alpha = primary["alpha"]
    decision = slope["model_based"]["p_value_two_sided_nominal"] < alpha
    intervals = {method: slope[method]["confidence_interval"] for method in ("model_based", "hc3")}
    report = {
        "stage": "Phase 4, Step 4", "status": "UDS-only and integrated models estimated; diagnostics and final Phase 4 deliverables remain pending.",
        "sample": "model_normal", "participants": len(integrated.response), "source_files_sha256": tracked,
        "environment": {"python": python_version(), "scipy": version("scipy"), "numpy_runtime": version("numpy")},
        "coding": {"explanatory_variables": plan["variables"]["explanatory"],
                   "reference_profile": plan["prediction"]["profile"]["values"],
                   "male_indicator": 0, "female_indicator": 1,
                   "response_unit": plan["variables"]["response"]["unit"],
                   "common_sample": "Identical 842 participants and response order in both multiple regressions and the completed MRI-only model; all eligible observations retained."},
        "models": models,
        "primary_test": {**primary, "parameter": plan["primary_estimand"]["parameter"],
                         "estimate": slope["estimate"], "model_based": slope["model_based"], "hc3_companion": slope["hc3"],
                         "reject_primary_null_at_alpha": decision,
                         "companion_nominal_reject_at_alpha": slope["hc3"]["p_value_two_sided_nominal"] < alpha,
                         "intervals_include_zero": {method: interval["lower"] <= 0 <= interval["upper"] for method, interval in intervals.items()},
                         "decision_rule": "Use the registered model-based test for the primary decision; report HC3 alongside it without selecting a covariance method by p-value.",
                         "interpretation": (f"The registered primary test {'rejects' if decision else 'does not reject'} a zero adjusted thickness coefficient at alpha={alpha}. "
                                            "Interpret the estimate and both intervals in correct responses per 0.1 mm. Non-rejection does not establish no association or practical equivalence; no equivalence margin has been supplied. The result remains conditional on the model and sampling assumptions, with residual and influence review pending.")},
        "methods": {"estimator": plan["estimation"]["method"], "numerics": plan["estimation"]["numerics"],
                    "model_based_covariance": plan["estimation"]["model_based_covariance"],
                    "coefficient_intervals": plan["estimation"]["coefficient_intervals"],
                    "hc3_companion": plan["estimation"]["hc3_companion"],
                    "assumptions": plan["estimation"]["assumptions"],
                    "interval_scope": "Individual coefficient intervals; not simultaneous intervals or new-observation prediction intervals.",
                    "p_value": "Under the relevant zero-coefficient null and test assumptions, the p-value is the tail probability of a t statistic at least as far from zero as observed. It is not the probability that the null is true.",
                    "numeric_tail": "A computed coefficient p-value of zero is flagged as numerical underflow, not interpreted as a literal zero probability.",
                    "metrics": {"residual_SD": "sqrt(SSE/(n-p))", "in_sample_RMSE": "sqrt(SSE/n)",
                                "R_squared": "1-SSE/SST", "adjusted_R_squared": "1-(SSE/(n-p))/(SST/(n-1))"}},
        "testing_policy": {"multiplicity": plan["testing_policy"]["multiplicity"],
                           "other_coefficients": plan["testing_policy"]["other_coefficients"],
                           "scope": "The integrated thickness test is primary. All other coefficient p-values and HC3 companions are nominal supporting summaries; no familywise guarantee across them is claimed."},
        "validation": {"passed": True, "checks": [
            "Frozen plan, validated common-sample designs and saved MRI-only results reproduce before fitting.",
            "All models retain the same participant membership and response order with unchanged predictor coding.",
            "UDS-only and integrated parameter counts and residual degrees of freedom match the registered specifications.",
            "Full-rank pivoted-QR estimates match an independent LAPACK gelsd solution and its residual sum of squares.",
            "OLS residual orthogonality, sum-of-squares decomposition and leverage trace checks pass in the shared engine.",
            "Both covariance methods use the same OLS estimates and the registered residual degrees of freedom.",
            "Input and analysis source fingerprints remain unchanged during estimation."],
                       "limit": "Computational agreement does not validate the conditional mean, residual distribution, participant independence or sampling frame."},
        "limitations": {"sampling_and_preparation": plan["sampling_model"],
                        "response_scale": plan["variables"]["response"]["limitations"],
                        "causality": plan["primary_estimand"]["interpretation_limit"],
                        "functional_form": plan["primary_estimand"]["nonlinear_limit"],
                        "uncertainty": "Classical t calibration requires its error assumptions. HC3 t calibration is approximate and does not correct dependence, nonlinearity, selection or source-imputation uncertainty.",
                        "fit": "All fit metrics are in-sample summaries, not externally validated predictive accuracy or a basis for choosing a causal model.",
                        "no_automatic_exclusions": True, "no_equivalence_margin": True},
        "privacy": {"persisted": "Aggregate coefficients, small covariance matrices, fit metrics, numerical checks and source hashes only.",
                    "not_exported": ["Participant IDs", "Individual design rows", "Response vector", "Fitted-value vector", "Residual vector", "Leverage vector"],
                    "in_memory_api": "analyze() returns the aggregate report and private OLSFit objects for these two models."},
        "progress": {"completed": ["Step 1 frozen specification", "Step 2 common-sample designs", "Step 3 MRI-only OLS and Pearson correlation", "Step 4 UDS-only and integrated OLS with both covariance methods and conditional coefficient interpretations"],
                     "remaining": ["Step 5 common-sample model/coefficient comparisons and collinearity review", "Step 6 residual and influence diagnostics", "Step 7 documented model retention/extension and profile intervals", "Step 8 final numerical review, LaTeX report and executed notebook"],
                     "historical_records": "The frozen plan and earlier step results retain their original stage status; this report records Step 4 progress."},
        "sources": {"qr": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html",
                    "independent_svd_solver": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html",
                    "hc3": "https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.HC3_se.html"},
    }
    require(tracked == {name: fingerprint(Path(name)) for name in tracked}, "Inputs or sources changed during multiple-regression estimation.")
    return report, fits


def build_report() -> dict:
    return analyze()[0]


def load_multiple_results() -> dict:
    """Reject stale or edited results rather than carrying them into later comparisons."""
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_report(), "The saved multiple-regression result differs from the current plan, design or code. Review and rebuild Step 4.")
    return saved


def main() -> None:
    output = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if REPORT_FILE.exists() and REPORT_FILE.read_bytes() == output:
        return
    write_output(REPORT_FILE, output)


if __name__ == "__main__":
    main()

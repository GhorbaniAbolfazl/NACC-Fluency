"""Record the Phase 4 model specification before fitting any regression."""

import csv
import json
from collections import Counter
from decimal import Decimal, localcontext
from hashlib import sha256
from io import StringIO
from pathlib import Path
from statistics import median

from .build_dataset import (
    LEFT_THICKNESS, MODEL_FIELDS, NORMAL_COGNITION, OUTPUT_COLUMNS,
    STATUS, build_dataset, write_output,
)
from .data_audit import PARTICIPANT_ID
from .descriptive_statistics import FLOW_FILE, PARTICIPANTS_FILE

REPORT_FILE = Path("reports/phase4/analysis_plan.json")
GUIDE_FILE = Path("docs/amarpr.pdf")
PHASE3_REVIEW = Path("reports/phase3/inference_review.json")
EXPECTED_GUIDE_HASH = "ede1706ab440f1ae737a71afaf0a37438112e7b010e6a5118b8f334b4c0d03e3"
EXPECTED_COUNTS = {"matched": 1237, "normal_before_filters": 847, "paired_normal": 843, "model_normal": 842}
CONTINUOUS_PREDICTORS = ("NACCAGE", "EDUC", LEFT_THICKNESS)
CONFIDENCE_LEVEL = 0.95
ALPHA = 0.05


def fingerprint(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_inputs() -> tuple[list[dict], dict, dict[str, str]]:
    """Rebuild eligibility in memory and verify the completed Phase 3 provenance."""
    paths = (
        Path("data/udslong_imputed.csv"), Path("data/mriscanclean.csv"),
        PARTICIPANTS_FILE, FLOW_FILE, GUIDE_FILE, PHASE3_REVIEW,
        Path("src/build_dataset.py"), Path("src/data_audit.py"),
        Path("src/descriptive_statistics.py"), Path("src/phase4_plan.py"),
    )
    before = {str(path): fingerprint(path) for path in paths}
    require(before[str(GUIDE_FILE)] == EXPECTED_GUIDE_HASH, "The course guide changed; review the Phase 4 requirements.")
    review = json.loads(PHASE3_REVIEW.read_text(encoding="utf-8"))
    require(review["passed_checks"] == 11 and len(review["checks"]) == 11
            and all(item["status"] == "pass" for item in review["checks"]),
            "A successful completed Phase 3 review is required.")
    require(all(fingerprint(Path(name)) == digest for name, digest in review["source_files_sha256"].items()),
            "The Phase 3 review is stale; resolve changed inputs before planning Phase 4.")

    rows, flow = build_dataset(Path("data"))
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    rebuilt = buffer.getvalue().encode("utf-8")
    require(rebuilt == PARTICIPANTS_FILE.read_bytes(), "The participant table differs from its raw-data rebuild.")
    flow["participants_csv_sha256"] = sha256(rebuilt).hexdigest()
    require(flow == json.loads(FLOW_FILE.read_text(encoding="utf-8")), "The sample flow does not reproduce.")
    observed = {
        "matched": len(rows), "normal_before_filters": flow["matched_normal_cognition_before_eligibility"],
        "paired_normal": flow["samples"]["paired_normal"]["included"],
        "model_normal": flow["samples"]["model_normal"]["included"],
    }
    require(observed == EXPECTED_COUNTS, "Sample counts changed; review the Phase 4 plan rather than silently changing its population.")
    ids = [row[PARTICIPANT_ID] for row in rows]
    require(len(ids) == len(set(ids)) and all(ids), "One nonblank participant ID per row is required.")
    members = [row for row in rows if row["eligible_model_normal"]]
    require(len(members) == observed["model_normal"] and all(
        row[STATUS] == NORMAL_COGNITION and all(row[column] is not None for column in MODEL_FIELDS)
        for row in members
    ), "The common modeling sample does not match the Phase 1 eligibility rules.")
    require(set(row["SEX"] for row in members) == {1, 2}, "Both recorded sex categories are required for the planned indicator.")
    require(before == {str(path): fingerprint(path) for path in paths}, "Inputs changed while the plan was being checked.")
    return members, flow, before


def reference_profile(rows: list[dict]) -> dict:
    """Choose a central observed predictor profile without inspecting its response."""
    with localcontext() as context:
        context.prec = 50
        values = {name: [Decimal(str(row[name])) for row in rows] for name in CONTINUOUS_PREDICTORS}
        medians = {name: median(items) for name, items in values.items()}
        means = {name: sum(items) / len(items) for name, items in values.items()}
        variances = {
            name: sum((value - means[name]) ** 2 for value in items) / (len(items) - 1)
            for name, items in values.items()
        }
        require(all(value > 0 for value in variances.values()), "Every continuous predictor must vary.")

        def distance(row: dict) -> Decimal:
            return sum((Decimal(str(row[name])) - medians[name]) ** 2 / variances[name]
                       for name in CONTINUOUS_PREDICTORS)

        chosen = min((row for row in rows if row["SEX"] == 1), key=lambda row: (
            distance(row), *(Decimal(str(row[name])) for name in CONTINUOUS_PREDICTORS)
        ))
        return {
            "selection_rule": (
                "Among participants with recorded SEX=1, minimize the sum of squared deviations "
                "from the common sample's predictor medians, divided by the sample variance of each "
                "predictor (denominator n-1). Use age, education and left thickness only. "
                "Break distance ties by ascending age, education, then thickness; identical predictor "
                "tuples define the same profile. Decimal arithmetic uses 50-digit precision."
            ),
            "response_used_for_selection": False,
            "jointly_observed_predictor_profile": True,
            "values": {"SEX": 1, "NACCAGE": float(chosen["NACCAGE"]),
                       "EDUC": chosen["EDUC"], LEFT_THICKNESS: chosen[LEFT_THICKNESS]},
            "sample_predictor_medians": {name: float(value) for name, value in medians.items()},
            "sample_predictor_variances": {name: float(value) for name, value in variances.items()},
            "standardized_squared_distance": float(distance(chosen)),
            "privacy": "Record predictor values only; omit participant ID, row index and observed score.",
        }


def build_plan() -> dict:
    rows, flow, fingerprints = validate_inputs()
    profile = reference_profile(rows)
    reference = profile["values"]
    n = len(rows)
    sex_counts = Counter(row["SEX"] for row in rows)
    model_specs = {
        "uds_only": {"role": "UDS baseline", "columns": ["intercept", "age_c", "female", "education_c"]},
        "mri_only": {"role": "Unadjusted simple MRI regression", "columns": ["intercept", "thickness_0_1mm"]},
        "integrated": {"role": "Primary adjusted association", "columns": ["intercept", "thickness_0_1mm", "age_c", "female", "education_c"]},
    }
    for model in model_specs.values():
        model.update({"response": "ANIMALS", "sample": "model_normal", "participants": n,
                      "formula": "ANIMALS ~ " + " + ".join("1" if term == "intercept" else term for term in model["columns"]),
                      "parameters_including_intercept": len(model["columns"]),
                      "residual_df_if_full_rank": n - len(model["columns"])})
    plan = {
        "plan_version": 1,
        "stage": "Phase 4, Step 1",
        "status": "Specification recorded before Phase 4 regression fitting; no fitted coefficients, association tests or predictions are produced by this step.",
        "analysis_character": "Exploratory secondary observational analysis informed by the Phase 1-3 descriptions and unadjusted comparisons; this is a pre-fit project plan, not an external preregistration or independent confirmation.",
        "inputs": {"source_files_sha256": fingerprints, "course_guide_pages": [4, 13, 14, 15, 18, 19]},
        "question": "Among participants with normal cognition, how does the mean number of correct animal names in 60 seconds differ per 0.1 mm greater left inferior parietal cortical thickness, after accounting for age, recorded sex and years of education?",
        "sampling_model": {
            "unit": "One participant with the single supplied UDS record and the single supplied MRI record, joined one-to-one on NACCID.",
            "linkage": flow["linkage"],
            "target": "The source population represented by this teaching extract, conditional on model_normal eligibility and the sampling/model assumptions; no claim of representativeness for all NACC participants.",
            "independence": "One supplied row per ID is verified; independence across people remains a modeling assumption.",
            "selection": "Random selection is reported, but its sampling frame, selection probabilities and source selection eligibility are undocumented.",
            "source_preparation": "The source preparation and imputation methods and value-level masks are unknown. Treat the supplied values as fixed inputs; this analysis cannot propagate additional imputation/preparation uncertainty.",
            "visit_and_scan_timing": {
                "record_rule": "Use the single supplied record per participant; no within-person visit selection is possible or performed.",
                "clinical_visit_date": "Unavailable", "clinical_visit_number": "Unavailable",
                "mri_scan_date": "Unavailable", "source_visit_count": "Unavailable",
                "temporal_proximity": "Unverifiable from the supplied files; do not assert baseline, first visit, nearest scan or contemporaneous measurements.",
                "track_c5_status": "Participant-level input structure is documented. Temporal proximity remains unverified and would require authoritative metadata or course-team clarification.",
            },
            "scan_quality": "Finite positive thickness establishes numerical usability, not scan-level quality control.",
        },
        "sample": {
            "name": "model_normal", "eligibility_flag": "eligible_model_normal", "common_to_all_models": True,
            "required_columns": list(MODEL_FIELDS), "uds_cognitive_status": {"column": STATUS, "code": 1},
            "counts": EXPECTED_COUNTS.copy(),
            "exclusions_within_normal_cognition": {"unusable_score": 4, "additional_unknown_education": 1},
            "phase_one_definition": flow["samples"]["model_normal"],
            "sex_counts": {"male_code_1": sex_counts[1], "female_code_2": sex_counts[2]},
            "membership_sha256": sha256(("\n".join(sorted(row[PARTICIPANT_ID] for row in rows)) + "\n").encode("utf-8")).hexdigest(),
            "membership_hash_rule": "Sorted eligible NACCID values, UTF-8, one per line with a final newline; IDs themselves are not exported.",
            "rules": ["Reuse Phase 1 usability rules and special-code handling; no additional imputation or rounding.",
                      "Use exactly the common 842-person sample even for models needing fewer predictors.",
                      "Retain IQR-flagged measurements and all valid model-eligible participants.",
                      "Do not use cognitive status as a predictor within the normal-cognition-only primary sample."],
        },
        "variables": {
            "response": {"column": "ANIMALS", "source": "UDS", "unit": "Correct animal names in 60 seconds",
                         "usable_range": [0, 77], "scale": "Discrete nonnegative count; larger scores indicate better performance on this task.",
                         "special_codes_excluded": [-4, 95, 96, 97, 98],
                         "limitations": "A task score is not a diagnosis or a calibrated clinical change. The count is analyzed with OLS as a conditional-mean approximation; inspect residuals and predictions against its scale."},
            "explanatory": {
                "age_c": {"source": "UDS", "column": "NACCAGE", "expression": f"NACCAGE - {reference['NACCAGE']:g}", "reference_years": reference["NACCAGE"], "unit": "One year", "rationale": "Prespecified demographic adjustment related to cognitive performance and the thickness association."},
                "female": {"source": "UDS", "column": "SEX", "code_map": {"1": 0, "2": 1}, "reference": "Recorded male, SEX=1", "unit": "Recorded female minus recorded male, holding other included predictors fixed", "rationale": "Prespecified demographic adjustment; use the recorded binary codes without interpreting them as a broader gender measure."},
                "education_c": {"source": "UDS", "column": "EDUC", "expression": f"EDUC - {reference['EDUC']}", "reference_years": reference["EDUC"], "unit": "One year", "rationale": "Prespecified adjustment for years of education in a cognitive-task analysis."},
                "thickness_0_1mm": {"source": "MRI", "column": LEFT_THICKNESS,
                                    "expression": f"({LEFT_THICKNESS} - {reference[LEFT_THICKNESS]}) / 0.1",
                                    "reference_mm": reference[LEFT_THICKNESS], "unit": "0.1 mm greater thickness",
                                    "rationale": "The single anatomical measure selected by the research question before Phase 4 fitting.",
                                    "continuous": True, "phase_three_cutoff_used": False},
            },
            "centering": "Center continuous predictors at the pre-fit observed reference profile to make intercepts interpretable. Centering changes the intercept but not slopes, residuals or R-squared for these intercept-containing additive models.",
            "not_selected": "No new MRI regions, adjacency-derived predictors, cognitive scores, interactions or covariates are selected by correlation, p-values or model fit. VEG, bilateral thickness and full-sample analyses remain candidates for a separately specified Phase 5 plan.",
        },
        "models": model_specs,
        "primary_estimand": {"model": "integrated", "term": "thickness_0_1mm",
                             "parameter": "Adjusted linear-model mean score difference per 0.1 mm greater left thickness at the same age, recorded sex and education.",
                             "interpretation_limit": "An observational association under the stated functional form and assumptions; neither a within-person change nor a causal effect.",
                             "nonlinear_limit": "If the additive linear mean is inadequate, qualify this constant-slope summary and report any diagnostic extension separately."},
        "estimation": {
            "method": "Ordinary least squares with an intercept, equal participant weights and no regularization.",
            "numerics": "Use a numerically stable QR/SVD least-squares solution; check full column rank before interpretation, rather than solving through an explicitly inverted normal equation.",
            "model_based_covariance": "s^2 (X'X)^(-1), with s^2 = SSE/(n-p) and p including the intercept.",
            "coefficient_intervals": "Individual two-sided 95% Student t intervals, estimate +/- t(0.975,n-p) SE. Exact t calibration requires the classical independent, homoskedastic normal-error linear model; otherwise qualifications are required.",
            "hc3_companion": {
                "planned": True, "models": list(model_specs),
                "role": "Prespecified companion coefficient uncertainty, reported alongside the model-based results regardless of diagnostic significance; never choose whichever gives a smaller p-value.",
                "covariance": "(X'X)^(-1) X' diag(e_i^2/(1-h_ii)^2) X (X'X)^(-1)",
                "reference_distribution": "Use t with n-p degrees of freedom as an explicit finite-sample convention; HC3 interval/test calibration is approximate, not exact normal-theory t inference.",
                "limit": "HC3 addresses coefficient covariance under heteroskedasticity; it does not correct nonlinearity, dependence, selection bias or source imputation uncertainty and does not supply a heteroskedastic individual prediction interval.",
            },
            "assumptions": ["Independent participant observations under the sampling model.", "An adequate additive linear conditional-mean specification over the observed predictor support.", "A full-rank design matrix and finite error variance.", "Constant conditional error variance for model-based SEs and classical individual prediction intervals; normal errors for exact t-based calibration."],
            "confidence_level": CONFIDENCE_LEVEL,
        },
        "testing_policy": {
            "primary": {"model": "integrated", "term": "thickness_0_1mm", "null": "beta_thickness = 0", "alternative": "beta_thickness != 0", "alpha": ALPHA, "method": "Two-sided model-based OLS t test with n-p residual degrees of freedom."},
            "secondary_simple_slope": {"model": "mri_only", "term": "thickness_0_1mm", "null": "Unadjusted beta_thickness = 0", "alternative": "Unadjusted beta_thickness != 0", "reporting": "Report the two-sided nominal t-test p-value and 95% interval for the course simple-regression exercise; this is secondary to the adjusted estimand."},
            "other_coefficients": "Report estimates, SEs, t statistics, nominal p-values and individual intervals as supporting model descriptions. Do not make a collection of unadjusted significance claims.",
            "multiplicity": "One primary Phase 4 coefficient test at alpha 0.05. Secondary and companion p-values are explicitly nominal; no familywise guarantee across coefficients, methods or phases is claimed. The separate completed Phase 3 Holm family remains unchanged.",
            "interpretation": "Report effect size and interval width in correct responses; p-values are conditional on the null and assumptions, not probabilities that hypotheses are true. Non-rejection does not establish equivalence. No practical-equivalence margin has been supplied.",
        },
        "model_comparison": {
            "common_sample_required": True, "metrics": ["participants", "parameters_including_intercept", "residual_df", "SSE", "R_squared", "adjusted_R_squared", "residual_SD", "in_sample_RMSE"],
            "residual_SD": "sqrt(SSE/(n-p))", "in_sample_RMSE": "sqrt(SSE/n); label its different denominator explicitly.",
            "coefficient_comparisons": "Compare the same thickness slope in mri_only versus integrated and the same demographic coefficients in uds_only versus integrated, including their uncertainty.",
            "scope": "Training-sample association/fit comparison, not validated predictive accuracy. No train/test optimization, best-model search, automatic selection or claim of external prediction performance.",
            "correlation": "Report Pearson r between score and continuous thickness on these same 842 participants. With an intercept, verify R_squared = r^2 and the corresponding simple-slope scaling identity; do not reuse a coefficient from the 843-person Phase 2 sample as if its sample were identical.",
        },
        "diagnostics": {
            "all_models": ["Residuals versus fitted values", "Normal Q-Q plot", "Scale-location plot", "Residual spread and fitted-value range"],
            "integrated_detail": ["Residuals versus each continuous predictor", "Internally studentized residuals", "Hat-matrix leverage", "Cook's distance", "Joint review of flagged observations and their influence on the thickness association"],
            "review_thresholds": {"internally_studentized_residual_abs_gt": 3, "leverage_gt": "2*p/n", "cooks_distance_gt": "4/n", "meaning": "Heuristic review flags, not tests, proof of error or automatic exclusion rules; p includes the intercept."},
            "collinearity": {"models": ["uds_only", "integrated"], "checks": ["Predictor correlation matrix including the sex indicator", "VIF for each non-intercept predictor", "Full-rank check and standardized-design condition diagnostics"], "vif_review_gt": 5, "response": "Describe uncertainty and coefficient stability; no automatic covariate removal. The MRI-only model has only one explanatory variable, so between-MRI-predictor collinearity is inapplicable."},
            "influence_followup": "Retain all eligible participants in the primary fit. If flags require follow-up, document their combined and individual effect on the thickness slope and uncertainty using clearly labeled leave-one-out diagnostic fits; report aggregate influence summaries without exporting participant IDs. Do not replace the primary sample with a cleaned sample.",
            "limits": "Graphs and diagnostics cannot prove independence or validate the unknown sampling design. Do not declare the model valid from R-squared or one normality-test p-value alone.",
        },
        "conditional_model_extension": {
            "required_only_if_justified": True,
            "trigger": "A clearly documented pattern of curvature or variance misspecification in the prespecified diagnostic review. Record the graphical evidence and intended action before fitting an extension; a smaller scientific p-value is not a trigger.",
            "curvature_option": "One exploratory integrated extension adding thickness_0_1mm squared, with the same sample, main effects and reference profile. Compare its mean curve, residuals, interval widths and fit with the retained linear reference; interpret thickness effects as depending on thickness in that extension.",
            "variance_option": "Use the prespecified HC3 companion to assess coefficient uncertainty. Any response transformation requires a separately recorded substantive and diagnostic rationale, zero/bound handling and back-transformation interpretation before fitting; no automatic log transform or transform search is planned.",
            "audit_rule": "Always retain and report the original three model specifications. Label extensions as diagnostic/data-informed and do not silently redefine the primary test or use them to manufacture significance. If no change is justified, record why the original form was retained.",
        },
        "prediction": {
            "model": "integrated", "profile": profile,
            "design_row": {"intercept": 1, "thickness_0_1mm": 0, "age_c": 0, "female": 0, "education_c": 0},
            "target": "The conditional mean at this prespecified profile and a new independent eligible participant with the same predictor values; the profile-selection participant's response is not used as a validation observation.",
            "mean_interval": "yhat0 +/- t(0.975,n-p) * s * sqrt(x0' (X'X)^(-1) x0)",
            "individual_prediction_interval": "yhat0 +/- t(0.975,n-p) * s * sqrt(1 + x0' (X'X)^(-1) x0)",
            "coverage": "Individual nominal 95% classical model intervals, conditional on an adequate mean form, independent errors and constant conditional residual variance; normal errors give exact calibration. A new individual's interval includes residual variation and is wider than the mean interval.",
            "scope": "An educational model-based prediction illustration, not external validation or a clinical prediction tool. HC3 mean uncertainty does not by itself validate the individual interval under heteroskedasticity.",
            "scale_rule": "Report the raw OLS limits without rounding scores to integers or clipping limits to 0-77. Note any limits or fitted values outside the usable score scale and qualify their interpretation.",
        },
        "requirement_map": {
            "general_1": "Question, response and prespecified explanatory variables",
            "general_2": "Step 3: continuous-thickness scatterplot, least-squares line, slope and centered intercept",
            "general_3": "Step 3: simple-regression coefficient interval/test and association-versus-causation interpretation",
            "general_4": "Step 4: uds_only and integrated OLS, conditional coefficient interpretations",
            "general_5": "Step 5: common-sample coefficient and fit comparisons",
            "general_6": "Recorded sex indicator: male=0, female=1; reference and conditional interpretation fixed",
            "general_7": "Step 6: prespecified residual and influence diagnostics",
            "general_8": "Step 7: justified extension only if diagnostic evidence requires it; otherwise document retention",
            "general_9": "Step 7: prespecified predictor profile with mean and new-observation intervals",
            "general_10": "Step 3: same-sample Pearson correlation and simple-regression identity",
            "general_11": "Step 5: VIF, rank and predictor dependence review",
            "track_c1": "ANIMALS response scale, direction, usability and limitations",
            "track_c2": "uds_only model fixed before fitting",
            "track_c3": "mri_only model with one prespecified continuous MRI measure",
            "track_c4": "integrated model and three-model common-sample comparisons",
            "track_c5": "One supplied matched record per ID verified; visit selection documented as not performed; temporal proximity remains unverified",
            "track_c6": "Integrated residual, variance, outlier and influence review without automatic deletion",
            "status": "Coverage specified, not yet implemented or passed; the temporal-proximity limitation is unresolved by the available data.",
        },
        "work_sequence": [
            {"step": 1, "task": "Register the pre-fit specification", "status": "Specified here; no regression fitted"},
            {"step": 2, "task": "Prepare and verify common-sample model matrices and transformations", "status": "Not started"},
            {"step": 3, "task": "Estimate Pearson correlation and MRI-only simple regression", "status": "Not started"},
            {"step": 4, "task": "Fit UDS-only and integrated regression", "status": "Not started"},
            {"step": 5, "task": "Compare models and assess collinearity", "status": "Not started"},
            {"step": 6, "task": "Review residuals and influential observations", "status": "Not started"},
            {"step": 7, "task": "Document justified model checks and profile intervals", "status": "Not started"},
            {"step": 8, "task": "Validate calculations, interpretations and course deliverables", "status": "Not started"},
        ],
        "delivery": {"report": "reports/phase4.pdf", "source_directory": "reports/phase4", "code_directory": "src", "notebook": "notebooks/phase4.ipynb", "figures": "figures/phase4", "tag": "phase-4", "status": "Only the analysis plan exists at Step 1; the report, notebook, figures and tag are later deliverables.", "rules": "Use LaTeX, grayscale figures, reproducible result tables and an English public report/notebook/documentation. Git publication, tagging and CW packaging are manual. Do not commit raw or derived participant-level CSVs."},
        "sources": {
            "course_guide": str(GUIDE_FILE),
            "model_diagnostics": "https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd44.htm",
            "hc3_definition": "https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.HC3_se.html",
            "prediction_mean_and_new_observation": "https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.get_prediction.html",
            "vif_definition": "https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.variance_inflation_factor.html",
            "implementation_note": "Method references do not add a statsmodels dependency at the planning stage; numerical implementation and dependency choices are checked during model-matrix preparation.",
        },
    }
    require(fingerprints == {path: fingerprint(Path(path)) for path in fingerprints}, "Inputs changed during plan construction.")
    return plan


def load_registered_plan() -> dict:
    """Reject a stale or edited specification before later Phase 4 work."""
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_plan(), "The saved Phase 4 specification differs from its inputs or code. Review and version the plan before continuing.")
    return saved


def main() -> None:
    plan = build_plan()
    if REPORT_FILE.exists():
        require(json.loads(REPORT_FILE.read_text(encoding="utf-8")) == plan,
                "An existing Phase 4 plan differs. Review and version it explicitly; this command will not overwrite a changed plan.")
        return
    write_output(REPORT_FILE, (json.dumps(plan, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()

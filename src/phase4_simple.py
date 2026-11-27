"""Estimate the prespecified MRI-only regression and same-sample Pearson r."""

import json
from decimal import Decimal
from importlib.metadata import version
from math import isclose, sqrt
from pathlib import Path
from platform import python_version
from statistics import stdev

from scipy.stats import pearsonr

from .build_dataset import LEFT_THICKNESS, write_output
from .phase4_design import ModelDesign, REPORT_FILE as DESIGN_FILE, canonical_hash, load_prepared_designs
from .phase4_plan import REPORT_FILE as PLAN_FILE, fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary

REPORT_FILE = Path("reports/phase4/simple_regression.json")
FIGURE_DIR = Path("figures/phase4")
FIGURE_NAME = "01_score_thickness_simple_regression.pdf"
STYLE_FILES = (Path("src/phase1_figures.py"),)


def analyze() -> tuple[dict, ModelDesign]:
    """Fit only MRI-only; demographic and integrated models remain for Step 4."""
    plan, designs = load_prepared_designs()
    design = designs["mri_only"]
    require(design.columns == ("intercept", "thickness_0_1mm"), "The simple model differs from its registered columns.")
    tracked = {**plan["inputs"]["source_files_sha256"], str(PLAN_FILE): fingerprint(PLAN_FILE),
               str(DESIGN_FILE): fingerprint(DESIGN_FILE)}
    for name in ("src/phase4_design.py", "src/phase4_regression.py", "src/phase4_simple.py", "requirements.txt"):
        tracked[name] = fingerprint(Path(name))
    tracked.update({str(path): fingerprint(path) for path in STYLE_FILES})
    fit = fit_ols(design)
    coefficients = coefficient_table(fit, plan["estimation"]["confidence_level"])
    x = tuple(row[1] for row in design.matrix)
    reference = plan["variables"]["explanatory"]["thickness_0_1mm"]["reference_mm"]
    # Reconstruct the displayed mm scale from the same registered design, without rejoining rows.
    thickness = tuple(float(Decimal(str(value)) * Decimal("0.1") + Decimal(str(reference))) for value in x)
    correlation = pearsonr(thickness, design.response, alternative="two-sided")
    r, correlation_p = float(correlation.statistic), float(correlation.pvalue)
    slope = coefficients["thickness_0_1mm"]
    expected_slope = r * stdev(design.response) / stdev(x)
    require(isclose(fit.r_squared, r * r, rel_tol=1e-9, abs_tol=1e-12), "R-squared does not match squared Pearson r.")
    require(isclose(slope["estimate"], expected_slope, rel_tol=1e-10, abs_tol=1e-12), "The simple slope does not match Pearson r * s_y / s_x.")
    expected_t = r * sqrt(fit.residual_df / (1 - r * r))
    require(isclose(slope["model_based"]["t_statistic"], expected_t, rel_tol=1e-10, abs_tol=1e-12), "The simple-slope t statistic does not match its correlation identity.")
    require(isclose(slope["model_based"]["p_value_two_sided_nominal"], correlation_p, rel_tol=1e-9, abs_tol=1e-12), "Classical slope and Pearson nominal p-values disagree.")
    require(tracked == {name: fingerprint(Path(name)) for name in tracked}, "Inputs changed during simple-regression estimation.")
    classical_interval = slope["model_based"]["confidence_interval"]
    hc3_interval = slope["hc3"]["confidence_interval"]
    classical_contains_zero = classical_interval["lower"] <= 0 <= classical_interval["upper"]
    hc3_contains_zero = hc3_interval["lower"] <= 0 <= hc3_interval["upper"]
    rejects = slope["model_based"]["p_value_two_sided_nominal"] < plan["testing_policy"]["primary"]["alpha"]
    report = {
        "stage": "Phase 4, Step 3", "status": "MRI-only simple regression and same-sample Pearson correlation estimated; adjusted primary model not yet fitted.",
        "sample": "model_normal", "participants": fit.participants, "model": "mri_only",
        "source_files_sha256": tracked,
        "environment": {"python": python_version(), "scipy": version("scipy"), "numpy_runtime": version("numpy"), "matplotlib": version("matplotlib")},
        "design": {"columns": list(design.columns), "matrix_sha256": canonical_hash(design.matrix),
                   "response_sha256": canonical_hash(design.response), "membership_sha256": plan["sample"]["membership_sha256"],
                   "thickness_reference_mm": reference, "expression": plan["variables"]["explanatory"]["thickness_0_1mm"]["expression"],
                   "response_unit": "Correct animal names in 60 seconds", "slope_unit": "Correct responses per 0.1 mm greater thickness"},
        "coefficients": coefficients,
        "covariance": {"column_order": list(fit.columns), "unscaled": [list(row) for row in fit.unscaled_covariance],
                       "model_based": [list(row) for row in fit.model_based_covariance],
                       "hc3": [list(row) for row in fit.hc3_covariance]},
        "fit": fit_summary(fit),
        "correlation": {"method": "scipy.stats.pearsonr, two-sided, analytic default", "participants": fit.participants,
                        "coefficient": r, "squared_coefficient": r * r, "p_value_two_sided_nominal": correlation_p,
                        "null": "Population Pearson correlation is zero", "alternative": "Population Pearson correlation is nonzero",
                        "calibration": "The analytic Pearson p-value has a bivariate-normal reference; non-normality and the undocumented sampling design limit its interpretation. The same-sample classical slope t test has the equivalent numeric p-value under its own linear-model assumptions.",
                        "scope": "Recomputed on the 842-person common modeling sample; no reuse of the earlier 843-person Phase 2 result."},
        "hypothesis_test": {**plan["testing_policy"]["secondary_simple_slope"], "nominal_alpha": plan["testing_policy"]["primary"]["alpha"],
                            "model_based_p": slope["model_based"]["p_value_two_sided_nominal"],
                            "reject_at_nominal_alpha": slope["model_based"]["p_value_two_sided_nominal"] < plan["testing_policy"]["primary"]["alpha"],
                            "multiplicity": plan["testing_policy"]["multiplicity"],
                            "relation_to_correlation": "The Pearson and classical slope p-values are algebraically linked for this one-predictor intercept model; they are not two independent confirmations."},
        "methods": {"estimator": "Equal-weight OLS using economic pivoted QR and triangular solves, with coefficients returned to registered order; no explicit inversion of X'X.",
                    "model_based_covariance": plan["estimation"]["model_based_covariance"],
                    "coefficient_intervals": plan["estimation"]["coefficient_intervals"],
                    "hc3_companion": plan["estimation"]["hc3_companion"],
                    "interval_scope": "Individual 95% intervals, not simultaneous or prediction intervals.",
                    "numeric_tail_note": "A zero computed coefficient p-value is flagged as numerical underflow, not interpreted as a literal zero probability.",
                    "fit_metrics": {"residual_SD": "sqrt(SSE/(n-p))", "in_sample_RMSE": "sqrt(SSE/n)"}},
        "numeric_checks": {"passed": True, "identities": ["R_squared equals Pearson r squared", "Slope equals r * sample_SD(score) / sample_SD(centered thickness per 0.1 mm)",
                              "Classical slope t equals r * sqrt((n-2)/(1-r^2))", "Classical slope p matches Pearson p",
                              "SST equals SSE plus explained SS", "Residuals sum to zero and are orthogonal to each model column", "The hat-matrix trace equals the model rank"],
                           "limits": "These checks establish computational consistency, not distributional validity or independent scientific replication."},
        "interpretation": {
            "slope": "The slope is an unadjusted between-person mean-score association per 0.1 mm greater left thickness within the eligible normal-cognition sample; it is not a within-person change or causal effect.",
            "intercept": f"The centered intercept is the fitted mean score at {reference} mm, within the observed thickness range. This MRI-only model does not condition on age, recorded sex or education.",
            "evidence": (f"The model-based slope interval {'includes' if classical_contains_zero else 'excludes'} zero; "
                         f"the HC3 companion interval {'includes' if hc3_contains_zero else 'excludes'} zero. "
                         f"The nominal model-based test {'rejects' if rejects else 'does not reject'} zero slope at the registered reference level. "
                         "Read the effect and interval widths in score units. Non-rejection does not establish absence of an association or practical equivalence."),
            "p_value": "Under the zero-slope null and the classical linear-model assumptions, the nominal p-value is the tail probability of a t statistic at least as far from zero as observed. It is not the probability that the null is true.",
            "fit": "R-squared describes in-sample squared score variation accounted for by this linear thickness model; it is not external predictive accuracy or a causal fraction.",
            "hc3": "The companion changes coefficient uncertainty, not the OLS line. Its t calibration is approximate; agreement does not remove sampling, functional-form or imputation limitations.",
            "remaining_work": "Fit UDS-only and integrated models in Step 4 to answer the adjusted primary question. Residual/influence diagnostics and profile prediction intervals remain later work.",
        },
        "limitations": {"sampling_and_preparation": plan["sampling_model"],
                        "scope": "Exploratory unadjusted inference. Model-based t calibration depends on the stated mean, independence, variance and error-distribution assumptions; count-scale coverage is approximate. HC3 does not establish a correct mean form or participant independence.",
                        "scale": "The response is a count. No raw scores, fitted values or intervals are rounded or clipped to force the usable range.",
                        "no_equivalence_margin": True, "no_automatic_exclusion": True},
        "figure": {"path": str(FIGURE_DIR / FIGURE_NAME), "purpose": "Observed score-thickness scatter and the MRI-only least-squares line over observed thickness support.",
                   "privacy": "No participant IDs or hover labels. All 842 observations are plotted without jitter or exclusions; overlapping observations may coincide.",
                   "intervals": "No prediction or mean-response ribbon is shown at this step."},
        "sources": {"qr": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html",
                    "triangular_solve": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_triangular.html",
                    "pearson": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html",
                    "hc3": "https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.HC3_se.html"},
    }
    return report, design


def build_report() -> dict:
    return analyze()[0]


def regression_figure(report: dict, design: ModelDesign):
    """Draw all supplied observations and the fitted line as a grayscale vector figure."""
    from .phase1_figures import GRAYSCALE, make_figure
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter
    reference = report["design"]["thickness_reference_mm"]
    thickness = [float(Decimal(str(row[1])) * Decimal("0.1") + Decimal(str(reference))) for row in design.matrix]
    intercept = report["coefficients"]["intercept"]["estimate"]
    slope = report["coefficients"]["thickness_0_1mm"]["estimate"]
    figure, axis = make_figure(
        "Animal-naming score and left cortical thickness",
        f"Normal cognition   |   n = {len(thickness):,}   |   Unadjusted MRI-only OLS\n"
        f"Slope / 0.1 mm = {slope:.4f}   |   Pearson r = {report['correlation']['coefficient']:.4f}   |   R-squared = {report['fit']['R_squared']:.6f}",
    )
    figure.set_size_inches(8.6, 5.3)
    axis.scatter(thickness, design.response, s=15, color=GRAYSCALE.accent, alpha=0.5,
                 edgecolors="none", label="Participants", rasterized=False)
    endpoints = [min(thickness), max(thickness)]
    axis.plot(endpoints, [intercept + slope * (value - reference) / 0.1 for value in endpoints],
              color=GRAYSCALE.ink, linewidth=1.7, label="Least-squares line")
    axis.set(xlabel="Left inferior parietal thickness (mm)", ylabel="Correct animal names in 60 seconds",
             ylim=(0, max(design.response) + 4), xlim=(min(thickness) - 0.035, max(thickness) + 0.035))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    axis.legend(loc="upper left", frameon=False, fontsize=9)
    return figure


def load_simple_results() -> dict:
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_report(), "The saved simple-regression result differs from the current plan, design or code. Review and rebuild Step 3.")
    return saved


def main() -> None:
    from .phase1_figures import GRAYSCALE, plt, save_pdf
    report, design = analyze()
    write_output(REPORT_FILE, (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8"))
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelcolor": GRAYSCALE.ink, "axes.labelpad": 9,
                         "pdf.fonttype": 42, "figure.facecolor": "white"}):
        figure = regression_figure(report, design)
        try:
            save_pdf(figure, FIGURE_NAME, FIGURE_DIR)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    main()

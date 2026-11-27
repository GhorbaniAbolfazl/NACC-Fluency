"""Review residuals and influence while retaining the registered full-sample fits."""

import json
from dataclasses import replace
from importlib.metadata import version
from math import fsum, sqrt
from pathlib import Path
from platform import python_version

import numpy as np
from scipy.stats import kurtosis, norm, skew

from .build_dataset import write_output
from .phase4_comparison import REPORT_FILE as COMPARISON_FILE, load_comparison_results
from .phase4_design import ModelDesign, load_prepared_designs
from .phase4_multiple import REPORT_FILE as MULTIPLE_FILE
from .phase4_plan import fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary
from .phase4_simple import REPORT_FILE as SIMPLE_FILE

REPORT_FILE = Path("reports/phase4/residual_diagnostics.json")
FIGURE_DIR = Path("figures/phase4")
FIGURE_NAMES = ("03_residual_diagnostics.pdf", "04_integrated_predictor_residuals.pdf", "05_integrated_influence.pdf")
FLAG_NAMES = ("studentized", "leverage", "cooks_distance")


def summary(values) -> dict:
    values = tuple(float(v) for v in values)
    require(len(values) > 0 and all(np.isfinite(values)), "A diagnostic summary requires finite nonempty values.")
    q = np.quantile(values, [0, .25, .5, .75, 1], method="linear")
    return {"n": len(values), "minimum": float(q[0]), "q1": float(q[1]), "median": float(q[2]),
            "q3": float(q[3]), "maximum": float(q[4]), "mean": fsum(values) / len(values)}


def binned_summary(x, residual, studentized, requested_bins: int = 8) -> list[dict]:
    """Quantile bins keep tied x values together; lines are descriptive, with no CI."""
    x, residual, studentized = (np.asarray(v, dtype=float) for v in (x, residual, studentized))
    require(len(x) == len(residual) == len(studentized) > 1 and requested_bins >= 2,
            "Bin inputs must align and request at least two bins.")
    require(all(np.all(np.isfinite(v)) for v in (x, residual, studentized)), "Bin inputs must be finite.")
    edges = np.unique(np.quantile(x, np.linspace(0, 1, requested_bins + 1), method="linear"))
    groups = np.searchsorted(edges[1:-1], x, side="right")
    result = []
    for group in range(max(1, len(edges) - 1)):
        selected = groups == group
        if not np.any(selected):
            continue
        e, z, values = residual[selected], studentized[selected], x[selected]
        result.append({"n": int(selected.sum()), "x_minimum": float(values.min()), "x_maximum": float(values.max()),
                       "x_mean": float(values.mean()), "residual_mean": float(e.mean()),
                       "residual_SD": float(e.std(ddof=1)) if len(e) > 1 else None,
                       "mean_sqrt_abs_studentized": float(np.sqrt(np.abs(z)).mean())})
    require(sum(v["n"] for v in result) == len(x), "Every observation must enter exactly one descriptive bin.")
    return result


def diagnostic_arrays(fit, thresholds: dict) -> dict:
    e, h = np.asarray(fit.residuals), np.asarray(fit.leverages)
    p, n = len(fit.columns), fit.participants
    studentized = e / (fit.residual_sd * np.sqrt(1 - h))
    cooks = studentized ** 2 * h / (p * (1 - h))
    require(np.all(np.isfinite(studentized)) and np.all(np.isfinite(cooks)) and np.all(cooks >= 0),
            "Residual influence diagnostics must be finite and nonnegative where required.")
    limits = {"studentized": float(thresholds["internally_studentized_residual_abs_gt"]),
              "leverage": 2 * p / n, "cooks_distance": 4 / n}
    masks = {"studentized": np.abs(studentized) > limits["studentized"],
             "leverage": h > limits["leverage"], "cooks_distance": cooks > limits["cooks_distance"]}
    union = np.logical_or.reduce(tuple(masks.values()))
    quantiles = norm.ppf((np.arange(n, dtype=float) + .5) / n)
    observed = np.sort(studentized)
    oq = np.quantile(observed, [.25, .75], method="linear")
    nq = norm.ppf([.25, .75])
    slope = float((oq[1] - oq[0]) / (nq[1] - nq[0]))
    intercept = float(oq[0] - slope * nq[0])
    return {"studentized": studentized, "cooks_distance": cooks, "limits": limits, "masks": masks,
            "union": union, "qq_theoretical": quantiles, "qq_observed": observed,
            "qq_line_slope": slope, "qq_line_intercept": intercept}


def flag_summary(arrays: dict) -> dict:
    masks = arrays["masks"]
    combinations = {}
    for bits in range(8):
        selected = np.ones(len(arrays["union"]), dtype=bool)
        names = []
        for j, name in enumerate(FLAG_NAMES):
            active = bool(bits & (1 << j))
            selected &= masks[name] if active else ~masks[name]
            if active:
                names.append(name)
        combinations[" + ".join(names) if names else "none"] = int(selected.sum())
    require(sum(combinations.values()) == len(arrays["union"]), "Flag combinations must partition the sample.")
    return {"thresholds_strict_gt": arrays["limits"], "counts": {name: int(mask.sum()) for name, mask in masks.items()},
            "union_count": int(arrays["union"].sum()), "mutually_exclusive_combinations": combinations,
            "meaning": "Review heuristics on the full fitted model, not tests, proof of data error or automatic exclusion rules."}


def uncertainty(fit) -> dict:
    q = coefficient_table(fit)["thickness_0_1mm"]
    return {"estimate": q["estimate"], **{method: {"SE": q[method]["standard_error"],
             "interval": q[method]["confidence_interval"]} for method in ("model_based", "hc3")}}


def subset(design: ModelDesign, keep: tuple[int, ...]) -> ModelDesign:
    require(len(keep) == len(set(keep)) and all(0 <= i < len(design.response) for i in keep), "Invalid diagnostic subset indices.")
    return replace(design, name="integrated_diagnostic_omission", participant_ids=tuple(design.participant_ids[i] for i in keep),
                   matrix=tuple(design.matrix[i] for i in keep), response=tuple(design.response[i] for i in keep))


def influence_followup(design: ModelDesign, fit, arrays: dict) -> tuple[dict, tuple[dict, ...]]:
    """Refit each flagged omission and their union as diagnostics, never new primary fits."""
    flagged = tuple(int(i) for i in np.flatnonzero(arrays["union"]))
    results, max_identity_error = [], 0.0
    g = np.asarray(fit.unscaled_covariance)
    beta = np.asarray(fit.coefficients)
    for index in flagged:
        keep = tuple(i for i in range(fit.participants) if i != index)
        omitted = fit_ols(subset(design, keep))
        expected = beta - (g @ np.asarray(design.matrix[index])) * fit.residuals[index] / (1 - fit.leverages[index])
        error = float(np.max(np.abs(np.asarray(omitted.coefficients) - expected)))
        require(np.allclose(omitted.coefficients, expected, rtol=1e-10, atol=1e-11), "Direct leave-one-out refit disagrees with the analytic deletion identity.")
        max_identity_error = max(max_identity_error, error)
        results.append(uncertainty(omitted))
    baseline = uncertainty(fit)

    def summarize_omissions(items: tuple[dict, ...]) -> dict:
        if not items:
            return {"omissions": 0}
        estimates = [q["estimate"] for q in items]
        return {"omissions": len(items), "estimate": summary(estimates),
                "maximum_absolute_estimate_change": max(abs(v - baseline["estimate"]) for v in estimates),
                "negative_estimates": sum(v < 0 for v in estimates), "positive_estimates": sum(v > 0 for v in estimates),
                "uncertainty": {method: {"SE": summary([q[method]["SE"] for q in items]),
                    "lower_endpoint": summary([q[method]["interval"]["lower"] for q in items]),
                    "upper_endpoint": summary([q[method]["interval"]["upper"] for q in items]),
                    "intervals_including_zero": sum(q[method]["interval"]["lower"] <= 0 <= q[method]["interval"]["upper"] for q in items)}
                    for method in ("model_based", "hc3")}}
    joint = None
    if flagged:
        excluded = set(flagged)
        keep = tuple(i for i in range(fit.participants) if i not in excluded)
        joint_fit = fit_ols(subset(design, keep))
        joint = {"omitted": len(flagged), "retained": len(keep), "parameters": len(fit.columns),
                 "residual_df": joint_fit.residual_df, "thickness": uncertainty(joint_fit),
                 "estimate_change_from_primary": joint_fit.coefficients[1] - baseline["estimate"]}
    by_flag = {name: summarize_omissions(tuple(q for index, q in zip(flagged, results) if arrays["masks"][name][index])) for name in FLAG_NAMES}
    return {"primary_full_sample": {"participants": fit.participants, "residual_df": fit.residual_df, "thickness": baseline},
            "flag_set": "Union of the three full-model heuristic flags, fixed before diagnostic omissions; not recalculated after deleting observations.",
            "individual_leave_one_out": {"participants_per_fit": fit.participants - 1, "residual_df_per_fit": fit.residual_df - 1,
                                          **summarize_omissions(tuple(results)), "by_original_flag": by_flag,
                                          "maximum_coefficient_identity_error": max_identity_error},
            "joint_union_omission": joint,
            "interpretation": "Outcome-dependent diagnostic omissions do not define a cleaned or independently validated sample. Recomputed model-based and HC3 intervals describe sensitivity using usual formulas; their nominal coverage is not justified after flag-based selection. Endpoint ranges are not simultaneous confidence bands. No new confirmatory p-values are reported. All 842 participants remain in the primary analysis.",
            "privacy": "Publish aggregate omission summaries only; no IDs, row indices, per-person omission results or ordered individual diagnostics."}, tuple(results)


def analyze() -> tuple[dict, dict]:
    comparison = load_comparison_results()
    plan, designs = load_prepared_designs()
    simple = json.loads(SIMPLE_FILE.read_text())
    multiple = json.loads(MULTIPLE_FILE.read_text())
    tracked = {**comparison["source_files_sha256"], str(COMPARISON_FILE): fingerprint(COMPARISON_FILE),
               "src/phase4_diagnostics.py": fingerprint(Path("src/phase4_diagnostics.py")),
               "src/phase4_diagnostic_figures.py": fingerprint(Path("src/phase4_diagnostic_figures.py"))}
    thresholds = plan["diagnostics"]["review_thresholds"]
    require(thresholds["leverage_gt"] == "2*p/n" and thresholds["cooks_distance_gt"] == "4/n",
            "Registered diagnostic thresholds changed.")
    fits, arrays, models = {}, {}, {}
    for name, design in designs.items():
        fit = fit_ols(design)
        previous = simple if name == "mri_only" else multiple["models"][name]
        require(fit_summary(fit) == previous["fit"] and coefficient_table(fit) == previous["coefficients"],
                "A full-sample diagnostic fit differs from the saved primary model.")
        a = diagnostic_arrays(fit, thresholds)
        fits[name], arrays[name] = fit, a
        model = {"participants": fit.participants, "parameters": len(fit.columns), "residual_df": fit.residual_df,
                 "fitted_values": summary(fit.fitted_values), "residuals": summary(fit.residuals),
                 "internally_studentized": summary(a["studentized"]), "leverage": summary(fit.leverages),
                 "cooks_distance": summary(a["cooks_distance"]), "review_flags": flag_summary(a),
                 "fitted_outside_usable_score_range": sum(v < 0 or v > 77 for v in fit.fitted_values),
                 "residual_skewness_bias_corrected": float(skew(fit.residuals, bias=False)),
                 "residual_excess_kurtosis_bias_corrected": float(kurtosis(fit.residuals, fisher=True, bias=False)),
                 "fitted_value_bins": binned_summary(fit.fitted_values, fit.residuals, a["studentized"]),
                 "normal_QQ_reference": {"positions": "Normal quantiles at (i-0.5)/n, i=1,...,n; sorted internally studentized residuals.",
                                         "line": "Reference through the observed 25th/75th percentiles and the corresponding normal quantiles; not a confidence band or a normality test.",
                                         "slope": a["qq_line_slope"], "intercept": a["qq_line_intercept"]}}
        models[name] = model
    integrated, fit, a = designs["integrated"], fits["integrated"], arrays["integrated"]
    predictor_bins = {}
    references = {"thickness_0_1mm": (2.344, .1), "age_c": (68, 1), "education_c": (16, 1)}
    for term, (reference, scale) in references.items():
        j = integrated.columns.index(term)
        raw_x = [row[j] * scale + reference for row in integrated.matrix]
        predictor_bins[term] = binned_summary(raw_x, fit.residuals, a["studentized"])
    sex_j = integrated.columns.index("female")
    sex_groups = {}
    for code, label in ((0, "recorded_male"), (1, "recorded_female")):
        indices = [i for i, row in enumerate(integrated.matrix) if row[sex_j] == code]
        e = [fit.residuals[i] for i in indices]
        sex_groups[label] = {"residuals": summary(e), "residual_SD": float(np.std(e, ddof=1))}
    followup, omissions = influence_followup(integrated, fit, a)
    report = {
        "stage": "Phase 4, Step 6", "status": "Residual/influence calculations, figures and graphical review completed; final model retention/extension and profile intervals remain Step 7.",
        "sample": "model_normal", "participants": fit.participants, "source_files_sha256": tracked,
        "environment": {"python": python_version(), "scipy": version("scipy"), "numpy_runtime": version("numpy"), "matplotlib": version("matplotlib")},
        "models": models, "integrated_predictor_bins": predictor_bins, "integrated_sex_residual_summaries": sex_groups,
        "integrated_influence_followup": followup, "primary_test_unchanged": comparison["primary_test"],
        "graphical_review": {"review_date": "2026-11-26", "scope": "Visual review of the three final grayscale figures for the frozen 842-person sample; descriptive evidence, not an automatic validity test.",
                             "mean_form": "Residual-versus-fitted bin means remain near zero relative to the residual spread in UDS-only and integrated models. Predictor plots do not show an obvious sustained broad curvature pattern in the observed central support; sparse thickness/age/education tails are less informative. This does not prove an adequate linear mean.",
                             "normal_reference": "Central Q-Q points are close to the quartile reference, while upper positive tails rise above it and lower tails flatten. Positive residual skewness and excess kurtosis support qualifying exact normal-error calibration rather than declaring normality.",
                             "variance": "Scale-location bin means and within-bin residual SDs vary, with somewhat greater spread at the upper fitted end; no dramatic funnel dominates the integrated plot. Constant conditional variance is not established. Retain both the classical and already planned HC3 summaries.",
                             "influence": "Residual, leverage and Cook flags overlap without identifying data errors. Individual flagged omissions and their joint omission change the thickness estimate; both uncertainty methods must be read alongside those changes. They do not authorize replacement of the primary sample.",
                             "action_for_step_7": "Use this evidence to record a justified retention/extension decision before profile intervals or any optional extension. No response transformation, quadratic term or alternative primary model is fitted here."},
        "methods": {"residual": "e_i = y_i - fitted_i", "internally_studentized": "r_i = e_i/(s*sqrt(1-h_ii)); s=sqrt(SSE/(n-p)) from the full model, not an externally studentized residual.",
                    "leverage": "h_ii is the sum of squared entries in row i of the full-rank economic Q matrix; trace(H)=p.",
                    "cooks_distance": "D_i = (r_i^2/p) * h_ii/(1-h_ii); a model-wide coefficient influence measure, not a thickness-specific test.",
                    "LOO_identity": "beta_(-i) = beta - (X'X)^(-1) x_i e_i/(1-h_ii); checked against direct QR refits for each flagged observation.",
                    "binning": "Eight requested empirical-quantile bins of the displayed x, using linear quantiles and unique edges; all tied x values stay together. Lines connect within-bin means, with no smoothing fit or confidence ribbon. Each point remains in the scatter.",
                    "scale_location": "sqrt(abs(internally studentized residual)) against fitted values; bin means are descriptive.",
                    "formal_tests": "No new normality, variance or coefficient-difference hypothesis test is introduced. Quantiles, skewness, excess kurtosis and bin summaries are descriptive diagnostics."},
        "interpretation_limits": {"calibration": "Normal Q-Q departures and spread patterns qualify exact classical t and prediction-interval calibration. Internally studentized residuals are dependent; the normal reference is diagnostic, not an exact joint reference distribution.",
                                  "HC3": "Retain the already planned HC3 coefficient companion regardless of the visual pattern. It does not correct nonlinearity, participant dependence, selection, source-imputation uncertainty or individual prediction variance.",
                                  "sampling": plan["sampling_model"], "no_automatic_deletion": plan["diagnostics"]["influence_followup"],
                                  "causality": "Neither graphical consistency nor stable diagnostic omissions proves a causal association, equivalence or independence.",
                                  "next_step": "Record the graphical evidence and a justified retention/extension decision in Step 7 before any exploratory extension or profile prediction intervals. No quadratic term, transformation or model search is fitted in this step."},
        "validation": {"passed": True, "checks": ["Steps 1-5 reproduce from the frozen inputs.", "All three full-sample refits exactly match their saved coefficients, uncertainty and fit metrics.",
                          "Finite residual, internally studentized, leverage and Cook diagnostics use the registered strict thresholds.", "Mutually exclusive flag combinations partition each common sample.",
                          "Each bin partitions all observations without splitting ties or excluding points.", "Direct flagged leave-one-out QR coefficients agree with the analytic deletion identity.",
                          "Primary sample, predictors, coefficient test and previous outputs remain unchanged.", "Inputs and source fingerprints remain unchanged during diagnostics."],
                       "limit": "Computational checks do not validate sampling, independence, mean form or constant variance."},
        "figures": [{"path": str(FIGURE_DIR / name), "page": 1} for name in FIGURE_NAMES],
        "privacy": "JSON contains aggregate summaries and counts only. Diagnostic plots show unlabeled coordinates or histograms without IDs, source row indices or hover labels; no participant-level diagnostic arrays or omission tables are exported.",
        "progress": {"completed_through_step": 6, "remaining": ["Step 7 documented model retention/extension and profile intervals", "Step 8 final numerical review, LaTeX report and executed notebook"],
                     "historical_records": "Earlier step status statements remain historical snapshots."},
        "sources": {"graphical_diagnostics": "https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd44.htm",
                    "influence": "https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.OLSInfluence.html",
                    "normal_quantiles": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html"},
    }
    require(tracked == {name: fingerprint(Path(name)) for name in tracked}, "Inputs or sources changed during residual diagnostics.")
    return report, {"designs": designs, "fits": fits, "arrays": arrays, "omissions": omissions}


def build_report() -> dict:
    return analyze()[0]


def load_diagnostic_results() -> dict:
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_report(), "Saved residual diagnostics differ from current inputs or code. Review and rebuild Step 6.")
    return saved


def main() -> None:
    from .phase1_figures import plt, save_pdf
    from .phase4_diagnostic_figures import build_figures
    report, context = analyze()
    output = (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != output:
        write_output(REPORT_FILE, output)
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42, "figure.facecolor": "white"}):
        figures = build_figures(report, context)
        for name, figure in zip(FIGURE_NAMES, figures):
            try:
                save_pdf(figure, name, FIGURE_DIR)
            finally:
                plt.close(figure)


if __name__ == "__main__":
    main()

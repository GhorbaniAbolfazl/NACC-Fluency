"""Compare the three registered models and review predictor collinearity."""

import json
from importlib.metadata import version
from math import fsum, isclose, isfinite, sqrt
from pathlib import Path
from platform import python_version

from scipy.linalg import lstsq

from .build_dataset import write_output
from .phase4_design import ModelDesign, canonical_hash, load_prepared_designs, matrix_diagnostics
from .phase4_multiple import REPORT_FILE as MULTIPLE_FILE, load_multiple_results
from .phase4_plan import fingerprint, require
from .phase4_simple import REPORT_FILE as SIMPLE_FILE, load_simple_results

REPORT_FILE = Path("reports/phase4/model_comparison.json")
FIGURE_DIR = Path("figures/phase4")
FIGURE_NAME = "02_coefficient_comparison.pdf"
MODEL_ORDER = ("uds_only", "mri_only", "integrated")
COMPARISONS = (("thickness_0_1mm", "mri_only"), ("age_c", "uds_only"),
               ("female", "uds_only"), ("education_c", "uds_only"))
METRICS = ("participants", "parameters_including_intercept", "residual_df", "SSE", "R_squared",
           "adjusted_R_squared", "residual_SD", "in_sample_RMSE")


def collinearity(design: ModelDesign, unscaled_covariance: list[list[float]], threshold: float) -> dict:
    """Regress each predictor on the others, including an intercept; omit intercept VIF."""
    diagnostics = matrix_diagnostics(design)
    n, p = len(design.response), len(design.columns)
    require(len(unscaled_covariance) == p and all(len(row) == p for row in unscaled_covariance),
            "Covariance dimensions must match the model columns.")
    require(isfinite(threshold) and threshold > 1, "The VIF review threshold must be finite and greater than one.")
    columns = tuple(zip(*design.matrix))
    centers = [fsum(values) / n for values in columns]
    centered = [tuple(value - mean for value in values) for values, mean in zip(columns, centers)]
    ss = [fsum(value * value for value in values) for values in centered]
    correlation = [[fsum(a * b for a, b in zip(centered[j], centered[k])) / sqrt(ss[j] * ss[k])
                    for k in range(1, p)] for j in range(1, p)]
    vif = {}
    for j in range(1, p):
        auxiliary = tuple(tuple(value for k, value in enumerate(row) if k != j) for row in design.matrix)
        beta, _, rank, _ = lstsq(auxiliary, columns[j], lapack_driver="gelsd", check_finite=True)
        require(rank == p - 1, "A VIF auxiliary regression is rank deficient.")
        sse = fsum((y - fsum(x * float(b) for x, b in zip(row, beta))) ** 2
                   for row, y in zip(auxiliary, columns[j]))
        require(isfinite(sse) and sse > 0, "VIF requires positive residual variation in each auxiliary regression.")
        value = ss[j] / sse
        require(isfinite(value) and value >= 1 - 1e-10, "A VIF must be finite and at least one within numeric tolerance.")
        value = max(1.0, value)  # Remove only a possible rounding-sized deficit below one.
        gram_identity = ss[j] * unscaled_covariance[j][j]
        require(isclose(value, gram_identity, rel_tol=1e-10, abs_tol=1e-11),
                "Auxiliary-SVD VIF disagrees with the earlier QR inverse-Gram identity.")
        vif[design.columns[j]] = {"VIF": value, "tolerance": 1 / value, "auxiliary_R_squared": 1 - 1 / value,
                                 "auxiliary_SSE": sse, "predictor_SST": ss[j],
                                 "model_based_SE_inflation_factor": sqrt(value),
                                 "review_flag": value > threshold, "inverse_gram_check": gram_identity}
    flagged = [term for term, value in vif.items() if value["review_flag"]]
    return {"predictor_order": list(design.columns[1:]), "correlation_matrix": correlation,
            "correlation_scope": "Pearson correlations among explanatory columns on the common sample; includes the recorded-sex indicator. No response or intercept is included.",
            "vif": vif, "review_threshold_gt": threshold, "flagged_terms": flagged,
            "full_column_rank": diagnostics["full_column_rank"], "rank": diagnostics["raw_design"]["svd_rank"],
            "rank_tolerance": diagnostics["raw_design"]["rank_tolerance"],
            "raw_design_condition_number": diagnostics["raw_design"]["condition_number_2norm"],
            "standardized_design_condition_number": diagnostics["standardized_design"]["condition_number_2norm"],
            "standardization": diagnostics["standardization"],
            "interpretation": (f"{len(flagged)} explanatory terms exceed the registered VIF review threshold. "
                               "Review rank, VIF and conditioning together. VIF describes classical coefficient-variance inflation from predictor dependence, holding residual variance and predictor spread fixed; its square root is the corresponding SE factor. It is not the observed HC3-to-model-based SE ratio. A low VIF does not establish model validity, causality or a precisely estimated thickness association."),
            "action": "Retain all registered covariates. VIF is a review heuristic, not an automatic selection or deletion rule."}


def coefficient_change(before: dict, after: dict, source_model: str) -> dict:
    methods = {}
    for method in ("model_based", "hc3"):
        a, b = before[method], after[method]
        aw = a["confidence_interval"]["upper"] - a["confidence_interval"]["lower"]
        bw = b["confidence_interval"]["upper"] - b["confidence_interval"]["lower"]
        methods[method] = {"before_SE": a["standard_error"], "after_SE": b["standard_error"],
                           "SE_ratio_after_over_before": b["standard_error"] / a["standard_error"],
                           "before_interval": a["confidence_interval"], "after_interval": b["confidence_interval"],
                           "before_interval_width": aw, "after_interval_width": bw,
                           "interval_width_ratio_after_over_before": bw / aw}
    return {"from_model": source_model, "to_model": "integrated", "before_estimate": before["estimate"],
            "after_estimate": after["estimate"], "change_integrated_minus_baseline": after["estimate"] - before["estimate"],
            "uncertainty": methods,
            "scope": "Descriptive change between models on the same participants. The estimates are dependent and target different conditional associations; this is not an independent-samples difference test. No percent change is reported for a near-zero baseline coefficient."}


def thickness_decomposition(design: ModelDesign, models: dict) -> dict:
    """Check the exact marginal/conditional slope identity, without a causal claim."""
    j = design.columns.index("thickness_0_1mm")
    x = tuple(row[j] for row in design.matrix)
    mean_x = fsum(x) / len(x)
    sxx = fsum((value - mean_x) ** 2 for value in x)
    pieces = {}
    for term in ("age_c", "female", "education_c"):
        k = design.columns.index(term)
        z = tuple(row[k] for row in design.matrix)
        mean_z = fsum(z) / len(z)
        slope_z_on_x = fsum((a - mean_x) * (b - mean_z) for a, b in zip(x, z)) / sxx
        adjusted_coefficient = models["integrated"]["coefficients"][term]["estimate"]
        pieces[term] = {"predictor_on_thickness_slope": slope_z_on_x,
                        "integrated_score_coefficient": adjusted_coefficient,
                        "contribution_to_marginal_minus_adjusted": slope_z_on_x * adjusted_coefficient}
    adjusted = models["integrated"]["coefficients"]["thickness_0_1mm"]["estimate"]
    marginal = models["mri_only"]["coefficients"]["thickness_0_1mm"]["estimate"]
    recovered = adjusted + fsum(value["contribution_to_marginal_minus_adjusted"] for value in pieces.values())
    require(isclose(marginal, recovered, rel_tol=1e-10, abs_tol=1e-11), "Marginal/conditional thickness slope identity failed.")
    return {"marginal_slope": marginal, "adjusted_slope": adjusted, "components": pieces,
            "reconstructed_marginal_slope": recovered, "absolute_identity_error": abs(marginal - recovered),
            "identity": "MRI-only thickness slope = integrated thickness slope + sum(integrated covariate coefficient * slope of that covariate regressed on thickness alone with an intercept).",
            "interpretation": "An algebraic description of the observed coefficient change, not a causal decomposition or proof of confounding. Covariate correlations and score coefficients jointly affect the marginal slope; sign change alone does not establish a protective/harmful effect or severe multicollinearity. Both thickness intervals still need their stated uncertainty and assumptions."}


def build_report() -> dict:
    simple, multiple = load_simple_results(), load_multiple_results()
    plan, designs = load_prepared_designs()
    tracked = {**multiple["source_files_sha256"], str(SIMPLE_FILE): fingerprint(SIMPLE_FILE),
               str(MULTIPLE_FILE): fingerprint(MULTIPLE_FILE),
               "src/phase4_comparison.py": fingerprint(Path("src/phase4_comparison.py")),
               "src/phase1_figures.py": fingerprint(Path("src/phase1_figures.py"))}
    models = {name: ({"role": plan["models"][name]["role"], "fit": simple["fit"],
                      "coefficients": simple["coefficients"], "covariance": simple["covariance"],
                      "matrix_sha256": simple["design"]["matrix_sha256"],
                      "response_sha256": simple["design"]["response_sha256"],
                      "membership_sha256": simple["design"]["membership_sha256"]}
                     if name == "mri_only" else multiple["models"][name]) for name in MODEL_ORDER}
    first = designs["integrated"]
    for name, d in designs.items():
        saved = models[name]
        require(d.participant_ids == first.participant_ids and d.response == first.response
                and saved["membership_sha256"] == plan["sample"]["membership_sha256"]
                and saved["response_sha256"] == canonical_hash(d.response)
                and saved["matrix_sha256"] == canonical_hash(d.matrix), "A comparison sample or matrix differs from its validated fitted model.")
        require(saved["fit"]["participants"] == len(d.response) == 842
                and saved["fit"]["parameters_including_intercept"] == len(d.columns)
                and saved["fit"]["residual_df"] == 842 - len(d.columns), "Fit dimensions do not match the common design.")
    require(all(models[name]["fit"]["SST"] == models["integrated"]["fit"]["SST"] for name in MODEL_ORDER),
            "Common-sample models must use identical total response sum of squares.")
    fit_table = {name: {metric: models[name]["fit"][metric] for metric in METRICS} for name in MODEL_ORDER}
    increments = {}
    for baseline in ("uds_only", "mri_only"):
        before, after = fit_table[baseline], fit_table["integrated"]
        indices = [first.columns.index(term) for term in designs[baseline].columns]
        require(designs[baseline].matrix == tuple(tuple(row[j] for j in indices) for row in first.matrix),
                "The claimed baseline model is not nested in the integrated design.")
        require(after["SSE"] <= before["SSE"] + 1e-9 and after["R_squared"] >= before["R_squared"] - 1e-12,
                "Nested-model OLS must not worsen SSE or R-squared on the same observations.")
        increments[baseline + "_to_integrated"] = {"added_terms": [term for term in first.columns if term not in designs[baseline].columns],
            "SSE_reduction": before["SSE"] - after["SSE"],
            "R_squared_increase": after["R_squared"] - before["R_squared"],
            "adjusted_R_squared_change": after["adjusted_R_squared"] - before["adjusted_R_squared"],
            "residual_SD_change": after["residual_SD"] - before["residual_SD"],
            "in_sample_RMSE_change": after["in_sample_RMSE"] - before["in_sample_RMSE"],
            "interpretation": "Descriptive nested-model change only. SSE and unadjusted R-squared cannot worsen when columns are added; adjusted R-squared and residual SD can worsen because residual degrees of freedom decrease. No extra hypothesis test or model-selection rule is introduced."}
    changes = {term: coefficient_change(models[baseline]["coefficients"][term], models["integrated"]["coefficients"][term], baseline)
               for term, baseline in COMPARISONS}
    threshold = plan["diagnostics"]["collinearity"]["vif_review_gt"]
    dependence = {name: collinearity(designs[name], models[name]["covariance"]["unscaled"], threshold)
                  for name in plan["diagnostics"]["collinearity"]["models"]}
    report = {
        "stage": "Phase 4, Step 5", "status": "Three-model comparison and collinearity review completed; residual/influence diagnostics remain pending.",
        "sample": "model_normal", "participants": len(first.response), "source_files_sha256": tracked,
        "environment": {"python": python_version(), "scipy": version("scipy"), "numpy_runtime": version("numpy"), "matplotlib": version("matplotlib")},
        "common_sample": {"membership_sha256": plan["sample"]["membership_sha256"], "response_sha256": canonical_hash(first.response),
                          "identical_participants_and_response_order": True, "identical_shared_predictor_values": True,
                          "no_exclusions_or_new_predictors": True},
        "model_order": list(MODEL_ORDER), "fit_comparison": fit_table,
        "metric_definitions": {"SSE": "Sum of squared residuals in squared score units", "R_squared": "1-SSE/SST",
                               "adjusted_R_squared": "1-(SSE/(n-p))/(SST/(n-1))", "residual_SD": "sqrt(SSE/(n-p)), in correct responses",
                               "in_sample_RMSE": "sqrt(SSE/n), in correct responses", "p": "Number of parameters including the intercept"},
        "nested_model_changes": increments, "coefficient_comparisons": changes,
        "thickness_change_identity": thickness_decomposition(first, models),
        "collinearity": dependence,
        "collinearity_methods": {"VIF": "SST_j/SSE_j = 1/(1-R_j^2), regressing explanatory column j on an intercept and every other explanatory column on the same sample.",
                                 "numerics": "Auxiliary OLS uses SVD-based scipy.linalg.lstsq(lapack_driver='gelsd'); VIF is cross-checked with predictor SST times the matching diagonal of the earlier pivoted-QR unscaled covariance.",
                                 "intercept": "Excluded from correlations and VIF; retained in all main and auxiliary models.",
                                 "single_MRI_predictor": "MRI-only has one explanatory term; dependence among MRI predictors is inapplicable, so no multi-predictor VIF table is reported for it.",
                                 "conditioning": "Full rank and raw/standardized-design condition numbers are checked again against the validated designs. Center/scale transformations are diagnostics only; original coefficient units remain unchanged.",
                                 "review_threshold": "VIF > 5 is the frozen review heuristic; no automatic covariate removal."},
        "testing_policy": multiple["testing_policy"], "primary_test": multiple["primary_test"],
        "interpretation": {"comparison": "All three models use the same 842 people; the fit/estimate changes therefore reflect added explanatory terms rather than sample replacement. These are in-sample association summaries, not out-of-sample predictive accuracy.",
                           "uncertainty": "Compare each coefficient with both of its individual 95% intervals. Intervals across models are dependent; overlap or a change in sign is not a test of the difference between coefficients. HC3 remains a companion and is not selected by its p-value.",
                           "thickness": "Read the unadjusted and adjusted thickness coefficients as different conditional summaries. Low VIF does not establish absence of an association, equivalence or causality. The primary test and uncertainty remain those in Step 4, pending residual and influence diagnostics.",
                           "selection": "No best-model search, automatic term removal, new contrast test, information criterion or cross-validation is performed. The registered integrated model remains the primary specification."},
        "validation": {"passed": True, "checks": ["Frozen plan, common-sample designs and both prior regression outputs reproduce.",
                         "Participant membership, response order, common predictor columns and SST agree across all three models.",
                         "Nested design relationships and expected SSE/R-squared monotonicity hold.",
                         "Each coefficient comparison preserves units and both individual uncertainty methods.",
                         "The marginal thickness slope reconstructs from the integrated coefficients and same-sample covariate-on-thickness slopes.",
                         "VIF from auxiliary SVD regressions agrees with the earlier QR inverse-Gram identity.",
                         "Both multiple designs retain full rank with recorded numerical conditioning.",
                         "Inputs and source fingerprints remain unchanged during comparison."],
                       "limit": "These checks establish numerical consistency; residual assumptions, participant independence and sampling design remain separate questions."},
        "privacy": "Only aggregate metrics, coefficient intervals, small predictor-correlation matrices, VIFs and source hashes are saved. No participant IDs, design rows, responses or individual fitted/residual/influence arrays are exported.",
        "figure": {"path": str(FIGURE_DIR / FIGURE_NAME), "purpose": "Compare the four common non-intercept coefficients in their own units, with model-based and HC3 individual 95% intervals.",
                   "scale": "Separate x-axis ranges by term; do not compare magnitudes across panels as standardized effects."},
        "limitations": multiple["limitations"],
        "progress": {"completed_through_step": 5, "remaining": ["Step 6 residual and influence diagnostics", "Step 7 justified model retention/extension and profile intervals", "Step 8 final numerical review, LaTeX report and executed notebook"],
                     "historical_records": "Earlier outputs and the frozen plan preserve their original stage status."},
        "sources": {"vif": "https://www.statsmodels.org/stable/generated/statsmodels.stats.outliers_influence.variance_inflation_factor.html",
                    "auxiliary_ols": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lstsq.html"},
    }
    require(tracked == {name: fingerprint(Path(name)) for name in tracked}, "Inputs or sources changed during model comparison.")
    return report


def comparison_figure(report: dict):
    """Draw four coefficient panels in grayscale with distinct interval styles."""
    from .phase1_figures import GRAYSCALE, plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    figure, axes = plt.subplots(2, 2, figsize=(10.4, 6.6), layout="constrained")
    figure.suptitle("Coefficient estimates across the registered models\n"
                    "Normal cognition | n = 842 | Individual 95% intervals", fontsize=15, fontweight="bold", color=GRAYSCALE.ink)
    labels = {"mri_only": "MRI-only", "uds_only": "UDS-only", "integrated": "UDS + MRI"}
    units = {"thickness_0_1mm": ("Left inferior parietal thickness", "Correct responses per 0.1 mm"),
             "age_c": ("Age", "Correct responses per year"), "female": ("Recorded sex", "Correct responses: female minus male"),
             "education_c": ("Education", "Correct responses per year")}
    for axis, (term, _) in zip(axes.flat, COMPARISONS):
        entry = report["coefficient_comparisons"][term]
        for y, prefix in ((1, "before"), (0, "after")):
            estimate = entry[prefix + "_estimate"]
            for method, offset, style, marker, color in (("model_based", 0.09, "-", "o", GRAYSCALE.ink),
                                                         ("hc3", -0.09, "--", "s", GRAYSCALE.muted)):
                interval = entry["uncertainty"][method][prefix + "_interval"]
                axis.plot([interval["lower"], interval["upper"]], [y + offset] * 2, linestyle=style, color=color, linewidth=1.6)
                axis.plot([estimate], [y + offset], marker=marker, markersize=4.8, color=color,
                          markerfacecolor=color if method == "model_based" else "white", linestyle="none")
        axis.axvline(0, color=GRAYSCALE.outline, linestyle=":", linewidth=1)
        axis.set_yticks([1, 0], [labels[entry["from_model"]], labels[entry["to_model"]]])
        axis.set_ylim(-0.4, 1.4)
        axis.set_title(units[term][0], loc="left", fontsize=11, fontweight="bold", color=GRAYSCALE.ink)
        axis.set_xlabel(units[term][1], color=GRAYSCALE.ink)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.spines["bottom"].set_color(GRAYSCALE.outline)
        axis.tick_params(axis="y", length=0, colors=GRAYSCALE.muted)
        axis.tick_params(axis="x", colors=GRAYSCALE.muted)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.grid(axis="x", color=GRAYSCALE.grid, linewidth=0.7)
        axis.set_axisbelow(True)
        axis.margins(x=0.14)
    legend = [Line2D([0], [0], color=GRAYSCALE.ink, marker="o", linewidth=1.6, label="Model-based"),
              Line2D([0], [0], color=GRAYSCALE.muted, marker="s", markerfacecolor="white", linestyle="--", linewidth=1.6, label="HC3 companion")]
    figure.legend(handles=legend, loc="outside lower center", ncols=2, frameon=False, fontsize=10)
    return figure


def load_comparison_results() -> dict:
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    require(saved == build_report(), "The saved model comparison differs from the current inputs or code. Review and rebuild Step 5.")
    return saved


def main() -> None:
    from .phase1_figures import plt, save_pdf
    report = build_report()
    output = (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != output:
        write_output(REPORT_FILE, output)
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "figure.facecolor": "white"}):
        figure = comparison_figure(report)
        try:
            save_pdf(figure, FIGURE_NAME, FIGURE_DIR)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    main()

"""Validate the five-phase evidence and generate the final tables and figures."""

import json
from math import isclose
from pathlib import Path

from . import phase5_bilateral, phase5_influence, phase5_samples, phase5_subgroups
from .build_dataset import LEFT_THICKNESS, write_output
from .data_audit import CATEGORY_LABELS
from .descriptive_statistics import REPORT_FILE as DESCRIPTION, describe, load_participants
from .phase1_figures import GRAYSCALE, plt, save_pdf
from .phase2_binomial_count import summarize_count
from .phase2_distributions import summarize_distributions
from .phase2_event_probabilities import summarize_events
from .phase2_independence_bayes import compare_probabilities
from .phase2_samples import load_phase2_samples
from .phase4_plan import fingerprint, require
from .phase4_tables import interval, p_tex, rows_tex
from .phase5_plan import load_registered_plan, read_json

DIRECTORY = Path("reports/final")
REPORT_FILE = DIRECTORY / "synthesis.json"
TABLE_FILE = DIRECTORY / "tables.tex"
FIGURE_DIR = Path("figures/phase5")
THICKNESS = "thickness_0_1mm"
SOURCES = {
    "description": DESCRIPTION,
    "events": Path("reports/phase2/event_probabilities.json"),
    "bayes": Path("reports/phase2/independence_bayes.json"),
    "counts": Path("reports/phase2/binomial_count.json"),
    "distributions": Path("reports/phase2/distributions.json"),
    "means": Path("reports/phase3/mean_intervals.json"),
    "means_bootstrap": Path("reports/phase3/sampling_distributions.json"),
    "contrast": Path("reports/phase3/score_difference_by_mri.json"),
    "contrast_bootstrap": Path("reports/phase3/score_difference_bootstrap.json"),
    "group_thickness": Path("reports/phase3/thickness_by_status.json"),
    "tests": Path("reports/phase3/hypothesis_tests.json"),
    "simple": Path("reports/phase4/simple_regression.json"),
    "multiple": Path("reports/phase4/multiple_regression.json"),
    "comparison": Path("reports/phase4/model_comparison.json"),
    "diagnostics": Path("reports/phase4/residual_diagnostics.json"),
    "prediction": Path("reports/phase4/profile_prediction.json"),
}
SENSITIVITY = {
    "samples": phase5_samples, "bilateral": phase5_bilateral,
    "influence": phase5_influence, "subgroups": phase5_subgroups,
}
MODEL_LABELS = {"uds_only": "UDS-only", "mri_only": "MRI-only", "integrated": "UDS + MRI"}


def compare(expected, observed, context: str) -> None:
    """Allow numerical round-off across supported environments, with exact structure."""
    if isinstance(expected, dict):
        require(isinstance(observed, dict) and expected.keys() == observed.keys(), f"Fields differ: {context}")
        for key, value in expected.items():
            compare(value, observed[key], f"{context}.{key}")
    elif isinstance(expected, list):
        require(isinstance(observed, list) and len(expected) == len(observed), f"Length differs: {context}")
        for i, (a, b) in enumerate(zip(expected, observed)):
            compare(a, b, f"{context}[{i}]")
    elif isinstance(expected, float):
        require(type(observed) in (int, float) and isclose(expected, observed, rel_tol=2e-10, abs_tol=2e-11),
                f"Value differs: {context}")
    else:
        require(type(expected) is type(observed) and expected == observed, f"Value differs: {context}")


def evidence() -> tuple[dict, dict, dict]:
    """Revalidate the reference, descriptive/probability evidence and sensitivity fits."""
    paths = {**SOURCES, **{name: module.REPORT_FILE for name, module in SENSITIVITY.items()}}
    tracked = (*paths.values(), Path("src/phase5_summary.py"), DIRECTORY / "analysis_plan.json",
               DIRECTORY / "literature_review.tex", Path("docs/amarpr.pdf"))
    before = {str(path): fingerprint(path) for path in tracked}
    plan = load_registered_plan()
    results = {name: read_json(path) for name, path in paths.items()}
    rows, flow = load_participants()
    compare(results["description"], describe(rows, flow), "Phase 1 descriptive results")
    samples, flow = load_phase2_samples()
    for name, data in (
        ("events", summarize_events(samples[results["events"]["sample"]])),
        ("bayes", compare_probabilities(summarize_events(samples[results["bayes"]["sample"]]))),
        ("counts", summarize_count(results["counts"]["reference_sample_size"],
                                   results["events"]["contingency_table"]["A"]["total"])),
        ("distributions", summarize_distributions(samples[results["distributions"]["sample"]])),
    ):
        compare(results[name], {"participants_csv_sha256": flow["participants_csv_sha256"],
                                "sample": results[name]["sample"], **data}, f"Phase 2 {name}")
    for name, module in SENSITIVITY.items():
        compare(results[name], module.build_report(), f"Phase 5 {name}")
    require(before == {str(path): fingerprint(path) for path in tracked},
            "A source changed during final validation.")
    return plan, results, before


def build_synthesis() -> dict:
    plan, results, hashes = evidence()
    primary = results["multiple"]["primary_test"]
    interval_normal = results["influence"]["reference"]["coefficients"][THICKNESS]
    intervals = results["influence"]["participant_bootstrap"]
    normal, broad = [results["samples"]["sample_definition"][name] for name in ("reference", "broader")]
    bilateral = results["bilateral"]["models"]["bilateral_mean_on_same_support"]
    omission = results["influence"]["influence"]
    subgroups = results["subgroups"]["groups"]
    models = results["multiple"]["models"]
    comparison = results["comparison"]
    simple = results["simple"]

    def slope(coefficient: dict, n: int) -> dict:
        return {"n": n, "estimate": coefficient["estimate"],
                **{method: coefficient[method]["confidence_interval"] for method in ("model_based", "hc3")}}

    sensitivity = {
        "reference": slope(interval_normal, normal["fit"]["participants"]),
        "broader": slope(broad["coefficients"][THICKNESS], broad["fit"]["participants"]),
        "bilateral": slope(bilateral["coefficients"][THICKNESS], bilateral["fit"]["participants"]),
        "cook_omission": slope(omission["temporary_fit"]["coefficients"][THICKNESS], omission["retained_n"]),
        "mci": slope(subgroups["mci"]["coefficients"][THICKNESS], subgroups["mci"]["participants"]),
    }
    coverage = {
        "general_1": "Related studies: three original articles with a structured comparison.",
        "general_2": "Sample sensitivity: the integrated model on the broader sample.",
        "general_3": "Definition sensitivity: left versus bilateral mean on identical membership.",
        "general_4": "Influence: the full fit and the temporary omission of the five largest flagged Cook distances.",
        "general_5": "Uncertainty: model-based, HC3 and participant-bootstrap intervals.",
        "general_6": "Exploratory subgroups: group-specific slopes and the direct interaction.",
        "general_7": "Final report: synthesis of descriptive, probability, inference and regression evidence.",
        "general_8": "Discussion and conclusion: scope, limitations and the primary answer.",
        "track_c1": "Common-sample model tables: coefficients, uncertainty and fit across the three models.",
        "track_c2": "MRI completeness: raw-UDS denominator and unchanged UDS-only fit on identical membership.",
        "track_c3": "Visit rule: applicability assessed; one supplied pair per person and missing timing fields.",
        "track_c4": "Normal cognition and MCI: one row per participant and an exploratory direct contrast.",
        "track_c5": "Discussion: incremental contribution and limits of combining UDS and MRI.",
    }
    require(coverage.keys() == plan["requirement_map"].keys(), "Final coverage differs from the registered requirements.")

    return {
        "stage": "Phase 5, Step 7",
        "question": plan["reference"]["question"],
        "primary": {"n": plan["reference"]["n"], "estimate": primary["estimate"], "model_based": primary["model_based"],
                    "hc3": primary["hc3_companion"]},
        "description_by_status": [
            {"status": name, "total": item["participants"],
             **{label: {key: item["numeric"][column][key]
                        for key in ("usable", "mean", "sample_sd", "median", "q1", "q3")}
                for label, column in (("score", "ANIMALS"), ("left_thickness", LEFT_THICKNESS))}}
            for name in CATEGORY_LABELS["NACCUDSD"].values()
            for item in (results["description"]["by_cognitive_status"][name],)
        ],
        "probability": {
            "cutoff_mm": results["events"]["events"]["B"]["median_mm"],
            "table": results["events"]["contingency_table"],
            "conditional": results["events"]["conditional_probabilities"],
            "independence": results["bayes"]["independence"], "bayes": results["bayes"]["bayes"],
            "draws": results["counts"]["number_of_draws"], "at_least": results["counts"]["tail_minimum"],
            "binomial_tail": results["counts"]["with_replacement"]["probability_at_least_tail_minimum"],
            "hypergeometric_tail": results["counts"]["without_replacement"]["probability_at_least_tail_minimum"],
            "score_mode": max(results["distributions"]["animal_score_frequency"]["levels"], key=lambda row: row["count"]),
            "thickness_lower_tail": results["distributions"]["thickness_normal_reference"]["quantile_checks"][0],
        },
        "inference": {
            "means": results["means"]["estimates"], "mean_bootstrap": results["means_bootstrap"]["estimates"],
            "score_contrast": results["contrast"]["estimate"],
            "score_contrast_bootstrap": results["contrast_bootstrap"]["estimate"]["percentile_interval"],
            "tests": results["tests"]["tests"],
            "thickness_eta_squared": results["group_thickness"]["descriptive_effect"]["eta_squared"],
        },
        "regression": {
            "models": {name: {"fit": report["fit"], "coefficients": report["coefficients"]}
                       for name, report in (("uds_only", models["uds_only"]), ("mri_only", simple), ("integrated", models["integrated"]))},
            "pearson": simple["correlation"]["coefficient"], "nested_changes": comparison["nested_model_changes"],
            "max_vif": max(item["VIF"] for item in comparison["collinearity"]["integrated"]["vif"].values()),
            "profile": results["prediction"]["intervals"],
        },
        "sensitivity": sensitivity,
        "bootstrap": {"valid_fits_n": intervals["valid_fits_n"], "interval": intervals["confidence_interval"],
                      "se": intervals["empirical_standard_error"]},
        "influence": {key: omission[key] for key in ("flagged_n", "omitted_n", "threshold_strict_gt")},
        "mri_coverage": results["samples"]["mri_completeness"],
        "visits": results["samples"]["visit_rule"],
        "mci_minus_normal": results["subgroups"]["direct_comparison"],
        "coverage": coverage,
        "checks": ["Phase 1 descriptions and Phase 2 event, distribution and count summaries recomputed.",
                   "Phase 3 and Phase 4 evidence checked by the frozen input fingerprints and reference refits.",
                   "Phase 5 sample, bilateral, influence, bootstrap and subgroup reports reproduced numerically."],
        "source_files_sha256": hashes,
    }


def tex_interval(values, digits=3):
    return interval(values["lower"], values["upper"], digits)


def render_tables(result: dict) -> bytes:
    primary, probability, inference = (result[key] for key in ("primary", "probability", "inference"))
    regression, sensitivity = result["regression"], result["sensitivity"]
    fit = regression["models"]
    contrast = result["mci_minus_normal"]
    classical = primary["model_based"]
    adjusted = regression["nested_changes"]["uds_only_to_integrated"]
    with_demographics = regression["nested_changes"]["mri_only_to_integrated"]
    boot = result["bootstrap"]
    scalars = {
        "FinalN": primary["n"], "FinalSlope": f"{primary['estimate']:.3f}",
        "FinalCI": tex_interval(classical["confidence_interval"]),
        "FinalHC": tex_interval(primary["hc3"]["confidence_interval"]),
        "FinalP": f"{classical['p_value_two_sided_nominal']:.3f}",
        "FinalR": f"{regression['pearson']:.3f}", "FinalVIF": f"{regression['max_vif']:.2f}",
        "MRIAddedR": f"{100*adjusted['R_squared_increase']:.3f}",
        "MRIAddedAdjusted": f"{100*adjusted['adjusted_R_squared_change']:.3f}",
        "MRIAddedRMSE": f"{adjusted['in_sample_RMSE_change']:.4f}",
        "UDSAddedR": f"{100*with_demographics['R_squared_increase']:.2f}",
        "CutoffMM": f"{probability['cutoff_mm']:.3f}",
        "ConditionalHi": f"{100*probability['conditional']['A_given_B']['probability']:.2f}",
        "ConditionalLo": f"{100*probability['conditional']['A_given_not_B']['probability']:.2f}",
        "JointObserved": f"{100*probability['independence']['observed_joint_A_and_B']['probability']:.2f}",
        "JointIndep": f"{100*probability['independence']['joint_if_independent']['probability']:.2f}",
        "ConditionalPrior": f"{100*probability['bayes']['prior_A']['probability']:.2f}",
        "BinTail": f"{100*probability['binomial_tail']:.2f}", "HypTail": f"{100*probability['hypergeometric_tail']:.2f}",
        "ModeScore": probability["score_mode"]["value"], "TailObserved": f"{probability['thickness_lower_tail']['empirical_mm']:.3f}",
        "TailFitted": f"{probability['thickness_lower_tail']['fitted_normal_mm']:.3f}",
        "ScoreContrast": f"{inference['score_contrast']['mean_difference']:.3f}",
        "ScoreContrastCI": tex_interval(inference["score_contrast"]["confidence_interval"]),
        "ScoreContrastBoot": tex_interval(inference["score_contrast_bootstrap"]),
        "ScoreContrastP": f"{inference['tests']['score_difference_by_mri']['p_value_holm']:.3f}",
        "StatusF": f"{inference['tests']['thickness_by_status']['statistic']:.2f}",
        "StatusP": p_tex(inference["tests"]["thickness_by_status"]["p_value_holm"], False),
        "StatusEta": f"{100*inference['thickness_eta_squared']:.2f}",
        "ProfileMean": f"{regression['profile']['point_estimate']:.2f}",
        "ProfileMeanCI": tex_interval(regression["profile"]["mean_confidence_interval"], 2),
        "ProfileIndividualCI": tex_interval(regression["profile"]["individual_prediction_interval"], 2),
        "BroaderSlope": f"{sensitivity['broader']['estimate']:.3f}",
        "FlaggedCook": result["influence"]["flagged_n"], "OmittedCook": result["influence"]["omitted_n"],
        "BootstrapN": f"{boot['valid_fits_n']:,}".replace(",", r"\,"),
        "BootstrapCI": tex_interval(boot["interval"]), "BootstrapSE": f"{boot['se']:.3f}",
        "InteractionSlope": f"{contrast['estimate']:.3f}",
        "InteractionCI": tex_interval(contrast["hc3"]["confidence_interval"]),
        "InteractionP": f"{contrast['hc3']['p_value_two_sided_nominal']:.4f}",
    }
    cohort_labels = (r"\NormalLabel", r"\ImpairedLabel", r"\MCILabel", r"\DementiaLabel")
    table_rows = {
        "FinalDescriptionRows": [
            (label, data["total"], data["score"]["usable"],
             rf"{data['score']['mean']:.2f} $\pm$ {data['score']['sample_sd']:.2f}",
             data["left_thickness"]["usable"],
             rf"{data['left_thickness']['mean']:.3f} $\pm$ {data['left_thickness']['sample_sd']:.3f}")
            for label, data in zip(cohort_labels, result["description_by_status"], strict=True)
        ],
        "FinalMeanRows": [
            (label, f"{item['mean']:.{digits}f}", tex_interval(item, digits),
             tex_interval(inference["mean_bootstrap"][name]["percentile_interval"], digits))
            for name, label, digits in (("mean_score", r"\ScoreLabel", 3),
                                        ("mean_thickness", r"\ThicknessLabel", 4))
            for item in (inference["means"][name],)
        ],
        "FinalModelRows": [
            (label, f"{item['R_squared']:.6f}", f"{item['adjusted_R_squared']:.6f}",
             f"{item['residual_SD']:.3f}", f"{item['in_sample_RMSE']:.3f}")
            for label, model in zip((r"\UDSLabel", r"\MRILabel", r"\IntegratedLabel"), fit.values())
            for item in (model["fit"],)
        ],
        "FinalCoefficientRows": [
            (label, *[
                r"\textemdash" if term not in model["coefficients"] else
                rf"${model['coefficients'][term]['estimate']:.3f}$ " + tex_interval(
                    model["coefficients"][term]["model_based"]["confidence_interval"])
                for model in fit.values()
            ])
            for term, label in ((THICKNESS, r"\ThicknessTermLabel"), ("age_c", r"\AgeLabel"),
                                ("female", r"\SexLabel"), ("education_c", r"\EducationLabel"))
        ],
        "FinalSensitivityRows": [
            (label, item["n"], rf"${item['estimate']:.3f}$", tex_interval(item["model_based"]), tex_interval(item["hc3"]))
            for label, item in zip((r"\ReferenceLabel", r"\BroaderLabel", r"\BilateralLabel", r"\CookLabel", r"\MCILabel"),
                                   sensitivity.values(), strict=True)
        ],
    }
    lines = ["% Generated by python3 -m src.phase5_summary. Numbers are shared by the final reports."]
    lines.extend(rf"\newcommand{{\{name}}}{{{value}}}" for name, value in scalars.items())
    lines.extend(rf"\newcommand{{\{name}}}{{%" + "\n" + rows_tex(rows) + "\n}" for name, rows in table_rows.items())
    return ("\n".join(lines) + "\n").encode("utf-8")


def build_figures(result: dict) -> list:
    """Two grayscale views of aggregate model evidence."""
    colors = GRAYSCALE
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42}):
        models = result["regression"]["models"]
        figure, ax = plt.subplots(figsize=(7.8, 2.75), layout="constrained")
        order = ("mri_only", "uds_only", "integrated")
        values = [100 * models[name]["fit"]["R_squared"] for name in order]
        bars = ax.barh([MODEL_LABELS[name] for name in order], values, color=colors.accent, height=.52)
        ax.bar_label(bars, [f"{value:.3f}%" for value in values], padding=5)
        ax.set(xlim=(0, 15.2), xlabel="In-sample variance explained (%)",
               title=f"Three models on the same {result['primary']['n']} participants")
        ax.invert_yaxis()

        forest, bx = plt.subplots(figsize=(7.8, 3.30), layout="constrained")
        values = result["sensitivity"]
        names = ("Reference", "Broader sample", "Bilateral mean", "Cook omission", "MCI subgroup")
        for i, (label, item) in enumerate(zip(names, values.values(), strict=True)):
            point, bounds = item["estimate"], item["hc3"]
            bx.errorbar(point, i, xerr=[[point-bounds["lower"]], [bounds["upper"]-point]],
                        marker="o", color=colors.ink, linestyle="none", capsize=3)
        bx.set_yticks(range(len(names)), [f"{label} (n={item['n']:,})" for label, item in zip(names, values.values())])
        bx.axvline(0, linestyle="--", color=colors.muted, linewidth=1)
        bx.set(xlabel="Difference in correct animal names per 0.1 mm",
               title="Thickness slopes across sensitivity analyses (95% HC3 intervals)")
        bx.invert_yaxis()
        for axis in (ax, bx):
            axis.spines[["top", "right", "left"]].set_visible(False)
            axis.set_axisbelow(True)
            axis.grid(axis="x", color=colors.grid)
    return [("01_model_comparison.pdf", figure), ("02_thickness_sensitivity.pdf", forest)]


def main() -> None:
    result = build_synthesis()
    write_output(REPORT_FILE, (json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode())
    write_output(TABLE_FILE, render_tables(result))
    for name, figure in build_figures(result):
        try:
            save_pdf(figure, name, FIGURE_DIR)
        finally:
            plt.close(figure)


if __name__ == "__main__":
    main()

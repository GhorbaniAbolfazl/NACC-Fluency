"""Build the shared numeric tables for both Phase 4 LaTeX reports."""

import json
from pathlib import Path

from .build_dataset import write_output
from .phase4_plan import fingerprint, require
from .phase4_prediction import load_prediction_results

DIRECTORY = Path("reports/phase4")
TABLE_FILE = DIRECTORY / "tables.tex"
RESULT_FILES = {
    "plan": "analysis_plan", "design": "design_matrices", "simple": "simple_regression",
    "multiple": "multiple_regression", "comparison": "model_comparison",
    "diagnostic": "residual_diagnostics", "prediction": "profile_prediction",
}
MODEL_LABELS = {"uds_only": "UDS-only", "mri_only": "MRI-only", "integrated": "UDS + MRI"}
TERM_LABELS = {"intercept": "Intercept", "thickness_0_1mm": "Thickness", "age_c": "Age",
               "female": "Recorded sex", "education_c": "Education"}


def load_results() -> dict:
    prediction = load_prediction_results()
    results = {key: json.loads((DIRECTORY / (name + ".json")).read_text(encoding="utf-8"))
               for key, name in RESULT_FILES.items() if key != "prediction"}
    for key, name in RESULT_FILES.items():
        if key != "prediction":
            path = DIRECTORY / (name + ".json")
            digest = prediction["source_files_sha256"].get(str(path))
            require(digest == fingerprint(path), f"Source changed during Phase 4 review: {path}")
    results["prediction"] = prediction
    return results


def coefficient_rows(model: dict, method: str):
    for term, coefficient in model["coefficients"].items():
        stats = coefficient[method]
        ci = stats["confidence_interval"]
        yield (TERM_LABELS[term], coefficient["estimate"], stats["standard_error"],
               stats["t_statistic"], stats["p_value_two_sided_nominal"], ci["lower"], ci["upper"],
               stats["p_value_underflow"])


def p_tex(value: float, underflow: bool) -> str:
    if underflow:
        return r"\textemdash"
    if value >= .0001:
        return f"{value:.5f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return rf"\ensuremath{{{mantissa}\times 10^{{{int(exponent)}}}}}"


def interval(lower: float, upper: float, digits: int = 4) -> str:
    return rf"\ensuremath{{[{lower:.{digits}f},\,{upper:.{digits}f}]}}"


def rows_tex(rows) -> str:
    return "\n".join(" & ".join(str(value) for value in row) + r" \\" for row in rows)


def render_tables(results: dict) -> bytes:
    plan, simple, multiple, comparison, diagnostic, prediction = (
        results[name] for name in ("plan", "simple", "multiple", "comparison", "diagnostic", "prediction"))
    sample = plan["sample"]["counts"]
    profile = prediction["profile"]["predictor_values"]
    intervals = prediction["intervals"]
    primary = multiple["primary_test"]
    slope = simple["coefficients"]["thickness_0_1mm"]
    integrated = diagnostic["models"]["integrated"]
    followup = diagnostic["integrated_influence_followup"]
    loo = followup["individual_leave_one_out"]
    joint = followup["joint_union_omission"]
    change = comparison["nested_model_changes"]["uds_only_to_integrated"]
    maximum_vif = max(item["VIF"] for item in comparison["collinearity"]["integrated"]["vif"].values())
    scalars = {
        "MatchedN": sample["matched"], "NormalN": sample["normal_before_filters"],
        "PairedN": sample["paired_normal"], "ModelN": sample["model_normal"],
        "AgeReference": f"{profile['NACCAGE']:.0f}", "EducationReference": profile["EDUC"],
        "ThicknessReference": f"{profile['LH_INFERIORPARIETAL_AVGTH']:.3f}",
        "SimpleIntercept": f"{simple['coefficients']['intercept']['estimate']:.4f}",
        "SimpleSlope": f"{slope['estimate']:.5f}", "PearsonR": f"{simple['correlation']['coefficient']:.5f}",
        "SimpleRtwo": f"{simple['fit']['R_squared']:.6f}",
        "SimpleProbability": f"{slope['model_based']['p_value_two_sided_nominal']:.5f}",
        "PrimarySlope": f"{primary['estimate']:.5f}",
        "PrimaryProbability": f"{primary['model_based']['p_value_two_sided_nominal']:.5f}",
        "PrimaryRobustProbability": f"{primary['hc3_companion']['p_value_two_sided_nominal']:.5f}",
        "PrimaryInterval": interval(**{k: primary["model_based"]["confidence_interval"][k] for k in ("lower", "upper")}, digits=5),
        "PrimaryRobustInterval": interval(**{k: primary["hc3_companion"]["confidence_interval"][k] for k in ("lower", "upper")}, digits=5),
        "DeltaRtwo": f"{change['R_squared_increase']:.6f}", "DeltaAdjustedRtwo": f"{change['adjusted_R_squared_change']:.6f}",
        "MaxVif": f"{maximum_vif:.5f}", "FlaggedN": integrated["review_flags"]["union_count"],
        "LooLow": f"{loo['estimate']['minimum']:.5f}", "LooHigh": f"{loo['estimate']['maximum']:.5f}",
        "JointN": joint["retained"], "JointSlope": f"{joint['thickness']['estimate']:.5f}",
        "PredictedMean": f"{intervals['point_estimate']:.5f}", "PredictionDF": intervals["residual_df"],
        "ResidualVariance": f"{intervals['residual_variance']:.5f}",
        "ProfileVarianceFactor": f"{intervals['profile_variance_factor']:.8f}",
        "PredictionCritical": f"{intervals['t_critical']:.5f}", "WidthRatio": f"{intervals['prediction_to_mean_width_ratio']:.2f}",
    }
    tables = {}
    for prefix, model in (("Simple", simple), ("Uds", multiple["models"]["uds_only"]),
                          ("Integrated", multiple["models"]["integrated"])):
        for method, suffix in (("model_based", "Classical"), ("hc3", "Robust")):
            tables[prefix + suffix] = rows_tex(
                (label, f"{beta:.4f}", f"{se:.4f}", f"{stat:.4f}", p_tex(p, underflow), interval(lo, hi))
                for label, beta, se, stat, p, lo, hi, underflow in coefficient_rows(model, method))
    tables["FitRows"] = rows_tex(
        (label, fit["parameters_including_intercept"], fit["residual_df"], f"{fit['SSE']:.2f}",
         f"{fit['R_squared']:.6f}", f"{fit['adjusted_R_squared']:.6f}",
         f"{fit['residual_SD']:.4f}", f"{fit['in_sample_RMSE']:.4f}")
        for name, label in MODEL_LABELS.items() for fit in [comparison["fit_comparison"][name]])
    tables["VifRows"] = rows_tex(
        (TERM_LABELS[term], f"{uds[term]['VIF']:.4f}" if term in uds else r"\textemdash", f"{data['VIF']:.4f}")
        for uds in [comparison["collinearity"]["uds_only"]["vif"]]
        for term, data in comparison["collinearity"]["integrated"]["vif"].items())
    tables["FlagRows"] = rows_tex(
        (label, flags["counts"]["studentized"], flags["counts"]["leverage"], flags["counts"]["cooks_distance"], flags["union_count"])
        for name, label in MODEL_LABELS.items() for flags in [diagnostic["models"][name]["review_flags"]])
    tables["PredictionRows"] = rows_tex(
        (label, f"{item['standard_error']:.5f}", interval(item["lower"], item["upper"], 5), f"{item['width']:.5f}")
        for label, item in ((r"\MeanTarget", intervals["mean_confidence_interval"]),
                            (r"\IndividualTarget", intervals["individual_prediction_interval"])))
    lines = ["% Generated by python3 -m src.phase4_tables. Shared across the Phase 4 reports."]
    lines.extend(rf"\newcommand{{\{name}}}{{{value}}}" for name, value in scalars.items())
    lines.extend(rf"\newcommand{{\{name}}}{{%" + "\n" + value + "\n}" for name, value in tables.items())
    return ("\n".join(lines) + "\n").encode("utf-8")


def main() -> None:
    content = render_tables(load_results())
    if not TABLE_FILE.exists() or TABLE_FILE.read_bytes() != content:
        write_output(TABLE_FILE, content)


if __name__ == "__main__":
    main()

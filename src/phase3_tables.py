"""Generate numeric macros for the Phase 3 LaTeX report."""

from pathlib import Path

from .build_dataset import write_output
from . import phase3_review as review

TABLE_FILE = Path("reports/phase3/tables.tex")


def load_results() -> dict:
    """Require a current successful integrated review before reading results."""
    audit = review.read_report(review.REPORT_FILE)
    if audit["passed_checks"] != 11 or any(check["status"] != "pass" for check in audit["checks"]):
        raise ValueError("The Phase 3 review is incomplete. Run python3 -m src.phase3_review.")
    fingerprints = audit["source_files_sha256"]
    if any(review.file_hash(Path(name)) != digest for name, digest in fingerprints.items()) or any(
        str(path) not in fingerprints for path in Path("src").glob("phase3_*.py")
    ):
        raise ValueError("The Phase 3 review is stale. Run python3 -m src.phase3_review before building the reports.")
    return {
        "plan": review.read_report(review.PLAN_FILE),
        "means": review.read_report(review.MEANS_FILE),
        "sampling": review.read_report(review.sampling.REPORT_FILE),
        "contrast": review.read_report(review.contrast.REPORT_FILE),
        "bootstrap": review.read_report(review.bootstrap.REPORT_FILE),
        "groups": review.read_report(review.comparison.REPORT_FILE),
        "tests": review.read_report(review.comparison.FAMILY_FILE),
        "review": audit,
    }


def scientific(value: float) -> str:
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"\ensuremath{{{mantissa}\times10^{{{int(exponent)}}}}}"


def table_values(results: dict) -> dict[str, str]:
    plan = results["plan"]
    score, thickness = (results["means"]["estimates"][name] for name in ("mean_score", "mean_thickness"))
    normal = plan["samples"]["paired_normal"]["before_score_and_thickness_eligibility"]
    c = results["contrast"]
    delta = c["estimate"]
    b = results["bootstrap"]["estimate"]
    g = results["groups"]
    anova = g["welch_anova"]
    tests = results["tests"]["tests"]
    values = {
        "PairedN": str(results["means"]["sample_size"]),
        "StatusN": str(g["sample_size"]),
        "NormalN": str(normal),
        "MissingScoreN": str(normal - results["means"]["sample_size"]),
        "ModelN": str(results["review"]["samples"]["planned_phase_four_model_normal"]),
        "ConfidencePct": f"{100*plan['confidence_level']:.0f}",
        "TestAlpha": f"{plan['testing_policy']['familywise_alpha']:.2f}",
        "ThicknessCut": f"{plan['mri_grouping']['cutoff_mm']:.3f}",
        "HigherN": str(c["groups"]["at_or_above_cutoff"]["participants"]),
        "LowerN": str(c["groups"]["below_cutoff"]["participants"]),
        "TiedN": str(c["mri_grouping"]["participants_equal_to_cutoff"]),
        "Replicates": str(plan["bootstrap"]["replicates"]),
        "MeansSeed": str(plan["bootstrap"]["whole_sample"]["seed"]),
        "ContrastSeed": str(plan["bootstrap"]["mri_contrast"]["seed"]),
        "MeanDF": str(score["degrees_of_freedom"]),
        "MeanCritical": f"{score['t_critical']:.6f}",
        "DeltaEstimate": f"{delta['mean_difference']:.4f}",
        "DeltaSE": f"{delta['standard_error']:.5f}",
        "DeltaT": f"{delta['t_statistic']:.4f}",
        "DeltaDF": f"{delta['degrees_of_freedom']:.4f}",
        "DeltaRawP": f"{delta['p_value_raw']:.5f}",
        "DeltaHolmP": f"{tests['score_difference_by_mri']['p_value_holm']:.5f}",
        "DeltaLow": f"{delta['confidence_interval']['lower']:.4f}",
        "DeltaHigh": f"{delta['confidence_interval']['upper']:.4f}",
        "DeltaWidth": f"{delta['confidence_interval']['interval_width']:.4f}",
        "DeltaBootMean": f"{b['bootstrap_mean_difference']:.4f}",
        "DeltaBootSE": f"{b['bootstrap_standard_error']:.5f}",
        "DeltaBootLow": f"{b['percentile_interval']['lower']:.4f}",
        "DeltaBootHigh": f"{b['percentile_interval']['upper']:.4f}",
        "DeltaBootWidth": f"{b['percentile_interval']['interval_width']:.4f}",
        "DeltaWidthChange": f"{100*(1-b['interval_comparison']['percentile_to_welch_width_ratio']):.2f}",
        "AnovaF": f"{anova['f_statistic']:.4f}",
        "AnovaDFOne": str(anova["df_numerator"]),
        "AnovaDFTwo": f"{anova['df_denominator']:.4f}",
        "AnovaRawP": scientific(anova["p_value_raw"]),
        "AnovaHolmP": scientific(tests["thickness_by_status"]["p_value_holm"]),
        "EtaSquared": f"{g['descriptive_effect']['eta_squared']:.4f}",
        "EtaPercent": f"{100*g['descriptive_effect']['eta_squared']:.2f}",
        "VarianceRatio": f"{g['diagnostics']['largest_to_smallest_sample_variance_ratio']:.3f}",
        "SizeRatio": f"{g['diagnostics']['largest_to_smallest_group_size_ratio']:.0f}",
        "ReviewCount": str(results["review"]["passed_checks"]),
    }
    for name, estimate, places in (("Score", score, 4), ("Thickness", thickness, 6)):
        for key, suffix in (("mean", "Mean"), ("sample_sd", "SD"), ("standard_error", "SE"), ("lower", "Low"), ("upper", "High")):
            precision = 5 if name == "Score" and key == "standard_error" else places
            values[name + suffix] = f"{estimate[key]:.{precision}f}"
        boot = results["sampling"]["estimates"]["mean_score" if name == "Score" else "mean_thickness"]
        values[name + "BootSE"] = f"{boot['bootstrap_standard_error']:.{5 if name=='Score' else 6}f}"
        values[name + "BootLow"] = f"{boot['percentile_interval']['lower']:.{places}f}"
        values[name + "BootHigh"] = f"{boot['percentile_interval']['upper']:.{places}f}"
        values[name + "WidthChange"] = f"{100*(1-boot['interval_comparison']['percentile_to_t_width_ratio']):.2f}"
    for key, prefix in (("at_or_above_cutoff", "Higher"), ("below_cutoff", "Lower")):
        s = c["groups"][key]["score_summary"]
        values[prefix + "ScoreMean"] = f"{s['mean']:.4f}"
        values[prefix + "ScoreSD"] = f"{s['sample_sd']:.4f}"
    for code, prefix in (("1", "Normal"), ("2", "Impaired"), ("3", "MCI"), ("4", "Dementia")):
        item = g["groups"][code]
        ci = item["mean_interval"]
        values[prefix + "GroupN"] = str(ci["participants"])
        for key, suffix in (("mean", "Mean"), ("sample_sd", "SD"), ("lower", "Low"), ("upper", "High")):
            values[prefix + suffix] = f"{ci[key]:.4f}"
        values[prefix + "Skew"] = f"{item['empirical_skewness']:.3f}"
        flags = item["summary"]["iqr_review"]
        values[prefix + "Flags"] = str(flags["below"] + flags["above"])
    return values


def main() -> None:
    lines = ["% Shared numeric values generated by python3 -m src.phase3_tables."]
    lines.extend(rf"\newcommand{{\{name}}}{{{value}}}" for name, value in table_values(load_results()).items())
    write_output(TABLE_FILE, ("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()

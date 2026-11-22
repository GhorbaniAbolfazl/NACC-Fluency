"""Review Phase 3 consistency and interpretation without rewriting earlier results."""

import csv
import json
from hashlib import sha256
from importlib.metadata import version
from io import BytesIO, StringIO
from math import isclose
from pathlib import Path
from platform import python_version
from statistics import quantiles, stdev

from .build_dataset import OUTPUT_COLUMNS, build_dataset, write_output
from .descriptive_statistics import FLOW_FILE, PARTICIPANTS_FILE
from .phase1_figures import GRAYSCALE
from .phase2_samples import load_phase2_samples
from .phase3_mean_intervals import PLAN_FILE, REPORT_FILE as MEANS_FILE, load_registered_plan, mean_interval
from . import phase3_sampling as sampling
from . import phase3_mri_contrast as contrast
from . import phase3_contrast_bootstrap as bootstrap
from . import phase3_group_comparison as comparison

from matplotlib import pyplot as plt
from matplotlib.colors import to_rgb

REPORT_FILE = Path("reports/phase3/inference_review.json")
REPORTS = (PLAN_FILE, MEANS_FILE, sampling.REPORT_FILE, contrast.REPORT_FILE, bootstrap.REPORT_FILE, comparison.REPORT_FILE, comparison.FAMILY_FILE)
FIGURES = tuple(sampling.FIGURE_DIR / name for name in sampling.FIGURES.values()) + (
    bootstrap.FIGURE_DIR / bootstrap.FIGURE_NAME,
    comparison.FIGURE_DIR / comparison.MEANS_FIGURE,
    comparison.FIGURE_DIR / comparison.DIAGNOSTICS_FIGURE,
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_report(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_bootstrap(values: list[float], standard_error: float, interval: dict, expected_count: int) -> None:
    """Cross-check retained draws with sample SD and inclusive linear quantiles."""
    require(len(values) == expected_count, "The bootstrap replicate count differs from the plan.")
    require(isclose(stdev(values), standard_error, rel_tol=1e-13), "The saved bootstrap SE differs from its draws.")
    cutpoints = quantiles(values, n=40, method="inclusive")
    require(
        isclose(cutpoints[0], interval["lower"], rel_tol=0, abs_tol=1e-12)
        and isclose(cutpoints[-1], interval["upper"], rel_tol=0, abs_tol=1e-12),
        "The saved percentile interval differs from the retained bootstrap draws.",
    )


def verify_figures(mean_bootstrap: dict, score_bootstrap: dict, group_comparison: dict) -> None:
    """Regenerate PDFs in memory using the same styles and deterministic metadata."""
    require(all(max(to_rgb(color)) == min(to_rgb(color)) for color in GRAYSCALE), "The shared figure palette contains a non-grayscale color.")
    _, groups, _ = comparison.load_status_groups()
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelcolor": GRAYSCALE.ink, "axes.labelpad": 9, "pdf.fonttype": 42, "figure.facecolor": "white"}):
        figures = []
        for name, filename in sampling.FIGURES.items():
            figure = sampling.sampling_figure(name, mean_bootstrap["estimates"][name], mean_bootstrap["bootstrap_means"][name], mean_bootstrap["bootstrap"])
            figures.append((sampling.FIGURE_DIR / filename, figure))
        figures.extend((
            (bootstrap.FIGURE_DIR / bootstrap.FIGURE_NAME, bootstrap.contrast_figure(score_bootstrap)),
            (comparison.FIGURE_DIR / comparison.MEANS_FIGURE, comparison.means_figure(group_comparison)),
            (comparison.FIGURE_DIR / comparison.DIAGNOSTICS_FIGURE, comparison.qq_figure(group_comparison, groups)),
        ))
        try:
            for path, figure in figures:
                with BytesIO() as buffer:
                    figure.savefig(buffer, format="pdf", metadata={"Title": path.name, "CreationDate": None, "ModDate": None})
                    require(buffer.getvalue() == path.read_bytes(), f"The saved figure does not reproduce: {path}")
        finally:
            for _, figure in figures:
                plt.close(figure)


def build_review() -> dict:
    tracked = set(REPORTS + FIGURES + (PARTICIPANTS_FILE, FLOW_FILE))
    tracked.update(Path(name) for name in ("data/udslong_imputed.csv", "data/mriscanclean.csv", "docs/amarpr.pdf", "reports/phase2/event_probabilities.json", "requirements.txt", "src/build_dataset.py", "src/data_audit.py", "src/descriptive_statistics.py", "src/phase1_figures.py", "src/phase2_samples.py"))
    tracked.update(Path("src").glob("phase3_*.py"))
    before = {str(path): file_hash(path) for path in sorted(tracked)}
    checks = []

    def passed(name: str, detail: str) -> None:
        checks.append({"check": name, "status": "pass", "detail": detail})

    rows, reconstructed_flow = build_dataset(Path("data"))
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    table_bytes = buffer.getvalue().encode("utf-8")
    require(table_bytes == PARTICIPANTS_FILE.read_bytes(), "The participant table differs from a raw-data rebuild.")
    reconstructed_flow["participants_csv_sha256"] = sha256(table_bytes).hexdigest()
    require(reconstructed_flow == read_report(FLOW_FILE), "The saved sample flow differs from a raw-data rebuild.")
    passed("raw_data_and_linkage", "The participant table and sample flow reproduce from the supplied CSVs with one matched row per ID.")

    plan = load_registered_plan()
    samples, flow = load_phase2_samples()
    paired = samples["paired_normal"]
    passed("registered_plan", "The saved plan, input hashes, eligibility rules, fixed cutoff and group counts match the current data.")

    means = read_report(MEANS_FILE)
    require(means["sample_size"] == len(paired) and means["sample"] == "paired_normal", "The mean-estimation sample differs from the paired-normal sample.")
    require(means["participants_csv_sha256"] == flow["participants_csv_sha256"] and means["analysis_plan_sha256"] == file_hash(PLAN_FILE), "The mean estimates have different input hashes.")
    for name in sampling.ESTIMANDS:
        specification = plan["estimands"][name]
        expected = {
            "column": specification["column"],
            "parameter": specification["parameter"],
            "unit": specification["unit"],
            **mean_interval([row[specification["column"]] for row in paired], plan["confidence_level"]),
        }
        require(means["estimates"][name] == expected, f"The saved mean estimate differs from the data: {name}")
    passed("mean_estimates", "Both means, sample SDs, SEs, t degrees of freedom and 95% interval endpoints reproduce on the same paired-normal rows.")

    mean_bootstrap = read_report(sampling.REPORT_FILE)
    require(mean_bootstrap == sampling.build_report(), "The whole-sample bootstrap report does not reproduce.")
    passed("whole_sample_bootstrap", "All 10,000 paired mean draws and their summaries reproduce with seed 314159.")
    score_contrast = read_report(contrast.REPORT_FILE)
    require(score_contrast == contrast.build_report(), "The saved Welch score contrast does not reproduce.")
    passed("score_contrast", "The higher-minus-lower contrast, fixed group sizes, Welch statistic, df, raw p-value and interval reproduce.")
    score_bootstrap = read_report(bootstrap.REPORT_FILE)
    require(score_bootstrap == bootstrap.build_report(), "The within-group contrast bootstrap does not reproduce.")
    passed("within_group_bootstrap", "All 10,000 stratified score differences and their summaries reproduce with seed 271828 and fixed group sizes.")
    group_comparison = read_report(comparison.REPORT_FILE)
    require(group_comparison == comparison.build_report(), "The cognitive-group comparison does not reproduce.")
    passed("cognitive_group_comparison", "All four group summaries, individual mean intervals, Welch ANOVA, diagnostics and descriptive eta-squared reproduce.")
    family = read_report(comparison.FAMILY_FILE)
    require(family == comparison.build_family_report(group_comparison), "The completed Holm table does not reproduce.")
    passed("completed_test_family", "Exactly the two registered tests appear with current raw p-values, Holm adjustment and familywise decisions.")

    for name in sampling.ESTIMANDS:
        estimate = mean_bootstrap["estimates"][name]
        verify_bootstrap(mean_bootstrap["bootstrap_means"][name], estimate["bootstrap_standard_error"], estimate["percentile_interval"], plan["bootstrap"]["replicates"])
    b = score_bootstrap["estimate"]
    verify_bootstrap(score_bootstrap["bootstrap_mean_differences"], b["bootstrap_standard_error"], b["percentile_interval"], plan["bootstrap"]["replicates"])
    c = score_contrast["estimate"]
    ci = c["confidence_interval"]
    alpha = family["policy"]["familywise_alpha"]
    zero_excluded = ci["lower"] > 0 or ci["upper"] < 0
    require(zero_excluded == (c["p_value_raw"] < 1 - ci["confidence_level"]), "The two-sided Welch interval and raw p-value give incompatible decisions.")
    h, l = (score_contrast["groups"][name]["score_summary"] for name in contrast.GROUP_ORDER)
    require(isclose(c["mean_difference"], h["mean"] - l["mean"], abs_tol=1e-13), "The contrast direction differs from its group means.")
    for test in family["tests"].values():
        require(test["p_value_holm"] >= test["p_value_raw"] and test["reject_holm"] == (test["p_value_holm"] <= alpha), "A Holm p-value or decision is inconsistent.")
    effect = group_comparison["descriptive_effect"]
    require(isclose(effect["between_group_sum_of_squares"] + effect["within_group_sum_of_squares"], effect["total_sum_of_squares"], rel_tol=1e-13), "The sums of squares do not decompose correctly.")
    require(isclose(effect["eta_squared"], effect["between_group_sum_of_squares"] / effect["total_sum_of_squares"], rel_tol=1e-13), "Descriptive eta-squared differs from its sums of squares.")
    passed("numerical_relationships", "Independent inclusive quantiles and SD checks agree with retained draws; contrast signs, CI/test decisions, Holm monotonicity and the SS decomposition agree.")
    verify_figures(mean_bootstrap, score_bootstrap, group_comparison)
    passed("figures", "All five vector PDFs reproduce byte for byte with the shared grayscale palette.")
    require(before == {str(path): file_hash(path) for path in sorted(tracked)}, "Inputs or analysis files changed during review.")
    passed("stable_inputs", "Tracked inputs, analysis code, prior results and figures remained unchanged throughout the review.")

    return {
        "stage": "Phase 3, Step 7",
        "review_status": "Numerical and structural consistency checks passed; interpretation remains conditional on the stated assumptions and limitations.",
        "passed_checks": len(checks),
        "checks": checks,
        "source_files_sha256": before,
        "environment": {"python": python_version(), "scipy": version("scipy"), "matplotlib": version("matplotlib")},
        "samples": {
            "paired_normal": len(paired),
            "mri_groups": {name: score_contrast["groups"][name]["participants"] for name in contrast.GROUP_ORDER},
            "status_and_mri": group_comparison["sample_size"],
            "cognitive_groups": {code: group["summary"]["usable"] for code, group in group_comparison["groups"].items()},
            "planned_phase_four_model_normal": flow["samples"]["model_normal"]["included"],
        },
        "findings": {
            "mean_estimates": means["estimates"],
            "within_normal_cognition": {
                "mean_difference": c["mean_difference"],
                "unit": "Correct animal names in 60 seconds, higher MRI group minus lower MRI group",
                "welch_interval": ci,
                "bootstrap_interval": b["percentile_interval"],
                "p_value_holm": family["tests"]["score_difference_by_mri"]["p_value_holm"],
                "interpretation": "The registered grouped comparison gives insufficient evidence of a score difference at familywise alpha 0.05. Its intervals leave uncertainty about the size and direction of the contrast; practical equivalence has not been assessed.",
            },
            "across_cognitive_groups": {
                "welch_anova": group_comparison["welch_anova"],
                "p_value_holm": family["tests"]["thickness_by_status"]["p_value_holm"],
                "descriptive_eta_squared": effect["eta_squared"],
                "interpretation": "The model-based omnibus test supports a difference in at least one mean thickness across the four cognitive groups. It does not identify a specific pair. Eta-squared describes this sample's SS partition and is sensitive to group composition.",
            },
            "bootstrap_comparisons": {
                "score_mean_se_ratio": mean_bootstrap["estimates"]["mean_score"]["bootstrap_to_analytic_se_ratio"],
                "thickness_mean_se_ratio": mean_bootstrap["estimates"]["mean_thickness"]["bootstrap_to_analytic_se_ratio"],
                "thickness_mean_interval_width_ratio": mean_bootstrap["estimates"]["mean_thickness"]["interval_comparison"]["percentile_to_t_width_ratio"],
                "score_contrast_se_ratio": b["bootstrap_to_welch_se_ratio"],
                "score_contrast_interval_width_ratio": b["interval_comparison"]["percentile_to_welch_width_ratio"],
                "interpretation": "Standard and percentile procedures give similar uncertainty estimates conditional on these observed data. Resampling the same extract provides an internal comparison of methods; it is not an independent replication of the study.",
            },
        },
        "interpretation_rules": {
            "different_questions": "The score contrast compares MRI-defined groups within normal cognition; the omnibus test compares MRI means across UDS cognitive groups. They address different parameters in different eligibility samples.",
            "confidence_intervals": "Interpret the 95% intervals through repeated-sampling coverage under the model assumptions. Each interval concerns its stated population parameter in its original units; it is an individual interval.",
            "p_values": "Interpret raw and adjusted p-values under the null and method assumptions, alongside effect sizes and uncertainty. They do not measure the probability that a hypothesis is true or establish practical importance.",
            "continuous_question": "The primary Phase 4 parameter is the adjusted score difference per 0.1 mm of thickness. A binary MRI-group contrast is a different parameter and leaves that continuous association to be estimated.",
            "observational_scope": "All Phase 3 comparisons are exploratory and unadjusted. The supplied single-record participant data support association descriptions under the stated eligibility rules.",
        },
        "material_limitations": {
            "sampling_frame": plan["sampling_model"]["scope"],
            "independence": "Unique participant IDs establish one supplied record per person. Independence across people remains a modeling assumption; the review cannot validate the sampling design.",
            "source_preparation": "Preparation and imputation details and value-level masks are unavailable. Estimates and resampling condition on the supplied values; additional source-preparation uncertainty cannot be quantified from these files.",
            "measurement_timing_and_qc": "Clinical and MRI timing and scan-level QC details are unavailable. Numerical usability checks do not establish measurement quality or temporal alignment.",
            "small_group": group_comparison["diagnostics"]["small_group"],
            "normality": group_comparison["diagnostics"]["normality"],
            "cutoff": "The pooled Phase 2 median is held fixed. Inference is conditional on that cutoff and, for the stratified bootstrap, the observed group sizes.",
            "covariates": "Age, sex and education will be incorporated in Phase 4 on the common eligible modeling sample.",
        },
        "next_step": "Assemble the numeric tables, LaTeX report and notebook, then check the Phase 3 deliverables against the course guide.",
    }


def main() -> None:
    write_output(REPORT_FILE, comparison.encode_report(build_review()))


if __name__ == "__main__":
    main()

"""Validate and record the analysis samples and inferential plan for Phase 3."""

import csv
import json
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from io import StringIO
from pathlib import Path
from statistics import median

from .build_dataset import LEFT_THICKNESS, PAIR_FIELDS, STATUS, index_participants, write_output
from .data_audit import CATEGORY_LABELS
from .descriptive_statistics import PARTICIPANTS_FILE
from .phase2_samples import load_phase2_samples

REPORT_FILE = Path("reports/phase3/analysis_plan.json")
EVENTS_FILE = Path("reports/phase2/event_probabilities.json")
GUIDE_FILE = Path("docs/amarpr.pdf")
THICKNESS_CUTOFF = Decimal("2.319")
CONFIDENCE_LEVEL = 0.95
ALPHA = 0.05
BOOTSTRAP_REPLICATES = 10000


def define_samples() -> tuple[dict, dict]:
    """Check participant identifiers, phase-one membership and the fixed phase-two cutoff."""
    samples, flow = load_phase2_samples()
    table_bytes = PARTICIPANTS_FILE.read_bytes()
    if sha256(table_bytes).hexdigest() != flow["participants_csv_sha256"]:
        raise ValueError("The participant table changed during validation. Rebuild and retry.")
    records = list(csv.DictReader(StringIO(table_bytes.decode("utf-8")), strict=True))
    index_participants(records, str(PARTICIPANTS_FILE))

    events = json.loads(EVENTS_FILE.read_text(encoding="utf-8"))
    reference = samples["status_and_mri"]
    observed_median = median(Decimal(str(row[LEFT_THICKNESS])) for row in reference)
    if (
        events["participants_csv_sha256"] != flow["participants_csv_sha256"]
        or events["sample"] != "status_and_mri"
        or events["sample_size"] != len(reference)
        or events["events"]["B"]["column"] != LEFT_THICKNESS
        or Decimal(str(events["events"]["B"]["median_mm"])) != THICKNESS_CUTOFF
        or observed_median != THICKNESS_CUTOFF
    ):
        raise ValueError("The phase-two reference differs from the fixed Phase 3 cutoff. Review the analysis plan.")

    paired = samples["paired_normal"]
    upper = [row for row in paired if Decimal(str(row[LEFT_THICKNESS])) >= THICKNESS_CUTOFF]
    lower = [row for row in paired if Decimal(str(row[LEFT_THICKNESS])) < THICKNESS_CUTOFF]
    group_counts = Counter(row[STATUS] for row in reference)
    if min(len(upper), len(lower), *(group_counts[code] for code in CATEGORY_LABELS[STATUS])) < 2:
        raise ValueError("Each planned comparison group must contain at least two participants.")

    profile = {
        "paired_normal": {
            "use": "Mean score, mean thickness, sampling distributions and the score contrast across MRI groups",
            **flow["samples"]["paired_normal"],
            "before_score_and_thickness_eligibility": flow["matched_normal_cognition_before_eligibility"],
            "mri_groups": {
                "below_cutoff": len(lower),
                "at_or_above_cutoff": len(upper),
            },
        },
        "status_and_mri": {
            "use": "Supplementary comparison of mean MRI thickness across all four UDS cognitive status groups",
            **flow["samples"]["status_and_mri"],
        },
    }
    return profile, flow


def build_plan() -> dict:
    profile, flow = define_samples()
    return {
        "plan_version": 1,
        "stage": "Phase 3, Step 1",
        "analysis_status": (
            "Specified before computing Phase 3 tests and confidence intervals. "
            "Exploratory analyses informed by Phase 1 and Phase 2 descriptions of the same extract."
        ),
        "inputs": {
            "participant_table": str(PARTICIPANTS_FILE),
            "participants_csv_sha256": flow["participants_csv_sha256"],
            "phase_two_events": str(EVENTS_FILE),
            "phase_two_events_sha256": sha256(EVENTS_FILE.read_bytes()).hexdigest(),
            "course_guide": str(GUIDE_FILE),
            "course_guide_sha256": sha256(GUIDE_FILE.read_bytes()).hexdigest(),
        },
        "sampling_model": {
            "unit": "One participant, using the single supplied UDS record and the single supplied MRI record",
            "linkage": flow["linkage"],
            "target": (
                "Means and contrasts in the source population represented by the teaching extract, "
                "conditional on each stated eligibility rule."
            ),
            "assumptions": "Participant observations are treated as independent under the sampling model.",
            "scope": (
                "The sampling frame and selection probabilities are undocumented. "
                "Confidence intervals depend on the sampling assumptions and the stated eligibility rules."
            ),
            "visit_rule": (
                "Use the supplied record for each participant. Clinical visit date and number, "
                "scan date and the interval between measurements are unavailable."
            ),
            "value_rules": (
                "Reuse the Phase 1 usability rules and sample flags. Retain measurements outside IQR review fences."
            ),
        },
        "samples": profile,
        "confidence_level": CONFIDENCE_LEVEL,
        "mri_grouping": {
            "column": LEFT_THICKNESS,
            "cutoff_mm": float(THICKNESS_CUTOFF),
            "source": "Pooled median in the Phase 2 status_and_mri sample",
            "comparison_order": "At or above cutoff minus below cutoff",
            "ties": "Include values equal to the cutoff in the at-or-above group",
            "bootstrap_rule": "Keep the cutoff fixed across bootstrap draws",
        },
        "estimands": {
            "mean_score": {
                "sample": "paired_normal",
                "parameter": "Population mean ANIMALS score under the paired_normal eligibility rules",
                "estimator": "Arithmetic mean",
                "column": PAIR_FIELDS[0],
                "unit": "Correct responses in 60 seconds",
                "standard_interval": "Two-sided Student t interval for a mean using sample SD / sqrt(n)",
            },
            "mean_thickness": {
                "sample": "paired_normal",
                "parameter": "Population mean left inferior parietal cortical thickness under the paired_normal rules",
                "estimator": "Arithmetic mean",
                "column": LEFT_THICKNESS,
                "unit": "mm",
                "standard_interval": "Two-sided Student t interval for a mean using sample SD / sqrt(n)",
            },
            "score_difference_by_mri": {
                "sample": "paired_normal",
                "parameter": "Mean score at or above the cutoff minus mean score below the cutoff",
                "estimator": "Difference between the two independent-group sample means",
                "unit": "Correct responses",
                "standard_interval": "Welch t interval using Satterthwaite degrees of freedom",
            },
            "thickness_by_status": {
                "sample": "status_and_mri",
                "parameter": "The four population mean thicknesses defined by UDS cognitive status",
                "estimator": "Arithmetic mean within each cognitive status group",
                "unit": "mm",
                "group_codes": list(CATEGORY_LABELS[STATUS]),
                "standard_interval": "Separate Student t intervals for the group means",
            },
        },
        "hypothesis_tests": {
            "score_difference_by_mri": {
                "role": "Main Phase 3 cross-modal contrast in the normal-cognition sample",
                "sample": "paired_normal",
                "null": "The mean score difference between the MRI groups is zero",
                "alternative": "The mean score difference between the MRI groups is nonzero",
                "method": "Two-sided independent-samples Welch t test",
                "effect": "Raw difference in mean score with its confidence interval",
            },
            "thickness_by_status": {
                "role": "Supplementary omnibus cross-modal comparison across cognitive groups",
                "sample": "status_and_mri",
                "null": "All four mean thicknesses are equal",
                "alternative": "At least one mean thickness differs",
                "method": "Welch one-way ANOVA",
                "effect": "Group means and descriptive sample eta-squared, computed from sums of squares",
                "diagnostics": (
                    "Review each group's size, distribution, spread and extreme values. "
                    "Keep the seven-person impairment-without-MCI group with an explicit precision limitation."
                ),
            },
        },
        "testing_policy": {
            "family": ["score_difference_by_mri", "thickness_by_status"],
            "familywise_alpha": ALPHA,
            "adjustment": "Holm adjustment across the two planned tests",
            "reporting": "Report the test statistic, degrees of freedom, raw p-value and Holm-adjusted p-value",
            "interval_scope": "The 95% intervals are individual intervals.",
        },
        "bootstrap": {
            "replicates": BOOTSTRAP_REPLICATES,
            "confidence_interval": "Percentile interval with linear interpolation at the 2.5th and 97.5th percentiles",
            "whole_sample": {
                "sample": "paired_normal",
                "seed": 314159,
                "draw_rule": "Resample n complete participant rows with replacement",
                "statistics": ["mean_score", "mean_thickness"],
                "purpose": (
                    "Show sampling distributions and empirical standard errors. "
                    "Compare the thickness mean's bootstrap interval with the t interval, "
                    "given the lower-tail asymmetry seen in Phase 2."
                ),
            },
            "mri_contrast": {
                "sample": "paired_normal",
                "seed": 271828,
                "draw_rule": (
                    "Resample complete participant rows with replacement separately within the two MRI groups, "
                    "preserving each group's observed size."
                ),
                "statistic": "score_difference_by_mri",
                "comparison": "Compare the bootstrap percentile interval with the Welch interval",
            },
        },
        "interpretation": {
            "confidence_interval": (
                "Under repeated sampling with the stated assumptions, a nominal 95% procedure "
                "covers the target population parameter in about 95% of repeated samples."
            ),
            "p_value": (
                "Under the null hypothesis and test assumptions, the probability of a test statistic "
                "at least as incompatible with the null as the observed statistic."
            ),
            "practical_importance": "Interpret the raw effect and interval width in their measurement units.",
            "phase_four": "Estimate the adjusted score difference per 0.1 mm using age, sex and education in Phase 4.",
        },
    }


def main() -> None:
    plan = build_plan()
    output = json.dumps(plan, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

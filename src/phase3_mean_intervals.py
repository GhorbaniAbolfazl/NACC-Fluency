"""Estimate both paired-normal means and their individual Student t intervals."""

import json
from hashlib import sha256
from math import isfinite, sqrt
from pathlib import Path
from statistics import fmean, stdev

from scipy.stats import t

from .build_dataset import write_output
from .phase2_samples import load_phase2_samples
from .phase3_plan import REPORT_FILE as PLAN_FILE, build_plan

REPORT_FILE = Path("reports/phase3/mean_intervals.json")
SAMPLE_NAME = "paired_normal"
ESTIMANDS = ("mean_score", "mean_thickness")


def load_registered_plan() -> dict:
    """Require the saved plan to agree with the current inputs and planning code."""
    plan = json.loads(PLAN_FILE.read_text(encoding="utf-8"))
    if plan != build_plan():
        raise ValueError("The saved Phase 3 plan differs from the current specification or inputs. Review the plan.")
    return plan


def mean_interval(values: list[float], confidence_level: float) -> dict:
    """Use the sample standard deviation and n - 1 degrees of freedom."""
    n = len(values)
    if n < 2 or not all(isfinite(value) for value in values):
        raise ValueError("A mean interval requires at least two finite observations.")
    if not 0 < confidence_level < 1:
        raise ValueError("The confidence level must be strictly between zero and one.")

    mean = fmean(values)
    sample_sd = stdev(values)
    standard_error = sample_sd / sqrt(n)
    degrees_of_freedom = n - 1
    critical_value = float(t.ppf((1 + confidence_level) / 2, degrees_of_freedom))
    margin = critical_value * standard_error
    return {
        "participants": n,
        "mean": mean,
        "sample_sd": sample_sd,
        "standard_error": standard_error,
        "degrees_of_freedom": degrees_of_freedom,
        "confidence_level": confidence_level,
        "t_critical": critical_value,
        "margin_of_error": margin,
        "lower": mean - margin,
        "upper": mean + margin,
        "interval_width": 2 * margin,
    }


def main() -> None:
    plan = load_registered_plan()
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    if (
        flow["participants_csv_sha256"] != plan["inputs"]["participants_csv_sha256"]
        or len(rows) != plan["samples"][SAMPLE_NAME]["included"]
    ):
        raise ValueError("The mean-interval sample differs from the registered Phase 3 plan.")

    estimates = {}
    for name in ESTIMANDS:
        specification = plan["estimands"][name]
        estimates[name] = {
            "column": specification["column"],
            "parameter": specification["parameter"],
            "unit": specification["unit"],
            **mean_interval([row[specification["column"]] for row in rows], plan["confidence_level"]),
        }

    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "analysis_plan_sha256": sha256(PLAN_FILE.read_bytes()).hexdigest(),
        "sample": SAMPLE_NAME,
        "sample_size": len(rows),
        "method": {
            "estimator": "Arithmetic mean",
            "standard_error": "Sample standard deviation divided by sqrt(n)",
            "interval": "Mean plus or minus t(1 - alpha / 2, n - 1) times the standard error",
            "interval_scope": plan["testing_policy"]["interval_scope"],
            "assumptions": (
                "Independent participants under the sampling model. The Student t interval is exact for "
                "normal observations and is used here as a large-sample approximation for the means."
            ),
            "value_rules": plan["sampling_model"]["value_rules"],
            "coverage_interpretation": plan["interpretation"]["confidence_interval"],
            "generalization": plan["sampling_model"]["scope"],
            "t_quantile": "scipy.stats.t.ppf",
        },
        "estimates": estimates,
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()

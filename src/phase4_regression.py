"""Fit full-rank participant-level OLS with model-based and HC3 uncertainty."""

from dataclasses import dataclass
from math import fsum, isclose, isfinite, sqrt
from sys import float_info

from scipy.linalg import qr, solve_triangular
from scipy.stats import t

from .phase4_design import ModelDesign, matrix_diagnostics
from .phase4_plan import require


@dataclass(frozen=True)
class OLSFit:
    """Private in-memory fit arrays; public reports export aggregate results only."""

    name: str
    columns: tuple[str, ...]
    participants: int
    residual_df: int
    coefficients: tuple[float, ...]
    unscaled_covariance: tuple[tuple[float, ...], ...]
    model_based_covariance: tuple[tuple[float, ...], ...]
    hc3_covariance: tuple[tuple[float, ...], ...]
    fitted_values: tuple[float, ...]
    residuals: tuple[float, ...]
    leverages: tuple[float, ...]
    sse: float
    sst: float
    explained_ss: float
    r_squared: float
    adjusted_r_squared: float
    residual_sd: float
    in_sample_rmse: float
    qr_pivot_order: tuple[int, ...]


def fit_ols(design: ModelDesign) -> OLSFit:
    """Solve pivoted QR without inverting X'X or changing the registered columns."""
    matrix_diagnostics(design)
    n, p = len(design.matrix), len(design.columns)
    require(len(set(design.response)) > 1, "A varying response is required for R-squared and the planned inference.")
    q, r, pivot = qr(design.matrix, mode="economic", pivoting=True, check_finite=True)
    qty = [fsum(float(q[i, j]) * design.response[i] for i in range(n)) for j in range(p)]
    permuted_beta = solve_triangular(r, qty, check_finite=True)
    pivot = tuple(int(index) for index in pivot)
    original_to_pivoted = tuple(pivot.index(index) for index in range(p))
    beta = tuple(float(permuted_beta[j]) for j in original_to_pivoted)
    fitted = tuple(fsum(value * coefficient for value, coefficient in zip(row, beta)) for row in design.matrix)
    residuals = tuple(y - fitted_y for y, fitted_y in zip(design.response, fitted))
    sse = fsum(value * value for value in residuals)
    mean_y = fsum(design.response) / n
    sst = fsum((value - mean_y) ** 2 for value in design.response)
    exact_fit_tolerance = (32 * float_info.epsilon) ** 2 * max(1.0, fsum(y * y for y in design.response))
    require(isfinite(sse) and sse > exact_fit_tolerance and sst > 0,
            "Residual variation must exceed machine-scale exact-fit error; do not report degenerate t or HC3 inference.")
    explained = fsum((value - mean_y) ** 2 for value in fitted)
    df = n - p
    variance = sse / df
    identity = [[float(i == j) for j in range(p)] for i in range(p)]
    r_inverse = solve_triangular(r, identity, check_finite=True)
    # Rows of R^-1 Q' are returned to the original, registered coefficient order.
    bread_rows = tuple(tuple(fsum(float(r_inverse[j, k]) * float(q[i, k]) for k in range(p))
                             for i in range(n)) for j in original_to_pivoted)
    gram_inverse = tuple(tuple(fsum(float(r_inverse[j, k]) * float(r_inverse[l, k]) for k in range(p))
                               for l in original_to_pivoted) for j in original_to_pivoted)
    model_cov = tuple(tuple(variance * value for value in row) for row in gram_inverse)
    leverage = tuple(fsum(float(q[i, j]) ** 2 for j in range(p)) for i in range(n))
    require(all(0 <= value < 1 - 16 * float_info.epsilon for value in leverage),
            "HC3 requires leverage strictly below one; review an unsupported high-leverage design.")
    hc3_scale = tuple((error / (1 - h)) ** 2 for error, h in zip(residuals, leverage))
    hc3_cov = tuple(tuple(fsum(bread_rows[j][i] * bread_rows[k][i] * hc3_scale[i] for i in range(n))
                          for k in range(p)) for j in range(p))
    for covariance in (model_cov, hc3_cov):
        require(all(isfinite(value) for row in covariance for value in row)
                and all(covariance[j][j] > 0 for j in range(p)), "Coefficient uncertainty must be finite and positive.")
    scale = max(1.0, sst)
    require(abs(sse + explained - sst) <= 1e-10 * scale, "OLS sum-of-squares decomposition failed.")
    require(abs(fsum(residuals)) <= 1e-9 * max(1.0, fsum(abs(y) for y in design.response)),
            "Residuals do not sum to zero for an intercept-containing fit.")
    for j in range(p):
        dot = fsum(row[j] * error for row, error in zip(design.matrix, residuals))
        tolerance = 1e-9 * max(1.0, fsum(abs(row[j] * y) for row, y in zip(design.matrix, design.response)))
        require(abs(dot) <= tolerance, "A residual is not numerically orthogonal to its model column.")
    require(isclose(fsum(leverage), p, rel_tol=1e-11, abs_tol=1e-11), "The hat-matrix trace differs from the design rank.")
    return OLSFit(design.name, design.columns, n, df, beta, gram_inverse, model_cov, hc3_cov,
                  fitted, residuals, leverage, sse, sst, explained, 1 - sse / sst,
                  1 - (sse / df) / (sst / (n - 1)), sqrt(variance), sqrt(sse / n), pivot)


def coefficient_table(fit: OLSFit, confidence_level: float = 0.95) -> dict:
    """Report both preplanned covariance methods with an explicit t reference."""
    require(0 < confidence_level < 1, "Confidence level must lie strictly between zero and one.")
    critical = float(t.ppf((1 + confidence_level) / 2, fit.residual_df))
    table = {}
    for j, name in enumerate(fit.columns):
        methods = {}
        for label, covariance in (("model_based", fit.model_based_covariance), ("hc3", fit.hc3_covariance)):
            se = sqrt(covariance[j][j])
            statistic = fit.coefficients[j] / se
            require(isfinite(statistic), "A coefficient t statistic must be finite.")
            probability = float(2 * t.sf(abs(statistic), fit.residual_df))
            methods[label] = {
                "standard_error": se, "t_statistic": statistic, "degrees_of_freedom": fit.residual_df,
                "p_value_two_sided_nominal": probability, "p_value_underflow": probability == 0,
                "confidence_interval": {"confidence_level": confidence_level,
                                        "lower": fit.coefficients[j] - critical * se,
                                        "upper": fit.coefficients[j] + critical * se},
            }
        table[name] = {"estimate": fit.coefficients[j], **methods}
    return table


def fit_summary(fit: OLSFit) -> dict:
    return {"participants": fit.participants, "parameters_including_intercept": len(fit.columns),
            "residual_df": fit.residual_df, "SSE": fit.sse, "SST": fit.sst,
            "explained_sum_of_squares": fit.explained_ss,
            "R_squared": fit.r_squared, "adjusted_R_squared": fit.adjusted_r_squared,
            "residual_SD": fit.residual_sd, "in_sample_RMSE": fit.in_sample_rmse}

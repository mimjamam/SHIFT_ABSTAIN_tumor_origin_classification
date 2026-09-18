"""Selective-prediction metrics: risk-coverage curves, AURC, selective
accuracy at fixed coverage, and expected calibration error (ECE)."""
import numpy as np


def risk_coverage_curve(scores: np.ndarray, correct: np.ndarray):
    """scores: higher = more confident. correct: bool array, same length.
    Returns (coverage, risk) arrays: for each k=1..n, coverage=k/n and
    risk=error rate among the k most-confident samples (ties broken by
    original order)."""
    order = np.argsort(-scores, kind="stable")
    correct_sorted = correct[order]
    n = len(scores)
    cum_correct = np.cumsum(correct_sorted)
    k = np.arange(1, n + 1)
    coverage = k / n
    risk = 1.0 - cum_correct / k
    return coverage, risk


def aurc(coverage: np.ndarray, risk: np.ndarray) -> float:
    """Area under the risk-coverage curve, as the mean risk over all n
    discrete coverage levels (Geifman & El-Yaniv 2017 definition)."""
    return float(np.mean(risk))


def selective_accuracy_at_coverage(scores: np.ndarray, correct: np.ndarray, coverage_level: float) -> float:
    n = len(scores)
    k = max(1, int(round(coverage_level * n)))
    order = np.argsort(-scores, kind="stable")
    return float(correct[order[:k]].mean())


def expected_calibration_error(confidences: np.ndarray, correct: np.ndarray, n_bins: int = 15) -> float:
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    n = len(confidences)
    ece = 0.0
    for i in range(n_bins):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        if i == 0:
            mask = (confidences >= lo) & (confidences <= hi)
        else:
            mask = (confidences > lo) & (confidences <= hi)
        count = mask.sum()
        if count == 0:
            continue
        acc_bin = correct[mask].mean()
        conf_bin = confidences[mask].mean()
        ece += (count / n) * abs(acc_bin - conf_bin)
    return float(ece)

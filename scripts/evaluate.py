"""
Phase 4: turn Phase 3's saved predictions into the actual evaluation the
project is about -- risk-coverage curves, AURC, selective accuracy at
80/90/95% coverage, and ECE, computed separately for the in-distribution
(exome) test set and the shifted (panel) test set.

Usage: python3 evaluate.py
Reads: results/predictions/{test_exome_id,test_panel_shifted}.npz
Writes: results/metrics_evaluation.json, results/accuracy_vs_coverage.csv,
        results/figures/*.png
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from data import RESULTS_DIR
from selective_metrics import aurc, expected_calibration_error, risk_coverage_curve, selective_accuracy_at_coverage
from uncertainty import ensemble_disagreement, predictive_entropy, softmax_response

PRED_DIR = os.path.join(RESULTS_DIR, "predictions")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
SPLITS = {"test_exome_id": "In-distribution (exome)", "test_panel_shifted": "Shifted (panel)"}
COVERAGE_LEVELS = [0.80, 0.90, 0.95]
COVERAGE_TABLE_LEVELS = [round(x, 1) for x in np.arange(0.1, 1.01, 0.1)]


def load_split(name):
    d = np.load(os.path.join(PRED_DIR, f"{name}.npz"))
    return d["member_probs"], d["mc_probs"], d["y_true"]


def scoring_methods(member_probs, mc_probs, y_true):
    """Returns {method_name: (scores, correct)} -- each method paired with
    the correctness array for the prediction it actually accompanies."""
    mean_probs = member_probs.mean(axis=0)
    mc_mean_probs = mc_probs.mean(axis=0)
    correct_ensemble = (mean_probs.argmax(axis=-1) == y_true)
    correct_mc = (mc_mean_probs.argmax(axis=-1) == y_true)

    return {
        "softmax_response": (softmax_response(mean_probs), correct_ensemble),
        "predictive_entropy": (predictive_entropy(mean_probs), correct_ensemble),
        "ensemble_disagreement": (ensemble_disagreement(member_probs), correct_ensemble),
        "mc_dropout_disagreement": (ensemble_disagreement(mc_probs), correct_mc),
    }


def calibration_variants(member_probs, mc_probs, y_true):
    """Returns {variant_name: (confidence, correct)} for ECE / reliability."""
    single_probs = member_probs[0]
    mean_probs = member_probs.mean(axis=0)
    mc_mean_probs = mc_probs.mean(axis=0)
    return {
        "single_model": (single_probs.max(axis=-1), single_probs.argmax(axis=-1) == y_true),
        "deep_ensemble": (mean_probs.max(axis=-1), mean_probs.argmax(axis=-1) == y_true),
        "mc_dropout": (mc_mean_probs.max(axis=-1), mc_mean_probs.argmax(axis=-1) == y_true),
    }


def reliability_bins(confidence, correct, n_bins=10):
    edges = np.linspace(0, 1, n_bins + 1)
    centers, accs, counts = [], [], []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        mask = (confidence > lo) & (confidence <= hi) if i > 0 else (confidence >= lo) & (confidence <= hi)
        if mask.sum() == 0:
            continue
        centers.append((lo + hi) / 2)
        accs.append(correct[mask].mean())
        counts.append(mask.sum())
    return np.array(centers), np.array(accs), np.array(counts)


if __name__ == "__main__":
    os.makedirs(FIG_DIR, exist_ok=True)
    metrics = {}
    table_rows = []

    fig_rc, axes_rc = plt.subplots(1, 2, figsize=(11, 4.5))
    fig_cal, axes_cal = plt.subplots(1, 2, figsize=(11, 4.5))

    for split_idx, (split_name, split_label) in enumerate(SPLITS.items()):
        member_probs, mc_probs, y_true = load_split(split_name)
        methods = scoring_methods(member_probs, mc_probs, y_true)
        split_metrics = {"n": len(y_true), "methods": {}}

        ax_rc = axes_rc[split_idx]
        for method_name, (scores, correct) in methods.items():
            coverage, risk = risk_coverage_curve(scores, correct)
            aurc_val = aurc(coverage, risk)
            sel_acc = {
                f"selective_accuracy_at_{int(c*100)}": selective_accuracy_at_coverage(scores, correct, c)
                for c in COVERAGE_LEVELS
            }
            split_metrics["methods"][method_name] = {"aurc": aurc_val, **sel_acc}
            ax_rc.plot(coverage, 1 - risk, label=method_name, linewidth=1.5)

            if method_name == "predictive_entropy":
                for c in COVERAGE_TABLE_LEVELS:
                    table_rows.append({
                        "split": split_name,
                        "method": method_name,
                        "coverage": c,
                        "selective_accuracy": selective_accuracy_at_coverage(scores, correct, c),
                    })

        full_coverage_acc = methods["softmax_response"][1].mean()
        ax_rc.axhline(full_coverage_acc, color="gray", linestyle=":", linewidth=1, label="full coverage (no abstention)")
        ax_rc.set_title(split_label)
        ax_rc.set_xlabel("Coverage")
        ax_rc.set_ylabel("Selective accuracy")
        ax_rc.set_ylim(0, 1)
        ax_rc.legend(fontsize=7, loc="lower left")

        cal_variants = calibration_variants(member_probs, mc_probs, y_true)
        split_metrics["calibration"] = {}
        ax_cal = axes_cal[split_idx]
        for variant_name, (confidence, correct) in cal_variants.items():
            ece_val = expected_calibration_error(confidence, correct)
            split_metrics["calibration"][variant_name] = {"ece": ece_val}
            if variant_name == "deep_ensemble":
                centers, accs, counts = reliability_bins(confidence, correct)
                ax_cal.plot(centers, accs, marker="o", label=f"deep_ensemble (ECE={ece_val:.3f})")

        ax_cal.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1, label="perfect calibration")
        ax_cal.set_title(split_label)
        ax_cal.set_xlabel("Confidence (max predicted prob.)")
        ax_cal.set_ylabel("Accuracy")
        ax_cal.set_xlim(0, 1)
        ax_cal.set_ylim(0, 1)
        ax_cal.legend(fontsize=8, loc="upper left")

        metrics[split_name] = split_metrics
        print(f"\n=== {split_label} (n={len(y_true)}) ===")
        for method_name, m in split_metrics["methods"].items():
            print(f"  {method_name}: AURC={m['aurc']:.4f}  "
                  f"acc@80={m['selective_accuracy_at_80']:.4f}  "
                  f"acc@90={m['selective_accuracy_at_90']:.4f}  "
                  f"acc@95={m['selective_accuracy_at_95']:.4f}")
        for variant_name, m in split_metrics["calibration"].items():
            print(f"  ECE ({variant_name}): {m['ece']:.4f}")

    fig_rc.suptitle("Risk-coverage curves (deep ensemble MLP)")
    fig_rc.tight_layout()
    fig_rc.savefig(os.path.join(FIG_DIR, "risk_coverage_curves.png"), dpi=150)

    fig_cal.suptitle("Calibration (reliability diagram, deep ensemble MLP)")
    fig_cal.tight_layout()
    fig_cal.savefig(os.path.join(FIG_DIR, "calibration.png"), dpi=150)

    # AURC comparison bar chart
    fig_bar, ax_bar = plt.subplots(figsize=(7, 4.5))
    method_names = list(next(iter(metrics.values()))["methods"].keys())
    x = np.arange(len(method_names))
    width = 0.35
    for i, (split_name, split_label) in enumerate(SPLITS.items()):
        vals = [metrics[split_name]["methods"][m]["aurc"] for m in method_names]
        ax_bar.bar(x + i * width, vals, width, label=split_label)
    ax_bar.set_xticks(x + width / 2)
    ax_bar.set_xticklabels(method_names, rotation=20, ha="right", fontsize=8)
    ax_bar.set_ylabel("AURC (lower is better)")
    ax_bar.set_title("AURC by abstention method")
    ax_bar.legend(fontsize=8)
    fig_bar.tight_layout()
    fig_bar.savefig(os.path.join(FIG_DIR, "aurc_comparison.png"), dpi=150)

    out_path = os.path.join(RESULTS_DIR, "metrics_evaluation.json")
    with open(out_path, "w") as f:
        json.dump(metrics, f, indent=2)

    import csv
    csv_path = os.path.join(RESULTS_DIR, "accuracy_vs_coverage.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["split", "method", "coverage", "selective_accuracy"])
        writer.writeheader()
        writer.writerows(table_rows)

    print(f"\nwrote {out_path}")
    print(f"wrote {csv_path}")
    print(f"wrote figures to {FIG_DIR}")

"""
Phase 2: baseline tumor-type classifiers, no uncertainty/abstention yet.

Trains on TCGA (exome) only, split into train/val/test. MSK-IMPACT (panel)
is never touched during training or model selection -- it's held out
entirely as the shifted evaluation set, matching a realistic deployment
scenario where you don't get labeled shifted data to tune on.

Usage: python3 train_baseline.py
Writes: results/metrics_baseline.json
"""
import json
import os
import sys

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
from data import RESULTS_DIR, SEED, load_data, panel_split, split_exome


def eval_split(name, model, X, y):
    pred = model.predict(X)
    acc = accuracy_score(y, pred)
    macro_f1 = f1_score(y, pred, average="macro")
    print(f"  {name}: n={len(y)}  accuracy={acc:.4f}  macro_f1={macro_f1:.4f}")
    return {"n": len(y), "accuracy": acc, "macro_f1": macro_f1}


if __name__ == "__main__":
    df, feature_cols = load_data()
    print(f"feature count: {len(feature_cols)}")

    (X_train, y_train), (X_val, y_val), (X_test, y_test) = split_exome(df, feature_cols)
    print(f"exome train/val/test sizes: {len(y_train)}/{len(y_val)}/{len(y_test)}")

    X_panel, y_panel = panel_split(df, feature_cols)
    print(f"panel (shifted, held out entirely from training): {len(y_panel)}")

    results = {}

    print("\n=== Gradient-boosted trees (HistGradientBoostingClassifier) ===")
    gbm = HistGradientBoostingClassifier(random_state=SEED)
    gbm.fit(X_train, y_train)
    results["gbm"] = {
        "train": eval_split("train", gbm, X_train, y_train),
        "val": eval_split("val", gbm, X_val, y_val),
        "test_exome_id": eval_split("test (exome, in-distribution)", gbm, X_test, y_test),
        "test_panel_shifted": eval_split("panel (shifted, out-of-distribution)", gbm, X_panel, y_panel),
    }

    print("\n=== 2-layer MLP ===")
    mlp = make_pipeline(
        StandardScaler(),
        MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=500, random_state=SEED, early_stopping=True),
    )
    mlp.fit(X_train, y_train)
    results["mlp"] = {
        "train": eval_split("train", mlp, X_train, y_train),
        "val": eval_split("val", mlp, X_val, y_val),
        "test_exome_id": eval_split("test (exome, in-distribution)", mlp, X_test, y_test),
        "test_panel_shifted": eval_split("panel (shifted, out-of-distribution)", mlp, X_panel, y_panel),
    }

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, "metrics_baseline.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n=== shift penalty (in-distribution accuracy - shifted accuracy) ===")
    for model_name in ["gbm", "mlp"]:
        id_acc = results[model_name]["test_exome_id"]["accuracy"]
        shift_acc = results[model_name]["test_panel_shifted"]["accuracy"]
        print(f"  {model_name}: {id_acc:.4f} -> {shift_acc:.4f}  (drop: {id_acc - shift_acc:.4f})")

    print(f"\nwrote {out_path}")

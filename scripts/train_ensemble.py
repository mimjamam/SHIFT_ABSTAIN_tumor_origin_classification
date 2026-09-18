"""
Phase 3: uncertainty via a deep ensemble (5 seeds) and MC-dropout, both
built on the same 2-layer MLP architecture as the Phase 2 baseline.

Deep ensemble: 5 independently-initialized copies of the MLP, trained on
the same data with different random init + minibatch shuffling (the
standard Lakshminarayanan et al. 2017 recipe). Ensemble disagreement
across members is one abstention signal.

MC-dropout: ensemble member 0's architecture, but dropout is left ACTIVE
at inference and sampled T times (Gal & Ghahramani 2016) to get a second,
cheaper (single-model) source of predictive uncertainty.

Usage: python3 train_ensemble.py
Writes: models/ (scaler + 5 MLP checkpoints), results/metrics_uncertainty.json,
        results/predictions/*.npz (per-split probability arrays for Phase 4)
"""
import json
import os
import sys

import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(__file__))
from data import CLASSES, MODELS_DIR, RESULTS_DIR, SEED, load_data, panel_split, split_exome
from mlp_model import TumorMLP

N_ENSEMBLE = 5
N_MC_SAMPLES = 30
MAX_EPOCHS = 200
PATIENCE = 20
BATCH_SIZE = 128
DEVICE = torch.device("cpu")

CLASS_TO_IDX = {c: i for i, c in enumerate(CLASSES)}


def encode_labels(y: np.ndarray) -> np.ndarray:
    return np.array([CLASS_TO_IDX[c] for c in y], dtype=np.int64)


def fit_scaler(X: np.ndarray):
    mean = X.mean(axis=0)
    scale = X.std(axis=0)
    scale[scale == 0] = 1.0
    return mean, scale


def apply_scaler(X, mean, scale):
    return (X - mean) / scale


def train_one_model(seed, X_train, y_train, X_val, y_val, n_features, n_classes):
    torch.manual_seed(seed)
    model = TumorMLP(n_features, n_classes).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    loss_fn = nn.CrossEntropyLoss()

    Xt = torch.tensor(X_train, dtype=torch.float32)
    yt = torch.tensor(y_train, dtype=torch.long)
    Xv = torch.tensor(X_val, dtype=torch.float32)
    yv = torch.tensor(y_val, dtype=torch.long)

    n = len(Xt)
    best_val_loss = float("inf")
    best_state = None
    epochs_since_improve = 0
    rng = np.random.default_rng(seed)

    for epoch in range(MAX_EPOCHS):
        model.train()
        perm = rng.permutation(n)
        for start in range(0, n, BATCH_SIZE):
            idx = perm[start:start + BATCH_SIZE]
            xb, yb = Xt[idx], yt[idx]
            opt.zero_grad()
            logits = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            opt.step()

        model.eval()
        with torch.no_grad():
            val_loss = loss_fn(model(Xv), yv).item()

        if val_loss < best_val_loss - 1e-4:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            epochs_since_improve = 0
        else:
            epochs_since_improve += 1
            if epochs_since_improve >= PATIENCE:
                break

    model.load_state_dict(best_state)
    return model, best_val_loss, epoch + 1


def predict_proba(model, X):
    model.eval()
    with torch.no_grad():
        logits = model(torch.tensor(X, dtype=torch.float32))
        return torch.softmax(logits, dim=-1).numpy()


def mc_dropout_proba(model, X, n_samples=N_MC_SAMPLES):
    """Keep dropout active (model.train()) but don't update weights."""
    model.train()
    Xt = torch.tensor(X, dtype=torch.float32)
    draws = []
    with torch.no_grad():
        for _ in range(n_samples):
            logits = model(Xt)
            draws.append(torch.softmax(logits, dim=-1).numpy())
    return np.stack(draws)  # (n_samples, n, n_classes)


def accuracy(probs, y_idx):
    pred = probs.argmax(axis=-1)
    return float((pred == y_idx).mean())


if __name__ == "__main__":
    df, feature_cols = load_data()
    (X_train, y_train_raw), (X_val, y_val_raw), (X_test, y_test_raw) = split_exome(df, feature_cols)
    X_panel, y_panel_raw = panel_split(df, feature_cols)

    mean, scale = fit_scaler(X_train)
    X_train_s = apply_scaler(X_train, mean, scale)
    X_val_s = apply_scaler(X_val, mean, scale)
    X_test_s = apply_scaler(X_test, mean, scale)
    X_panel_s = apply_scaler(X_panel, mean, scale)

    y_train = encode_labels(y_train_raw)
    y_val = encode_labels(y_val_raw)
    y_test = encode_labels(y_test_raw)
    y_panel = encode_labels(y_panel_raw)

    n_features = X_train_s.shape[1]
    n_classes = len(CLASSES)

    os.makedirs(MODELS_DIR, exist_ok=True)
    np.savez(os.path.join(MODELS_DIR, "scaler.npz"), mean=mean, scale=scale)
    with open(os.path.join(MODELS_DIR, "classes.json"), "w") as f:
        json.dump({"classes": CLASSES, "feature_cols": feature_cols}, f, indent=2)

    print(f"training {N_ENSEMBLE}-member deep ensemble (2-layer MLP, {n_features} features -> {n_classes} classes)")
    models = []
    for seed in range(N_ENSEMBLE):
        model, val_loss, n_epochs = train_one_model(seed, X_train_s, y_train, X_val_s, y_val, n_features, n_classes)
        torch.save(model.state_dict(), os.path.join(MODELS_DIR, f"mlp_seed{seed}.pt"))
        models.append(model)
        print(f"  seed {seed}: {n_epochs} epochs, best val loss {val_loss:.4f}")

    splits = {
        "val": (X_val_s, y_val),
        "test_exome_id": (X_test_s, y_test),
        "test_panel_shifted": (X_panel_s, y_panel),
    }

    os.makedirs(os.path.join(RESULTS_DIR, "predictions"), exist_ok=True)
    metrics = {"n_ensemble": N_ENSEMBLE, "n_mc_samples": N_MC_SAMPLES}

    for split_name, (X_s, y_idx) in splits.items():
        member_probs = np.stack([predict_proba(m, X_s) for m in models])  # (5, n, C)
        ensemble_mean_probs = member_probs.mean(axis=0)
        single_model_probs = member_probs[0]
        mc_probs = mc_dropout_proba(models[0], X_s)  # (T, n, C)
        mc_mean_probs = mc_probs.mean(axis=0)

        metrics[split_name] = {
            "n": len(y_idx),
            "single_model_accuracy": accuracy(single_model_probs, y_idx),
            "deep_ensemble_accuracy": accuracy(ensemble_mean_probs, y_idx),
            "mc_dropout_accuracy": accuracy(mc_mean_probs, y_idx),
        }
        print(f"{split_name}: single={metrics[split_name]['single_model_accuracy']:.4f}  "
              f"ensemble={metrics[split_name]['deep_ensemble_accuracy']:.4f}  "
              f"mc_dropout={metrics[split_name]['mc_dropout_accuracy']:.4f}")

        np.savez(
            os.path.join(RESULTS_DIR, "predictions", f"{split_name}.npz"),
            member_probs=member_probs,
            mc_probs=mc_probs,
            y_true=y_idx,
        )

    out_path = os.path.join(RESULTS_DIR, "metrics_uncertainty.json")
    with open(out_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nwrote {out_path}")
    print(f"wrote model artifacts to {MODELS_DIR}")

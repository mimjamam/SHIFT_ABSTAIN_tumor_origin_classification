"""
Shared data loading/splitting for every model script (baseline, ensemble,
MC-dropout). Keeping this in one place guarantees every model sees the
exact same train/val/test/panel split.
"""
import os
import pandas as pd
from sklearn.model_selection import train_test_split

ROOT = os.path.join(os.path.dirname(__file__), "..")
FEATURES_PATH = os.path.join(ROOT, "data", "processed", "features.parquet")
RESULTS_DIR = os.path.join(ROOT, "results")
MODELS_DIR = os.path.join(ROOT, "models")
SEED = 0

CLASSES = [
    "Bladder", "Breast", "Colorectal", "Endometrial", "Esophagogastric",
    "Glioma", "Head and Neck", "Melanoma", "NSCLC", "Ovarian", "Pancreatic",
    "Prostate", "Renal", "Thyroid",
]


def load_data():
    df = pd.read_parquet(FEATURES_PATH)
    feature_cols = [c for c in df.columns if c.startswith("sbs_") or c.startswith("driver_")] + ["tmb"]
    return df, feature_cols


def split_exome(df: pd.DataFrame, feature_cols: list[str]):
    """Same 70/15/15 stratified split every model script must use."""
    exome = df[df["cohort"] == "exome"]
    X = exome[feature_cols].to_numpy(dtype=float)
    y = exome["cancer_type"].to_numpy(dtype=str)

    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=0.30, stratify=y, random_state=SEED
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.50, stratify=y_temp, random_state=SEED
    )
    return (X_train, y_train), (X_val, y_val), (X_test, y_test)


def panel_split(df: pd.DataFrame, feature_cols: list[str]):
    panel = df[df["cohort"] == "panel"]
    X_panel = panel[feature_cols].to_numpy(dtype=float)
    y_panel = panel["cancer_type"].to_numpy(dtype=str)
    return X_panel, y_panel

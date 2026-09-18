"""
Core inference logic for the demo: load the trained deep ensemble once,
then turn an uploaded MAF file into predictions.

Reuses the exact same feature extraction (scripts/extract_features.py) and
normalization constants (scripts/build_dataset.py) as training, so a
prediction is computed the same way the model was evaluated in Phase 4 --
this is not a separate, simplified inference path.
"""
import json
import os
import sys

import numpy as np
import py2bit
import torch

ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPTS_DIR = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from build_dataset import DRIVER_GENES, EXOME_FOOTPRINT_MB, PANEL_FOOTPRINT_MB  # noqa: E402
from extract_features import extract_maf_features  # noqa: E402
from mlp_model import TumorMLP  # noqa: E402

MODELS_DIR = os.path.join(ROOT, "models")
REFERENCE_2BIT = os.path.join(ROOT, "data", "reference", "hg19.2bit")
N_ENSEMBLE = 5

# assay_type -> (footprint_mb, regime). "regime" selects which Phase-4 test
# set's threshold/calibration applies -- exome and panel are calibrated
# very differently (see docs/phase4_evaluation.md), so an IMPACT341 upload
# must be judged against the panel threshold, not the exome one.
ASSAY_TYPES = {
    "exome": {"footprint_mb": EXOME_FOOTPRINT_MB, "regime": "exome", "label": "Whole exome (38 Mb)"},
    "IMPACT341": {"footprint_mb": PANEL_FOOTPRINT_MB["IMPACT341"], "regime": "panel", "label": "MSK-IMPACT IMPACT341 (0.90 Mb)"},
    "IMPACT410": {"footprint_mb": PANEL_FOOTPRINT_MB["IMPACT410"], "regime": "panel", "label": "MSK-IMPACT IMPACT410 (1.02 Mb)"},
}


class Predictor:
    def __init__(self):
        with open(os.path.join(MODELS_DIR, "classes.json")) as f:
            meta = json.load(f)
        self.classes = meta["classes"]
        self.feature_cols = meta["feature_cols"]

        scaler = np.load(os.path.join(MODELS_DIR, "scaler.npz"))
        self.mean = scaler["mean"]
        self.scale = scaler["scale"]

        with open(os.path.join(MODELS_DIR, "thresholds.json")) as f:
            self.thresholds = json.load(f)

        n_features = len(self.feature_cols)
        n_classes = len(self.classes)
        self.models = []
        for seed in range(N_ENSEMBLE):
            model = TumorMLP(n_features, n_classes)
            state = torch.load(os.path.join(MODELS_DIR, f"mlp_seed{seed}.pt"), map_location="cpu")
            model.load_state_dict(state)
            model.eval()
            self.models.append(model)

        if not os.path.exists(REFERENCE_2BIT):
            raise FileNotFoundError(
                f"{REFERENCE_2BIT} not found -- run scripts/build_dataset.py's reference "
                "download step first (see Snakefile rule fetch_reference)."
            )
        self.tb = py2bit.open(REFERENCE_2BIT)

    def predict_maf(self, maf_path: str, assay_type: str, coverage: int = 90):
        if assay_type not in ASSAY_TYPES:
            raise ValueError(f"unknown assay_type {assay_type!r}, must be one of {list(ASSAY_TYPES)}")
        assay = ASSAY_TYPES[assay_type]
        footprint_mb = assay["footprint_mb"]
        regime = assay["regime"]
        threshold = self.thresholds[regime][str(coverage)]

        raw = extract_maf_features(maf_path, self.tb, DRIVER_GENES)
        if len(raw) == 0:
            return []

        sbs_cols = [c for c in raw.columns if c.startswith("sbs_")]
        raw[sbs_cols] = raw[sbs_cols].div(footprint_mb)
        raw["tmb"] = raw["n_nonsyn"] / footprint_mb

        X = raw[self.feature_cols].to_numpy(dtype=float)
        X_scaled = (X - self.mean) / self.scale
        Xt = torch.tensor(X_scaled, dtype=torch.float32)

        with torch.no_grad():
            member_probs = np.stack([
                torch.softmax(model(Xt), dim=-1).numpy() for model in self.models
            ])
        mean_probs = member_probs.mean(axis=0)

        results = []
        for i, sample_id in enumerate(raw.index):
            probs = mean_probs[i]
            order = np.argsort(-probs)
            top3 = [{"cancer_type": self.classes[j], "probability": float(probs[j])} for j in order[:3]]
            confidence = float(probs[order[0]])
            decision = "PREDICT" if confidence >= threshold else "ABSTAIN"
            results.append({
                "sample_id": sample_id,
                "assay_type": assay_type,
                "regime": regime,
                "footprint_mb": footprint_mb,
                "n_mutations_used": int(raw.loc[sample_id, "n_snv_total"]),
                "tmb": float(raw.loc[sample_id, "tmb"]),
                "predicted_type": top3[0]["cancer_type"],
                "confidence": confidence,
                "decision": decision,
                "coverage_operating_point": coverage,
                "threshold_used": threshold,
                "top3": top3,
            })
        return results


_predictor = None


def get_predictor() -> Predictor:
    global _predictor
    if _predictor is None:
        _predictor = Predictor()
    return _predictor

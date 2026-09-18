"""
Precompute the softmax-response threshold for each (regime, coverage level)
pair, from the exact same predictions already evaluated in Phase 4 --
softmax_response was the best-ranking abstention score there (see
docs/phase4_evaluation.md), so it's what the demo API uses for its
PREDICT/ABSTAIN decision.

"regime" is "exome" (uses the in-distribution test set's score
distribution) or "panel" (uses the shifted test set's) -- calibration
differs sharply between them (see the ECE numbers in Phase 4), so a
panel-assay upload is judged against the panel threshold, not the exome
one.

Usage: python3 compute_thresholds.py
Writes: models/thresholds.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from data import MODELS_DIR, RESULTS_DIR
from uncertainty import softmax_response

PRED_DIR = os.path.join(RESULTS_DIR, "predictions")
COVERAGE_LEVELS = [0.80, 0.90, 0.95]
REGIME_TO_SPLIT = {"exome": "test_exome_id", "panel": "test_panel_shifted"}

if __name__ == "__main__":
    thresholds = {}
    for regime, split_name in REGIME_TO_SPLIT.items():
        d = np.load(os.path.join(PRED_DIR, f"{split_name}.npz"))
        mean_probs = d["member_probs"].mean(axis=0)
        scores = softmax_response(mean_probs)
        thresholds[regime] = {}
        for c in COVERAGE_LEVELS:
            # keep the top-c fraction most confident => threshold is the
            # (1-c) quantile of the score distribution.
            t = float(np.quantile(scores, 1 - c))
            thresholds[regime][str(int(c * 100))] = t
            print(f"{regime} @ {int(c*100)}% coverage: threshold={t:.4f}")

    out_path = os.path.join(MODELS_DIR, "thresholds.json")
    with open(out_path, "w") as f:
        json.dump(thresholds, f, indent=2)
    print(f"\nwrote {out_path}")

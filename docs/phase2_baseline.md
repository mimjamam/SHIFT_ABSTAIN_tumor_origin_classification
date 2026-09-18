# Phase 2 — baseline models

## Setup

`scripts/train_baseline.py`. Trains on TCGA (exome) only:
5,750 train / 1,232 val / 1,233 test, stratified by cancer_type
(70/15/15 split, seed 0). MSK-IMPACT (7,803 panel samples) is **never**
touched during training, validation, or model selection — it's held out
entirely as the shifted evaluation set, matching a realistic deployment
scenario where labeled shifted data usually isn't available to tune on.

127 features: 96 `sbs_*` per-Mb rates, `tmb`, 30 `driver_*` flags. No
uncertainty/abstention yet — that's Phase 3. This phase just answers "does
the shift actually hurt, and by how much."

Two models, both intentionally small per the project scope:
- **Gradient-boosted trees**: `sklearn.HistGradientBoostingClassifier`,
  default hyperparameters.
- **2-layer MLP**: `sklearn.MLPClassifier`, hidden layers (64, 32),
  standardized inputs, early stopping.

## Results

| Model | Train acc | Val acc | Test (exome, ID) acc | Panel (shifted) acc | Drop |
|---|---|---|---|---|---|
| GBM | 1.0000 | 0.7094 | 0.7380 | 0.3176 | **0.4205** |
| MLP | 0.7351 | 0.6648 | 0.6504 | 0.3077 | **0.3427** |

(macro-F1 tells the same story: GBM 0.7267 -> 0.2910, MLP 0.6431 -> 0.2481 —
so this isn't just a majority-class artifact from class imbalance.)

Full numbers in `results/metrics_baseline.json`.

## Reading

- Both models generalize reasonably within-distribution (65-74% accuracy on
  a 14-class problem where chance is ~7%, and where several classes are
  genuinely hard to separate by mutation profile alone, e.g. Esophagogastric
  vs. Colorectal, or Renal subtypes).
- Both models lose roughly **half** their accuracy under the panel shift —
  a real, substantial, non-trivial degradation. This is the effect the
  whole project exists to characterize and then mitigate with abstention.
- GBM trains to 100% on the training set (expected for an unregularized
  tree ensemble on 127 features with ~5.7k samples) but still generalizes
  better than the MLP both in- and out-of-distribution, so it's not pure
  overfitting noise.
- Neither model collapses to chance-level (~7%) under shift — they retain
  partial signal (~31-32%), which is exactly the regime where a calibrated
  abstention policy has room to help: confidently-wrong panel predictions
  can in principle be separated from confidently-right ones and deferred.

## What's next (Phase 3)

Deep ensemble (5 seeds) + MC-dropout uncertainty, and abstention baselines
(softmax response, predictive entropy, ensemble disagreement) — evaluated
via risk-coverage curves in Phase 4.

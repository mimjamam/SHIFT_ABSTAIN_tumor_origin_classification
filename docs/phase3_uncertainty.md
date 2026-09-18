# Phase 3 — uncertainty and abstention signals

## What was built

- `scripts/mlp_model.py` — the 2-layer MLP (same shape as Phase 2's
  `sklearn` MLP: 64 -> 32 hidden units), now in PyTorch with dropout
  (p=0.3) so it can serve double duty as both an ensemble member and an
  MC-dropout model.
- `scripts/train_ensemble.py` — trains a **5-seed deep ensemble**
  (Lakshminarayanan et al. 2017: same architecture and data, different
  random init + minibatch order) and runs **MC-dropout** (Gal &
  Ghahramani 2016: dropout left active at inference, 30 stochastic
  forward passes) on ensemble member 0. Saves model weights to `models/`
  (committed — small, ~46KB per member, and needed later for the API
  without retraining) and per-split probability arrays to
  `results/predictions/` (gitignored — regenerable from fixed seeds).
- `scripts/uncertainty.py` — three abstention score functions, all
  oriented so **higher = more confident**:
  - `softmax_response`: max predicted probability.
  - `predictive_entropy`: negative entropy of the predictive distribution.
  - `ensemble_disagreement`: negative mutual information between the
    prediction and ensemble membership (the "BALD" score) — isolates
    genuine disagreement between models from a single model just being
    unsure between two classes.
- `scripts/data.py` — split logic factored out of `train_baseline.py` so
  every model (GBM, sklearn MLP, deep ensemble, MC-dropout) trains and
  evaluates on the *exact same* train/val/test/panel split.

## Results

| Split | Single model | Deep ensemble (mean) | MC-dropout (mean) |
|---|---|---|---|
| val | 0.6826 | 0.6778 | 0.6834 |
| test (exome, ID) | 0.6853 | 0.6869 | 0.6837 |
| panel (shifted) | 0.3095 | **0.3351** | 0.3099 |

The deep ensemble gives a small but real accuracy bump specifically on the
*shifted* set (31.0% -> 33.5%) — consistent with ensembling helping most
exactly where a single model is least reliable. MC-dropout tracks the
single model closely everywhere, as expected (it's smoothing the same
model's own uncertainty, not combining independently-trained models).

## Sanity check: do the abstention scores actually track correctness?

On the shifted (panel) test set, correlation between each score and
whether the ensemble-mean prediction was correct:

| Score | Correlation with correctness |
|---|---|
| softmax_response | 0.260 |
| predictive_entropy | 0.234 |
| ensemble_disagreement | 0.156 |
| MC-dropout disagreement | 0.161 |

All four are positive, which is the precondition for abstention to work at
all — thresholding on any of them and declining the lowest-confidence
predictions should raise accuracy on what's left. Softmax response and
entropy (which are closely related) carry the strongest signal here;
ensemble/MC-dropout disagreement is weaker alone but measures something
different (epistemic uncertainty specifically) and may combine usefully
with the others. Phase 4 turns this into actual risk-coverage curves and
AURC rather than a single correlation number.

## What's next (Phase 4 — stop for review after)

Risk-coverage curves, AURC, selective accuracy at 80/90/95% coverage, and
ECE, computed separately for in-distribution and shifted test sets, for
each abstention method and each model.

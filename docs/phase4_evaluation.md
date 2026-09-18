# Phase 4 — evaluation

## Setup

`scripts/evaluate.py`, `scripts/selective_metrics.py`. All metrics are
computed on the deep-ensemble MLP's predictions from Phase 3, using four
abstention scores (`softmax_response`, `predictive_entropy`,
`ensemble_disagreement`, `mc_dropout_disagreement`), separately for the
in-distribution (exome) test set (n=1,233) and the shifted (panel) test
set (n=7,803) — the panel set was never touched during training or
threshold selection. Full numbers: `results/metrics_evaluation.json`,
`results/accuracy_vs_coverage.csv`. Figures: `results/figures/`.

## Headline result: abstention helps, but does not close the shift gap

| Coverage | In-distribution accuracy | Shifted accuracy |
|---|---|---|
| 100% (no abstention) | 68.7% | 33.5% |
| 90% | 72.5% | 35.0% |
| 80% | 77.6% | 37.1% |
| 50% | 91.1% | 43.3% |
| 10% | 99.2% | 56.4% |

(`predictive_entropy` scoring; `results/accuracy_vs_coverage.csv` has the
full 10%-step table for both methods and both splits.)

Abstaining on the least-confident cases works well in-distribution — even
giving up half the coverage gets accuracy above 90%, and the most
confident 10% of predictions are essentially always right (99.2%).

Under shift, abstention still helps — accuracy does climb as coverage
drops — but nowhere near enough to close the gap. Even discarding 90% of
panel samples and keeping only the 10% the model is most confident about,
accuracy reaches just 56.4%: still well below the *unfiltered* exome
accuracy of 68.7%. **The honest reading is that a calibrated abstention
policy recovers only a modest fraction of the shift-induced loss, not
most of it** — this is a negative result relative to the project's
original hypothesis, reported as-is per the project's own ground rules.

## Why: the model's confidence itself is unreliable under shift

Expected Calibration Error (10-bin, max-predicted-probability confidence):

| Variant | ID ECE | Shifted ECE |
|---|---|---|
| Single model | 0.033 | 0.404 |
| Deep ensemble | 0.035 | 0.299 |
| MC-dropout | 0.080 | 0.236 |

In-distribution, the model is close to well-calibrated (ECE ~0.03-0.08).
Under shift, it becomes badly **overconfident** — the reliability diagram
(`results/figures/calibration.png`) shows the shifted curve sitting well
below the diagonal across the whole confidence range (e.g. predictions
made at ~0.95 confidence are only right about 51% of the time). This is
the direct mechanism behind the limited abstention recovery above: a
confidence score that has decoupled from actual correctness can't rank
samples for abstention as well as one that hasn't. Interestingly,
MC-dropout has the *worst* ID calibration but the *best* shifted
calibration of the three — its stochastic averaging seems to specifically
soften shifted-set overconfidence, even though (see below) it's also the
worst-ranking method for selective accuracy. Calibration quality and
selective-ranking quality are related but not the same thing.

## AURC by method (lower is better; `results/figures/aurc_comparison.png`)

| Method | ID AURC | Shifted AURC |
|---|---|---|
| softmax_response | 0.114 | 0.547 |
| predictive_entropy | 0.117 | 0.554 |
| ensemble_disagreement | 0.163 | 0.560 |
| mc_dropout_disagreement | 0.198 | 0.599 |

`softmax_response` and `predictive_entropy` (closely related, and closely
tied in practice) are the best-ranking abstention scores on both splits.
The two disagreement-based scores, which isolate epistemic
(model-uncertainty) signal specifically, rank worse here — on this
problem, "how confident is the model" is a better abstention signal than
"how much do independent copies of the model disagree." This is a
legitimate, useful negative result about which of the three specified
abstention baselines actually earns its complexity: the simplest one
(softmax response) wins.

## Selective accuracy at fixed coverage (`results/metrics_evaluation.json`)

| Method | ID acc@80/90/95 | Shifted acc@80/90/95 |
|---|---|---|
| softmax_response | 0.774 / 0.729 / 0.712 | 0.372 / 0.352 / 0.345 |
| predictive_entropy | 0.776 / 0.725 / 0.712 | 0.371 / 0.350 / 0.343 |
| ensemble_disagreement | 0.729 / 0.706 / 0.695 | 0.357 / 0.349 / 0.342 |
| mc_dropout_disagreement | 0.714 / 0.692 / 0.692 | 0.336 / 0.324 / 0.317 |

## Summary for the README

- Shift penalty at full coverage: **68.7% -> 33.5%**, a 35.2-point drop.
- Best abstention method (softmax response), shifted set: **33.5% -> 37.1%**
  at 80% coverage — a real but modest recovery (+3.6 points, keeping 80%
  of samples), nowhere close to erasing the shift penalty.
- Root cause, evidenced directly: the model becomes markedly overconfident
  under shift (ECE 0.033 -> 0.404), so its own confidence is a
  meaningfully weaker signal for "should I trust this" than it is
  in-distribution.

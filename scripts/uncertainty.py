"""
Abstention score functions. Each takes predictive probabilities and returns
a per-sample score where HIGHER = MORE confident (so thresholding "abstain
if score < t" behaves the same way for every method). Risk-coverage curves
sweep t in Phase 4.
"""
import numpy as np

EPS = 1e-12


def softmax_response(probs: np.ndarray) -> np.ndarray:
    """probs: (n_samples, n_classes) from a single model (or an ensemble
    mean). Confidence = max predicted probability."""
    return probs.max(axis=-1)


def predictive_entropy(probs: np.ndarray) -> np.ndarray:
    """probs: (n_samples, n_classes). Confidence = negative entropy (so
    higher is more confident, consistent with the other scores here)."""
    p = np.clip(probs, EPS, 1.0)
    entropy = -(p * np.log(p)).sum(axis=-1)
    return -entropy


def ensemble_disagreement(member_probs: np.ndarray) -> np.ndarray:
    """member_probs: (n_members, n_samples, n_classes) -- softmax output of
    each ensemble member (or each MC-dropout draw). Uses mutual information
    between the prediction and the ensemble/dropout mask (the "BALD" score
    from Houlston et al. / Gal & Ghahramani): the entropy of the mean
    prediction minus the mean of each member's entropy. This isolates
    epistemic (disagreement) uncertainty from aleatoric uncertainty -- two
    members that individually assign 50/50 to two classes contribute little
    here if they agree on *which* two classes, but a lot if they disagree.
    Confidence = negative mutual information (higher = more confident,
    i.e. the ensemble agrees)."""
    mean_probs = member_probs.mean(axis=0)
    mean_p = np.clip(mean_probs, EPS, 1.0)
    entropy_of_mean = -(mean_p * np.log(mean_p)).sum(axis=-1)

    member_p = np.clip(member_probs, EPS, 1.0)
    mean_of_entropy = (-(member_p * np.log(member_p)).sum(axis=-1)).mean(axis=0)

    mutual_info = entropy_of_mean - mean_of_entropy
    return -mutual_info


SCORE_FUNCS = {
    "softmax_response": softmax_response,
    "predictive_entropy": predictive_entropy,
}

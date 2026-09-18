"""2-layer MLP with dropout, used both for the deep ensemble (dropout as a
plain regularizer, disabled at inference) and for MC-dropout (same
architecture, dropout left active at inference and sampled multiple times)."""
import torch.nn as nn


class TumorMLP(nn.Module):
    def __init__(self, n_features: int, n_classes: int, hidden=(64, 32), dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, hidden[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden[0], hidden[1]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden[1], n_classes),
        )

    def forward(self, x):
        return self.net(x)  # logits

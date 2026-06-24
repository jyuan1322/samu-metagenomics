import torch
import torch.nn as nn
from skorch import NeuralNetClassifier

class FC_NN(nn.Module):
    def __init__(self, L, h=[18, 12], p=0.2):
        super().__init__()
        hnew = []
        for i, hh in enumerate(h):
            if hh == 0:
                if i == 0:
                    hnew.append(int((3*L)**0.5 + 2*(L/3)**0.5))
                else:
                    hnew.append(int((L/3)**0.5))
            else:
                hnew.append(hh)
        h = hnew
        # set final layer to 2 outputs instead of 1
        self.hidden_sizes = [L] + h + [2]
        layers = []
        for k in range(len(self.hidden_sizes) - 1):
            layers.append(nn.Linear(self.hidden_sizes[k], self.hidden_sizes[k + 1]))
            if k <= len(h) - 2: # only hidden layers
                layers.append(nn.GELU())
                layers.append(nn.Dropout(p))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        x = self.net(x)
        return torch.log_softmax(x, dim=-1)  # skorch expects log-probs


def make_ffnn_classifier(input_dim, h=[18, 12], p=0.2, lr=0.001, weight_decay=0.0, max_epochs=50):
    return NeuralNetClassifier(
        module=FC_NN,
        module__L=input_dim,
        module__h=h,
        module__p=p,
        optimizer=torch.optim.Adam,
        lr=lr,
        optimizer__weight_decay=weight_decay,  # <── L2 regularization
        max_epochs=max_epochs,
        batch_size=32,
        iterator_train__shuffle=True,
        train_split=None,  # disable internal validation; outer CV handles it
        verbose=0,
        device='cpu',       # or 'cuda' when GPU ready
        criterion=torch.nn.CrossEntropyLoss
    )

# If running neural net with skorch and getting nan training loss,
# switch criterion to torch.nn.CrossEntropyLoss
# https://stackoverflow.com/a/70258008
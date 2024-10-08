import torch
import torch.nn as nn
import torch.nn.functional as F

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device('cpu')
DEVICE = device

class BiologicalModule(nn.Module):
    def __init__(self, input_dim, pathway_dim, output_dim, nonzero_ind=None, activation='tanh', dropout = 0.4, use_bias=True):
        super(BiologicalModule, self).__init__()

        self.input_dim = input_dim
        self.pathway_dim = pathway_dim
        self.output_dim = output_dim
        self.activation_fn = self.get_activation_fn(activation)

        # Initialize the linear transformation (input_dim -> output_dim)
        self.linear1 = nn.Linear(input_dim, pathway_dim)
        self.linear2 = nn.Linear(pathway_dim, output_dim)
        self.dropout = nn.Dropout(dropout)
        self.sigmoid = nn.Sigmoid()

        # If nonzero indices are provided, we apply masking
        # import pdb; pdb.set_trace()
        self.nonzero_ind = nonzero_ind

    def forward(self, inputs):
        # If nonzero_ind is provided, mask the weight matrix
        if self.nonzero_ind is not None:
            weight = self.linear1.weight.data.to(device)
            mask = torch.zeros_like(weight).to(device)
            mask[self.nonzero_ind[1], self.nonzero_ind[0]] = 1
            self.linear1.weight.data = weight * mask

        output = self.linear1(inputs)
        output = self.dropout(output)
        output = self.sigmoid(output)
        output = self.linear2(output)

        if self.activation_fn is not None:
            output = self.activation_fn(output)

        return output

    def get_activation_fn(self, activation):
        if activation == 'tanh':
            return torch.tanh
        elif activation == 'relu':
            return F.relu
        elif activation == 'sigmoid':
            return torch.sigmoid
        else:
            return None

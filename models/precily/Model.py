import torch
import random
import numpy as np
from optuna.samplers import TPESampler
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.utils import to_dense_batch
from BM import *

from Model_MHA import MultiHeadAttentionLayer
from Model_DAE import ldim

features_dim_gene = ldim
features_dim_bionic = 512

class AttentionLayer_DIPK(nn.Module):
    def __init__(self, heads):
        super(AttentionLayer_DIPK, self).__init__()
        self.fc_layer = nn.Linear(features_dim_gene, 768)
        self.attention = MultiHeadAttentionLayer(hid_dim=768, n_heads=1, dropout=0.3, device=DEVICE)

    def forward(self, x, g, gene, bionic):
        gene = F.relu(self.fc_layer(gene))
        bionic = F.relu(self.fc_layer(bionic))
        x = to_dense_batch(x, g.batch)
        query_0 = torch.unsqueeze(gene, 1)
        query_1 = torch.unsqueeze(bionic, 1)
        key = x[0]
        value = x[0]
        mask = torch.unsqueeze(torch.unsqueeze(x[1], 1), 1)
        x_att = self.attention(query_0, key, value, mask)
        x = torch.squeeze(x_att[0])
        x_att = self.attention(query_1, key, value, mask)
        x += torch.squeeze(x_att[0])
        return x


class DenseLayers_DIPK(nn.Module):
    def __init__(self, heads, fc_layer_num, fc_layer_dim, dropout_rate):
        super(DenseLayers_DIPK, self).__init__()
        self.fc_layer_num = fc_layer_num
        self.fc_layer = nn.Linear(512, 512)
        self.fc_input = nn.Linear(768 + 512, 768 + 512)
        self.fc_layers = torch.nn.Sequential(
            nn.Linear(768 + 512, 512),
            nn.Linear(512, fc_layer_dim[0]),
            nn.Linear(fc_layer_dim[0], fc_layer_dim[1]),
            nn.Linear(fc_layer_dim[1], fc_layer_dim[2]),
            nn.Linear(fc_layer_dim[2], fc_layer_dim[3]),
            nn.Linear(fc_layer_dim[3], fc_layer_dim[4]),
            nn.Linear(fc_layer_dim[4], fc_layer_dim[5])
        )
        self.dropout_layers = torch.nn.ModuleList(
            [nn.Dropout(p=dropout_rate) for _ in range(fc_layer_num)]
        )
        self.fc_output = nn.Linear(fc_layer_dim[fc_layer_num - 2], 1)

    def forward(self, x, gene, bionic):
        gene = F.relu(self.fc_layer(gene))
        bionic = F.relu(self.fc_layer(bionic))
        f = torch.cat((x, gene + bionic), 1)
        f = F.relu(self.fc_input(f))
        for layer_index in range(self.fc_layer_num):
            f = F.relu(self.fc_layers[layer_index](f))
            f = self.dropout_layers[layer_index](f)
        f = self.fc_output(f)
        return f


class Predictor_DIPK(nn.Module):
    def __init__(self, embedding_dim, heads, fc_layer_num, fc_layer_dim, dropout_rate):
        super(Predictor_DIPK, self).__init__()
        # self.graph_encoder = GraphEncoder(embedding_dim, heads)
        self.attention_layer = AttentionLayer_DIPK(heads)
        self.dense_layers = DenseLayers_DIPK(heads, fc_layer_num, fc_layer_dim, dropout_rate)

    def forward(self, x, g, gene, bionic):
        # x = self.graph_encoder(g)
        x = self.attention_layer(x, g, gene, bionic)
        f = self.dense_layers(x, gene, bionic)
        return f

class Net_DIPK(nn.Module):
    def __init__(self, encoder, model):
        super(Net_DIPK, self).__init__()
        self.encoder = encoder
        self.model = model

    def forward(self, features):
        EXPR = features['EXPR']
        BNF = features['BNF']
        GRAPH = features['GRAPH']

        EXPR = self.encoder(EXPR)
        prediction = self.model(GRAPH.x, GRAPH, EXPR, BNF)
        return torch.squeeze(prediction)

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    sampler = TPESampler(seed=seed)
    return sampler

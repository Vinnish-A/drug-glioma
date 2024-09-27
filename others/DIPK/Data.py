# import os; os.chdir('F:/0Local/project/drug-glioma/code')
# import sys; sys.path.append(os.getcwd())

from abc import ABC
import pandas as pd
import joblib
import torch
from torch_geometric.data import Batch
from torch_geometric.data import Data
from torch_geometric.data import Dataset
from BM import *

gene_add_num = 256
debug = 0

def Zscore(vector):
    return (vector - torch.mean(vector)) / (torch.std(vector))



class MyDataSet(Dataset, ABC):
    def __init__(self, graphs):
        self._graphs = graphs

    def __getitem__(self, idx):
        graph = self._graphs[idx]
        return graph

    def __len__(self):
        return len(self._graphs)

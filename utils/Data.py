# import os; os.chdir('F:/0Local/project/drug-glioma/code')
# import sys; sys.path.append(os.getcwd())

from abc import ABC
import pandas as pd
import joblib
import torch
from torch_geometric.data import Batch
from torch_geometric.data import Data
from torch_geometric.data import Dataset
from BM import DEVICE

gene_add_num = 256
debug = 0

df_kegg = joblib.load('Dataset/df_kegg.pkl').to_numpy()
GRAPH_dict = joblib.load('Dataset/GRAPH_dict.pkl')
RNA_dict = joblib.load('Dataset/RNAseq_Dict.pkl')
MolGNet_dict = joblib.load('Dataset/MolGNet_dict.pkl')
BIONIC_dict = joblib.load('Dataset/BIONIC_dict.pkl')
CNV_dict = joblib.load('Dataset/CNV_Dict.pkl')

GEF = []
Cells = list(RNA_dict.keys())
for each in Cells:
    GEF.append(RNA_dict[each])
GEF = torch.tensor(GEF, dtype=torch.float32)
values, indices = torch.sort(GEF, descending=True)
BNF_dict = dict()
for i in range(len(Cells)):
    k = 0
    feature = torch.tensor([0.0] * 512, dtype=torch.float32)
    for j in range(gene_add_num):
        if int(indices[i, j]) in BIONIC_dict:
            k += 1
            feature += torch.tensor(BIONIC_dict[int(indices[i, j])], dtype=torch.float32)
    feature = (feature / k).tolist()
    BNF_dict[Cells[i]] = feature

def Zscore(vector):
    return (vector - torch.mean(vector)) / (torch.std(vector))


def GetData(dataset):
    Cell_ = dataset.iloc[:, 0].tolist()
    Drug_ = dataset.iloc[:, 1].tolist()
    IC50_ = dataset.iloc[:, 2].tolist()

    Cell = []
    Drug = []
    IC50 = []
    for ii in range(len(Cell_)):
        if Cell_[ii] in Cells:
            Cell.append(Cell_[ii])
            Drug.append(Drug_[ii])
            IC50.append(IC50_[ii])

    Graph = []
    for ii in range(len(Cell)):
        graph = GRAPH_dict[Drug[ii]]
        x, edge_index, edge_attr = MolGNet_dict[Drug[ii]], graph.edge_index, graph.edge_attr
        graph = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, cell=Cell[ii], drug=Drug[ii],
                     GEF=Zscore(torch.tensor((RNA_dict[Cell[ii]]), dtype=torch.float32)),
                     CNV=torch.tensor((CNV_dict[Cell[ii]]), dtype=torch.float32),
                     ic50=torch.tensor([IC50[ii]], dtype=torch.float32))
        Graph.append(graph)

    print('Prepared! ')

    return Graph


def GetData_DIPK(dataset):
    Cell_ = dataset.iloc[:, 0].tolist()
    Drug_ = dataset.iloc[:, 1].tolist()
    IC50_ = dataset.iloc[:, 2].tolist()

    Cell = []
    Drug = []
    IC50 = []
    for ii in range(len(Cell_)):
        if Cell_[ii] in Cells:
            Cell.append(Cell_[ii])
            Drug.append(Drug_[ii])
            IC50.append(IC50_[ii])

    Graph = []
    for ii in range(len(Cell)):
        graph = GRAPH_dict[Drug[ii]]
        x, edge_index, edge_attr = MolGNet_dict[Drug[ii]], graph.edge_index, graph.edge_attr
        graph = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, cell=Cell[ii], drug=Drug[ii],
                     GEF=Zscore(torch.tensor((RNA_dict[Cell[ii]]), dtype=torch.float32)),
                     BNF=torch.tensor(BNF_dict[Cell[ii]], dtype=torch.float32),
                     ic50=torch.tensor([IC50[ii]], dtype=torch.float32))
        Graph.append(graph)

    return Graph

class CollateFn:
    def __init__(self, test=False, follow_batch=None, exclude_keys=None):
        self.test = test
        self.follow_batch = follow_batch
        self.exclude_keys = exclude_keys

    def __call__(self, batch):
        pyg_list = [Data(x=g.x, edge_index=g.edge_index, edge_attr=g.edge_attr, ic50=g.ic50) for g in batch]
        pyg_batch = Batch.from_data_list(pyg_list, self.follow_batch, self.exclude_keys).to(DEVICE)
        GeneFt = torch.stack([g.GEF for g in batch]).to(DEVICE)
        CNVFt = torch.stack([g.CNV for g in batch]).to(DEVICE)
        BionicFt = GeneFt

        features = dict(zip(['GRAPH', 'EXPR', 'PATHWAY', 'CNV'], [pyg_batch, GeneFt, BionicFt, CNVFt]))
        labels = pyg_batch.ic50.to(DEVICE)

        if self.test:
            samples = [g.cell + '\t' + g.drug for g in batch]
            return features, labels, samples

        return features, labels


class CollateFn_DIPK:
    def __init__(self, test = False, follow_batch=None, exclude_keys=None):
        self.test = test
        self.follow_batch = follow_batch
        self.exclude_keys = exclude_keys

    def __call__(self, batch):
        pyg_list = [Data(x=g.x, edge_index=g.edge_index, edge_attr=g.edge_attr, ic50=g.ic50) for g in batch]
        pyg_batch = Batch.from_data_list(pyg_list, self.follow_batch, self.exclude_keys).to(DEVICE)
        GeneFt = torch.stack([g.GEF for g in batch]).to(DEVICE)
        BionicFt = torch.stack([g.BNF for g in batch]).to(DEVICE)
        features = dict(zip(['GRAPH', 'EXPR', 'BNF'], [pyg_batch, GeneFt, BionicFt]))
        labels = pyg_batch.ic50

        if self.test:
            samples = [g.cell + '\t' + g.drug for g in batch]
            return features, labels, samples

        return features, labels


class MyDataSet(Dataset, ABC):
    def __init__(self, graphs):
        self._graphs = graphs

    def __getitem__(self, idx):
        graph = self._graphs[idx]
        return graph

    def __len__(self):
        return len(self._graphs)

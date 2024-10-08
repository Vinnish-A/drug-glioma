import os
import sys

import joblib
import torch

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./models')
sys.path.append('./models/Mine_eval')

from torch.utils.data import DataLoader
from torchmetrics import MeanSquaredError

from utils.Data import *
from utils.Train import *

from models.model_Mine_eval import net_mine_eval
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem

def get_color(score):
    return (1.0, 1.0 - score, 1.0 - score)  # (R, G, B)

def vis_mol(smiles, node_score, path):
    mol = Chem.MolFromSmiles(smiles)

    min_score, max_score = min(node_score.values()), max(node_score.values())
    normalized_scores = {k: (v - min_score) / (max_score - min_score) for k, v in node_score.items()}
    highlightAtomColors = {idx: get_color(score) for idx, score in normalized_scores.items()}

    d = rdMolDraw2D.MolDraw2DSVG(500, 500)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=list(node_score.keys()),
                                       highlightAtomColors=highlightAtomColors)
    d.FinishDrawing()

    with open(path, "w") as f:
        f.write(d.GetDrawingText())


SMILES_dict = joblib.load('Dataset/SMILES_dict.pkl')

part1 = pd.read_csv('Dataset/Train.csv', sep='\t')
part2 = pd.read_csv('Dataset/Test.csv', sep='\t')
tcga = pd.read_csv('Dataset/TCGA.csv', sep='\t')
sample_all = pd.concat([part1, part2], axis=0)
sample_att = sample_all.loc[sample_all.Drug.isin(['Temozolomide'])]
dl_att  = DataLoader(MyDataSet(GetData(sample_att)), batch_size=sample_att.shape[0], shuffle=True, collate_fn=CollateFn(True))

net_mine_eval.load_state_dict(torch.load('checkpoint/mine_trained.pt'))
with torch.no_grad():
    for it, (features, labels, samples) in enumerate(dl_att):
        prediction, lst_att = net_mine_eval(features)

score_att = dict(zip(range(14), (torch.sum(lst_att[0][1].squeeze([1,2]), dim = 0)/14).tolist()))

vis_mol(SMILES_dict['Temozolomide'], score_att, "scratch/output.svg")


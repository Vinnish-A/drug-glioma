import os
import sys

import joblib
import pandas as pd
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
from PIL.ImageColor import getrgb

SMILES_dict = joblib.load('Dataset/SMILES_dict.pkl')
net_mine_eval.load_state_dict(torch.load('checkpoint/mine_trained.pt'))

def get_color(score, HEX = '#94697a'):

    color_start = np.array([255, 255, 255])
    color_end = np.array(list(getrgb(HEX)))
    color = ((1 - score) * color_start + score * color_end)/256

    return tuple(color)

def vis_mol(smiles, node_score, path):
    mol = Chem.MolFromSmiles(smiles)

    min_score, max_score = min(node_score.values()), max(node_score.values())
    normalized_scores = {k: (v - min_score) / (max_score - min_score) for k, v in node_score.items()}
    highlightAtomColors = {idx: get_color(score) for idx, score in normalized_scores.items()}

    d = rdMolDraw2D.MolDraw2DSVG(400, 400)
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=list(node_score.keys()),
                                       highlightAtomColors=highlightAtomColors)
    d.drawOptions().clearBackground = False
    d.FinishDrawing()

    with open(path, "w") as f:
        f.write(d.GetDrawingText())

def save_score(score, name):
    pd.DataFrame({'key': score.keys(), 'value': score.values()}).to_csv(f'result/att_{name}.csv', index=None)

def flow(path, name):
    sample_att = pd.read_csv(f'Dataset/sliced/{path}')
    dl_att = DataLoader(MyDataSet(GetData(sample_att)), batch_size=sample_att.shape[0], shuffle=True, collate_fn=CollateFn(True))

    with torch.no_grad():
        for it, (features, labels, samples) in enumerate(dl_att):
            prediction, lst_att = net_mine_eval(features)
    score_att = dict(zip(range(14), (torch.sum(lst_att[0][1].squeeze([1, 2]), dim=0) / 14).tolist()))
    vis_mol(SMILES_dict['Temozolomide'], score_att, f"result/fig/att_{name}.svg")
    save_score(score_att, name)

    return(score_att)


score_ccle_all = flow('CCLE_glioma_temo_all.csv', 'CCLE_glioma_temo_all')
score_ccle_sensitive = flow('CCLE_glioma_temo_sensitive.csv', 'CCLE_glioma_temo_sensitive')

score_gdsc_all = flow('GDSC_glioma_temo_all.csv', 'GDSC_glioma_temo_all')
score_gdsc_sensitive = flow('GDSC_glioma_temo_sensitive.csv', 'GDSC_glioma_temo_sensitive')

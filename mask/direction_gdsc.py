import os
import sys

import joblib
import pandas as pd

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./models')
sys.path.append('./models/Mine')
sys.path.append('./models/DIPK')

from torch.utils.data import DataLoader
from torchmetrics.regression import PearsonCorrCoef

from utils.Data import *
from utils.Train import *

from models.model_Mine import net_mine, optimizer_mine
from models.model_DIPK import net_DIPK, optimizer_DIPK
from utils.TrainConfig import *
from copy import deepcopy

from scipy.stats import pearsonr


cancer = 'LGG'
drug = 'Temozolomide'

hallmark = pd.read_csv('Dataset/reference/hallmark_hsa.gmt', sep='\t', header=None)

pathway = {}
for i in range(hallmark.shape[0]):
    tmp = hallmark.iloc[i, :].to_list()
    pathway[tmp[0]] = tmp[2:]
    pathway[tmp[0]] = [ele for ele in pathway[tmp[0]] if not pd.isna(ele)]

genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()
gdsc = pd.read_csv('Dataset/sample/GDSC2.csv')
meta = pd.read_csv('Dataset/sample/GDSC_meta.txt', sep='\t')
meta['COSMIC_ID'] = meta['COSMIC_ID'].map(lambda x: 'COSMIC_' + str(x))
sample_sliced = set(meta['COSMIC_ID'].loc[meta['TCGA_DESC'].isin([cancer])])
gdsc_sliced = gdsc.loc[gdsc.Drug.isin([drug])].loc[gdsc.Cell.isin(sample_sliced)]

toKeep = set(gdsc_sliced.Cell.to_list())
RNA_dict_keep = {key: value for key, value in RNA_dict.items() if key in toKeep}

# mine
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

list(pathway.keys())

def mask_pathway(id, alias):
    RNA_dict_TCGA2 = deepcopy(RNA_dict_keep)
    genes_delete = list(set(pathway[id]).intersection(set(genes)))
    for ind in [i for i, gene in enumerate(genes) if gene in genes_delete]:
        for key in RNA_dict_TCGA2.keys():
            RNA_dict_TCGA2[key][ind] = 0
    dl_mine_tcga = DataLoader(MyDataSet(GetData(gdsc_sliced, RNA_dict=RNA_dict_TCGA2)), batch_size=batch_size, shuffle=True,
                              collate_fn=CollateFn(True))
    res_mine_tcga = predict_model(net_mine, dl_mine_tcga)
    pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_mine_tcga))).to_csv(f'result/direction_gdsc_LGG_Temozolomide_{alias}.csv', index=None)


mask_pathway('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT')
mask_pathway('HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'UPR')
mask_pathway('HALLMARK_HEDGEHOG_SIGNALING', 'HEDGEHOG')
mask_pathway('HALLMARK_TGF_BETA_SIGNALING', 'TGF-Beta')

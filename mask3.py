import os
import sys

import joblib
import pandas as pd
import torch

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./models')
sys.path.append('./models/Mine')
sys.path.append('./models/DIPK')

from torch.utils.data import DataLoader

from utils.Data import *
from utils.Train import *
from utils.TrainConfig import *

from models.model_Mine import net_mine, optimizer_mine

from torchmetrics.regression import PearsonCorrCoef
from scipy.stats import pearsonr

tcga = pd.read_csv('Dataset/sample/TCGA.csv')
meta = pd.read_csv('Dataset/sample/TCGA_meta.txt', sep='\t', encoding='GB2312')
meta['patient.arr'] = list(map(lambda x: x + '-01', meta['patient.arr']))
sample_LGG = set(meta.loc[meta.cancers.isin(['LGG'])]['patient.arr'].to_list())
glioma = tcga.loc[tcga.Cell.isin(sample_LGG).to_list() and tcga.Drug.isin(['Temozolomide']).to_list()]
toKeep = list(map(lambda x: x not in ['TCGA-DD-A4NS-01', 'TCGA-QT-A5XJ-01'], glioma.Cell))
toKeep = set(glioma.loc[toKeep].Cell.to_list())
RNA_dict_TCGA = joblib.load('Dataset/RNAseq_Dict_TCGA.pkl')
RNA_dict_keep = {key: value for key, value in RNA_dict_TCGA.items() if key in toKeep}
genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

# mine
dl_base = DataLoader(MyDataSet(GetData(glioma.loc[glioma.Cell.isin(toKeep)])), batch_size=batch_size, shuffle=True,collate_fn=CollateFn(True))
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
result_base = predict_model(net_mine, dl_base)
loss_base = list(pearsonr(result_base[1], result_base[0]))[0]
loss_base

lst_pred = joblib.load('result/lst_pred.pkl')

lst_loss = []
for i in range(len(lst_pred)):
    df = pd.DataFrame(dict(zip(['pred', 'label', 'sample'], lst_pred[i])))
    df['sample'] = df['sample'].map(lambda x: x.split('\t')[0])
    df = df.loc[df['sample'].isin(toKeep)]
    lst_loss.append(list(pearsonr(df.label, df.pred))[0])

pd.DataFrame({'gene': genes, 'delta': lst_loss - loss_base}).to_csv('result/mask2.csv', index=None)

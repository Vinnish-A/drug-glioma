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

from utils.Data import *
from utils.Train import *
from utils.TrainConfig import *

from models.model_Mine import net_mine, optimizer_mine

from scipy.stats import pearsonr

tcga = pd.read_csv('Dataset/sample/TCGA.csv')
meta = pd.read_csv('Dataset/sample/TCGA_meta.txt', sep='\t', encoding='GB2312')
sample_LGG = meta.loc[meta.cancers.isin(['LGG'])]['patient.arr'].to_list()
glioma = tcga.loc[tcga.Cell.isin(sample_LGG).to_list() and tcga.Drug.isin(['Temozolomide']).to_list()]

toKeep = set(glioma.Cell.to_list())
RNA_dict_TCGA = joblib.load('Dataset/RNAseq_Dict_TCGA.pkl')
RNA_dict_keep = {key: value for key, value in RNA_dict_TCGA.items() if key in toKeep}
genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

# mine
dl_base = DataLoader(MyDataSet(GetData(glioma)), batch_size=batch_size, shuffle=True,collate_fn=CollateFn(True))
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
result_base = predict_model(net_mine, dl_base)
loss_base = list(pearsonr(result_base[1], result_base[0]))[0]
loss_base

lst_pred = joblib.load('result/lst_pred.pkl')

lst_pred[0][2]
lst_pred[1][2]

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], lst_pred[0]))).sort_values(by=['sample'])
pd.DataFrame(dict(zip(['pred', 'label', 'sample'], lst_pred[1]))).sort_values(by=['sample'])

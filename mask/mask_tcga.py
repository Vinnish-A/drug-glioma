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

tcga = pd.read_csv('Dataset/sample/TCGA.csv')
meta = pd.read_csv('Dataset/sample/TCGA_meta.txt', sep='\t', encoding='GB2312')
sample_LGG = meta.loc[meta.cancers.isin(['LGG'])]['patient.arr'].to_list()
glioma = tcga.loc[tcga.Cell.isin(sample_LGG).to_list() and tcga.Drug.isin(['Temozolomide']).to_list()]

toKeep = set(glioma.Cell.to_list())
RNA_dict_TCGA = joblib.load('Dataset/RNAseq_Dict_TCGA.pkl')
RNA_dict_keep = {key: value for key, value in RNA_dict_TCGA.items() if key in toKeep}

# mine
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

def mask(ind):
    RNA_dict_TCGA2 = deepcopy(RNA_dict_keep)
    for key in RNA_dict_TCGA2.keys():
        RNA_dict_TCGA2[key][ind] = 0
    dl_mine_tcga = DataLoader(MyDataSet(GetData(glioma, RNA_dict=RNA_dict_TCGA2)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
    res_mine_tcga = predict_model(net_mine, dl_mine_tcga)
    return(res_mine_tcga)

def run_pool():  # main process
    from multiprocessing import Pool

    cpu_worker_num = 4
    process_args = [(i) for i in range(len(genes))]

    print(f'| inputs:  {process_args}')
    start_time = time.time()
    with Pool(cpu_worker_num) as p:
        lst_pred = p.map(mask, process_args)
    joblib.dump(lst_pred, f'result/lst_mask_Temozolomide.pkl')

if __name__ == '__main__':
    run_pool()

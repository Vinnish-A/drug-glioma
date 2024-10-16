import os
import sys

import joblib
import pandas as pd

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./models')
sys.path.append('./models/Mine')

from torch.utils.data import DataLoader

from utils.Data import *
from utils.Train import *

from models.model_Mine import net_mine, optimizer_mine
from utils.TrainConfig import *
from copy import deepcopy

from scipy.stats import pearsonr

drug = 'Temozolomide'

gdsc = pd.read_csv('Dataset/sample/GDSC2.csv')
gdsc_sliced = gdsc.loc[gdsc.Drug.isin([drug])]

toKeep = set(gdsc_sliced.Cell.to_list())
RNA_dict_keep = {key: value for key, value in RNA_dict.items() if key in toKeep}

net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

def mask(ind):
    RNA_dict_TCGA2 = deepcopy(RNA_dict_keep)
    for key in RNA_dict_TCGA2.keys():
        RNA_dict_TCGA2[key][ind] = 0
    dl_mine_tcga = DataLoader(MyDataSet(GetData(gdsc_sliced, RNA_dict=RNA_dict_TCGA2)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
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
    joblib.dump(lst_pred, f'result/lst_mask_{drug}.pkl')

if __name__ == '__main__':
    run_pool()

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

ccle = pd.read_csv('Dataset/sample/CCLE.csv')
tcga = pd.read_csv('Dataset/sample/TCGA.csv')
gdsc = pd.read_csv('Dataset/sample/GDSC2.csv')

RNA_dict_TCGA = joblib.load('Dataset/RNAseq_Dict_TCGA.pkl')

loss_func = nn.MSELoss()
metrics_dict = {'MSE': PearsonCorrCoef().to(DEVICE)}

# mine
dl_mine_tcga = DataLoader(MyDataSet(GetData(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))

net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
res_mine_tcga = predict_model(net_mine, dl_mine_tcga)

genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()

# lst_mask = []
# for i in range(len(genes)):
#     RNA_dict_TCGA2 = deepcopy(RNA_dict_TCGA)
#     for key in RNA_dict_TCGA2.keys():
#         RNA_dict_TCGA2[key][i] = 0
#     print(i)
#     dl_mine_tcga = DataLoader(MyDataSet(GetData(tcga, RNA_dict = RNA_dict_TCGA2)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
#     res_mine_tcga = predict_model(net_mine, dl_mine_tcga)
#     lst_mask.append(res_mine_tcga)
# joblib.dump(lst_mask, 'result/lst_mask.pkl')
#
#
# loss_base = list(pearsonr(res_mine_tcga[1], res_mine_tcga[0]))[0]
# lst_loss = []
# for i in range(len(lst_mask)):
#     loss_mask = list(pearsonr(lst_mask[i][1], lst_mask[i][0]))[0]
#     loss_delta = loss_mask - loss_base
#     lst_loss.append(loss_delta.tolist())


# pd.DataFrame({'gene': genes, 'delta': lst_loss}).to_csv('result/mask2.csv', index=None)

def mask(ind):
    RNA_dict_TCGA2 = deepcopy(RNA_dict_TCGA)
    for key in RNA_dict_TCGA2.keys():
        RNA_dict_TCGA2[key][ind] = 0
    dl_mine_tcga = DataLoader(MyDataSet(GetData(tcga, RNA_dict=RNA_dict_TCGA2)), batch_size=batch_size, shuffle=True,
                              collate_fn=CollateFn(True))
    res_mine_tcga = predict_model(net_mine, dl_mine_tcga)
    return(res_mine_tcga)

def run_pool():  # main process
    from multiprocessing import Pool

    cpu_worker_num = 6
    process_args = [(i) for i in range(len(genes))]

    print(f'| inputs:  {process_args}')
    start_time = time.time()
    with Pool(cpu_worker_num) as p:
        lst_pred = p.map(mask, process_args)
    joblib.dump(lst_pred, f'result/lst_pred.pkl')

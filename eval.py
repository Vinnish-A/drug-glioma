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

from models.model_Mine import net_mine
from models.model_DIPK import net_DIPK
from utils.TrainConfig import *

part1 = pd.read_csv('Dataset/Train.csv', sep='\t')
part2 = pd.read_csv('Dataset/Test.csv', sep='\t')
tcga = pd.read_csv('Dataset/TCGA.csv', sep='\t')
sample_all = pd.concat([part1, part2], axis=0)

loss_func = nn.MSELoss()
metrics_dict = {'MSE': PearsonCorrCoef().to(DEVICE)}

# mine
dl_test_mine = DataLoader(MyDataSet(GetData(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))

net_mine.load_state_dict(torch.load('checkpoint/mine_trained.pt'))
res_tcga_mine = predict_model(net_mine, dl_test_mine)

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_tcga_mine))).to_csv('result/TCGA_mine.csv', index=None)

# DIPK
dl_test_DIPK = DataLoader(MyDataSet(GetData_DIPK(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn_DIPK(True))

net_DIPK.load_state_dict(torch.load('checkpoint/DIPK_cv4.pt'))
res_tcga_DIPK = predict_model(net_DIPK, dl_test_DIPK)

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_tcga_DIPK))).to_csv('result/TCGA_DIPK.csv', index=None)

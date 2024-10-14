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

ccle = pd.read_csv('Dataset/sample/CCLE.csv')
tcga = pd.read_csv('Dataset/sample/TCGA.csv')
gdsc = pd.read_csv('Dataset/sample/GDSC2.csv')

loss_func = nn.MSELoss()
metrics_dict = {'MSE': PearsonCorrCoef().to(DEVICE)}

# mine
dl_mine_tcga = DataLoader(MyDataSet(GetData(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
dl_mine_gdsc = DataLoader(MyDataSet(GetData(gdsc)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))

net_mine.load_state_dict(torch.load('checkpoint/mine_trained.pt'))
res_mine_tcga = predict_model(net_mine, dl_mine_tcga)
res_mine_gdsc = predict_model(net_mine, dl_mine_gdsc)

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_mine_tcga))).to_csv('result/TCGA_mine.csv', index=None)
pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_mine_gdsc))).to_csv('result/GDSC_mine.csv', index=None)

# DIPK
dl_DIPK_tcga = DataLoader(MyDataSet(GetData_DIPK(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn_DIPK(True))
dl_DIPK_gdsc = DataLoader(MyDataSet(GetData_DIPK(gdsc)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn_DIPK(True))


net_DIPK.load_state_dict(torch.load('checkpoint/DIPK_cv4.pt'))
res_DIPK_tcga = predict_model(net_DIPK, dl_DIPK_tcga)
res_DIPK_gdsc = predict_model(net_DIPK, dl_DIPK_gdsc)

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_DIPK_tcga))).to_csv('result/TCGA_DIPK.csv', index=None)
pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_DIPK_gdsc))).to_csv('result/GDSC_DIPK.csv', index=None)

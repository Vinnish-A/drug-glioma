import
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
from torchmetrics.regression import PearsonCorrCoef

from utils.Data import *
from utils.Train import *

from models.model_Mine import net_mine, optimizer_mine
from utils.TrainConfig import *

loss_func = nn.MSELoss()
metrics_dict = {'PCC': PearsonCorrCoef().to(DEVICE)}
net_mine.load_state_dict(torch.load('checkpoint/mine_trained.pt'))


# ccle
ccle = pd.read_csv('Dataset/generated/CCLE_generated.csv')
dl_mine_ccle = DataLoader(MyDataSet(GetData(ccle)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
res_mine_ccle = predict_model(net_mine, dl_mine_ccle)
pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_mine_ccle))).to_csv('result/screen_CCLE.csv', index=None)

# gdsc
gdsc = pd.read_csv('Dataset/generated/GDSC_generated.csv')
dl_mine_gdsc = DataLoader(MyDataSet(GetData(gdsc)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
res_mine_gdsc = predict_model(net_mine, dl_mine_gdsc)
pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_mine_gdsc))).to_csv('result/screen_GDSC.csv', index=None)

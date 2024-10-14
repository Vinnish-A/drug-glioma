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
from utils.TrainConfig import *

ccle = pd.read_csv('Dataset/sample/CCLE.csv')
tcga = pd.read_csv('Dataset/sample/TCGA.csv')

loss_func = nn.MSELoss()
metrics_dict = {'MSE': PearsonCorrCoef().to(DEVICE)}

dl_train_mine = DataLoader(MyDataSet(GetData(ccle)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
dl_test_mine  = DataLoader(MyDataSet(GetData(tcga)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))

dfhistory_mine = train_model(
    net_mine,
    optimizer_mine,
    loss_func,
    metrics_dict,
    train_data=dl_train_mine,
    val_data=dl_test_mine,
    epochs=EPOCHS,
    patience=EPOCHS,
    monitor="val_MSE",
    mode="min",
    ckpt_path=f'checkpoint/mine_trained.pt'
)
net_mine.load_state_dict(torch.load('checkpoint/mine_trained.pt'))
res_pair_tcga = predict_model(net_mine, dl_test_mine)

pd.DataFrame(dict(zip(['pred', 'label', 'sample'], res_pair_tcga))).to_csv('result/TCGA_mine.csv')

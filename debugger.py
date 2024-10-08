import os
import sys

import joblib

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

part1 = pd.read_csv('Dataset/Train.csv', sep='\t')
part2 = pd.read_csv('Dataset/Test.csv', sep='\t')
tcga = pd.read_csv('Dataset/TCGA.csv', sep='\t')
sample_all = pd.concat([part1, part2], axis=0)

loss_func = nn.MSELoss()
metrics_dict = {'MSE': PearsonCorrCoef().to(DEVICE)}

dl_train_mine = DataLoader(MyDataSet(GetData(sample_all)), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
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
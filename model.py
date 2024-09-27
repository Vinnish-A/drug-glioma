import os; os.chdir('F:/0Local/project/drug-glioma')
import sys; sys.path.append(os.getcwd()); sys.path.append('F:/0Local/project/drug-glioma/utils')

import pandas as pd
import torch.optim as optim
from torch.utils.data import DataLoader
import time

import importlib
from utils.Model import *
from utils.Data import *
from utils.TrainConfig import *
from utils.BM import *
from utils.Train import *
from torchmetrics import MeanSquaredError

sampler = setup_seed(seed)
test_pre_rec = dict()

# create model
encoder, _ = joblib.load('Dataset/PreTrain.pkl')
encoder = encoder.to(DEVICE)
biology = BiologicalModule(df_kegg.shape[1], 297, 512, df_kegg.T.nonzero()).to(DEVICE)
model = Predictor(embedding_dim, heads, fc_layer_num, fc_layer_dim, dropout_rate).to(DEVICE)
loss_func = nn.MSELoss()
params = [
    {'params': encoder.parameters(), 'lr': 0.5 * lr, 'weight_decay': 0.0001},
    {'params': biology.parameters(), 'lr': 0.5 * lr, 'weight_decay': 0.0001},
    {'params': model.parameters(), 'weight_decay': 0.0001}
]
optimizer = optim.Adam(params, lr=lr)
net = Net(encoder, biology, model)

# load data
my_collate = CollateFn()
train_loader = DataLoader(MyDataSet(GetData(pd.read_csv('Dataset/Train.csv', sep='\t'))), batch_size=batch_size, shuffle=True, collate_fn=my_collate)
test_loader = DataLoader(MyDataSet(GetData(pd.read_csv('Dataset/Test.csv', sep='\t'))), batch_size=batch_size, shuffle=False, collate_fn=my_collate)

metrics_dict = {"MSE": MeanSquaredError().to(DEVICE)}

dfhistory = train_model(
    net,
    optimizer,
    loss_func,
    metrics_dict,
    train_data = train_loader,
    val_data= test_loader,
    epochs=EPOCHS,
    patience=EPOCHS,
    monitor="val_MSE",
    mode="min",
    ckpt_path='checkpoint/checkpoint.pt'
)

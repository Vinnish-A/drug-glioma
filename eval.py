import os
import sys

import joblib
import pandas as pd

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./others')
sys.path.append('./others/DIPK')

from torch.utils.data import DataLoader
from torchmetrics import MeanSquaredError

from utils.Data import *
from utils.Train import *

from others.model_DIPK import net_DIPK, optimizer_DIPK, CollateFn_DIPK
from utils.Mine import net_mine, optimizer_mine, CollateFn
from utils.Model import setup_seed
from utils.TrainConfig import *

from scipy.stats import pearsonr

sampler = setup_seed(seed)


# mine

collate_mine = CollateFn(True)
tcga_loader = DataLoader(MyDataSet(GetData(pd.read_csv('Dataset/TCGA.csv', sep='\t'))), batch_size=batch_size, shuffle=True, collate_fn=collate_mine)
net_mine.load_state_dict(torch.load('checkpoint/200epoch.pt'))
res_mine = predict_model(net_mine, tcga_loader)

pearsonr(res_mine[1], res_mine[0])

# pd.DataFrame({'sample': samples, 'drug': drugs, 'response': true, 'pred': pred}).to_csv('result/pred_cnv.csv', index=0)

# DIPK

collate_DIPK = CollateFn_DIPK(True)
tcga_loader_DIPK = DataLoader(MyDataSet(GetData_DIPK(pd.read_csv('Dataset/TCGA.csv', sep='\t'))), batch_size=batch_size, shuffle=True, collate_fn=collate_DIPK)
net_DIPK.load_state_dict(torch.load('checkpoint/DIPK_cv4.pt'))
res_DIPK = predict_model(net_DIPK, tcga_loader_DIPK)

pearsonr(res_DIPK[1], res_DIPK[0])

tcga_loader_DIPK = DataLoader(MyDataSet(GetData_DIPK('Dataset/TCGA.csv')), batch_size=batch_size, shuffle=True, collate_fn=collate_DIPK)

epoch_DIPK, model_DIPK, encoder_DIPK, test_pre_rec_DIPK = joblib.load('result/Train_DIPK.pkl')
model_DIPK.eval()
encoder_DIPK.eval()
pred_DIPK = list()
true_DIPK = list()
samples = list()
drugs = list()
with torch.no_grad():
    for it, (pyg_batch, GeneFt, BionicFt, sample, drug) in enumerate(tcga_loader_DIPK):
        pyg_batch, GeneFt, BionicFt = pyg_batch.to(DEVICE), GeneFt.to(DEVICE), BionicFt.to(DEVICE)
        GeneFt = encoder_DIPK(GeneFt)
        prediction = model_DIPK(pyg_batch.x, pyg_batch, GeneFt, BionicFt)
        pred_DIPK += torch.squeeze(prediction).tolist()
        true_DIPK += pyg_batch.ic50.tolist()
        samples += sample
        drugs +=drug

pearsonr(true_DIPK, pred_DIPK)

pd.DataFrame({'sample': samples, 'drug': drugs, 'response': true_DIPK, 'pred': pred_DIPK}).to_csv('result/pred_DIPK.csv', index=0)

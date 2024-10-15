import os;

import pandas as pd
import torch

os.chdir('F:/0Local/project/drug-glioma')
import sys; sys.path.append(os.getcwd()); sys.path.append('F:/0Local/project/drug-glioma/utils')

import torch.optim as optim
from torch.utils.data import DataLoader
import time

import importlib
from utils.Model import *
from utils.Data import *
from utils.TrainConfig import *
from utils.BM import *
from scipy.stats import pearsonr
from collections import Counter
import shap
import random

sampler = setup_seed(seed)

epoch, model, encoder, biology, test_pre_rec = joblib.load('result/Train.pkl')

phen_tcga = pd.read_csv('Dataset/TCGA.csv', header=0, sep='\t')
data_response = pd.read_csv('DataPreprocess/TPM/drug_response.txt', sep='\t', encoding='GB2312')
data_response['patient.arr'] += '-01'

toKeep = data_response.loc[data_response.cancers.isin(['GBM', 'LGG'])].loc[data_response['drug.name'] == 'Temozolomide']['patient.arr'].tolist()

df = phen_tcga.loc[phen_tcga['patient.arr'].isin(toKeep), :].loc[phen_tcga['drug.name'] == 'Temozolomide']
sample_train = pd.read_csv('Dataset/Train.csv', sep='\t', header=0)
df_background = sample_train.loc[sample_train.Drug.isin(['Temozolomide'])]


# def specify(df, model = model, encoder = encoder, biology = biology):
#     my_collate = CollateFnEval()
#
#     explain_loader = DataLoader(MyDataSet(TransData(df)), batch_size=len(df) + 1, shuffle=True, collate_fn=my_collate)
#     for it, (pyg_batch, GeneFt, BionicFt, sample, drug) in enumerate(explain_loader):
#         pyg_batch, rna1, rna2 = pyg_batch.to(DEVICE), GeneFt.to(DEVICE), BionicFt.to(DEVICE)
#
#     class Wraped(nn.Module):
#         def __init__(self, model, encoder, biology):
#             super(Wraped, self).__init__()
#             self.model = model
#             self.encoder = encoder
#             self.biology = biology
#
#         def forward(self, rna, graph):
#
#             GeneFt = encoder(rna)
#             BionicFt = biology(rna)
#             prediction = model(graph.x, graph, GeneFt, BionicFt)
#
#             return prediction
#
#     wraped = Wraped(model, encoder, biology)
#
#     return rna1, pyg_batch, wraped
#
#
# RNA, batch, wraped = specify(df)
# RNA_background, batch_background, wraped = specify(df_background.iloc[random.sample(range(len(df_background)), 104)])
#
# explainer = shap.DeepExplainer(wraped, [RNA_background, batch_background])
# shap_values = explainer.shap_values([RNA, batch])

my_collate = CollateFnEval()
explain_loader = DataLoader(MyDataSet(GetData(df)), batch_size=len(df) + 1, shuffle=True, collate_fn=my_collate)
model.eval()
encoder.eval()
biology.eval()
pred = list()
true = list()
samples = list()
drugs = list()
with torch.no_grad():
    for it, (pyg_batch, RNA, BionicFt, sample, drug) in enumerate(explain_loader):
        pyg_batch, RNA, BionicFt = pyg_batch.to(DEVICE), RNA.to(DEVICE), BionicFt.to(DEVICE)
        GeneFt = encoder(RNA)
        BionicFt = biology(RNA)
        prediction = model(pyg_batch.x, pyg_batch, GeneFt, BionicFt)
        pred += torch.squeeze(prediction).tolist()
        true += pyg_batch.ic50.tolist()
        samples += sample
        drugs +=drug

start = pearsonr(true, pred)[0]

loss = []
for i in range(5863):
    model.eval()
    encoder.eval()
    biology.eval()
    pred = list()
    true = list()
    samples = list()
    drugs = list()
    with torch.no_grad():
        for it, (pyg_batch, RNA, BionicFt, sample, drug) in enumerate(explain_loader):
            pyg_batch, RNA, BionicFt = pyg_batch.to(DEVICE), RNA.to(DEVICE), BionicFt.to(DEVICE)
            RNA[:, i] = 0
            GeneFt = encoder(RNA)
            BionicFt = biology(RNA)
            prediction = model(pyg_batch.x, pyg_batch, GeneFt, BionicFt)
            pred += torch.squeeze(prediction).tolist()
            true += pyg_batch.ic50.tolist()
            samples += sample
            drugs += drug
    print(i)
    loss.append(start - pearsonr(true, pred)[0])

loss_func = nn.MSELoss()
my_collate = CollateFnEval()
explain_loader = DataLoader(MyDataSet(GetData(df_background)), batch_size=64, shuffle=True, collate_fn=my_collate)
model.eval()
encoder.eval()
biology.eval()
pred = list()
true = list()
samples = list()
drugs = list()
epoch_loss = 0
with torch.no_grad():
    for it, (pyg_batch, RNA, BionicFt, sample, drug) in enumerate(explain_loader):
        pyg_batch, RNA, BionicFt = pyg_batch.to(DEVICE), RNA.to(DEVICE), BionicFt.to(DEVICE)
        GeneFt = encoder(RNA)
        BionicFt = biology(RNA)
        prediction = model(pyg_batch.x, pyg_batch, GeneFt, BionicFt)
        loss = loss_func(torch.squeeze(prediction), pyg_batch.ic50)
        epoch_loss += loss

start = epoch_loss

lst = []
for i in range(5863):
    model.eval()
    encoder.eval()
    biology.eval()
    pred = list()
    true = list()
    samples = list()
    drugs = list()
    epoch_loss = 0
    with torch.no_grad():
        for it, (pyg_batch, RNA, BionicFt, sample, drug) in enumerate(explain_loader):
            pyg_batch, RNA, BionicFt = pyg_batch.to(DEVICE), RNA.to(DEVICE), BionicFt.to(DEVICE)
            GeneFt = encoder(RNA)
            BionicFt = biology(RNA)
            prediction = model(pyg_batch.x, pyg_batch, GeneFt, BionicFt)
            loss = loss_func(torch.squeeze(prediction), pyg_batch.ic50)
            epoch_loss += loss
    print(i)
    lst.append((start - epoch_loss).tolist())

lst_mask = joblib.load('result/loss_mask.pkl')
lst_gene = pd.read_csv('DataPreprocess/TPM/gene_list_sel.txt', header=None).iloc[:,0].to_list()

pd.DataFrame({'gene': lst_gene, 'mask': lst_mask}).to_csv('result/mask.csv', index=0)

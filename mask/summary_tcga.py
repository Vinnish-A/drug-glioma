import os
import sys

os.chdir('F:/0Local/project/drug-glioma')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./models')
sys.path.append('./models/Mine')
sys.path.append('./models/DIPK')

from torch.utils.data import DataLoader
from collections import Counter
from scipy.stats import pearsonr

from utils.Data import *
from utils.Train import *
from utils.TrainConfig import *

from models.model_Mine import net_mine, optimizer_mine

drug = 'Temozolomide'
cancer = 'LGG'
strangers = ['TCGA-DD-A4NS-01', 'TCGA-QT-A5XJ-01']

genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()
tcga = pd.read_csv('Dataset/sample/TCGA.csv')
meta = pd.read_csv('Dataset/sample/TCGA_meta.txt', sep='\t', encoding='GB2312')
meta['patient.arr'] = meta['patient.arr'].map(lambda x: x + '-01')

sample_sliced = set(meta.loc[meta.cancers.isin([cancer])]['patient.arr'].to_list())
tcga_sliced = tcga.loc[tcga.Cell.isin(sample_sliced)].loc[tcga.Drug.isin([drug])]

# Temozolomide
dl_base = DataLoader(MyDataSet(GetData(tcga_sliced)), batch_size=batch_size, shuffle=True,collate_fn=CollateFn(True))
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
result_base = predict_model(net_mine, dl_base)
loss_base = list(pearsonr(result_base[1], result_base[0]))[0]

lst_pred = joblib.load('result/lst_pred.pkl')

lst_loss = []
for i in range(len(lst_pred)):
    df = pd.DataFrame(dict(zip(['pred', 'label', 'sample'], lst_pred[i])))
    df['sample'] = df['sample'].map(lambda x: x.split('\t')[0])
    df = df.loc[df['sample'].isin(sample_sliced)]
    lst_loss.append(list(pearsonr(df.label, df.pred))[0])

pd.DataFrame({'gene': genes, 'delta': lst_loss - loss_base}).to_csv(f'result/mask_{cancer}_{drug}.csv', index=None)



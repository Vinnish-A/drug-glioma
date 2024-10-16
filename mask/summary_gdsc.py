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

cancer = 'LGG'
drug = 'Temozolomide'

genes = pd.read_csv('Dataset/symbols.txt', header=None).loc[:, 0].to_list()
gdsc = pd.read_csv('Dataset/sample/GDSC2.csv')
meta = pd.read_csv('Dataset/sample/GDSC_meta.txt', sep='\t')
meta['COSMIC_ID'] = meta['COSMIC_ID'].map(lambda x: 'COSMIC_' + str(x))
sample_sliced = set(meta['COSMIC_ID'].loc[meta['TCGA_DESC'].isin([cancer])])
gdsc_sliced = gdsc.loc[gdsc.Drug.isin([drug])].loc[gdsc.Cell.isin(sample_sliced)]

# Temozolomide
dl_base = DataLoader(MyDataSet(GetData(gdsc_sliced)), batch_size=batch_size, shuffle=True,collate_fn=CollateFn(True))
net_mine.load_state_dict(torch.load('checkpoint/tcga_trained.pt'))
result_base = predict_model(net_mine, dl_base)
loss_base = list(pearsonr(result_base[1], result_base[0]))[0]

lst_pred = joblib.load(f'result/lst_mask_gdsc_all_{drug}.pkl')

lst_loss = []
for i in range(len(lst_pred)):
    df = pd.DataFrame(dict(zip(['pred', 'label', 'sample'], lst_pred[i])))
    df['sample'] = df['sample'].map(lambda x: x.split('\t')[0])
    df = df.loc[df['sample'].isin(sample_sliced)]
    lst_loss.append(list(pearsonr(df.label, df.pred))[0])

pd.DataFrame({'gene': genes, 'delta': lst_loss - loss_base}).to_csv(f'result/mask_GDSC_{cancer}_{drug}.csv', index=None)


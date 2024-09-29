import os
import sys

import joblib

os.chdir('.')
sys.path.append(os.getcwd())
sys.path.append('./utils')
sys.path.append('./others')
sys.path.append('./others/DIPK')

from torch.utils.data import DataLoader
from torchmetrics import MeanSquaredError

from utils.Data import *
from utils.Train import *

from others.model_DIPK import net_DIPK, optimizer_DIPK
from utils.Mine import net_mine, optimizer_mine
from utils.Model import setup_seed
from utils.TrainConfig import *

sample_all = pd.read_csv('Dataset/GDSC2.csv', sep='\t')

sampler = setup_seed(seed)
# splits = split(sample_all.shape[0])
# joblib.dump(splits, 'checkpoint/splits_GDSC.pkl')
splits = joblib.load('checkpoint/splits_GDSC.pkl')
loss_func = nn.MSELoss()
metrics_dict = {'MSE': MeanSquaredError().to(DEVICE)}

res_history_mine = []
res_pairs_mine = []
for i in range(len(splits)):
    net_init_mine = net_mine
    net_init_mine.load_state_dict(torch.load('checkpoint/origin_mine.pt'))
    ind_train = splits[i]['train']
    ind_test  = splits[i]['test']
    dl_train_mine = DataLoader(MyDataSet(GetData(sample_all.iloc[ind_train])), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))
    dl_test_mine  = DataLoader(MyDataSet(GetData(sample_all.iloc[ind_test])), batch_size=batch_size, shuffle=True, collate_fn=CollateFn(True))

    dfhistory_mine = train_model(
        net_init_mine,
        optimizer_mine,
        loss_func,
        metrics_dict,
        train_data=dl_train_mine,
        val_data=dl_test_mine,
        epochs=EPOCHS,
        patience=EPOCHS,
        monitor="val_MSE",
        mode="min",
        ckpt_path=f'checkpoint/mine_cv_GDSC{i+1}.pt'
    )
    res_pair_mine = predict_model(net_init_mine, dl_test_mine)

    res_history_mine.append(dfhistory_mine)
    res_pairs_mine.append(res_pair_mine)
    if i + 1 < len(splits):
        joblib.dump({'history': res_history_mine, 'pair': res_pairs_mine}, f'result/cv_mine_GDSC{i + 1}.pkl')
    else:
        joblib.dump({'history': res_history_mine, 'pair': res_pairs_mine}, 'result/cv_mine_mine.pkl')
        files_to_delete = [f'result/cv_mine_GDSC{i + 1}.pkl' for i in range(len(splits)-1)]
        for file_to_delete in files_to_delete:
            os.remove(file_to_delete)

joblib.dump({'history': res_history_mine, 'pair': res_pairs_mine}, 'result/cv_mine_GDSC.pkl')

EPOCHS = 80
batch_size = 128
res_history_DIPK = []
res_pairs_DIPK = []
for i in range(len(splits)):
    net_init_DIPK = net_DIPK
    net_init_DIPK.load_state_dict(torch.load('checkpoint/origin_DIPK.pt'))
    ind_train = splits[i]['train']
    ind_test  = splits[i]['test']
    dl_train_DIPK = DataLoader(MyDataSet(GetData_DIPK(sample_all.iloc[ind_train])), batch_size=batch_size, shuffle=True, collate_fn=CollateFn_DIPK(True))
    dl_test_DIPK  = DataLoader(MyDataSet(GetData_DIPK(sample_all.iloc[ind_test])), batch_size=batch_size, shuffle=True, collate_fn=CollateFn_DIPK(True))

    dfhistory_DIPK = train_model(
        net_init_DIPK,
        optimizer_DIPK,
        loss_func,
        metrics_dict,
        train_data=dl_train_DIPK,
        val_data=dl_test_DIPK,
        epochs=EPOCHS,
        patience=EPOCHS,
        monitor="val_MSE",
        mode="min",
        ckpt_path=f'checkpoint/DIPK_cv_GDSC{i+1}.pt'
    )
    res_pair_DIPK = predict_model(net_init_DIPK, dl_test_DIPK)

    res_history_DIPK.append(dfhistory_DIPK)
    res_pairs_DIPK.append(res_pair_DIPK)
    if i+1<len(splits):
        joblib.dump({'history': res_history_DIPK, 'pair': res_pairs_DIPK}, f'result/cv_DIPK_GDSC{i+1}.pkl')
    else:
        joblib.dump({'history': res_history_DIPK, 'pair': res_pairs_DIPK}, 'result/cv_DIPK_GDSC.pkl')
        files_to_delete = [f'result/cv_DIPK_GDSC{i+1}.pkl' for i in range(len(splits)-1)]
        for file_to_delete in files_to_delete:
            os.remove(file_to_delete)

joblib.dump({'history': res_history_DIPK, 'pair': res_pairs_DIPK}, 'result/cv_DIPK_GDSC.pkl')

# summary

cv_mine = joblib.load('result/cv_mine.pkl')
cv_DIPK = joblib.load('result/cv_DIPK.pkl')

MSE = MeanSquaredError()

pair_mine = cv_mine.get('pair')
pair_DIPK = cv_DIPK.get('pair')

def cal_cv(fun, pair):

    res = []
    for i in range(len(pair)):
        preds = torch.tensor([float(_) for _ in pair[i][0]])
        trues = torch.tensor([float(_) for _ in pair[i][1]])
        res.append(fun(preds, trues).tolist())

    return res

sum(cal_cv(MSE, pair_mine))/5
sum(cal_cv(MSE, pair_DIPK))/5


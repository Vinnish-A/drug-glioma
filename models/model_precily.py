import torch.optim as optim
import joblib
from DIPK.Model import *
from DIPK.Data import CollateFn_DIPK
from DIPK.TrainConfig import *
from DIPK.BM import DEVICE

sampler = setup_seed(seed)
test_pre_rec = dict()

# create model
encoder_DIPK, _ = joblib.load('Dataset/PreTrain.pkl')
encoder_DIPK = encoder_DIPK.to(DEVICE)

model_DIPK = Predictor_DIPK(embedding_dim, heads, fc_layer_num, fc_layer_dim, dropout_rate).to(DEVICE)
net_DIPK = Net_DIPK(encoder_DIPK, model_DIPK)

params = [
    {'params': encoder_DIPK.parameters(), 'lr': 0.5 * lr, 'weight_decay': 0.0001},
    {'params': model_DIPK.parameters(), 'weight_decay': 0.0001}
]
optimizer_DIPK = optim.Adam(params, lr=lr)


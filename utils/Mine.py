import joblib
import torch.optim as optim
from Model import *
from TrainConfig import *
from BM import *

sampler = setup_seed(seed)
test_pre_rec = dict()

# create model
df_kegg = joblib.load('Dataset/df_kegg.pkl').to_numpy()
encoder, _ = joblib.load('Dataset/PreTrain.pkl')
encoder = encoder.to(DEVICE)
biology = BiologicalModule(df_kegg.shape[1], 297, 512, df_kegg.T.nonzero()).to(DEVICE)
model = Predictor(embedding_dim, heads, fc_layer_num, fc_layer_dim, dropout_rate).to(DEVICE)
params = [
    {'params': encoder.parameters(), 'lr': 0.5 * lr, 'weight_decay': 0.0001},
    {'params': biology.parameters(), 'lr': 0.5 * lr, 'weight_decay': 0.0001},
    {'params': model.parameters(), 'weight_decay': 0.0001}
]

optimizer_mine = optim.Adam(params, lr=lr)
net_mine = Net(encoder, biology, model)
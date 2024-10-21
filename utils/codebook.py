import pandas as pd
import joblib
import torch
from rdkit.Chem import Descriptors, MolFromSmiles

RNA_dict = joblib.load('Dataset/RNAseq_Dict.pkl')
MolGNet_dict = joblib.load('Dataset/MolGNet_dict.pkl')
CNV_dict = joblib.load('Dataset/CNV_Dict.pkl')
SMILES_dict = joblib.load('Dataset/SMILES_dict.pkl')

pd.DataFrame({'sample': list(RNA_dict.keys())}).to_csv('Dataset/codebook/RNA_sample.txt', index=None, sep='\t')
pd.DataFrame({'drug': list(MolGNet_dict.keys())}).to_csv('Dataset/codebook/drug_item.txt', index=None, sep='\t')
pd.DataFrame({'gene': list(CNV_dict.keys())}).to_csv('Dataset/codebook/CNV_item.txt', index=None, sep='\t')

mol = MolFromSmiles(SMILES_dict['Mirin'])
logP = Descriptors.MolLogP(mol)
mol_weight = Descriptors.MolWt(mol)

from rdkit import Chem
from rdkit.Chem import Descriptors


def esol_predict(mol):
    # 使用 MolLogP, TPSA, 分子量 和 环结构预测水溶性
    logP = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    mol_weight = Descriptors.MolWt(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    # 简单的 ESOL 预测公式 (单位: logS)
    solubility = 0.16 - 0.63 * logP - 0.0062 * mol_weight + 0.066 * num_rotatable_bonds - 0.74

    return solubility


# 计算分子水溶性 (例如乙醇)
smiles = 'CCO'
mol = Chem.MolFromSmiles(smiles)

solubility = esol_predict(mol)

print(f"估算的水溶性 (logS): {solubility}")

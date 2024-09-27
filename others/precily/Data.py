import pandas as pd
import numpy as np
import copy

# genes-pathways annotation

path = 'scratch/20230205_kegg_hsa.gmt'

files = open(path,encoding='utf-8')

files = files.readlines()

paways_genes_dict = {}
for i in files:
    paways_genes_dict[i.split('\t')[0].split('_')[0]] = i.replace('\n','').split('\t')[2:]

#mirna-pathways annotation
path = 'scratch/kegg_anano.txt'

files = open(path,encoding='utf-8')

files = files.readlines()

paways_mirna_dict = {}
for i in files:
     keys = i.split(',')[0].split('|')[1]
     values1 = i.split(',')[1:-1]
     values2 =  i.split(',')[-1].replace('\n','')
     values1.append(values2)
     values1 =list(set(values1))
     paways_mirna_dict[keys] = values1

union_kegg = list(set(paways_genes_dict.keys()).intersection(set(paways_mirna_dict.keys())))

paways_genes_dicts = {}
paways_mirna_dicts = {}

for i in union_kegg:
    paways_genes_dicts[i] = paways_genes_dict[i]

for i in union_kegg:
    paways_mirna_dicts[i] = paways_mirna_dict[i]


genes_existed_pathway = []

mirna_existed_pathway = []

for index in paways_genes_dicts.keys():
    genes_existed_pathway = genes_existed_pathway+ list(paways_genes_dicts[index])
genes_existed_pathway = set(genes_existed_pathway)


for index in paways_mirna_dicts.keys():
    mirna_existed_pathway = mirna_existed_pathway+ list(paways_mirna_dicts[index])
mirna_existed_pathway = set(mirna_existed_pathway)

union_gene_snv = list(snv_data.columns)
union_gene_miRNA = list(miRNA_data.columns)
union_gene_mRNA = list(mRNA_data.columns)

pathway_union = list(paways_genes_dicts.keys())

mask_list = [union_gene_snv, union_gene_mRNA]

gene_pathway_bp_dfs = []

for i in range(len(mask_list)):
    pathways_genes = np.zeros((len(pathway_union), len(mask_list[i])))
    for p in pathway_union:
        gs = paways_genes_dicts[p]
        g_inds = [mask_list[i].index(g) for g in gs if g in mask_list[i]]
        p_ind = pathway_union.index(p)
        pathways_genes[p_ind, g_inds] = 1
    gene_pathway_bp = pd.DataFrame(pathways_genes, index=pathway_union, columns=mask_list[i])

    #     gene_pathway_bp = gene_pathway_bp.loc[:, (gene_pathway_bp != 0).any(axis=0)]
    gene_pathway_bp_dfs.append(gene_pathway_bp)

pathways_genes = np.zeros((len(pathway_union), len(union_gene_miRNA)))
for p in pathway_union:
    gs = paways_mirna_dicts[p]
    g_inds = [union_gene_miRNA.index(g) for g in gs if g in union_gene_miRNA]
    p_ind = pathway_union.index(p)
    pathways_genes[p_ind, g_inds] = 1
gene_pathway_bp = pd.DataFrame(pathways_genes, index=pathway_union, columns=union_gene_miRNA)

#     gene_pathway_bp = gene_pathway_bp.loc[:, (gene_pathway_bp != 0).any(axis=0)]
gene_pathway_bp_dfs.append(gene_pathway_bp)



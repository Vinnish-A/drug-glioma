
# make --------------------------------------------------------------------


library(tidyverse)

ccle = read_csv('Dataset/sample/CCLE.csv')
gdsc = read_csv('Dataset/sample/GDSC2.csv')
ccle_meta = read_tsv('Dataset/sample/CCLE_meta.txt')
gdsc_meta = read_tsv('Dataset/sample/GDSC_meta.txt')

## attention ----

# ccle
cell_glioma_temo_ccle_all = ccle_meta |> 
  filter(Histology == 'glioma') |> 
  pull(CCLE_ID) |> 
  str_split('_', 2) |> 
  map_chr(`[[`, 1) |> 
  unique()

cell_glioma_temo_ccle_sensitive = ccle_meta |> 
  filter(Histology == 'glioma') |> 
  mutate(Cell = CCLE_ID |> str_split('_', 2) |> map_chr(`[[`, 1)) |> 
  inner_join(ccle, by = 'Cell') |> 
  filter(Drug == 'Temozolomide') |> 
  slice_min(IC50, prop = 0.5) |> 
  pull(Cell)

ccle |> 
  filter(Cell %in% cell_glioma_temo_ccle_all) |> 
  filter(Drug == 'Temozolomide') |> 
  write_csv('Dataset/sliced/CCLE_glioma_temo_all.csv')

ccle |> 
  filter(Cell %in% cell_glioma_temo_ccle_sensitive) |> 
  filter(Drug == 'Temozolomide') |> 
  write_csv('Dataset/sliced/CCLE_glioma_temo_sensitive.csv')

# gdsc
cell_glioma_temo_gdsc_all = gdsc_meta |> 
  filter(TCGA_DESC == 'GBM') |> 
  filter(DRUG_NAME == 'Temozolomide') |> 
  mutate(COSMIC_ID = paste0('COSMIC_', COSMIC_ID)) |> 
  pull(COSMIC_ID) |> 
  unique()

cell_glioma_temo_gdsc_sensitive = gdsc_meta |> 
  filter(TCGA_DESC == 'GBM') |> 
  filter(DRUG_NAME == 'Temozolomide') |> 
  slice_min(LN_IC50, prop = 0.5) |> 
  mutate(COSMIC_ID = paste0('COSMIC_', COSMIC_ID)) |> 
  pull(COSMIC_ID) |> 
  unique()

gdsc |> 
  filter(Cell %in% cell_glioma_temo_gdsc_all) |> 
  filter(Drug == 'Temozolomide') |> 
  write_csv('Dataset/sliced/GDSC_glioma_temo_all.csv')

gdsc |> 
  filter(Cell %in% cell_glioma_temo_gdsc_sensitive) |> 
  filter(Drug == 'Temozolomide') |> 
  write_csv('Dataset/sliced/GDSC_glioma_temo_sensitive.csv')

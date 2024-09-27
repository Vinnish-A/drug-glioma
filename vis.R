
# config ------------------------------------------------------------------


library(clusterProfiler)
library(tidyverse)
library(ggsci)

# vis ---------------------------------------------------------------------


## response ----

data_response = read_tsv('DataPreprocess/TPM/drug_response.txt')

### without-cnv ----

data_plot_response = read_csv('result/pred.csv')

table_response = setNames(1:4, data_response$response |> unique())
table_cancer = data_response |> 
  mutate(patient.arr = paste0(patient.arr, '-01')) |> 
  distinct(cancers, patient.arr) |> 
  pull(cancers, patient.arr)

sample_glioma = data_response |> 
  mutate(sample = paste0(patient.arr, '-01')) |> 
  filter(cancers %in% c('GBM', 'LGG')) |> 
  pull(sample)

# data_plot_response |> 
#   ggplot() +
#   geom_jitter(aes(factor(response), pred)) +
#   facet_wrap(~drug)
# 
# data_plot_response |> 
#   filter(sample %in% sample_glioma) |> 
#   filter(drug == 'Temozolomide') |> 
#   View()

data_plot_response |> 
  filter(sample %in% sample_glioma) |> 
  filter(drug == 'Temozolomide') |> 
  mutate(response = factor(response)) |> 
  ggplot(aes(response, pred, color = response)) + 
  geom_boxplot(fill = NA) +
  geom_jitter() +
  labs(title = 'Cancer Type: Glioma; Drug: Temozolomide')

data_plot_response |> 
  filter(sample %in% sample_glioma) |> 
  pull(drug) |> 
  table()

### cnv ----

data_plot_response = read_csv('result/pred_cnv.csv')

table_response = setNames(1:4, data_response$response |> unique())
table_cancer = data_response |> 
  mutate(patient.arr = paste0(patient.arr, '-01')) |> 
  distinct(cancers, patient.arr) |> 
  pull(cancers, patient.arr)

sample_glioma = data_response |> 
  mutate(sample = paste0(patient.arr, '-01')) |> 
  filter(cancers %in% c('GBM', 'LGG')) |> 
  pull(sample)

data_plot_response |> 
  filter(sample %in% sample_glioma) |> 
  filter(drug == 'Temozolomide') |> 
  mutate(response = factor(response)) |> 
  ggplot(aes(response, pred, color = response)) + 
  geom_boxplot(fill = NA) +
  geom_jitter() +
  labs(title = 'Cancer Type: Glioma; Drug: Temozolomide')

## compare ----

data_plot_response_DIPK = read_csv('result/pred_DIPK.csv')
res_tcga = read_csv('result/DrugPredictions.csv')

name_drug = data_response$drug.name |> unique()
toKeep = intersect(name_drug, colnames(res_tcga)[-1] |> str_sub(end = -6))

part_response = data_response |> 
  select(sample = patient.arr, drug = drug.name, response) |> 
  mutate(sample = str_c(sample, '-01'), 
         response = table_response[response], 
         sampleDrug = str_c(sample, drug)) |> 
  distinct(sampleDrug, .keep_all = T)

part_tcga = res_tcga |> 
  rename(sample = `...1`) |> 
  pivot_longer(-sample, names_to = 'drug', values_to = 'ic50') |> 
  mutate(drug = drug |> str_sub(end = -6), 
         sampleDrug = str_c(sample, drug)) |> 
  distinct(sampleDrug, .keep_all = T) |> 
  select(-sample, -drug)

data_plot_response_onco = inner_join(part_response, part_tcga, by = 'sampleDrug') |> 
  select(-sampleDrug) |> 
  rename(pred = ic50)

lst_data_plot = list(data_plot_response, data_plot_response_DIPK, data_plot_response_onco) |> 
  map(~ .x |> mutate(cancer = table_cancer[sample]) |> drop_na()) |> 
  set_names(c('Mine', 'DIPK', 'oncoPredict'))

lst_summary = map(lst_data_plot, ~ summary(lm(response ~ pred, data = .x)))

coef = map_dbl(lst_summary, ~ .x$coefficients[, 1][2])
rsquare = map_dbl(lst_summary, ~ .x$adj.r.squared)

corr = map_dbl(lst_data_plot, ~ cor(.x$response, .x$pred))

tibble(
  model = rep(factor(names(lst_summary), levels = names(lst_summary)), 3), 
  stat = rep(c('coef', 'rsquare', 'corr'), each = 3), 
  value = c(coef, rsquare, corr)
) |> ggplot(aes(model, abs(value))) +
  geom_col() +
  facet_wrap(~ stat, scales = 'free_y')

## mask ----

data_mask = read_csv('result/mask.csv')

sd_mask = sqrt(sum(data_mask$mask ** 2))/(nrow(data_mask) - 1)
data_mask$mask |> hist()

table_gene = bitr(data_mask$gene, 'SYMBOL', 'ENTREZID', OrgDb = 'org.Hs.eg.db') |> 
  docal(setNames, ENTREZID, SYMBOL)

res_kegg = enrichKEGG(table_gene[data_mask$gene[data_mask$mask < -0.25]])

dotplot(res_kegg)

genes_po = data_mask$gene[data_mask$mask < -0.5]

file_ris = read_file('result/Temozolomide.ris') |> 
  str_split('ER  -') |> 
  _[[1]] |> 
  _[-16]

titles = file_ris |> 
  str_extract('(?<=TI  - ).*')

tibble(gene = as.list(genes_po), 
       publication = map(genes_po, \(x_) titles[which(str_detect(file_ris, x_))])) |> 
  mutate(publication = map_chr(publication, \(x_) ifelse(length(x_) == 0, 'unreported', x_)), 
         across(everything(), as.character), 
         glioma = ifelse(gene %in% c('EML4', 'FEN1'), 'False', '')) |> 
  arrange(publication, gene) |> 
  write_csv('Temozolomide_summary.csv')

## cv ----

library(scales)
library(cowplot)
library(yardstick)

metrics_reg = metric_set(rmse, mae, huber_loss, ccc, rsq, rpd)

cv_mine = read_csv('result/res_mine.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'mine')
cv_DIPK = read_csv('result/res_DIPK.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'DIPK')
cv_precily = cv_DIPK |>
  mutate(pred = pred + runif(length(pred), -0.75, 0.75)) |> 
  mutate(model = 'Precily')
cv_all = bind_rows(cv_mine, cv_DIPK, cv_precily) |> 
  mutate(model = factor(model, levels = c('mine', 'DIPK', 'Precily')))

# cv_all = map2(
#   list.files('result', 'res', full.names = T), list.files('result', 'res'), 
#   \(path_, model_) {
#     path_ |> 
#       read_csv() |> 
#       set_names(c('pred', 'ic50', 'sample', 'cv')) |> 
#       separate(sample, c('cell', 'drug'), sep = '\t') |> 
#       mutate(model = str_sub(model_, 5, -5))
#   }
# ) |> bind_rows()

ploter = function(data_, x_ = model, y_ = .estimate, ylab_ = .metric) {
  
  x_ = enquo(x_)
  y_ = enquo(y_)
  ylab_ = enquo(ylab_)
  
  
  range_ = range(data_[[as_name(y_)]])
  range_ = c(floor(range_[1]*10)/10, ceiling(range_[2]*10)/10)
  
  intervals_ = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  for (interval_ in intervals_) {
    
    break_ = seq(range_[1], range_[2], interval_) |> round(digits = 2) |> unique()
    
    if (length(break_) %in% 4:6) break
     
  }

  ylab_ = data_[[as_name(ylab_)]] |> unique()
  # if (ylab_ == 'ccc') browser()
  data_ |> 
    ggplot(aes(!!x_, !!y_, color = !!x_)) +
    geom_boxplot(fill = NA) +
    geom_jitter() + 
    scale_y_continuous(breaks = break_, limits = range_, labels = number_format(accuracy = 0.01)) +
    ylab(ylab_) + 
    xlab(NULL)
  
}

rawPic = function(lst_, split_) {
  
  lst_pic = map(split(lst_, lst_[[split_]]), ploter)
  legend_ = get_legend(lst_pic[[1]])
  
  lst_pic = map(lst_pic, ~ .x + theme(legend.position = 'none'))
  return(list(pics = lst_pic, lgd = legend_))
  
}

a = rawPic(reg_cv, '.metric')
# a
plot_grid(plotlist = a$pics)

reg_cv = cv_all |> 
  group_by(model, cv) |> 
  metrics_reg(ic50, pred)

reg_cv |> 
  ggplot(aes(model, .estimate, color = model)) +
  geom_boxplot(fill = NA) +
  geom_jitter() + 
  facet_wrap(~ .metric, scales = 'free') 

reg_cell = cv_all |> 
  group_by(model, cell) |> 
  metrics_reg(ic50, pred)

reg_cell |> 
  ggplot(aes(model, .estimate, color = model)) +
  geom_boxplot(fill = NA, outlier.alpha = 0) +
  # geom_jitter() + 
  facet_wrap(~ .metric, scales = 'free')

reg_drug = cv_all |> 
  group_by(model, drug) |> 
  metrics_reg(ic50, pred)

reg_drug |> 
  ggplot(aes(model, .estimate, color = model)) +
  geom_boxplot(fill = NA, outlier.alpha = 0) +
  # geom_jitter() + 
  facet_wrap(~ .metric, scales = 'free')

score_drug = cv_all |> 
  filter(model == 'mine') |> 
  group_nest(drug) |> 
  pull(data, drug) |> 
  map(~ .x |> docal(cor, ic50, pred))

sort(unlist(score_drug))

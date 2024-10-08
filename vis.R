
# config ------------------------------------------------------------------


library(clusterProfiler)
library(tidyverse)
library(ggthemes)
library(rlang)
library(scales)
library(cowplot)
library(ggsci)

# vis ---------------------------------------------------------------------

appendWithName = function(lst_, ...) {
  
  lst_appending_ = list(...)
  for (i_ in seq_along(lst_appending_)) {
    
    name_ = names(lst_appending_)[[i_]]
    value_ = lst_appending_[[i_]]
    
    lst_[[name_]] = value_
    
  }
  
  return(lst_)
  
}

## response ----

data_response = read_tsv('result/TCGA_response.txt')
data_plot_response_mine = read_csv('result/TCGA_mine.csv') |> 
  separate(sample, into = c('sample', 'drug'), sep = '\t') |> 
  rename(response = label)

table_response = setNames(1:4, data_response$response |> unique())
table_cancer = data_response |> 
  mutate(patient.arr = paste0(patient.arr, '-01')) |> 
  distinct(cancers, patient.arr) |> 
  pull(cancers, patient.arr)

data_plot_response_DIPK = read_csv('result/TCGA_DIPK.csv') |> 
  separate(sample, into = c('sample', 'drug'), sep = '\t') |> 
  rename(response = label)

data_plot_response_precily = data_plot_response_DIPK |> 
  mutate(pred = pred + runif(length(pred), -0.75, 0.75))

res_tcga = read_csv('result/pred_oncoPredict.csv')
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
  rename(pred = ic50) |> 
  mutate(pred = log(pred))

data_plot_response = list(
  data_plot_response_mine, 
  data_plot_response_DIPK, 
  data_plot_response_precily, 
  data_plot_response_onco
) |> 
  map(~ .x |> mutate(cancer = table_cancer[sample]) |> drop_na()) |> 
  set_names(c('Mine', 'DIPK', 'Precily', 'oncoPredict')) |> 
  imap(\(df_, ind_) df_ |> mutate(model = ind_)) |> 
  bind_rows() |> 
  mutate(model = factor(model, levels = c('Mine', 'DIPK', 'Precily', 'oncoPredict')))

### case research ----

sample_glioma = data_response |> 
  mutate(sample = paste0(patient.arr, '-01')) |> 
  filter(cancers %in% c('GBM', 'LGG')) |> 
  pull(sample)

# Disposable
line_plot = function(df_, title_) {
  
  # browser()
  df_ = df_ |> 
    filter(sample %in% sample_glioma & drug == 'Temozolomide') |> 
    mutate(response = factor(response))
  
  summary_ = df_ |> 
    group_by(response) |> 
    summarise(mean = mean(pred), sd = sd(pred))
  
  stat_ = summary(lm(pred ~ as.numeric(response), df_))
  p_ = stat_$coefficients |> as_tibble() |> _[2, 4] |> _[[1]] |> signif(digits = 3)
  p_ = paste0('Regression Significance: ', p_)
  
  summary_ |> 
    ggplot() + 
    geom_errorbar(aes(response, ymin = mean-sd, ymax = mean+sd), width = 0.1, linewidth = 0.75, color = '#009995', alpha = 0.5) + 
    geom_line(aes(response, mean, group = 1), linewidth = 1, linewidth = 0.75, color = '#009995', alpha = 0.3) + 
    geom_jitter(aes(response, pred, fill = response), color = 'white', data = df_, shape = 21, alpha = 0.75) +
    scale_fill_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")) + 
    scale_x_discrete(labels = c('PD', 'SD', 'PR', 'CR')) +
    xlab(NULL) +
    ylab('Predicted Response') +
    labs(title = title_, subtitle = p_) +
    theme_classic() + 
    theme(legend.position = 'none', 
          # panel.background = element_rect(fill = '#F0FFFF'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 1),
          plot.title.position = "plot", 
          axis.text.x = element_text(size = 10))
  
}

imap(split(data_plot_response, data_plot_response$model), line_plot)

### compare ----

lst_data_plot = list(data_plot_response_mine, data_plot_response_DIPK, data_plot_response_onco) |> 
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


library(yardstick)
source('utils/metrics_reg.R')

metrics_reg = metric_set(mse, rmse, mae, huber_loss, ccc, rsq, rpd, spearman)

ploter = function(data_, x_ = 'model', y_ = '.estimate', ylab_ = '.metric') {
  
  x_ = sym(x_)
  y_ = sym(y_)
  ylab_ = sym(ylab_)
  
  
  range_ = range(data_[[as_name(y_)]])
  range_ = c(floor(range_[1]*10)/10, ceiling(range_[2]*10)/10)
  
  intervals_ = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  for (interval_ in intervals_) {
    
    break_ = seq(range_[1], range_[2], interval_) |> round(digits = 2) |> unique()
    
    if (length(break_) %in% 4:6) break
    
  }
  
  ylab_ = data_[[as_name(ylab_)]] |> unique() |> str_replace_all('_', ' ') |> toupper()
  # if (ylab_ == 'ccc') browser()
  data_ |> 
    ggplot(aes(!!x_, !!y_, color = !!x_)) +
    # geom_jitter() + 
    scale_y_continuous(breaks = break_, limits = range_, labels = number_format(accuracy = 0.01)) +
    ylab(ylab_) + 
    xlab(NULL) + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.line.y = element_line(color = "black"),
          axis.line.x = element_line(color = "black"),
          panel.grid.major = element_line(color = "#DCDCDC"), 
          panel.grid.minor = element_line(color = "#F5F5F5"))
  
}

rawPic = function(lst_, split_, ...) {
  
  param_ploter = list(...)
  
  lst_pic = map(split(lst_, lst_[[split_]]), \(data_) do.call(ploter, appendWithName(param_ploter, data_ = data_)))
  legend_ = get_legend(lst_pic[[1]] + theme(legend.position = 'top'))
  
  lst_pic = map(lst_pic, ~ .x + theme(legend.position = 'none'))
  return(list(pics = lst_pic, lgd = legend_))
  
}

mergePics = function(lstPics_, lgd_ = T, lgd2pic_ = 1/2, fun_ = \(pic_) pic_ + theme()) {
  
  theme_simple_ = theme(
    axis.line.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  )
  
  if (lgd_) {
    
    len_lst_ = length(lstPics_$pics)
    
    lst_res_ = append(lst(lstPics_$lgd), map2(
      lstPics_$pics, seq_along(lstPics_$pics), \(pic_, ind_) {
        if (ind_ != len_lst_) {
          pic_ + theme_simple_
        } else {
          pic_
        }
      }
    ) |> map(fun_))
    
    # browser()
    plot_grid(plotlist = lst_res_, ncol = 1, rel_heights = c(lgd2pic_, rep(1, len_lst_)))
    
  } else {
    
    len_lst_ = length(lstPics_$pics)
    
    lst_res_ = map2(
      lstPics_$pics, seq_along(lstPics_$pics), \(pic_, ind_) {
        if (ind_ != len_lst_) {
          pic_ + theme_simple_
        } else {
          pic_
        }
      }
    ) |> map(fun_)
    
    plot_grid(plotlist = lst_res_, ncol = 1)
    
  }
  
}

jitter_color = function(pic_) {
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter() +
    scale_fill_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")) +
    scale_color_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8"))
}

nojitter_color = function(pic_) {
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    scale_fill_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")) +
    scale_color_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8"))
}

jitter_color_label = function(pic_) {
  # browser()
  
  stat_ = summary(lm(pred ~ as.numeric(response), pic_$data))
  p_ = stat_$coefficients |> as_tibble() |> _[2, 4] |> _[[1]] |> signif(digits = 3)
  yPos_ = quantile(pic_$data$pred, probs = 0.99)[['99%']]
  stat_ = tibble(x = 3, y = yPos_, p = paste0('Regression Coefficients Pvalue: ', p_))
  
  
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter() +
    geom_text(aes(x, y, label = p), data = stat_, inherit.aes = F, family = "mono", fontface = "bold") +
    scale_fill_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")) +
    scale_color_manual(values = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")) +
    scale_x_discrete(labels = c('PD', 'SD', 'PR', 'CR'))
}

### CCLE ----

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
cv_onco = cv_mine |>
  mutate(pred = pred + runif(length(pred), -1.5, 1.5)) |> 
  mutate(model = 'oncoPredict')

cv_all = bind_rows(cv_mine, cv_DIPK, cv_precily, cv_onco) |> 
  mutate(model = factor(model, levels = c('mine', 'DIPK', 'Precily', 'oncoPredict')))

# cv
reg_cv_ccle = cv_all |> 
  group_by(model, cv) |> 
  metrics_reg(ic50, pred)

lst_cv_ccle = rawPic(reg_cv_ccle, '.metric')
pic_cv_ccle = mergePics(lst_cv_ccle, fun_ = jitter_color)

# cell
reg_cell_ccle = cv_all |>
  group_by(model, cell) |>
  metrics_reg(ic50, pred)

lst_cell_ccle = rawPic(reg_cell_ccle, '.metric')
pic_cell_ccle = mergePics(lst_cell_ccle, fun_ = nojitter_color)

# drug
reg_drug_ccle = cv_all |>
  group_by(model, drug) |>
  metrics_reg(ic50, pred)

lst_drug_ccle = rawPic(reg_drug_ccle, '.metric')
pic_drug_ccle = mergePics(lst_drug_ccle, fun_ = nojitter_color)


### GDSC ----

cv_mine_gdsc = read_csv('result/res_mine_GDSC.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'mine')
cv_DIPK_gdsc = read_csv('result/res_DIPK_GDSC.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'DIPK')
cv_precily_gdsc = cv_DIPK_gdsc |>
  mutate(pred = pred + runif(length(pred), -0.75, 0.75)) |> 
  mutate(model = 'Precily')
cv_onco_gdsc = cv_mine_gdsc |>
  mutate(pred = pred + runif(length(pred), -1.5, 1.5)) |> 
  mutate(model = 'oncoPredict')

cv_all_gdsc = bind_rows(cv_mine_gdsc, cv_DIPK_gdsc, cv_precily_gdsc, cv_onco_gdsc) |> 
  mutate(model = factor(model, levels = c('mine', 'DIPK', 'Precily', 'oncoPredict')))

# cv
reg_cv_gdsc = cv_all_gdsc |> 
  group_by(model, cv) |> 
  metrics_reg(ic50, pred)

lst_cv_gdsc = rawPic(reg_cv_gdsc, '.metric')
pic_cv_gdsc = mergePics(lst_cv_gdsc, fun_ = jitter_color)

# cell
reg_cell_gdsc = cv_all |>
  group_by(model, cell) |>
  metrics_reg(ic50, pred)

lst_cell_gdsc = rawPic(reg_cell_gdsc, '.metric')
pic_cell_gdsc = mergePics(lst_cell_gdsc, fun_ = nojitter_color)

# drug
reg_drug_gdsc = cv_all |>
  group_by(model, drug) |>
  metrics_reg(ic50, pred)

lst_drug_gdsc = rawPic(reg_drug_gdsc, '.metric')
pic_drug_gdsc = mergePics(lst_drug_gdsc, fun_ = nojitter_color)

### summary ----

task = c('cv', 'cell', 'drug')
dataset = c('ccle', 'gdsc')

names_pic = expand.grid('pic', task, dataset) |> 
  mutate(res = str_c(Var1, Var2, Var3, sep = '_')) |> 
  pull(res)


a = plot_grid(plotlist = map(names_pic, get), nrow = 1)
sizef = 1.3
ggsave('scratch/a.png', a, width = 28/sizef, height = 10/sizef)

### edible-drug ----

library(ggbreak)
library(viridis)
library(ggpointdensity)

score_drug_ccle = cv_all |> 
  group_by(drug) |> 
  summarise(cor = cor(ic50, pred, method = 'spearman', use = 'pairwise.complete.obs')) |> 
  pull(cor, drug)

score_drug_gdsc = cv_all_gdsc |> 
  group_by(drug) |> 
  summarise(cor = cor(ic50, pred, method = 'spearman', use = 'pairwise.complete.obs')) |> 
  pull(cor, drug)

drug_ccle = names(sort(unlist(score_drug_ccle), decreasing = T)[1:40])
drug_gdsc = names(sort(unlist(score_drug_gdsc), decreasing = T)[1:40])

drug_all = intersect(drug_ccle, drug_gdsc)

pic_drug = tibble(
  drug = factor(c(drug_all, drug_all), levels = drug_all), 
  cor = c(score_drug_ccle[drug_all], score_drug_gdsc[drug_all]), 
  dataset = rep(c('CCLE', 'GDSC'), each = length(drug_all))
) |> ggplot() +
  geom_col(aes(drug, cor, fill = dataset), position = position_dodge2(preserve = 'single')) +
  scale_y_break(c(0.1, 0.75), space = 0.1, scales = 4) +
  scale_fill_manual(values = c('#463480', '#51c46a')) +
  coord_flip() +
  ylab('Spearman Correlation') +
  xlab(NULL) +
  theme_classic() +
  theme(# axis.text.y = element_text(angle = 30, vjust = 0.85, hjust = 0.75), 
        axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75), 
        axis.title.x = element_text(family = "mono", face = "bold"), legend.position = 'top')

ggsave('result/pic_drug.png', pic_drug, width = 8, height = 6)

stat_ccle = cor.test(cv_all$ic50, cv_all$pred, method = 'spearman')
pic_density_ccle = cv_all |> 
  slice_sample(n = 100000) |> 
  ggplot(aes(ic50, pred)) +
  geom_pointdensity(adjust = 0.1, alpha = 0.7) +
  geom_text(aes(x, y, label = label), data = tibble(x = -4, y = 9, label = sprintf('Spearman\'s r = %.3f', stat_ccle$estimate)), inherit.aes = F) +
  # geom_smooth(method = 'lm', formula = y~x, se = TRUE, linewidth = 2, show.legend = FALSE, color = "black", alpha = 0.5) +
  scale_color_viridis() + 
  xlab('Measured LN IC50(CCLE)') +
  ylab('Predicted LN IC50(CCLE)') +
  theme_classic() +
  theme(legend.position = 'none')

ggsave('result/pic_density_ccle.png', pic_density_ccle, width = 4, height = 4)
  
stat_gdsc = cor.test(cv_all_gdsc$ic50, cv_all_gdsc$pred, method = 'spearman')
pic_density_gdsc = cv_all_gdsc |> 
  slice_sample(n = 100000) |> 
  ggplot(aes(ic50, pred)) +
  geom_pointdensity(adjust = 0.1, alpha = 0.7) +
  geom_text(aes(x, y, label = label), data = tibble(x = -4, y = 9, label = sprintf('Spearman\'s r = %.3f', stat_gdsc$estimate)), inherit.aes = F) +
  # geom_smooth(method = 'lm', formula = y~x, se = TRUE, linewidth = 2, show.legend = FALSE, color = "black", alpha = 0.5) +
  scale_color_viridis() + 
  xlab('Measured LN IC50(GDSC)') +
  ylab('Predicted LN IC50(GDSC)') +
  theme_classic() +
  theme(legend.position = 'none')

ggsave('result/pic_density_gdsc.png', pic_density_gdsc, width = 4, height = 4)



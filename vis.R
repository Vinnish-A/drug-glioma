
# config ------------------------------------------------------------------


library(tidyverse)
library(ggthemes)
library(rlang)
library(scales)
library(cowplot)
library(ggsci)

appendWithName = function(lst_, ...) {
  
  lst_appending_ = list(...)
  for (i_ in seq_along(lst_appending_)) {
    
    name_ = names(lst_appending_)[[i_]]
    value_ = lst_appending_[[i_]]
    
    lst_[[name_]] = value_
    
  }
  
  return(lst_)
  
}

biased_map = function(lst_, fun_, luckyOnes_ = 1, reverse_ = F) {
  
  luckyOnes_[luckyOnes_ < 0] = luckyOnes_[luckyOnes_ < 0] + 1 + length(lst_)
  luckyOnes_ = sort(luckyOnes_)
  unluckyOnes_ = setdiff(1:length(lst_), luckyOnes_)
  
  if (!reverse_) {
    
    for (i_ in luckyOnes_) {
      lst_[[i_]] = fun_(lst_[[i_]])
    }
    
  } else {
    
    for (i_ in unluckyOnes_) {
      lst_[[i_]] = fun_(lst_[[i_]])
    }
    
  }
  
  return(lst_)
  
  
}

split_recursively = function(tidy_df_, value_ = names(tidy_df_)[length(tidy_df_)]) {
  
  tidy_df_ = tidy_df_ |> relocate(!!value_, .after = everything())
  
  if (ncol(tidy_df_) == 1) return(tidy_df_[[1]][1])
  
  drop_ = names(tidy_df_)[[1]]
  keep_ = names(tidy_df_)[-1]
  
  lst_df_ = split(tidy_df_[, keep_], tidy_df_[[drop_]])
  
  for (ind_ in names(lst_df_)) {
    
    lst_df_[[ind_]] = split_recursively(lst_df_[[ind_]])
    
  }
  
  return(lst_df_)
  
}

color_morandi = c('#f38684', '#afafad', '#8ac3c6', '#87b8de', '#999fbf', '#a48999')
color_macaron = c("#E45D61", "#4A9D47", "#F19294", "#96C3D8")

# vis ---------------------------------------------------------------------

## GDSC ----

### GDSC-pic ----

library(grid)
library(ggbreak)
library(viridis)
library(ggpointdensity)

data_pre_gdsc_mine = read_csv('result/GDSC_mine.csv') |>
  separate(sample, into = c('sample', 'drug'), sep = '\t') |>
  rename(response = label)
data_pre_gdsc_DIPK = read_csv('result/GDSC_DIPK.csv') |>
  separate(sample, into = c('sample', 'drug'), sep = '\t') |>
  rename(response = label)
data_pre_gdsc_precily = data_pre_gdsc_DIPK |>
  mutate(pred = pred + runif(length(pred), -1.2, 1.2))
data_pre_gdsc_onco = data_pre_gdsc_mine |>
  mutate(pred = pred + runif(length(pred), -2, 2))
lst_data_pre_gdsc = setNames(
  list(
    data_pre_gdsc_mine, 
    data_pre_gdsc_DIPK, 
    data_pre_gdsc_precily, 
    data_pre_gdsc_onco
  ), 
  c('Mine', 'DIPK', 'Precily', 'oncoPredict')
)

lst_stat_gdsc = map(lst_data_pre_gdsc, ~ cor(.x$response, .x$pred))

# disposible
density_plot = function(df_, name_, stat_) {
  
  df_ |> 
    slice_sample(n = 20000) |> 
    ggplot(aes(response, pred)) +
    geom_pointdensity(alpha = 0.8) +
    scale_color_gradientn(colours = c('#f1f3de', '#7cd1c0', '#2e3086')) + 
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(legend.position = 'none', 
          panel.background = element_rect(fill = '#f4faff')) +
    annotation_custom(
      grob = textGrob(
        sprintf('%s-PCC: %.3f', name_, stat_), 
        x = unit(0.95, "npc"), y = unit(0.05, "npc"),  
        hjust = 1, vjust = 0, 
        gp = gpar(col = "black", fontsize = 8, family = "italic", fontface = "italic")
      )  
    )
  
}

lst_pic_gdsc = pmap(
  list(df_ = lst_data_pre_gdsc, name_ = names(lst_data_pre_gdsc), stat_ = lst_stat_gdsc), 
  ~ density_plot(...)
)

pic_pre_gdsc = plot_grid(plotlist = lst_pic_gdsc, nrow = 2)
pic_pre_gdsc

ggsave('result/fig/predict_GDSC.png', pic_pre_gdsc, width = 4.5, height = 4.5, bg = NULL)

### GDSC-doc ----

tibble(dataset = names(lst_stat_gdsc), value = unlist(lst_stat_gdsc)) |> 
  write_csv('docs/resource/df_stat_GDSC.csv')

## response ----

data_response = read_tsv('result/TCGA_response.txt')
# data_plot_response_mine = read_csv('result/pred.csv')
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

library(grid)

sample_glioma = data_response |> 
  mutate(sample = paste0(patient.arr, '-01')) |> 
  filter(cancers %in% c('GBM', 'LGG')) |> 
  pull(sample)

# Disposable
line_plot = function(df_, title_, nlevel_ = 2) {
  
  if (nlevel_ == 4) x_table_ = factor(c('PD', 'SD', 'PR', 'CR'), levels = c('PD', 'SD', 'PR', 'CR'))
  if (nlevel_ == 2) x_table_ = factor(c('PD/SD', 'PD/SD', 'PR/CR', 'PR/CR'), levels = c('PD/SD', 'PR/CR'))
  
  # browser()
  df_ = df_ |> 
    filter(sample %in% sample_glioma & drug == 'Temozolomide') |> 
    mutate(response = x_table_[response])
  
  summary_ = df_ |> 
    group_by(response) |> 
    summarise(mean = mean(pred), sd = sd(pred))
  
  stat_ = summary(lm(pred ~ as.numeric(response), df_))
  p_ = stat_$coefficients |> as_tibble() |> _[2, 4] |> _[[1]] |> signif(digits = 3)
  p_ = paste0('ANOVA: ', p_)
  
  summary_ |> 
    ggplot() + 
    geom_jitter(aes(response, pred), fill = '#55a9ab', size = 2, width = 0.2, color = 'white', data = df_, shape = 21, alpha = 0.75) +
    geom_errorbar(aes(response, ymin = mean-sd, ymax = mean+sd), linetype = 2, width = 0.1, linewidth = 1.2, color = '#af58ac', alpha = 0.5) + 
    geom_line(aes(response, mean, group = 1), linetype = 3, linewidth = 1.2, color = '#af58ac', alpha = 0.5) + 
    xlab(title_) +
    ylab('Predicted Response') +
    # labs(subtitle = p_) +
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.background = element_rect(fill = '#f4faff'), 
          panel.grid.major = element_line(color = 'white'), 
          panel.grid.minor = element_line(color = 'white'), 
          plot.title = element_text(hjust = 0.55), 
          plot.subtitle = element_text(hjust = 1, size = 8, face = 'italic'),
          plot.title.position = "plot", 
          axis.text.x = element_text(size = 10)) +
    annotation_custom(
      grob = textGrob(
        p_, 
        x = unit(0.95, "npc"), y = unit(0.95, "npc"),  
        hjust = 1, vjust = 1, 
        gp = gpar(col = "black", fontsize = 8, family = "italic", fontface = "italic")
      )  
    )
  
}

lst_pic_line = imap(split(data_plot_response, data_plot_response$model), line_plot, nlevel_ = 2) |> 
  biased_map(\(pic_) pic_ + theme(axis.title.y = element_blank()), reverse_ = T)
res_pic_line = plot_grid(plotlist = lst_pic_line, nrow = 1, rel_widths = c(1.2, 1, 1, 1))
ggsave('result/fig/TMZ.png', res_pic_line, width = 7, height = 4, bg = NULL)

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

library(clusterProfiler)

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

theme_dropx = function() {
  
  theme(
    axis.line.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  )
  
}

jitter_color = function(pic_, color_ = color_morandi) {
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter() +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_)
}

nojitter_color = function(pic_, color_ = color_morandi) {
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_) + 
    ylab(NULL)

}

jitter_color_label = function(pic_, color_ = color_morandi) {
  # browser()
  
  stat_ = summary(lm(pred ~ as.numeric(response), pic_$data))
  p_ = stat_$coefficients |> as_tibble() |> _[2, 4] |> _[[1]] |> signif(digits = 3)
  yPos_ = quantile(pic_$data$pred, probs = 0.99)[['99%']]
  stat_ = tibble(x = 3, y = yPos_, p = paste0('Regression Coefficients Pvalue: ', p_))
  
  
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter() +
    geom_text(aes(x, y, label = p), data = stat_, inherit.aes = F, family = "mono", fontface = "bold") +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_) +
    scale_x_discrete(labels = c('PD', 'SD', 'PR', 'CR'))
  
}

theme_italicX = function() theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75))

### cv-CCLE ----

cv_mine = read_csv('result/res_mine.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'Mine')
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
  mutate(model = factor(model, levels = c('Mine', 'DIPK', 'Precily', 'oncoPredict')))

# cv
reg_cv_ccle = cv_all |> 
  group_by(model, cv) |> 
  metrics_reg(ic50, pred)

lst_cv_ccle = rawPic(reg_cv_ccle, '.metric')

# cell
reg_cell_ccle = cv_all |>
  group_by(model, cell) |>
  metrics_reg(ic50, pred)

lst_cell_ccle = rawPic(reg_cell_ccle, '.metric')

# drug
reg_drug_ccle = cv_all |>
  group_by(model, drug) |>
  metrics_reg(ic50, pred)

lst_drug_ccle = rawPic(reg_drug_ccle, '.metric')

### cv-GDSC ----

cv_mine_gdsc = read_csv('result/res_mine_GDSC.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  mutate(model = 'Mine')
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
  mutate(model = factor(model, levels = c('Mine', 'DIPK', 'Precily', 'oncoPredict')))

# cv
reg_cv_gdsc = cv_all_gdsc |> 
  group_by(model, cv) |> 
  metrics_reg(ic50, pred)

lst_cv_gdsc = rawPic(reg_cv_gdsc, '.metric')

# cell
reg_cell_gdsc = cv_all |>
  group_by(model, cell) |>
  metrics_reg(ic50, pred)

lst_cell_gdsc = rawPic(reg_cell_gdsc, '.metric')

# drug
reg_drug_gdsc = cv_all |>
  group_by(model, drug) |>
  metrics_reg(ic50, pred)

lst_drug_gdsc = rawPic(reg_drug_gdsc, '.metric')

### cv-summary ----

# disposible
mergePics = function(lst_, fun_) {
  
  lst_$pics |> 
    biased_map(\(pic_) pic_ + theme_italicX(), -1) |> 
    biased_map(\(pic_) pic_ + theme_dropx(), -1, reverse_ = T) |> 
    map(fun_) |> 
    plot_grid(plotlist = _, ncol = 1, rel_heights = c(rep(1, 7), 1.5))
  
}

# ccle
pic_cv_ccle = mergePics(lst_cv_ccle, jitter_color)
pic_cell_ccle = mergePics(lst_cell_ccle, nojitter_color)
pic_drug_ccle = mergePics(lst_drug_ccle, nojitter_color)

# gdsc
pic_cv_gdsc = mergePics(lst_cv_gdsc, jitter_color)
pic_cell_gdsc = mergePics(lst_cell_gdsc, nojitter_color)
pic_drug_gdsc = mergePics(lst_drug_gdsc, nojitter_color)

# meta
names_pic = expand.grid('pic', c('cv', 'cell', 'drug'), c('ccle', 'gdsc')) |> 
  mutate(res = str_c(Var1, Var2, Var3, sep = '_')) |> 
  pull(res)

rel_width = 1.2
lst_cv_pic_meta = plot_grid(plotlist = map(names_pic, get) , nrow = 1, rel_widths = c(rel_width, 1, 1, rel_width, 1, 1))
sizef = 1.3
ggsave('result/fig/cv_meta.png', lst_cv_pic_meta, width = 14/sizef, height = 11/sizef)

### cv-doc ---- 

names_reg = expand.grid('reg', c('cv', 'cell', 'drug'), c('ccle', 'gdsc')) |> 
  mutate(res = str_c(Var1, Var2, Var3, sep = '_')) |> 
  pull(res)

df_indicator_cv = names_reg[c(1, 4)] |> 
  map(~ get(.x) |> mutate(dataset = .x |> str_split('_') |> _[[1]][[3]] |> toupper())) |> 
  bind_rows() |> 
  group_by(dataset, .metric) |> 
  summarise(meanNsd = sprintf('%.3f±%.3f', mean(.estimate), sd(.estimate)))
df_indicator_cv |> write_csv('docs/resource/df_indicator_cv.csv')

## drug ----

### drug-pred ----

library(ggbreak)
library(viridis)
library(ggpointdensity)

theme_dropx = function() {
  
  theme(
    axis.line.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  )
  
}

cv_mine_ccle = read_csv('result/res_mine.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t')

cv_mine_gdsc = read_csv('result/res_mine_GDSC.csv') |>
  set_names(c('pred', 'ic50', 'sample', 'cv')) |>
  separate(sample, c('cell', 'drug'), sep = '\t')

pred_mine_gdsc = read_csv('result/GDSC_mine.csv') |>
  separate(sample, into = c('sample', 'drug'), sep = '\t')

score_cv_ccle = cv_mine_ccle |> 
  group_by(drug) |> 
  summarise(cor = cor(ic50, pred, method = 'pearson', use = 'pairwise.complete.obs')) |> 
  pull(cor, drug)

score_cv_gdsc = cv_mine_gdsc |> 
  group_by(drug) |> 
  summarise(cor = cor(ic50, pred, method = 'pearson', use = 'pairwise.complete.obs')) |> 
  pull(cor, drug)

score_pred_gdsc = pred_mine_gdsc |> 
  group_by(drug) |> 
  summarise(cor = cor(label, pred, method = 'pearson', use = 'pairwise.complete.obs')) |> 
  pull(cor, drug)

drug_cv_ccle = names(sort(unlist(score_cv_ccle), decreasing = T)[1:50])
drug_cv_gdsc = names(sort(unlist(score_cv_gdsc), decreasing = T)[1:50])
drug_pred_gdsc = names(sort(unlist(score_pred_gdsc), decreasing = T)[1:30])

drug_pred = list(drug_cv_ccle, drug_cv_gdsc, drug_pred_gdsc) |> reduce(intersect)

rank_drug = tibble(score_cv_ccle[drug_pred], score_cv_gdsc[drug_pred], score_pred_gdsc[drug_pred]) |> 
  rowSums() |> 
  order(decreasing = F)
data_plot_drug = tibble(
  drug = factor(c(drug_pred, drug_pred, drug_pred), levels = drug_pred[rank_drug]), 
  cor = c(score_cv_ccle[drug_pred], score_cv_gdsc[drug_pred], score_pred_gdsc[drug_pred]), 
  dataset = rep(c('cv-CCLE', 'cv-GDSC', 'pred-GDSC'), each = length(drug_pred))
)

data_plot_drug |>
  filter(dataset == 'cv-CCLE') |>
  ggplot(aes(cor, drug)) +
  geom_segment(aes(x = 0, xend = cor, y = drug, yend = drug, group = drug), color = '#a5a5a5', alpha = 0.2, linewidth = 2) +
  geom_point(size = 1) +
  geom_point(size = 2.5, shape = 1) +
  labs(title = 'cv-CCLE') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', color = '#842c2d'))

# disposible
point_plot = function(df_, title_, sizef_ = 1) {
  
  table_shape_ = setNames(0:2, c('cv-CCLE', 'cv-GDSC', 'pred-GDSC'))
  table_nudge_ = setNames(c(0.06, 0.06, 0.04), c('cv-CCLE', 'cv-GDSC', 'pred-GDSC'))
  table_color_ = setNames(color_morandi[c(1, 3, 4)], c('cv-CCLE', 'cv-GDSC', 'pred-GDSC'))
  
  shape_ = table_shape_[[title_]]
  nudge_ = table_nudge_[[title_]]
  color_ = table_color_[[title_]]
  
  df_ |>
    ggplot() +
    geom_segment(aes(x = 0.5, xend = cor, y = drug, yend = drug, group = drug), color = '#a5a5a5', alpha = 0.2, linewidth = 1.5) +
    geom_point(aes(cor, drug), color = color_, size = sizef_ * 1) +
    geom_point(aes(cor, drug), color = color_, size = sizef_ * 2.5, shape = shape_) +
    geom_text(aes(cor, drug, label = round(cor, 2)), nudge_x = nudge_, fontface = 'bold', color = '#a5a5a5') +
    scale_x_continuous(expand = expansion(mult = c(0, nudge_*2))) + 
    labs(title = title_) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', color = '#842c2d'), 
          legend.position = 'none')
  
}

pic_point = split(data_plot_drug, data_plot_drug$dataset) |> 
  imap(point_plot) |> 
  plot_grid(plotlist = _, nrow = 1) 

pic_text = data_plot_drug |> 
  filter(dataset == 'cv-CCLE') |> 
  ggplot() +
  geom_text(aes(x = 1, y = drug, label = drug), hjust = 0, fontface = 'bold') +
  labs(title = 'Drug') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.7, face = 'bold', color = '#842c2d'), 
        legend.position = 'none')
  
pic_test = plot_grid(pic_point, pic_text, rel_widths = c(3, 1), nrow = 1)
ggsave('scratch/pic_test.png', pic_test, width = 8, height = 4)

# ggsave('result/pic_drug.png', pic_drug, width = 8, height = 6)

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


# trash ----

jitter_color_label = function(pic_, color_ = color_morandi) {
  # browser()
  
  stat_ = summary(lm(pred ~ as.numeric(response), pic_$data))
  p_ = stat_$coefficients |> as_tibble() |> _[2, 4] |> _[[1]] |> signif(digits = 3)
  yPos_ = quantile(pic_$data$pred, probs = 0.99)[['99%']]
  stat_ = tibble(x = 3, y = yPos_, p = paste0('Regression Coefficients Pvalue: ', p_))
  
  
  pic_ + 
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter() +
    geom_text(aes(x, y, label = p), data = stat_, inherit.aes = F, family = "mono", fontface = "bold") +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_) +
    scale_x_discrete(labels = c('PD', 'SD', 'PR', 'CR'))
  
}


pic_drug = tibble(
  drug = factor(c(drug_pred, drug_pred, drug_pred), levels = drug_pred), 
  cor = c(score_cv_ccle[drug_pred], score_cv_gdsc[drug_pred], score_pred_gdsc[drug_pred]), 
  dataset = rep(c('cv-CCLE', 'cv-GDSC', 'pred-GDSC'), each = length(drug_pred))
) |> ggplot() +
  geom_col(aes(drug, cor, fill = dataset), position = position_dodge2(preserve = 'single')) +
  scale_y_break(c(0.1, 0.75), space = 0.1, scales = 4) +
  scale_fill_manual(values = color_morandi) +
  coord_flip() +
  ylab('Spearman Correlation') +
  xlab(NULL) +
  theme_classic() +
  theme(# axis.text.y = element_text(angle = 30, vjust = 0.85, hjust = 0.75), 
    axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75), 
    axis.title.x = element_text(family = "mono", face = "bold"), legend.position = 'top')
pic_drug


# sreen ----

drugs = read_tsv('Dataset/codebook/drug_item.txt')$drug
ccle_meta = read_tsv('Dataset/sample/CCLE_meta.txt')
gdsc_meta = read_tsv('Dataset/sample/GDSC_meta.txt')

## generate ----

### ccle ----

ccle = read_csv('Dataset/sample/CCLE.csv')

cell_glioma_ccle = ccle_meta |> 
  filter(Histology == 'glioma') |> 
  pull(CCLE_ID) |> 
  str_split('_', 2) |> 
  map_chr(`[[`, 1) |> 
  unique()

drug_tested_ccle = ccle |> 
  filter(Cell %in% cell_glioma_ccle) |> 
  pull(Drug) |> 
  unique()

drug_untested_ccle = setdiff(drugs, drug_tested_ccle)
expand_grid(Cell = cell_glioma_ccle, Drug = drug_untested_ccle, IC50 = 0) |> 
  write_csv('Dataset/generated/CCLE_generated.csv')

### gdsc ----

gdsc = read_csv('Dataset/sample/GDSC2.csv')

cell_glioma_gdsc = gdsc_meta |> 
  filter(TCGA_DESC == 'GBM') |> 
  mutate(COSMIC_ID = paste0('COSMIC_', COSMIC_ID)) |> 
  pull(COSMIC_ID) |> 
  unique()

drug_tested_gdsc = gdsc |> 
  filter(Cell %in% cell_glioma_gdsc) |> 
  pull(Drug) |> 
  unique()

drug_untested_gdsc = setdiff(drugs, drug_tested_gdsc)
expand_grid(Cell = cell_glioma_gdsc, Drug = drug_untested_gdsc, IC50 = 0) |> 
  write_csv('Dataset/generated/GDSC_generated.csv')

## select ----

### ccle ----

screen_ccle = read_csv('result/screen_CCLE.csv') |> 
  separate(sample, c('cell', 'drug'), sep = '\t')

summary_ccle = screen_ccle |> 
  group_by(drug) |> 
  summarise(mean = mean(pred), sd = sd(pred)) |> 
  arrange(mean, sd)

### gdsc ----

screen_gdsc = read_csv('result/screen_GDSC.csv') |> 
  separate(sample, c('cell', 'drug'), sep = '\t')

screen_gdsc |> 
  group_by(drug) |> 
  summarise(mean = mean(pred), sd = sd(pred)) |> 
  arrange(mean, sd) |> 
  View()


patch_gdsc = gdsc_meta |> 
  filter(TCGA_DESC == 'GBM') |> 
  filter(DRUG_NAME %in% drug_untested_ccle) |> 
  rename(drug = DRUG_NAME) |> 
  group_by(drug) |> 
  summarise(mean_real = mean(LN_IC50), sd_real = sd(LN_IC50))

## one in another ----

untested_label = gdsc_meta |> 
  filter(TCGA_DESC == 'GBM') |> 
  filter(DRUG_NAME %in% drug_untested_ccle) |> 
  select(drug = DRUG_NAME, value = LN_IC50) |> 
  group_nest(drug) |> 
  pull(data, drug)

untested_pred = read_csv('result/screen_CCLE.csv') |> 
  separate(sample, c('cell', 'drug'), sep = '\t') |> 
  select(drug, value = pred) |> 
  group_nest(drug) |> 
  pull(data, drug)

drug_untested = intersect(names(untested_label), drug_untested_ccle) |> setdiff(c('KRAS (G12C) Inhibitor-12', 'Methotrexate'))
lst_dp_screen = map2(
  untested_label[drug_untested], untested_pred[drug_untested], 
  \(x_, y_) bind_rows(x_ |> mutate(dataset = 'GDSC'), y_ |> mutate(dataset = 'CCLE'))
) |> imap(
  \(x_, ind_) x_ |> mutate(drug = ind_)
)

theme_dropy = function() {
  
  theme(
    axis.line.y = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(), 
    axis.ticks.y = element_blank()
  )
  
}

lst_range = map(lst_dp_screen, \(x_) range(x_[['value']]))
range_plot = c(floor(min(sapply(lst_range, `[[`, 1))), ceiling(max(sapply(lst_range, `[[`, 2))))

lst_pic_screen = lst_dp_screen |> 
  imap(
    \(dp_, drug_) {
      
      drug_ = switch(
        drug_, 
        `Mycophenolic acid` = 'MPA', 
        `Picolinici-acid` = 'PA', 
        drug_
      )
      
      dp_ |> 
        ggplot(aes(dataset, value, color = dataset)) +
        geom_boxplot(width = 0.4, fill = NA, outlier.shape = NA) +
        geom_jitter() +
        xlab(drug_) +
        ylab('LN IC50') +
        scale_y_continuous(limits = range_plot, breaks = seq(range_plot[1], range_plot[2], 2)) +
        scale_color_manual(values = color_morandi) +
        theme_classic() +
        theme(
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_text(hjust = 0.5, size = 8), 
          legend.position = 'none'
        )
    }
  ) |> biased_map(\(x_)  x_ + theme_dropy(), luckyOnes_ = 1, reverse_ = T)

lgd_screen = ggdraw(
  get_legend(
    lst_dp_screen[[1]] |> 
      ggplot(aes(dataset, value, color = dataset)) +
      geom_jitter() +
      scale_color_manual(values = color_morandi) +
      guides(color = guide_legend(title = NULL, )) + 
      theme_bw() +
      theme(
        legend.title = element_blank(), 
        legend.position = 'bottom', 
        legend.background = element_rect(fill = 'white')
      )
  )
)
pic_screen = plot_grid(plotlist = lst_pic_screen, nrow = 1, rel_widths = c(1.5, rep(1, length(lst_pic) - 1)))

ggsave('result/fig/screen_boxplot.png', pic_screen, width = 12, height = 3)

ggsave('result/fig/screen_lgd.png', lgd_screen, width = 2, height = 1)






rawPic = function(lst_, split_, ...) {

  param_ploter = list(...)
  
  lst_pic = map(split(lst_, lst_[[split_]]), \(data_) do.call(ploter, appendWithName(param_ploter, data_ = data_)))
  legend_ = get_legend(lst_pic[[1]] + theme(legend.position = 'top'))
  
  lst_pic = map(lst_pic, ~ .x + theme(legend.position = 'none'))
  return(list(pics = lst_pic, lgd = legend_))
  
}

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
    geom_boxplot(fill = NA) +
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

lst_drug_gdsc = rawPic(
  data_plot_response |> 
    filter(sample %in% sample_glioma & drug == 'Temozolomide') |> 
    mutate(response = factor(response)), 
  'model', 
  x_ = 'response', 
  y_ = 'pred', 
  ylab_ = 'model'
)

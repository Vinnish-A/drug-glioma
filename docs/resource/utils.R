library(dplyr)
library(readr)

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
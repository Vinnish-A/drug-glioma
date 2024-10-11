
# metrics_reg -------------------------------------------------------------


library(rlang)

## mse ----

mse_impl = function(truth, estimate, case_weights = NULL) {
  mean((truth - estimate) ^ 2)
}

mse_vec = function(truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  check_numeric_metric(truth, estimate, case_weights)
  
  if (na_rm) {
    result = yardstick_remove_missing(truth, estimate, case_weights)
    
    truth = result$truth
    estimate = result$estimate
    case_weights = result$case_weights
  } else if (yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  
  mse_impl(truth, estimate, case_weights = case_weights)
}

mse = function(data, ...) {
  UseMethod("mse")
}

mse = new_numeric_metric(mse, direction = "minimize")

mse.data.frame = function(data, truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  
  numeric_metric_summarizer(
    name = "mse",
    fn = mse_vec,
    data = data,
    truth = !!enquo(truth),
    estimate = !!enquo(estimate),
    na_rm = na_rm,
    case_weights = !!enquo(case_weights)
  )
}

## spearman ----

spearman_impl = function(truth, estimate, case_weights = NULL) {
  cor(truth, estimate, method = 'spearman')
}

spearman_vec = function(truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  check_numeric_metric(truth, estimate, case_weights)
  
  if (na_rm) {
    result = yardstick_remove_missing(truth, estimate, case_weights)
    
    truth = result$truth
    estimate = result$estimate
    case_weights = result$case_weights
  } else if (yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  
  spearman_impl(truth, estimate, case_weights = case_weights)
}

spearman = function(data, ...) {
  UseMethod("spearman")
}

spearman = new_numeric_metric(spearman, direction = "maximize")

spearman.data.frame = function(data, truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  
  numeric_metric_summarizer(
    name = "spearman",
    fn = spearman_vec,
    data = data,
    truth = !!enquo(truth),
    estimate = !!enquo(estimate),
    na_rm = na_rm,
    case_weights = !!enquo(case_weights)
  )
}

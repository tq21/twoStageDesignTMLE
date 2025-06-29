library(purrr)
library(dplyr)
library(devtools)
library(survey)
library(mice)
library(data.table)
library(marginaleffects)
load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")
source("../williamson_et_al/src_williamson_et_al/00_utils.R")
source("../williamson_et_al/src_williamson_et_al/01_generate_data.R")
source("../williamson_et_al/src_williamson_et_al/02_methods_design.R")
source("../williamson_et_al/src_williamson_et_al/02_methods_tmle.R")
source("../williamson_et_al/src_williamson_et_al/03_estimate.R")
set.seed(123)

run <- function(Y_type,
                miss_type,
                true_Pi,
                truth,
                seed) {
  set.seed(seed)
  B <- 100#500
  n_seq <- 2000#seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      cat("n: " %+% .n %+% ", b: " %+% .b %+% "... ")
      data_obj <- sim_data(n = .n, Y_type = Y_type, miss_type = miss_type)
      data <- data_obj$data

      # raking -----------------------------------------------------------------
      res_rak <- run_raking_lr(data = data,
                               formula = "Y ~ A + W1 + W2 + W3 + W4",
                               miss_formula = "Delta ~ W1 + W2 + A + Y",
                               NimpRaking = 20,
                               calOption = 1,
                               fam = "binomial",
                               coefficient_of_interest = "A",
                               missing_indicator = "Delta",
                               start_from_ipw = FALSE,
                               rake_on_y = FALSE)

      cur_df <- data.frame(n = .n,
                           b = .b,
                           est_name = "raking",
                           psi = res_rak$results[res_rak$results$estimand == "RD", "est"],
                           lower = res_rak$results[res_rak$results$estimand == "RD", "est"]+qnorm(0.025)*res_rak$results[res_rak$results$estimand == "RD", "SE"],
                           upper = res_rak$results[res_rak$results$estimand == "RD", "est"]+qnorm(0.975)*res_rak$results[res_rak$results$estimand == "RD", "SE"])

      cur_mse <- mean((cur_df$psi-truth)^2)
      cur_cover <- mean(cur_df$lower <= truth & cur_df$upper >= truth)
      cat("cur MSE: " %+% cur_mse %+% ", coverage: " %+% cur_cover %+% "\n")

      return(cur_df)
    })
  })

  return(res_df)
}

.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
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
library(future)
ncores <- 10
plan(multisession, workers = ncores)
set.seed(123)

run <- function(n,
                B,
                miss_type,
                Q_g_method,
                pi.SL.library,
                Q.SL.library,
                g.SL.library,
                DFullReg_sl_lib) {

  res_df <- map_dfr(seq(B), function(.b) {
    print("n: " %+% n %+% ", b: " %+% .b %+% "...")
    data_obj <- sim_data(n = n, miss_type = miss_type)
    data <- data_obj$data

    # TMLE
    res_tmle <- imp_plugin_logistic(Y = data$Y,
                                    A = data$A,
                                    W = data[, c("W1", "W2", "W3"), drop = FALSE],
                                    W.stage2 = data[complete.cases(data), c("W4"), drop = FALSE],
                                    Delta.W = data$Delta,
                                    condSetNames = c("W", "A", "Y"),
                                    Q.family = "binomial",
                                    Nimp = 20,
                                    Q_g_method = Q_g_method,
                                    pi.SL.library = pi.SL.library,
                                    Q.SL.library = Q.SL.library,
                                    g.SL.library = g.SL.library,
                                    DFullReg_sl_lib = DFullReg_sl_lib,
                                    augmentW = FALSE,
                                    verbose = FALSE,
                                    browse = FALSE)

    # raking -----------------------------------------------------------------
    res_rak <- run_raking_lr(data = data,
                             formula = "Y ~ A + W1 + W2 + W3 + W4",
                             miss_formula = "Delta ~ A + W1 + W2 + W3 + Y",
                             NimpRaking = 20,
                             calOption = 1,
                             fam = "binomial",
                             pi = as.numeric(res_tmle$pi),
                             coefficient_of_interest = "A",
                             missing_indicator = "Delta",
                             start_from_ipw = FALSE,
                             rake_on_y = FALSE)

    return(data.frame(n = n,
                      b = .b,
                      est_name = c("IPCW-TMLE",
                                   "IPCW-TMLE-imp-plugin-MI",
                                   "A-IPCW",
                                   "raking"),
                      psi = c(res_tmle$tmle$estimates$ATE$psi,
                              res_tmle$psi,
                              res_tmle$psi_aipcw,
                              res_rak$results[res_rak$results$estimand == "RD", "est"]),
                      lower = c(res_tmle$tmle$estimates$ATE$CI[1],
                                res_tmle$lower,
                                res_tmle$lower_aipcw,
                                res_rak$results[res_rak$results$estimand == "RD", "est"]+qnorm(0.025)*res_rak$results[res_rak$results$estimand == "RD", "SE"]),
                      upper = c(res_tmle$tmle$estimates$ATE$CI[2],
                                res_tmle$upper,
                                res_tmle$upper_aipcw,
                                res_rak$results[res_rak$results$estimand == "RD", "est"]+qnorm(0.975)*res_rak$results[res_rak$results$estimand == "RD", "SE"])))
  })

  return(res_df)
}

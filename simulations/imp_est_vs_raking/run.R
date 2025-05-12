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
                Delta_type,
                true_Pi,
                Q.family) {
  B <- 50#500
  n_seq <- 500#seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
      data_obj <- sim_data(.n, Y_type, Delta_type)
      data <- data_obj$data
      Pi <- data_obj$Pi

      args_tmle <- list(
        Y = data$Y,
        A = data$A,
        W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
        W.stage2 = data[complete.cases(data), c("W5"), drop = FALSE],
        Delta.W = data$Delta,
        condSetNames = c("W", "A", "Y"),
        pi.SL.library = "SL.glm", V.pi = 10,
        Q.family = Q.family,
        Q.SL.library = "SL.glm", V.Q = 10,
        g.SL.library = "SL.glm", V.g = 10,
        augmentW = FALSE,
        verbose = FALSE,
        browse = FALSE
      )
      if (true_Pi) {
        args_tmle$pi_oracle <- Pi
      }
      res_tmle <- do.call("twoStageTMLE_plugin_tmle_MI", args_tmle)

      # raking -----------------------------------------------------------------
      args_rak <- list(
        data = data,
        formula = "Y ~ A + W1 + W2 + W3 + W4 + W5",
        miss_formula = "Delta ~ A + W1 + W2 + W3 + W4",
        NimpRaking = 20,
        calOption = 1,
        fam = Q.family,
        coefficient_of_interest = "A",
        missing_indicator = "Delta",
        start_from_ipw = FALSE,
        rake_on_y = FALSE
      )
      res_rak <- do.call("run_raking_lr", args_rak)

      return(data.frame(n = .n,
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
  })

  return(res_df)
}

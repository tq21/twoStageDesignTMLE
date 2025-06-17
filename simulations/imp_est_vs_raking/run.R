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
                seed,
                truth) {

  set.seed(seed)
  B <- 100#500
  n_seq <- 2000#seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    ipcw_tmle_psi <- ipcw_tmle_imp_psi <- aipcw_tmle_imp_psi <- rak_psi <-
      ipcw_tmle_lower <- ipcw_tmle_imp_lower <- aipcw_tmle_imp_lower <- rak_lower <-
      ipcw_tmle_upper <- ipcw_tmle_imp_upper <- aipcw_tmle_imp_upper <- rak_upper <- rep(NA, B)

    walk(seq(B), function(.b) {
      cat("n: " %+% .n %+% ", b: " %+% .b %+% "...\n")
      data_obj <- sim_data(n = .n, Y_type = Y_type, miss_type = miss_type)
      data <- data_obj$data
      Pi <- data_obj$Pi

      args_tmle <- list(
        Y = data$Y,
        A = data$A,
        W = data[, c("W1", "W2"), drop = FALSE],
        W.stage2 = data[complete.cases(data), c("W3", "W4"), drop = FALSE],
        Delta.W = data$Delta,
        condSetNames = c("W", "A", "Y"),
        Q_g_method = "ipcw",
        pi.SL.library = c("SL.glm", "SL.gam", "SL.dbarts"), V.pi = 10,
        Q.family = "binomial",
        Q.SL.library = c("SL.glm", "SL.gam", "SL.dbarts"), V.Q = 10,
        g.SL.library = c("SL.glm", "SL.gam", "SL.dbarts"), V.g = 10,
        augmentW = FALSE,
        verbose = FALSE,
        browse = FALSE
      )
      if (true_Pi) {
        args_tmle$pi_oracle <- Pi
      }
      res_tmle <- do.call("imp_plugin_logistic", args_tmle)

      # raking -----------------------------------------------------------------
      args_rak <- list(
        data = data,
        formula = "Y ~ A + W1 + W2 + W3 + W4",
        miss_formula = "Delta ~ W1 + W2 + A + Y",
        NimpRaking = 20,
        calOption = 1,
        fam = "binomial",
        coefficient_of_interest = "A",
        missing_indicator = "Delta",
        start_from_ipw = FALSE,
        rake_on_y = FALSE
      )
      res_rak <- do.call("run_raking_lr", args_rak)

      # store results
      ipcw_tmle_psi[.b] <<- res_tmle$tmle$estimates$ATE$psi
      ipcw_tmle_lower[.b] <<- res_tmle$tmle$estimates$ATE$CI[1]
      ipcw_tmle_upper[.b] <<- res_tmle$tmle$estimates$ATE$CI[2]
      ipcw_tmle_imp_psi[.b] <<- res_tmle$psi
      ipcw_tmle_imp_lower[.b] <<- res_tmle$lower
      ipcw_tmle_imp_upper[.b] <<- res_tmle$upper
      aipcw_tmle_imp_psi[.b] <<- res_tmle$psi_aipcw
      aipcw_tmle_imp_lower[.b] <<- res_tmle$lower_aipcw
      aipcw_tmle_imp_upper[.b] <<- res_tmle$upper_aipcw
      rak_psi[.b] <<- res_rak$results[res_rak$results$estimand == "RD", "est"]
      rak_lower[.b] <<- res_rak$results[res_rak$results$estimand == "RD", "est"] + qnorm(0.025) * res_rak$results[res_rak$results$estimand == "RD", "SE"]
      rak_upper[.b] <<- res_rak$results[res_rak$results$estimand == "RD", "est"] + qnorm(0.975) * res_rak$results[res_rak$results$estimand == "RD", "SE"]

      # rak fails
      if (is.na(rak_psi[.b]) | is.na(rak_lower[.b]) | is.na(rak_upper[.b])) {
        rak_psi[.b] <- NA; rak_lower[.b] <- NA; rak_upper[.b] <- NA
      }

      # print running results
      cur_mse_rak <- mean((rak_psi-truth)^2, na.rm = TRUE)
      cur_cover_rak <- mean(rak_lower <= truth & rak_upper >= truth, na.rm = TRUE)
      cur_mse_ipcw <- mean((ipcw_tmle_imp_psi-truth)^2, na.rm = TRUE)
      cur_cover_ipcw <- mean(ipcw_tmle_imp_lower <= truth & ipcw_tmle_imp_upper >= truth, na.rm = TRUE)
      cat("raking MSE: " %+% cur_mse_rak %+% ", coverage: " %+% round(cur_cover_rak, 2) %+% "\n")
      cat("IPCW-TMLE-imp MSE: " %+% cur_mse_ipcw %+% ", coverage: " %+% round(cur_cover_ipcw, 2) %+% "\n")
    })

    return(rbind(data.frame(n = .n, b = seq_len(B),
                            est_name = "ipcw_tmle",
                            psi = ipcw_tmle_psi,
                            lower = ipcw_tmle_lower, upper = ipcw_tmle_upper),
                 data.frame(n = .n, b = seq_len(B),
                            est_name = "ipcw_tmle_imp",
                            psi = ipcw_tmle_imp_psi,
                            lower = ipcw_tmle_imp_lower, upper = ipcw_tmle_imp_upper),
                 data.frame(n = .n, b = seq_len(B),
                            est_name = "aipcw_tmle_imp",
                            psi = aipcw_tmle_imp_psi,
                            lower = aipcw_tmle_imp_lower, upper = aipcw_tmle_imp_upper),
                 data.frame(n = .n, b = seq_len(B),
                            est_name = "raking",
                            psi = rak_psi, lower = rak_lower, upper = rak_upper)))
  })

  return(res_df)
}

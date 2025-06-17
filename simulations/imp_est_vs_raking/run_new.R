library(purrr)
library(dplyr)
library(devtools)
library(survey)
library(mice)
library(data.table)
library(marginaleffects)
library(sl3)
library(doMC)
registerDoMC(cores = 5)
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

  B <- 100#500
  n_seq <- 2000#seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    ipcw_tmle_psi <- ipcw_tmle_imp_psi <- aipcw_tmle_imp_psi <- rak_psi <-
      ipcw_tmle_lower <- ipcw_tmle_imp_lower <- aipcw_tmle_imp_lower <- rak_lower <-
      ipcw_tmle_upper <- ipcw_tmle_imp_upper <- aipcw_tmle_imp_upper <- rak_upper <- rep(NA, B)

    walk(seq(B), function(.b) {
      cat("n: " %+% .n %+% ", b: " %+% .b %+% "...\n")
      set.seed(seed+.b)
      data_obj <- sim_data(n = .n, Y_type = Y_type, miss_type = miss_type)
      data <- data_obj$data
      Pi <- data_obj$Pi

      # IPCW-TMLE + efficient IPCW-TMLE
      res_tmle <- ipcw_tmle_imp(data = data,
                                W_name = c("W1", "W2", "W3", "W4"),
                                A_name = "A",
                                Y_name = "Y",
                                Delta_name = "Delta",
                                V_name = c("W1", "W2", "A", "Y"),
                                nfolds = 5,
                                Pi_method = "HAL",
                                Q_method = "HAL",
                                g_method = "HAL",
                                m_method = "HAL",
                                family = "binomial",
                                enumerate_basis_args = list(max_degree = 3,
                                                            smoothness_orders = 1,
                                                            num_knots = 5),
                                browse = FALSE)

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

      # store results
      ipcw_tmle_psi[.b] <<- res_tmle$psi_ipcw_tmle
      ipcw_tmle_lower[.b] <<- res_tmle$lower_ipcw_tmle
      ipcw_tmle_upper[.b] <<- res_tmle$upper_ipcw_tmle
      ipcw_tmle_imp_psi[.b] <<- res_tmle$psi_tmle
      ipcw_tmle_imp_lower[.b] <<- res_tmle$lower_tmle
      ipcw_tmle_imp_upper[.b] <<- res_tmle$upper_tmle
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
      cur_mse_ipcw_tmle <- mean((ipcw_tmle_psi-truth)^2, na.rm = TRUE)
      cur_cover_ipcw_tmle <- mean(ipcw_tmle_lower <= truth & ipcw_tmle_upper >= truth, na.rm = TRUE)
      cur_mse_tmle <- mean((ipcw_tmle_imp_psi-truth)^2, na.rm = TRUE)
      cur_cover_tmle <- mean(ipcw_tmle_imp_lower <= truth & ipcw_tmle_imp_upper >= truth, na.rm = TRUE)
      cat("raking MSE:              " %+% cur_mse_rak %+% ", coverage: " %+% round(cur_cover_rak, 2) %+% "\n")
      cat("IPCW-TMLE MSE:           " %+% cur_mse_ipcw_tmle %+% ", coverage: " %+% round(cur_cover_ipcw_tmle, 2) %+% "\n")
      cat("Efficient IPCW-TMLE MSE: " %+% cur_mse_tmle %+% ", coverage: " %+% round(cur_cover_tmle, 2) %+% "\n")
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

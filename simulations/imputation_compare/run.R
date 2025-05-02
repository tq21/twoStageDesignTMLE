library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")
set.seed(123)

run <- function(Y_type,
                Delta_type,
                true_Pi,
                Q.family) {
  B <- 50#100
  n_seq <- 500#seq(1000, 4000, 1000)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
      data_obj <- sim_data(.n, Y_type, Delta_type)
      data <- data_obj$data
      Pi <- data_obj$Pi

      call_args <- list(Y = data$Y,
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
                        browse = FALSE)
      if (true_Pi) {
        call_args$pi_oracle <- Pi
      }
      res <- do.call("twoStageTMLE_plugin_DFullReg_compare", call_args)

      return(data.frame(n = .n,
                        b = .b,
                        est_name = c("IPCW-TMLE",
                                     "SL-AIPCW",
                                     "target-DFullNCReg-SL",
                                     "target-DFullNCReg-MI",
                                     "target-DFullNCReg-MI-direct"),
                        psi = c(res$tmle$estimates$ATE$psi,
                                res$psi_aipcw,
                                res$psi,
                                res$psi_MI,
                                res$psi_MI_direct),
                        lower = c(res$tmle$estimates$ATE$CI[1],
                                  res$lower_aipcw,
                                  res$lower,
                                  res$lower_MI,
                                  res$lower_MI_direct),
                        upper = c(res$tmle$estimates$ATE$CI[2],
                                  res$upper_aipcw,
                                  res$upper,
                                  res$upper_MI,
                                  res$upper_MI_direct)))
    })
  })

  return(res_df)
}

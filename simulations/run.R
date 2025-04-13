library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")
set.seed(123)

run <- function(Y_type,
                Delta_type) {
  B <- 500
  n_seq <- seq(2000, 5000, 1000)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
      data <- sim_data(.n)
      res <- twoStageTMLE(
        Y = data$Y,
        A = data$A,
        W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
        W.stage2 = data[complete.cases(data), c("W5"), drop = FALSE],
        Delta.W = data$Delta,
        target_pi_max_iter=5,
        condSetNames = c("W", "A", "Y"),
        pi.SL.library = "SL.glm", V.pi = 10,
        Q.family = "gaussian",
        Q.SL.library = "SL.glm", V.Q = 10,
        g.SL.library = "SL.glm", V.g = 10,
        augmentW = FALSE,
        verbose = FALSE,
        browse = FALSE
      )

      return(data.frame(n = .n,
                        b = .b,
                        est_name = c("IPCW-TMLE",
                                     "IPCW-TMLE-target-Pi-1",
                                     "IPCW-TMLE-target-Pi-relearn-1",
                                     "IPCW-TMLE-target-Pi-2",
                                     "IPCW-TMLE-target-Pi-relearn-2",
                                     "IPCW-TMLE-target-Pi-3",
                                     "IPCW-TMLE-target-Pi-relearn-3",
                                     "IPCW-TMLE-target-Pi-4",
                                     "IPCW-TMLE-target-Pi-relearn-4",
                                     "IPCW-TMLE-target-Pi-5",
                                     "IPCW-TMLE-target-Pi-relearn-5",
                                     "IPCW-TMLE-imputation"),
                        psi = c(res$tmle_0$estimates$ATE$psi,
                                res$tmle_target_Pi_1$psi,
                                res$tmle_1$estimates$ATE$psi,
                                res$tmle_target_Pi_2$psi,
                                res$tmle_2$estimates$ATE$psi,
                                res$tmle_target_Pi_3$psi,
                                res$tmle_3$estimates$ATE$psi,
                                res$tmle_target_Pi_4$psi,
                                res$tmle_4$estimates$ATE$psi,
                                res$tmle_target_Pi_5$psi,
                                res$tmle_5$estimates$ATE$psi,
                                res$ATE),
                        lower = c(res$tmle_0$estimates$ATE$CI[1],
                                  res$tmle_target_Pi_1$lower,
                                  res$tmle_1$estimates$ATE$CI[1],
                                  res$tmle_target_Pi_2$lower,
                                  res$tmle_2$estimates$ATE$CI[1],
                                  res$tmle_target_Pi_3$lower,
                                  res$tmle_3$estimates$ATE$CI[1],
                                  res$tmle_target_Pi_4$lower,
                                  res$tmle_4$estimates$ATE$CI[1],
                                  res$tmle_target_Pi_5$lower,
                                  res$tmle_5$estimates$ATE$CI[1],
                                  res$CI[1]),
                        upper = c(res$tmle_0$estimates$ATE$CI[2],
                                  res$tmle_target_Pi_1$upper,
                                  res$tmle_1$estimates$ATE$CI[2],
                                  res$tmle_target_Pi_2$upper,
                                  res$tmle_2$estimates$ATE$CI[2],
                                  res$tmle_target_Pi_3$upper,
                                  res$tmle_3$estimates$ATE$CI[2],
                                  res$tmle_target_Pi_4$upper,
                                  res$tmle_4$estimates$ATE$CI[2],
                                  res$tmle_target_Pi_5$upper,
                                  res$tmle_5$estimates$ATE$CI[2],
                                  res$CI[2])))
    })
  })

  return(res_df)
}

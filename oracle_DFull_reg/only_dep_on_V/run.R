library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")
set.seed(123)

run <- function() {
  B <- 500
  n_seq <- seq(1000, 5000, 1000)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
      data_obj <- sim_data(.n)
      data <- data_obj$data

      # extract true values to compute true E(D^F|V,DELTA=1)
      g1W <- data_obj$g1W; Pi <- data_obj$Pi; psi <- data_obj$psi
      Q1W <- data_obj$Q1W; Q0W <- data_obj$Q0W; QAW <- data_obj$QAW
      DFullReg_true <- (data$A/g1W-(1-data$A)/(1-g1W))*(data$Y-QAW)+Q1W-Q0W-psi

      res <- twoStageTMLE(
        Y = data$Y,
        A = data$A,
        W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
        W.stage2 = data[complete.cases(data), c("W5"), drop = FALSE],
        Delta.W = data$Delta,
        condSetNames = c("W", "A", "Y"),
        pi_oracle = Pi,
        DFullReg_true = DFullReg_true,
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
                                     "IPCW-TMLE-target-Pi-true-Pi",
                                     "IPCW-TMLE-target-Pi-true-Pi-DFullReg"),
                        psi = c(res$tmle$estimates$ATE$psi,
                                res$tmle_Pi$psi,
                                res$tmle_Pi_DFullReg_true$psi),
                        lower = c(res$tmle$estimates$ATE$CI[1],
                                  res$tmle_Pi$lower,
                                  res$tmle_Pi_DFullReg_true$lower),
                        upper = c(res$tmle$estimates$ATE$CI[2],
                                  res$tmle_Pi$upper,
                                  res$tmle_Pi_DFullReg_true$upper)))
    })
  })

  return(res_df)
}

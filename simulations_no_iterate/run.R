library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)
source("sim_data.R")
set.seed(123)

run <- function(Y_type,
                Delta_type) {
  B <- 100
  n_seq <- 2000#seq(1000, 5000, 1000)

  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
      data_obj <- sim_data(.n, Y_type, Delta_type)
      data <- data_obj$data
      Pi <- data_obj$Pi
      res <- twoStageTMLE(
        Y = data$Y,
        A = data$A,
        W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
        W.stage2 = data[complete.cases(data), c("W5"), drop = FALSE],
        Delta.W = data$Delta,
        pi_oracle = Pi,
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
                                     "IPCW-TMLE-imputation"),
                        psi = c(res$tmle$estimates$ATE$psi,
                                res$impute$psi),
                        lower = c(res$tmle$estimates$ATE$CI[1],
                                  res$impute$lower),
                        upper = c(res$tmle$estimates$ATE$CI[2],
                                  res$impute$upper)))
    })
  })

  return(res_df)
}

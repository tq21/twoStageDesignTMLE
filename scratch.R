library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)

sim_data <- function(n,
                     A_counter = -1) {
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  W4 <- runif(n, -1, 1)
  W5 <- runif(n, -1, 1)
  UY <- rnorm(n, 0, 1)
  # UY <- rnorm(n, 0, 2)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-W1-W2+W3+W4))
  } else {
    A <- rep(A_counter, n)
  }
  # Y <- 2.1+2*A-W1-1.5*W2-1.5*W5+UY
  Y <- 0.2+0.4*A-W1-0.8*W2-0.8*W5+UY # good
  #Delta <- rbinom(n, 1, plogis(-W1^3+0.4*W2*A*W1)) # good
  Delta <- rbinom(n, 1, plogis(-6*W1+6*W2))
  #Delta <- rbinom(n, 1, plogis(1.3*W1-0.4*as.numeric(W2 > 0.5)+0.1*W3^2))
  W5[Delta == 0] <- NA

  return(data.frame(W1 = W1,
                    W2 = W2,
                    W3 = W3,
                    W4 = W4,
                    W5 = W5,
                    A = A,
                    Y = Y,
                    Delta = Delta))
}


set.seed(123)
truth <- 0.4

B <- 500
n_seq <- 2000#seq(500, 2000, 500)

res_df <- map_dfr(n_seq, function(.n) {
  map_dfr(seq(B), function(.b) {
    print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
    data <- sim_data(.n)

    # Susan's
    res <- twoStageTMLE(
      Y = data$Y,
      A = data$A,
      W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
      W.stage2 = data[complete.cases(data), c("W5"), drop = FALSE],
      Delta.W = data$Delta,
      #pi_oracle = rep(0.7, .n) + rnorm(.n, 0.05, 0.01),
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

res_df %>%
  summarize(abs_bias = abs(mean(psi - truth)),
            se = sd(psi),
            mse = mean((psi - truth)^2),
            bias_se = abs_bias/se,
            coverage = mean((lower <= truth) & (truth <= upper)),
            power = mean((lower >= 0) & (0 <= upper)),
            .by = c("n", "est_name"))

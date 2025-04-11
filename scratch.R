library(purrr)
library(dplyr)
library(devtools)
load_all()
`%+%` <- function(a, b) paste0(a, b)

sim_data <- function(n,
                     A_counter = -1) {
  W1 <- runif(n, 0, 1)
  W2 <- runif(n, 0, 1)
  W3 <- runif(n, 0, 1)
  W4 <- runif(n, 0, 1)
  W5 <- runif(n, 0, 1)
  UY <- rnorm(n, 0, 0.5)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-W1-W2+W3+W4))
  } else {
    A <- rep(A_counter, n)
  }
  #Y <- rbinom(n, 1, plogis(A-4*W1+A*W1-1.5*W2+W5))
  #Y <- rbinom(n, 1, plogis(A-4*W1+A*W1-1.5*W2+sin(W5)))
  Y <- A-4*W1-1.5*W2+sin(W5)+UY
  #Delta <- rbinom(n, 1, plogis(1.3*W1-0.4*as.numeric(W2 > 0.5)+0.1*W3^2))
  Delta <- rbinom(n, 1, plogis(0.5))
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

get_truth <- function() {
  Y1 <- mean(sim_data(1e6, 1)$Y)
  Y0 <- mean(sim_data(1e6, 0)$Y)
  return(mean(Y1-Y0))
}

set.seed(123)
truth <- 1

B <- 100
n_seq <- 500#seq(500, 2000, 500)

res_df <- map_dfr(n_seq, function(.n) {
  map_dfr(seq(B), function(.b) {
    print("n: " %+% .n %+% ", b: " %+% .b %+% "...")
    data <- sim_data(.n)

    # Susan's
    res_susan <- twoStageTMLE(
      Y = data$Y,
      A = data$A,
      W = data[, c("W1", "W2", "W3", "W4"), drop = FALSE],
      W.stage2 = data[complete.cases(data), c("W1", "W2", "W3", "W4", "W5"), drop = FALSE],
      Delta.W = data$Delta,
      condSetNames = c("W", "A", "Y"),
      pi.SL.library = "SL.glm", V.pi = 10,
      Q.family = "gaussian",
      Q.SL.library = "SL.xgboost", V.Q = 5,
      g.SL.library = "SL.glm", V.g = 10,
      augmentW = FALSE,
      verbose = FALSE
    )

    return(data.frame(n = .n,
                      b = .b,
                      est_name = c("IPCW-TMLE", "IPCW-TMLE-target-Pi", "IPCW-TMLE-imputation"),
                      psi = c(res_susan$tmle$estimates$ATE$psi, res_susan$tmle_Pi_star$estimates$ATE$psi, res_susan$ATE),
                      lower = c(res_susan$tmle$estimates$ATE$CI[1], res_susan$tmle_Pi_star$estimates$ATE$CI[1], res_susan$CI[1]),
                      upper = c(res_susan$tmle$estimates$ATE$CI[2], res_susan$tmle_Pi_star$estimates$ATE$CI[2], res_susan$CI[2])))
  })
})

res_df %>%
  summarize(abs_bias = abs(mean(psi - truth)),
            se = sd(psi),
            mse = mean((psi - truth)^2),
            coverage = mean((lower <= truth) & (truth <= upper)),
            .by = c("n", "est_name"))

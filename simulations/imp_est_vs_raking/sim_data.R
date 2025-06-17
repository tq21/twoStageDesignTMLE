sim_data <- function(n,
                     Y_type,
                     miss_type,
                     A_counter = -1) {
  # latent
  Z1 <- rnorm(n, 0, 1)
  Z2 <- rnorm(n, 0, 1)
  Z3 <- rnorm(n, 0, 1)
  Z4 <- rnorm(n, 0, 1)

  # observed
  W1 <- exp(Z1/2)
  W2 <- Z2/(1+exp(Z1))+10
  W3 <- (Z4*Z3/25+0.6)^3
  W4 <- (Z3+Z4+20)^2

  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-0.2*Z1-0.6*Z2+0.9*Z4))
  } else {
    A <- rep(A_counter, n)
  }

  # outcome regression
  if (Y_type == "complex") {
    lp_Q <- -1+0.6*Z1-0.4*Z2+0.2*Z3-0.5*Z4+A*(1.2+0.9*Z1-0.7*Z2)
    Y <- rbinom(n, 1, plogis(lp_Q))
  } else if (Y_type == "rare") {
    lp_Q <- -4.7+0.6*Z1-0.4*Z2+0.2*Z3-0.5*Z4+A*(1.2+0.9*Z1-0.7*Z2)
    Y <- rbinom(n, 1, plogis(lp_Q))
  }

  # missing mechanism
  if (miss_type == "strong") {
    # covariates strongly predictive of missingness
    Pi <- plogis(-1.5*Z1+1.5*Z2)
    Delta <- rbinom(n, 1, Pi)
  }
  W3[Delta == 0] <- NA
  W4[Delta == 0] <- NA

  return(list(data = data.frame(W1 = W1,
                                W2 = W2,
                                W3 = W3,
                                W4 = W4,
                                A = A,
                                Y = Y,
                                Delta = Delta),
              Pi = Pi))
}

get_truth <- function(Y_type) {
  data_A1 <- sim_data(1e7, Y_type = Y_type, miss_type = "strong", A_counter = 1)$data
  data_A0 <- sim_data(1e7, Y_type = Y_type, miss_type = "strong", A_counter = 0)$data
  return(mean(data_A1$Y-data_A0$Y))
}

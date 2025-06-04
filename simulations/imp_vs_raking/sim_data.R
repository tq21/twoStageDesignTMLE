sim_data <- function(n,
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
  W3 <- (Z1*Z3/25+0.6)^3
  W4 <- (Z2+Z4+20)^2

  # treatment mechanism
  if (A_counter == -1) {
    g1W <- plogis(-0.2*Z1-0.6*Z2+0.9*Z4)
    A <- rbinom(n, 1, g1W)
  } else {
    g1W <- A_counter
    A <- rep(A_counter, n)
  }

  # outcome regression
  Y_form <- function(A) {
    plogis(-2.25+0.45*A-0.2*Z1-0.4*Z2-0.4*Z3+0.50*Z4)
  }
  QAW <- Y_form(A)
  Q1W <- Y_form(1)
  Q0W <- Y_form(0)
  Y <- rbinom(n, 1, QAW)

  # missing mechanism
  if (miss_type == "med") {
    # ~50% missing
    Pi <- plogis(-0.2*Z1+0.2*Z2)
    Delta <- rbinom(n, 1, Pi)
  } else if (miss_type == "med_strong") {
    # ~50% missing, strongly predictive
    Pi <- plogis(-1.5*Z1+1.5*Z2)
    Delta <- rbinom(n, 1, Pi)
  } else if (miss_type == "large") {
    # ~80% missing
    Pi <- plogis(-1.4-0.2*Z1+0.2*Z2)
    Delta <- rbinom(n, 1, Pi)
  }
  W4[Delta == 0] <- NA

  return(list(data = data.frame(W1 = W1,
                                W2 = W2,
                                W3 = W3,
                                W4 = W4,
                                A = A,
                                Y = Y,
                                Delta = Delta),
              Pi = Pi,
              QAW = QAW,
              Q1W = Q1W,
              Q0W = Q0W,
              g1W = g1W))
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, miss_type = "med", A_counter = 1)$data
  data_A0 <- sim_data(1e7, miss_type = "med", A_counter = 0)$data
  return(mean(data_A1$Y-data_A0$Y))
}

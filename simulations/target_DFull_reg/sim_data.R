sim_data <- function(n,
                     Y_type,
                     Delta_type,
                     A_counter = -1) {
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  W4 <- runif(n, -1, 1)
  W5 <- runif(n, -1, 1)
  UY <- rnorm(n, 0, 1)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-W1-W2+W3+W4))
  } else {
    A <- rep(A_counter, n)
  }

  # outcome regression
  if (Y_type == "simple_gaussian") {
    # gaussian, linear
    Y <- 0.4*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == "complex_gaussian") {
    # gaussian, non-linear
    Y <- 0.1*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == "simple_binomial") {
    Y <- rbinom(n, 1, plogis(A-4*W1-1.5*W2+W5))
  } else if (Y_type == "complex_binomial") {
    Y <- rbinom(n, 1, plogis(A-4*W1+A*W1-1.5*W2*W1+sin(W5)))
  } else if (Y_type == "rare") {
    Y <- rbinom(n, 1, plogis(-A-4*abs(W1)-abs(W2)-abs(W5)))
  }

  # missing mechanism
  if (Delta_type == "med") {
    # ~50% missing
    Pi <- plogis(-0.5*W1+0.5*W2)
    Delta <- rbinom(n, 1, Pi)
  } else if (Delta_type == "med_strong") {
    # ~50% missing, strongly predictive
    Pi <- plogis(-6*W1+6*W2)
    Delta <- rbinom(n, 1, Pi)
  } else if (Delta_type == "large") {
    # ~80% missing
    Pi <- plogis(-2.3*abs(W1)-0.5*abs(W2))
    Delta <- rbinom(n, 1, Pi)
  }
  W5[Delta == 0] <- NA

  return(list(data = data.frame(W1 = W1,
                                W2 = W2,
                                W3 = W3,
                                W4 = W4,
                                W5 = W5,
                                A = A,
                                Y = Y,
                                Delta = Delta),
              Pi = Pi))
}

get_truth <- function(Y_type) {
  if (Y_type == "simple_gaussian") {
    return(0.4)
  } else if (Y_type == "complex_gaussian") {
    return(0.1)
  } else if (Y_type == "simple_binomial") {
    data_A1 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  } else if (Y_type == "complex_binomial") {
    data_A1 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  } else if (Y_type == "rare") {
    data_A1 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e5, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  }
}

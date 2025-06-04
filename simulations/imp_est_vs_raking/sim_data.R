sim_data <- function(n,
                     Y_type,
                     Delta_type,
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

  UY <- rnorm(n, 0, 1)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-0.2*Z1-0.6*Z2+0.9*Z4))
  } else {
    A <- rep(A_counter, n)
  }

  # hinge <- function(X, u) as.numeric(X >= u)*(X-u)

  # outcome regression
  if (Y_type == "simple_gaussian") {
    # gaussian, linear
    Y <- 0.4*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == "complex_gaussian") {
    # gaussian, non-linear
    Y <- 0.1*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == "simple_binomial") {
    Y <- rbinom(n, 1, plogis(-2.25+0.45*A-0.2*Z1-0.4*Z2-0.4*Z3+0.50*Z4))
  } else if (Y_type == "complex_binomial") {
    #Y_logit <- 0.2*A-4*hinge(Z1,0.5)+0.4*hinge(Z2,0.2)*hinge(Z4,-0.5)+sin(2*pi*Z1)+cos(2*pi*Z2)
    #Y <- rbinom(n, 1, plogis(Y_logit))
  } else if (Y_type == "rare") {
    Y <- rbinom(n, 1, plogis(-A-4*abs(W1)-abs(W2)-abs(W5)))
  }

  # missing mechanism
  if (Delta_type == "med") {
    # ~50% missing
    Pi <- plogis(-0.5*Z1+0.5*Z2)
    Delta <- rbinom(n, 1, Pi)
  } else if (Delta_type == "med_strong") {
    # ~50% missing, strongly predictive
    #Pi <- plogis(-6*hinge(Z1,0.4)*hinge(Z2,-0.4)+6*Z2^2-6*Z1-6*abs(Z1))
    Pi <- plogis(-1.5*Z1+1.5*Z2)
    Delta <- rbinom(n, 1, Pi)
  } else if (Delta_type == "large") {
    # ~80% missing
    Pi <- plogis(-2.3*abs(W1)-0.5*abs(W2))
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
              Pi = Pi))
}

get_truth <- function(Y_type) {
  if (Y_type == "simple_gaussian") {
    return(0.4)
  } else if (Y_type == "complex_gaussian") {
    return(0.1)
  } else if (Y_type == "simple_binomial") {
    data_A1 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  } else if (Y_type == "complex_binomial") {
    data_A1 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  } else if (Y_type == "rare") {
    data_A1 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 1)$data
    data_A0 <- sim_data(1e7, Y_type = Y_type, Delta_type = "med", A_counter = 0)$data
    return(mean(data_A1$Y-data_A0$Y))
  }
}

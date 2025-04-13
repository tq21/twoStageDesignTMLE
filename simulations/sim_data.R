sim_data <- function(n,
                     Y_type = 1.0,
                     Delta_type = 1.0,
                     A_counter = -1) {
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  W4 <- runif(n, -1, 1)
  W5 <- runif(n, -1, 1)
  UY <- rnorm(n, 0, 2)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-W1-W2+W3+W4))
  } else {
    A <- rep(A_counter, n)
  }

  if (Y_type == 1.0) {
    # gaussian, linear, canonical
    Y <- 0.4*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == 1.1) {
    # gaussian, linear, small effect size
    Y <- 0.1*A-W1-0.8*W2-0.8*W5+UY
  } else if (Y_type == 2.0) {
    # binomial, logit linear, canonical
    Y <- rbinom(n, 1, plogis(A-4*W1-1.5*W2+W5))
  } else if (Y_type == 2.1) {
    # binomial, logit linear, small effect size
    Y <- rbinom(n, 1, plogis(0.1*A-4*W1-1.5*W2+W5))
  }

  if (Delta_type == 1.0) {
    # 50% missing, logit linear
    Delta <- rbinom(n, 1, plogis(-0.5*W1+0.5*W2))
  } else if (Delta_type == 1.1) {
    # 50% missing, logit linear, strongly predictive
    Delta <- rbinom(n, 1, plogis(-6*W1+6*W2))
  }

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

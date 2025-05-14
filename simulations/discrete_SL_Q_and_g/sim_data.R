sim_data <- function(n,
                     A_counter = -1) {
  W1 <- round(runif(n, -1, 1), 3)
  W2 <- round(runif(n, -1, 1), 3)
  W3 <- rbinom(n, 1, 0.5)
  if (A_counter == -1) {
    A <- rbinom(n, 1, plogis(-W1-W2+W3))
  } else {
    A <- rep(A_counter, n)
  }

  Y <- rbinom(n, 1, plogis(A-4*W1-1.5*W2+0.5*W3+2*sin(3*W1*W2)-1.2*cos(W3^2)+0.8*A*W1-0.5*A*W2^2))
  Delta <- rbinom(n, 1, plogis(-0.5*W1+0.5*W2))
  W3[Delta == 0] <- NA

  return(data.frame(W1 = W1,
                    W2 = W2,
                    W3 = W3,
                    A = A,
                    Y = Y,
                    Delta = Delta))
}

get_truth <- function() {
  data_A1 <- sim_data(1e7, A_counter = 1)
  data_A0 <- sim_data(1e7, A_counter = 0)
  return(mean(data_A1$Y-data_A0$Y))
}

sim_data <- function(n,
                     A_counter = -1) {
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  W4 <- runif(n, -1, 1)
  W5 <- runif(n, -1, 1)
  UY <- rnorm(n, 0, 1)
  if (A_counter == -1) {
    # A independent of W5
    g1W <- plogis(-W1-W2+W3+W4)
    A <- rbinom(n, 1, g1W)
  } else {
    A <- rep(A_counter, n)
  }

  # Y independent of W5
  psi <- 0.4
  Y <- psi*A-W1-0.8*W2+UY
  Q1W <- psi*1-W1-0.8*W2
  Q0W <- psi*0-W1-0.8*W2
  QAW <- psi*A-W1-0.8*W2

  Pi <- plogis(-0.5*W1+0.5*W2)
  Delta <- rbinom(n, 1, Pi)
  W5[Delta == 0] <- NA

  return(list(data = data.frame(W1 = W1,
                                W2 = W2,
                                W3 = W3,
                                W4 = W4,
                                W5 = W5,
                                A = A,
                                Y = Y,
                                Delta = Delta),
              g1W = g1W,
              Pi = Pi,
              Q1W = Q1W,
              Q0W = Q0W,
              QAW = QAW,
              psi = psi))
}

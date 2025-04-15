sim_data <- function(n,
                     A_counter = -1) {
  W1 <- runif(n, -1, 1)
  W2 <- runif(n, -1, 1)
  W3 <- runif(n, -1, 1)
  W4 <- runif(n, -1, 1)
  W5 <- runif(n, -1, 1)
  UY <- rnorm(n, 0, 1)
  if (A_counter == -1) {
    g1W <- plogis(-W1-W2+W3+W4)
    A <- rbinom(n, 1, g1W)
  } else {
    A <- rep(A_counter, n)
  }

  psi <- 0.4
  Y <- psi*A-W1-0.8*W2-0.8*W5+UY
  Q1W <- psi*1-W1-0.8*W2-0.8*W5
  Q0W <- psi*0-W1-0.8*W2-0.8*W5
  QAW <- psi*A-W1-0.8*W2-0.8*W5

  DFullReg <- map_dbl(seq(n), function(i) {
    W5_grid <- seq(-1, 1, 0.001)
    Q1W_i <- psi*1-W1[i]-0.8*W2[i]-0.8*W5_grid
    Q0W_i <- psi*0-W1[i]-0.8*W2[i]-0.8*W5_grid
    QAW_i <- psi*A[i]-W1[i]-0.8*W2[i]-0.8*W5_grid
    DFullReg_i <- (A[i]/g1W[i]-(1-A[i])/(1-g1W[i]))*(Y[i]-QAW_i)+Q1W_i-Q0W_i-psi
    return(mean(DFullReg_i))
  })

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
              psi = psi,
              DFullReg = DFullReg))
}

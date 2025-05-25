library(purrr)
library(devtools)
load_all()
source("sim_data.R")

B <- 2
n_seq <- seq(500, 2000, 500)

res_df <- map_dfr(n_seq, function(.n) {
  map_dfr(1:B, function(.b) {
    print(paste0("n = ", .n, ", b = ", .b))
    data <- sim_data(.n)
    W1 <- data$W3
    W2 <- data[, c("W1", "W2")]
    Delta <- data$Delta
    A <- data$A
    Y <- data$Y

    # estimate Q and g
    Q_and_g <- learn_lik(W1 = W1,
                         W2 = as.matrix(W2),
                         Delta = Delta,
                         A = A,
                         Y = Y)
    QAW <- Q_and_g$QAW
    g1W <- Q_and_g$g1W

    loss_Q <- mean(log_lik_loss(QAW, Y))
    loss_g <- mean(log_lik_loss(g1W, A))

    return(data.frame(n = .n,
                      b = .b,
                      loss_Q = loss_Q,
                      loss_g = loss_g))
  })
})




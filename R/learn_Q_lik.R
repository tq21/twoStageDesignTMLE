#' Likelihood-based approach to estimate Q(A,W)
#' Assume W1 is one-dimensional binary for now, assume Y is binary
library(hal9001)
library(glmnet)
set.seed(123)
data <- sim_data(500)
W2 <- as.matrix(data[, c("W1", "W2")])
Delta <- data$Delta
W1 <- data$W3[Delta==1]
A <- data$A
Y <- data$Y

learn_lik <- function(W1,
                      W2,
                      Delta,
                      A,
                      Y) {

  W2_A_Y <- as.matrix(cbind(W2, A=A, Y=Y)[Delta==1,])
  W2_A0_Y <- as.matrix(cbind(W2, A=0, Y=Y))
  W2_A0_Y0 <- as.matrix(cbind(W2, A=0, Y=0))
  W2_A0_Y1 <- as.matrix(cbind(W2, A=0, Y=1))
  W2_A1_Y <- as.matrix(cbind(W2, A=1, Y=Y))
  W2_A1_Y0 <- as.matrix(cbind(W2, A=1, Y=0))
  W2_A1_Y1 <- as.matrix(cbind(W2, A=1, Y=1))

  # estimate p(W1|W2,A,Y)
  W2_A_Y_basis_list <- enumerate_basis(x = W2_A_Y, smoothness_orders = 1, max_degree = 3)
  W2_A_Y_design <- make_design_matrix(X = W2_A_Y, blist = W2_A_Y_basis_list)
  W2_A0_Y_design <- make_design_matrix(X = W2_A0_Y, blist = W2_A_Y_basis_list)
  W2_A0_Y0_design <- make_design_matrix(X = W2_A0_Y0, blist = W2_A_Y_basis_list)
  W2_A0_Y1_design <- make_design_matrix(X = W2_A0_Y1, blist = W2_A_Y_basis_list)
  W2_A1_Y_design <- make_design_matrix(X = W2_A1_Y, blist = W2_A_Y_basis_list)
  W2_A1_Y0_design <- make_design_matrix(X = W2_A1_Y0, blist = W2_A_Y_basis_list)
  W2_A1_Y1_design <- make_design_matrix(X = W2_A1_Y1, blist = W2_A_Y_basis_list)
  fit_W1 <- cv.glmnet(x = W2_A_Y_design,
                      y = W1,
                      alpha = 1,
                      family = "binomial",
                      keep = TRUE)
  pred_W1_W2_A0_Y <- as.numeric(predict(fit_W1, newx = W2_A0_Y_design, s = "lambda.min", type = "response"))
  pred_W1_W2_A0_Y0 <- as.numeric(predict(fit_W1, newx = W2_A0_Y0_design, s = "lambda.min", type = "response"))
  pred_W1_W2_A0_Y1 <- as.numeric(predict(fit_W1, newx = W2_A0_Y1_design, s = "lambda.min", type = "response"))
  pred_W1_W2_A1_Y <- as.numeric(predict(fit_W1, newx = W2_A1_Y_design, s = "lambda.min", type = "response"))
  pred_W1_W2_A1_Y0 <- as.numeric(predict(fit_W1, newx = W2_A1_Y0_design, s = "lambda.min", type = "response"))
  pred_W1_W2_A1_Y1 <- as.numeric(predict(fit_W1, newx = W2_A1_Y1_design, s = "lambda.min", type = "response"))

  # estimate p(Y|W2,A)
  W2_A_basis_list <- enumerate_basis(x = as.matrix(cbind(W2,A=A)), smoothness_orders = 1, max_degree = 3)
  W2_A_design <- make_design_matrix(X = as.matrix(cbind(W2,A=A)), blist = W2_A_basis_list)
  W2_A0_design <- make_design_matrix(X = as.matrix(cbind(W2,A=0)), blist = W2_A_basis_list)
  W2_A1_design <- make_design_matrix(X = as.matrix(cbind(W2,A=1)), blist = W2_A_basis_list)
  fit_Y <- cv.glmnet(x = W2_A_design,
                     y = Y,
                     alpha = 1,
                     family = "binomial",
                     keep = TRUE)
  pred_Y_W2_A0 <- as.numeric(predict(fit_Y, newx = W2_A0_design, s = "lambda.min", type = "response"))
  pred_Y_W2_A1 <- as.numeric(predict(fit_Y, newx = W2_A1_design, s = "lambda.min", type = "response"))

  # estimate p(A|W2)
  W2_basis_list <- enumerate_basis(x = W2, smoothness_orders = 1, max_degree = 3)
  W2_design <- make_design_matrix(X = W2, blist = W2_basis_list)
  fit_A <- cv.glmnet(x = W2_design,
                     y = A,
                     alpha = 1,
                     family = "binomial",
                     keep = TRUE)
  pred_A <- as.numeric(predict(fit_A, newx = W2_design, s = "lambda.min", type = "response"))

  # compute p(Y|W1,W2,A)
  pred_W1_W2_A1 <- pred_W1_W2_A1_Y1*pred_Y_W2_A1+pred_W1_W2_A1_Y0*(1-pred_Y_W2_A1)
  pred_W1_W2_A0 <- pred_W1_W2_A0_Y1*pred_Y_W2_A0+pred_W1_W2_A0_Y0*(1-pred_Y_W2_A0)
  Q1W <- (pred_W1_W2_A1_Y*pred_Y_W2_A1)/pred_W1_W2_A1
  Q0W <- (pred_W1_W2_A0_Y*pred_Y_W2_A0)/pred_W1_W2_A0
  QAW <- A*Q1W+(1-A)*Q0W

  # compute p(A|W)
  g1W <- (pred_W1_W2_A1*pred_A)/(pred_W1_W2_A1*pred_A+pred_W1_W2_A0*(1-pred_A))

  return(list(Q1W = Q1W,
              Q0W = Q0W,
              QAW = QAW,
              g1W = g1W))
}

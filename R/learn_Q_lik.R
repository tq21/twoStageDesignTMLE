#' Likelihood-based approach to estimate Q(A,W)
#' Assume W1 is one-dimensional binary for now, assume Y is binary
library(hal9001)
library(glmnet)

set.seed(123)
data <- sim_data(500)
W2 <- data[, c("W1", "W2")]
Delta <- data$Delta
W1 <- data$W3[Delta==1]
A <- data$A
Y <- data$Y

learn_Q_lik <- function(W1,
                        W2,
                        Delta,
                        A,
                        Y) {

  V <- as.matrix(cbind(W2, A=A, Y=Y)[Delta==1,])
  V_A0 <- as.matrix(cbind(W2, A=0, Y=Y))
  V_A0_Y1 <- as.matrix(cbind(W2, A=0, Y=1))
  V_A0_Y0 <- as.matrix(cbind(W2, A=0, Y=0))
  V_A1 <- as.matrix(cbind(W2, A=1, Y=Y))
  V_A1_Y0 <- as.matrix(cbind(W2, A=1, Y=0))
  V_A1_Y1 <- as.matrix(cbind(W2, A=1, Y=1))

  # estimate p(W1|V)
  W1_V_basis_list <- enumerate_basis(x = V, smoothness_orders = 1, max_degree = 3)
  W1_V_design <- make_design_matrix(X = V, blist = W1_V_basis_list)
  W1_V_A0_design <- make_design_matrix(X = V_A0, blist = W1_V_basis_list)
  W1_V_A0_Y0_design <- make_design_matrix(X = V_A0_Y0, blist = W1_V_basis_list)
  W1_V_A0_Y1_design <- make_design_matrix(X = V_A0_Y1, blist = W1_V_basis_list)
  W1_V_A1_design <- make_design_matrix(X = V_A1, blist = W1_V_basis_list)
  W1_V_A1_Y0_design <- make_design_matrix(X = V_A1_Y0, blist = W1_V_basis_list)
  W1_V_A1_Y1_design <- make_design_matrix(X = V_A1_Y1, blist = W1_V_basis_list)

  fit_W1_V <- cv.glmnet(x = W1_V_design,
                        y = W1,
                        alpha = 1,
                        family = "binomial")
  pred_W1_V_A0 <- as.numeric(predict(fit_W1_V, newx = W1_V_A0_design, s = "lambda.min", type = "response"))
  pred_W1_V_A0_Y0 <- as.numeric(predict(fit_W1_V, newx = W1_V_A0_Y0_design, s = "lambda.min", type = "response"))
  pred_W1_V_A0_Y1 <- as.numeric(predict(fit_W1_V, newx = W1_V_A0_Y1_design, s = "lambda.min", type = "response"))
  pred_W1_V_A1 <- as.numeric(predict(fit_W1_V, newx = W1_V_A1_design, s = "lambda.min", type = "response"))
  pred_W1_V_A1_Y0 <- as.numeric(predict(fit_W1_V, newx = W1_V_A1_Y0_design, s = "lambda.min", type = "response"))
  pred_W1_V_A1_Y1 <- as.numeric(predict(fit_W1_V, newx = W1_V_A1_Y1_design, s = "lambda.min", type = "response"))

  # estimate p(Y|W2,A)
  W2_A_basis_list <- enumerate_basis(x = as.matrix(cbind(W2,A=A)), smoothness_orders = 1, max_degree = 3)
  W2_A_design <- make_design_matrix(X = as.matrix(cbind(W2,A=A)), blist = W2_A_basis_list)
  W2_A0_design <- make_design_matrix(X = as.matrix(cbind(W2,A=0)), blist = W2_A_basis_list)
  W2_A1_design <- make_design_matrix(X = as.matrix(cbind(W2,A=1)), blist = W2_A_basis_list)

  fit_W2_A <- cv.glmnet(x = W2_A_design,
                        y = Y,
                        alpha = 1,
                        family = "binomial")
  pred_W2_A0 <- as.numeric(predict(fit_W2_A, newx = W2_A0_design, s = "lambda.min", type = "response"))
  pred_W2_A1 <- as.numeric(predict(fit_W2_A, newx = W2_A1_design, s = "lambda.min", type = "response"))

  # compute p(Y|W1,W2,A)
  Q1W <- (pred_W1_V_A1*pred_W2_A1)/(pred_W1_V_A1_Y0*(1-pred_W2_A1)+pred_W1_V_A1_Y1*pred_W2_A1)
  Q0W <- (pred_W1_V_A0*pred_W2_A0)/(pred_W1_V_A0_Y0*(1-pred_W2_A0)+pred_W1_V_A0_Y1*pred_W2_A0)
  QAW <- A*Q1W+(1-A)*Q0W

  # compute log-lik
  mean(Y*log(QAW)+(1-Y)*(1-QAW))

}

#' 3 default learners:
#' - glm inversely weighted by Delta/Pi
#' - MICE + glm
#' - likelihood-based HAL
library(origami)
library(furrr)
library(purrr)
estimate_risk <- function(nfolds = 5) {


  folds <- make_folds(n = length(Y), V = nfolds)

  V <- cbind(W2, A=A, Y=Y)
  AW <- cbind(A=A, W)
  A1W <- cbind(A=1, W)
  A0W <- cbind(A=0, W)

  map(folds, function(.fold) {
    train_idx <- .fold$training_set
    train_obs_idx <- train_idx[Delta[train_idx] == 1]
    valid_idx <- .fold$validation_set

    if (Q.family == "gaussian") {
      loss_fn <- mse_loss
    } else if (Q.family == "binomial") {
      loss_fn <- loglik_loss
    }

    # glm inversely weighted by Delta/Pi ---------------------------------------
    # estimate Pi
    fit_Pi <- glm(Delta[train_idx] ~ ., data = V[train_idx,], family = "binomial")
    Pi <- as.numeric(predict(fit_Pi, newdata = V[valid_idx,], type = "response"))
    wt <- Delta[train_idx_obs]/Pi[train_idx_obs]

    # estimate Q
    fit_Q_glm <- glm(Y[train_obs_idx] ~ ., data = W_A[train_obs_idx, ], weights = wt, family = Q.family)
    QAW_glm <- as.numeric(predict(fit_Q_glm, newdata = W_A[train_idx,], type = "response"))
    Q1W_glm <- as.numeric(predict(fit_Q_glm, newdata = W_A1[train_idx,], type = "response"))
    Q0W_glm <- as.numeric(predict(fit_Q_glm, newdata = W_A0[train_idx,], type = "response"))

    # estimate g
    fit_g_glm <- glm(A[train_obs_idx] ~ ., data = W[train_obs_idx, ], weights = wt, family = "binomial")
    g1W_glm <- as.numeric(predict(fit_g_glm, newdata = W[valid_idx,], type = "response"))

    # estimate risk and its regression
    Q_loss_glm <- loss_fn(QAW_glm, Y[train_idx])
    g_loss_glm <- loglik_loss(g1W_glm, A[train_idx])
    basis_list_V <- enumerate_basis(x = V[train_idx,], max_degree = 3, smoothness_orders = 1)
    design_V <- make_design_matrix(X = V[train_idx,], blist = basis_list_V)
    fit_Q_loss_glm <- cv.glmnet(x = design_V,
                                y = Q_loss_glm,
                                family = "gaussian",
                                alpha = 1)
    Q_loss_glm_reg <- as.numeric(predict(fit_Q_loss_glm, newx = design_V_valid, s = "lambda.min"))

    # MICE + glm ---------------------------------------------------------------
    Nimp <- 50
    W.stage2_aug <- as.data.frame(matrix(NA, nrow = length(Y), ncol = ncol(W.stage2)))
    W.stage2_aug[Delta.W == 1, ] <- W.stage2
    names(W.stage2_aug) <- names(W.stage2)
    df_mi <- data.frame(Y = Y,
                        A = A,
                        Delta.W = Delta.W,
                        W,
                        W.stage2_aug)
    init <- mice::mice(df_mi, maxit = 0)
    predM <- init$predictorMatrix
    imp <- mice::mice(df_mi, predictorMatrix = predM, m = Nimp, maxit = 20, print = FALSE)
    Q1W_all <- Q0W_all <- QAW_all <- g1W_all <- matrix(0, nrow = length(Y), ncol = Nimp)

    for (m in 1:Nimp) {
      comp <- mice::complete(imp, m)
      tmle_m <- tmle::tmle(Y = comp$Y,
                           A = comp$A,
                           W = cbind(comp[, colnames(W), drop = FALSE],
                                     comp[, colnames(W.stage2), drop = FALSE]),
                           family = Q.family,
                           g.SL.library = c("SL.glm"),
                           Q.SL.library = c("SL.glm"),
                           verbose = FALSE)
      Q1W_m <- tmle_m$Qinit$Q[, "Q1W"]
      Q0W_m <- tmle_m$Qinit$Q[, "Q0W"]
      QAW_m <- Q1W_m*A+Q0W_m*(1-A)
      g1W_m <- tmle_m$g$g1W

      Q1W_all[, m] <- Q1W_m; Q0W_all[, m] <- Q0W_m; QAW_all[, m] <- QAW_m
      g1W_all[, m] <- g1W_m
    }
    Q1W_MI <- rowMeans(Q1W_all)
    Q0W_MI <- rowMeans(Q0W_all)
    QAW_MI <- rowMeans(QAW_all)
    g1W_MI <- rowMeans(g1W_all)

    # estimate risk and its regression
    Q_loss_MI <- loss_fn(QAW_MI, Y[train_idx])
    g_loss_MI <- loglik_loss(g1W_MI, A[train_idx])
    basis_list_Delta_V <- enumerate_basis(x = Delta_V[train_idx,], max_degree = 3, smoothness_orders = 1)
    design_Delta_V <- make_design_matrix(X = Delta_V[train_idx,], blist = basis_list_V)
    fit_Q_loss_MI <- cv.glmnet(x = design_Delta_V,
                               y = Q_loss_MI,
                               family = "gaussian",
                               alpha = 1)
    fit_g_loss_MI <- cv.glmnet(x = design_Delta_V,
                               y = g_loss_MI,
                               family = "gaussian",
                               alpha = 1)
    Q_loss_MI_reg <- as.numeric(predict(fit_Q_loss_MI, newx = design_Delta1_V_valid, s = "lambda.min"))
    g_loss_MI_reg <- as.numeric(predict(fit_g_loss_MI, newx = design_Delta1_V_valid, s = "lambda.min"))

    # likelihood-based HAL -----------------------------------------------------
    lik_list <- learn_lik(W1,
                          W2,
                          Delta,
                          A,
                          Y)
    Q1W_lik <- lik_list$Q1W
    Q0W_lik <- lik_list$Q0W
    QAW_lik <- lik_list$QAW
    g1W_lik <- lik_list$g1W

    QAW <- cbind(QAW_glm, QAW_MI, QAW_lik)
    Q1W <- cbind(Q1W_glm, Q1W_MI, Q1W_lik)
    Q0W <- cbind(Q0W_glm, Q0W_MI, Q0W_lik)
    g1W <- cbind(g1W_glm, g1W_MI, g1W_lik)

    # estimate risk and its regression
    Q_loss_lik <- loss_fn(QAW_lik, Y[train_idx])
    g_loss_lik <- loglik_loss(g1W_lik, A[train_idx])
    fit_Q_loss_lik <- cv.glmnet(x = design_Delta_V,
                                y = Q_loss_lik,
                                family = "gaussian",
                                alpha = 1)
    fit_g_loss_lik <- cv.glmnet(x = design_Delta_V,
                                y = g_loss_lik,
                                family = "gaussian",
                                alpha = 1)
    Q_loss_lik_reg <- as.numeric(predict(fit_Q_loss_lik, newx = design_Delta1_V_valid, s = "lambda.min"))
    g_loss_lik_reg <- as.numeric(predict(fit_g_loss_lik, newx = design_Delta1_V_valid, s = "lambda.min"))

    return(data.frame(Q_loss_glm,
                      Q_loss_glm_reg,
                      Q_loss_MI,
                      Q_loss_MI_reg,
                      Q_loss_lik,
                      Q_loss_lik_reg,
                      g_loss_glm,
                      g_loss_glm_reg,
                      g_loss_MI,
                      g_loss_MI_reg,
                      g_loss_lik,
                      g_loss_lik_reg))
  })
}

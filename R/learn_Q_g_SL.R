#' 3 default learners:
#' - glm inversely weighted by Delta/Pi
#' - MICE + glm
#' - likelihood-based HAL
library(origami)
library(purrr)
sq_err_loss <- function(pred, truth) return((pred-truth)^2)
log_lik_loss <- function(pred, truth) return(-(truth*log(pred)+(1-truth)*log(1-pred)))
learn_Q_g_SL <- function(W1,
                         W2,
                         W,
                         A,
                         Y,
                         Delta,
                         Q.family,
                         nfolds = 5) {

  folds <- make_folds(n = length(Y), V = nfolds, strata_ids = Delta)
  W_names <- c(names(W1), names(W2))

  V <- cbind(W2, A=A, Y=Y)
  W_A <- cbind(W, A=A)
  W_A0 <- cbind(W, A=0)
  W_A1 <- cbind(W, A=1)

  df <- map_dfr(folds, function(.fold) {
    print(.fold$v)
    train_idx <- .fold$training_set
    train_obs_idx <- train_idx[Delta[train_idx] == 1]
    valid_idx <- .fold$validation_set
    valid_obs_idx <- valid_idx[Delta[valid_idx] == 1]

    if (Q.family == "gaussian") {
      loss_fn <- sq_err_loss
    } else if (Q.family == "binomial") {
      loss_fn <- log_lik_loss
    }

    # glm inversely weighted by Delta/Pi ---------------------------------------
    # estimate Pi
    fit_Pi <- glm(Delta[train_idx] ~ ., data = as.data.frame(V[train_idx,]), family = "binomial")
    Pi <- as.numeric(predict(fit_Pi, type = "response"))
    Pi_valid <- as.numeric(predict(fit_Pi, newdata = as.data.frame(V[valid_idx,]), type = "response"))
    Pi_valid <- Pi_valid[Delta[valid_idx] == 1]
    wt <- 1/Pi[Delta[train_idx] == 1]

    # estimate Q
    fit_Q_glm <- suppressWarnings(glm(Y[train_obs_idx] ~ ., data = as.data.frame(W_A[train_obs_idx,]), weights = wt, family = Q.family))
    QAW_train_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A[train_obs_idx,]), type = "response"))
    QAW_valid_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A[valid_obs_idx,]), type = "response"))
    Q0W_train_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A0[train_obs_idx,]), type = "response"))
    Q0W_valid_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A0[valid_obs_idx,]), type = "response"))
    Q1W_train_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A1[train_obs_idx,]), type = "response"))
    Q1W_valid_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A1[valid_obs_idx,]), type = "response"))

    # estimate g
    fit_g_glm <- suppressWarnings(glm(A[train_obs_idx] ~ ., data = as.data.frame(W[train_obs_idx,]), weights = wt, family = "binomial"))
    g1W_train_glm <- as.numeric(predict(fit_g_glm, newdata = as.data.frame(W[train_obs_idx,]), type = "response"))
    g1W_valid_glm <- as.numeric(predict(fit_g_glm, newdata = as.data.frame(W[valid_obs_idx,]), type = "response"))

    # estimate risk and its regression
    Q_loss_train_glm <- loss_fn(QAW_train_glm, Y[train_obs_idx])
    Q_loss_valid_glm <- loss_fn(QAW_valid_glm, Y[valid_obs_idx])
    g_loss_train_glm <- log_lik_loss(g1W_train_glm, A[train_obs_idx])
    g_loss_valid_glm <- log_lik_loss(g1W_valid_glm, A[valid_obs_idx])
    V_basis_list <- enumerate_basis(x = V[train_obs_idx,], max_degree = 3, smoothness_orders = 1)
    V_design_train <- make_design_matrix(X = V[train_obs_idx,], blist = V_basis_list)
    V_design_valid <- make_design_matrix(X = V[valid_obs_idx,], blist = V_basis_list)
    fit_Q_loss_glm <- cv.glmnet(x = V_design_train,
                                y = Q_loss_train_glm,
                                family = "gaussian",
                                alpha = 1)
    fit_g_loss_glm <- cv.glmnet(x = V_design_train,
                                y = g_loss_train_glm,
                                family = "gaussian",
                                alpha = 1)
    Q_loss_reg_valid_glm <- as.numeric(predict(fit_Q_loss_glm, newx = V_design_valid, s = "lambda.min"))
    g_loss_reg_valid_glm <- as.numeric(predict(fit_g_loss_glm, newx = V_design_valid, s = "lambda.min"))

    # MICE + glm ---------------------------------------------------------------
    Nimp <- 50
    df_mi <- data.frame(W, A=A, Y=Y, Delta=Delta)[train_idx,]
    init <- mice::mice(df_mi, maxit = 0)
    predM <- init$predictorMatrix
    imp <- suppressWarnings(mice::mice(df_mi, predictorMatrix = predM, m = Nimp, maxit = 20, print = FALSE))
    Q1W_train_MI_all <- Q0W_train_MI_all <- QAW_train_MI_all <- g1W_train_MI_all <- matrix(0, nrow = length(train_idx), ncol = Nimp)
    Q1W_valid_MI_all <- Q0W_valid_MI_all <- QAW_valid_MI_all <- g1W_valid_MI_all <- matrix(0, nrow = length(valid_obs_idx), ncol = Nimp)

    for (m in 1:Nimp) {
      comp <- mice::complete(imp, m)
      W_A_m <- comp[, c(W_names, "A")]
      W_A0_m <- cbind(comp[, W_names], A=0)
      W_A1_m <- cbind(comp[, W_names], A=1)

      # estimate Q
      fit_Q_m <- glm(comp$Y ~ ., data = as.data.frame(W_A_m), family = Q.family)
      QAW_train_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A_m), type = "response"))
      Q0W_train_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A0_m), type = "response"))
      Q1W_train_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A1_m), type = "response"))
      QAW_valid_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A[valid_obs_idx,]), type = "response"))
      Q0W_valid_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A0[valid_obs_idx,]), type = "response"))
      Q1W_valid_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A1[valid_obs_idx,]), type = "response"))

      # estimate g
      W_m <- comp[, W_names, drop=FALSE]
      fit_g_m <- glm(comp$A ~ ., data = as.data.frame(W_m), family = "binomial")
      g1W_train_m <- as.numeric(predict(fit_g_m, newdata = as.data.frame(W_m), type = "response"))
      g1W_valid_m <- as.numeric(predict(fit_g_m, newdata = as.data.frame(W[valid_obs_idx,,drop=FALSE]), type = "response"))

      QAW_train_MI_all[, m] <- QAW_train_m
      Q0W_train_MI_all[, m] <- Q0W_train_m
      Q1W_train_MI_all[, m] <- Q1W_train_m
      QAW_valid_MI_all[, m] <- QAW_valid_m
      Q0W_valid_MI_all[, m] <- Q0W_valid_m
      Q1W_valid_MI_all[, m] <- Q1W_valid_m
      g1W_train_MI_all[, m] <- g1W_train_m
      g1W_valid_MI_all[, m] <- g1W_valid_m
    }
    QAW_train_MI <- rowMeans(QAW_train_MI_all)
    Q0W_train_MI <- rowMeans(Q0W_train_MI_all)
    Q1W_train_MI <- rowMeans(Q1W_train_MI_all)
    QAW_valid_MI <- rowMeans(QAW_valid_MI_all)
    Q0W_valid_MI <- rowMeans(Q0W_valid_MI_all)
    Q1W_valid_MI <- rowMeans(Q1W_valid_MI_all)
    g1W_train_MI <- rowMeans(g1W_train_MI_all)
    g1W_valid_MI <- rowMeans(g1W_valid_MI_all)

    # estimate risk and its regression
    Q_loss_train_MI <- loss_fn(QAW_train_MI, Y[train_idx])
    Q_loss_valid_MI <- loss_fn(QAW_valid_MI, Y[valid_obs_idx])
    g_loss_train_MI <- log_lik_loss(g1W_train_MI, A[train_idx])
    g_loss_valid_MI <- log_lik_loss(g1W_valid_MI, A[valid_obs_idx])
    V_Delta_basis_list <- enumerate_basis(x = as.matrix(cbind(V, Delta=Delta)[train_idx,]), max_degree = 3, smoothness_orders = 1)
    V_Delta_design_train <- make_design_matrix(X = as.matrix(cbind(V, Delta=Delta)[train_idx,]), blist = V_Delta_basis_list)
    V_Delta1_design_train <- make_design_matrix(X = as.matrix(cbind(V, Delta=1)[train_idx,]), blist = V_Delta_basis_list)
    V_Delta_design_valid <- make_design_matrix(X = as.matrix(cbind(V, Delta=Delta)[valid_idx,]), blist = V_Delta_basis_list)
    V_Delta1_design_valid <- make_design_matrix(X = as.matrix(cbind(V, Delta=1)[valid_idx,]), blist = V_Delta_basis_list)
    fit_Q_loss_MI <- cv.glmnet(x = V_Delta_design_train,
                               y = Q_loss_train_MI,
                               family = "gaussian",
                               alpha = 1)
    fit_g_loss_MI <- cv.glmnet(x = V_Delta_design_train,
                               y = g_loss_train_MI,
                               family = "gaussian",
                               alpha = 1)
    Q_loss_reg_valid_MI <- as.numeric(predict(fit_Q_loss_MI, newx = V_Delta1_design_valid, s = "lambda.min"))
    Q_loss_reg_valid_MI <- Q_loss_reg_valid_MI[Delta[valid_idx] == 1]
    g_loss_reg_valid_MI <- as.numeric(predict(fit_g_loss_MI, newx = V_Delta1_design_valid, s = "lambda.min"))
    g_loss_reg_valid_MI <- g_loss_reg_valid_MI[Delta[valid_idx] == 1]

    # likelihood-based HAL -----------------------------------------------------
    lik_list <- learn_lik(W1 = W1[train_idx,,drop=FALSE],
                          W1_valid = W1[valid_idx,,drop=FALSE],
                          W2 = W2[train_idx,,drop=FALSE],
                          W2_valid = W2[valid_idx,,drop=FALSE],
                          Delta = Delta[train_idx],
                          Delta_valid = Delta[valid_idx],
                          A = A[train_idx],
                          A_valid = A[valid_idx],
                          Y = Y[train_idx],
                          Y_valid = Y[valid_idx])
    Q1W_train_lik <- lik_list$Q1W
    Q0W_train_lik <- lik_list$Q0W
    QAW_train_lik <- lik_list$QAW
    g1W_train_lik <- lik_list$g1W
    Q1W_valid_lik <- lik_list$Q1W_valid
    Q0W_valid_lik <- lik_list$Q0W_valid
    QAW_valid_lik <- lik_list$QAW_valid
    g1W_valid_lik <- lik_list$g1W_valid

    # estimate risk and its regression
    Q_loss_train_lik <- loss_fn(QAW_train_lik, Y[train_idx])
    Q_loss_valid_lik <- loss_fn(QAW_valid_lik, Y[valid_idx])
    Q_loss_valid_lik <- Q_loss_valid_lik[Delta[valid_idx] == 1]
    g_loss_train_lik <- log_lik_loss(g1W_train_lik, A[train_idx])
    g_loss_valid_lik <- log_lik_loss(g1W_valid_lik, A[valid_idx])
    g_loss_valid_lik <- g_loss_valid_lik[Delta[valid_idx] == 1]
    fit_Q_loss_lik <- cv.glmnet(x = V_Delta_design_train,
                                y = Q_loss_train_lik,
                                family = "gaussian",
                                alpha = 1)
    fit_g_loss_lik <- cv.glmnet(x = V_Delta_design_train,
                                y = g_loss_train_lik,
                                family = "gaussian",
                                alpha = 1)
    Q_loss_reg_valid_lik <- as.numeric(predict(fit_Q_loss_lik, newx = V_Delta1_design_valid, s = "lambda.min"))
    Q_loss_reg_valid_lik <- Q_loss_reg_valid_lik[Delta[valid_idx] == 1]
    g_loss_reg_valid_lik <- as.numeric(predict(fit_g_loss_lik, newx = V_Delta1_design_valid, s = "lambda.min"))
    g_loss_reg_valid_lik <- g_loss_reg_valid_lik[Delta[valid_idx] == 1]

    return(data.frame(Pi_valid,
                      Q_loss_valid_glm,
                      Q_loss_reg_valid_glm,
                      Q_loss_valid_MI,
                      Q_loss_reg_valid_MI,
                      Q_loss_valid_lik,
                      Q_loss_reg_valid_lik,
                      g_loss_valid_glm,
                      g_loss_reg_valid_glm,
                      g_loss_valid_MI,
                      g_loss_reg_valid_MI,
                      g_loss_valid_lik,
                      g_loss_reg_valid_lik))
  })

  # target En(L(Qn)|Delta=1,V)
  clever_cov <- 1/df$Pi_valid
  eps_Q_glm <- coef(glm(df$Q_loss_valid_glm ~ -1+offset(df$Q_loss_reg_valid_glm)+clever_cov, family = "gaussian"))
  eps_Q_glm[is.na(eps_Q_glm)] <- 0
  df$Q_loss_reg_valid_glm <- df$Q_loss_reg_valid_glm+eps_Q_glm*clever_cov

  eps_Q_MI <- coef(glm(df$Q_loss_valid_MI ~ -1+offset(df$Q_loss_reg_valid_MI)+clever_cov, family = "gaussian"))
  eps_Q_MI[is.na(eps_Q_MI)] <- 0
  df$Q_loss_reg_valid_MI <- df$Q_loss_reg_valid_MI+eps_Q_MI*clever_cov

  eps_Q_lik <- coef(glm(df$Q_loss_valid_lik ~ -1+offset(df$Q_loss_reg_valid_lik)+clever_cov, family = "gaussian"))
  eps_Q_lik[is.na(eps_Q_lik)] <- 0
  df$Q_loss_reg_valid_lik <- df$Q_loss_reg_valid_lik+eps_Q_lik*clever_cov

  eps_g_glm <- coef(glm(df$g_loss_valid_glm ~ -1+offset(df$g_loss_reg_valid_glm)+clever_cov, family = "gaussian"))
  eps_g_glm[is.na(eps_g_glm)] <- 0
  df$g_loss_reg_valid_glm <- df$g_loss_reg_valid_glm+eps_g_glm*clever_cov

  eps_g_MI <- coef(glm(df$g_loss_valid_MI ~ -1+offset(df$g_loss_reg_valid_MI)+clever_cov, family = "gaussian"))
  eps_g_MI[is.na(eps_g_MI)] <- 0
  df$g_loss_reg_valid_MI <- df$g_loss_reg_valid_MI+eps_g_MI*clever_cov

  eps_g_lik <- coef(glm(df$g_loss_valid_lik ~ -1+offset(df$g_loss_reg_valid_lik)+clever_cov, family = "gaussian"))
  eps_g_lik[is.na(eps_g_lik)] <- 0
  df$g_loss_reg_valid_lik <- df$g_loss_reg_valid_lik+eps_g_lik*clever_cov

  Q_loss_reg_valid <- cbind(df$Q_loss_reg_valid_glm, df$Q_loss_reg_valid_MI, df$Q_loss_reg_valid_lik)
  g_loss_reg_valid <- cbind(df$g_loss_reg_valid_glm, df$g_loss_reg_valid_MI, df$g_loss_reg_valid_lik)
  Q_choice <- which.min(colMeans(Q_loss_reg_valid))
  g_choice <- which.min(colMeans(g_loss_reg_valid))

  # estimate Q using all folds
  # glm inversely weighted by Delta/Pi -----------------------------------------
  fit_Pi <- glm(Delta ~ ., data = as.data.frame(V), family = "binomial")
  Pi <- as.numeric(predict(fit_Pi, type = "response"))
  wt <- 1/Pi[Delta == 1]
  fit_Q_glm <- suppressWarnings(glm(Y[Delta == 1] ~ ., data = as.data.frame(W_A[Delta == 1,]), weights = wt, family = Q.family))
  fit_g_glm <- suppressWarnings(glm(A[Delta == 1] ~ ., data = as.data.frame(W[Delta == 1,]), weights = wt, family = "binomial"))
  QAW_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A[Delta == 1,]), type = "response"))
  Q0W_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A0[Delta == 1,]), type = "response"))
  Q1W_glm <- as.numeric(predict(fit_Q_glm, newdata = as.data.frame(W_A1[Delta == 1,]), type = "response"))
  g1W_glm <- as.numeric(predict(fit_g_glm, newdata = as.data.frame(W[Delta == 1,]), type = "response"))

  # MICE + glm -----------------------------------------------------------------
  Nimp <- 50
  df_mi <- data.frame(W, A=A, Y=Y, Delta=Delta)
  init <- mice::mice(df_mi, maxit = 0)
  predM <- init$predictorMatrix
  imp <- suppressWarnings(mice::mice(df_mi, predictorMatrix = predM, m = Nimp, maxit = 20, print = FALSE))
  Q1W_all <- Q0W_all <- QAW_all <- g1W_all <- matrix(0, nrow = length(Y), ncol = Nimp)

  for (m in 1:Nimp) {
    comp <- mice::complete(imp, m)
    W_A_m <- comp[, c(W_names, "A")]
    W_A0_m <- cbind(comp[, W_names], A=0)
    W_A1_m <- cbind(comp[, W_names], A=1)

    # estimate Q
    fit_Q_m <- glm(comp$Y ~ ., data = as.data.frame(W_A_m), family = Q.family)
    QAW_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A_m), type = "response"))
    Q0W_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A0_m), type = "response"))
    Q1W_m <- as.numeric(predict(fit_Q_m, newdata = as.data.frame(W_A1_m), type = "response"))

    # estimate g
    W_m <- comp[, W_names, drop=FALSE]
    fit_g_m <- glm(comp$A ~ ., data = as.data.frame(W_m), family = "binomial")
    g1W_m <- as.numeric(predict(fit_g_m, newdata = as.data.frame(W_m), type = "response"))

    QAW_all[, m] <- QAW_m
    Q0W_all[, m] <- Q0W_m
    Q1W_all[, m] <- Q1W_m
    g1W_all[, m] <- g1W_m
  }
  QAW_MI <- rowMeans(QAW_all)[Delta == 1]
  Q0W_MI <- rowMeans(Q0W_all)[Delta == 1]
  Q1W_MI <- rowMeans(Q1W_all)[Delta == 1]
  g1W_MI <- rowMeans(g1W_all)[Delta == 1]

  # likelihood-based HAL -------------------------------------------------------
  lik_list <- learn_lik(W1 = W1,
                        W2 = W2,
                        Delta = Delta,
                        A = A,
                        Y = Y)
  QAW_lik <- lik_list$QAW[Delta == 1]
  Q0W_lik <- lik_list$Q0W[Delta == 1]
  Q1W_lik <- lik_list$Q1W[Delta == 1]
  g1W_lik <- lik_list$g1W[Delta == 1]

  if (Q_choice == 1) {
    QAW <- QAW_glm
    Q0W <- Q0W_glm
    Q1W <- Q1W_glm
  } else if (Q_choice == 2) {
    QAW <- QAW_MI
    Q0W <- Q0W_MI
    Q1W <- Q1W_MI
  } else if (Q_choice == 3) {
    QAW <- QAW_lik
    Q0W <- Q0W_lik
    Q1W <- Q1W_lik
  }

  if (g_choice == 1) {
    g1W <- g1W_glm
  } else if (g_choice == 2) {
    g1W <- g1W_MI
  } else if (g_choice == 3) {
    g1W <- g1W_lik
  }

  return(list(QAW = QAW,
              Q0W = Q0W,
              Q1W = Q1W,
              g1W = g1W))
}

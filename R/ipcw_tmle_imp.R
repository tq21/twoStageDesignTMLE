library(origami)
library(hal9001)
library(glmnet)
ipcw_tmle_imp <- function(data,
                          W_name,
                          A_name,
                          Y_name,
                          Delta_name,
                          V_name,
                          Pi = NULL,
                          Q = NULL,
                          g1W = NULL,
                          nfolds = 10,
                          Pi_method = "HAL",
                          Q_method = "HAL",
                          g_method = "HAL",
                          m_method = "HAL",
                          g_bounds = NULL,
                          family = "binomial", # FOR NOW.
                          DFullNCReg_pool = FALSE,
                          enumerate_basis_args = list(max_degree = 2,
                                                      smoothness_orders = 1,
                                                      num_knots = 5),
                          parallel = TRUE,
                          browse = FALSE) {
  if (browse) browser()

  # data pre-processing --------------------------------------------------------
  # extract variables
  A <- data[[A_name]]
  Y <- data[[Y_name]]
  Delta <- data[[Delta_name]]
  data_A0 <- data; data_A0[[A_name]] <- 0
  data_A1 <- data; data_A1[[A_name]] <- 1

  # cross-validation scheme
  folds <- make_folds(n = nrow(data), V = nfolds,
                      strata_ids = paste0(Delta, "-", A, "-", Y))
  foldid <- folds2foldvec(folds)
  phase_2_idx <- which(Delta == 1)
  folds_obs <- make_folds(n = sum(Delta), V = nfolds,
                          strata_ids = paste0(A[Delta == 1], "-", Y[Delta == 1]))
  foldid_obs <- folds2foldvec(folds_obs)

  # estimate Pi(V)=P(Delta=1|V)
  if (is.null(Pi)) {
    if (is.character(Pi_method) && Pi_method == "sl3") {
      Pi_method <- get_default_sl3_learners("binomial")
    }
    if (is.list(Pi_method)) {
      Pi_task <- sl3_Task$new(data = data, covariates = V_name,
                              outcome = Delta_name, folds = folds,
                              outcome_type = "binomial")
      Pi_lrnr_stack <- Stack$new(Pi_SL_lib)
      Pi_lrnr <- make_learner(Pipeline, Lrnr_cv$new(Pi_lrnr_stack),
                              Lrnr_cv_selector$new(loss_loglik_binomial))
      suppressMessages(Pi_fit <- Pi_lrnr$train(Pi_task))
      Pi <- Pi_fit$predict(Pi_task)
    } else if (Pi_method == "HAL") {
      basis_list <- enumerate_basis(
        x = as.matrix(data[, V_name, drop=FALSE]),
        max_degree = enumerate_basis_args$max_degree,
        smoothness_orders = enumerate_basis_args$smoothness_orders,
        num_knots = enumerate_basis_args$num_knots)
      hal_design <- make_design_matrix(
        X = as.matrix(data[, V_name, drop=FALSE]),
        blist = basis_list)
      Pi_fit <- cv.glmnet(x = hal_design, y = Delta, family = "binomial",
                          alpha = 1, foldid = foldid, parallel = parallel)
      Pi <- as.numeric(predict(Pi_fit, newx = hal_design, s = "lambda.min",
                               type = "response"))
    }
  }

  data$wt <- Delta/Pi; data_A0$wt <- Delta/Pi; data_A1$wt <- Delta/Pi

  # estimate full data Q(A,W)=E(Y|A,W), and g(W)=P(A=1|W)
  if (is.null(Q)) {
    if (is.character(Q_method) && Q_method == "sl3") {
      Q_method <- get_default_sl3_learners(family)
    }
    if (is.list(Q_method)) {
      Q_task <- sl3_Task$new(data = data[phase_2_idx, , drop=FALSE],
                             covariates = c(W_name, A_name), outcome = Y_name,
                             weights = "wt", folds = folds_obs,
                             outcome_type = family)
      Q0_task <- sl3_Task$new(data = data_A0[phase_2_idx, , drop=FALSE],
                              covariates = c(W_name, A_name), outcome = Y_name,
                              weights = "wt", folds = folds_obs,
                              outcome_type = family)
      Q1_task <- sl3_Task$new(data = data_A1[phase_2_idx, , drop=FALSE],
                              covariates = c(W_name, A_name), outcome = Y_name,
                              weights = "wt", folds = folds_obs,
                              outcome_type = family)
      Q_lrnr_stack <- Stack$new(Q_SL_lib)
      if (family == "binomial") {
        Q_lrnr <- make_learner(Pipeline, Lrnr_cv$new(Q_lrnr_stack),
                               Lrnr_cv_selector$new(loss_loglik_binomial))
      } else if (family == "gaussian") {
        Q_lrnr <- make_learner(Pipeline, Lrnr_cv$new(Q_lrnr_stack),
                               Lrnr_cv_selector$new(loss_squared_error))
      }
      suppressMessages(Q_fit <- Q_lrnr$train(Q_task))
      QAW <- Q0W <- Q1W <- numeric(nrow(data))
      QAW[Delta == 1] <- Q_fit$predict(Q_task)
      Q0W[Delta == 1] <- Q_fit$predict(Q0_task)
      Q1W[Delta == 1] <- Q_fit$predict(Q1_task)
      Q <- data.frame(QAW, Q0W, Q1W)
    } else if (Q_method == "HAL") {
      basis_list <- enumerate_basis(
        x = as.matrix(data[phase_2_idx, c(W_name, A_name), drop=FALSE]),
        max_degree = enumerate_basis_args$max_degree,
        smoothness_orders = enumerate_basis_args$smoothness_orders,
        num_knots = enumerate_basis_args$num_knots)
      hal_design <- make_design_matrix(
        X = as.matrix(data[phase_2_idx, c(W_name, A_name), drop=FALSE]),
        blist = basis_list)
      hal_design_A0 <- make_design_matrix(
        X = as.matrix(data_A0[phase_2_idx, c(W_name, A_name), drop=FALSE]),
        blist = basis_list)
      hal_design_A1 <- make_design_matrix(
        X = as.matrix(data_A1[phase_2_idx, c(W_name, A_name), drop=FALSE]),
        blist = basis_list)
      Q_fit <- cv.glmnet(x = hal_design, y = Y[phase_2_idx], family = family,
                         alpha = 1, foldid = foldid_obs, parallel = parallel)
      QAW <- Q0W <- Q1W <- numeric(nrow(data))
      QAW[Delta == 1] <- as.numeric(predict(Q_fit, newx = hal_design,
                                            s = "lambda.min",
                                            type = "response"))
      Q0W[Delta == 1] <- as.numeric(predict(Q_fit, newx = hal_design_A0,
                                            s = "lambda.min",
                                            type = "response"))
      Q1W[Delta == 1] <- as.numeric(predict(Q_fit, newx = hal_design_A1,
                                            s = "lambda.min",
                                            type = "response"))
      Q <- data.frame(QAW, Q0W, Q1W)
    }
  }

  if (is.null(g1W)) {
    if (is.character(g_method) && g_method == "sl3") {
      g_method <- get_default_sl3_learners("binomial")
    }
    if (is.list(g_method)) {
      g1W_task <- sl3_Task$new(data = data[phase_2_idx, , drop=FALSE],
                               covariates = W_name, outcome = A_name,
                               weights = "wt", folds = folds_obs,
                               outcome_type = "binomial")
      g1W_lrnr_stack <- Stack$new(g1W_SL_lib)
      g1W_lrnr <- make_learner(Pipeline, Lrnr_cv$new(g1W_lrnr_stack),
                               Lrnr_cv_selector$new(loss_loglik_binomial))
      suppressMessages(g1W_fit <- g1W_lrnr$train(g1W_task))
      g1W <- numeric(nrow(data))
      g1W[Delta == 1] <- g1W_fit$predict(g1W_task)
    } else if (g_method == "HAL") {
      basis_list <- enumerate_basis(
        x = as.matrix(data[phase_2_idx, W_name, drop=FALSE]),
        max_degree = enumerate_basis_args$max_degree,
        smoothness_orders = enumerate_basis_args$smoothness_orders,
        num_knots = enumerate_basis_args$num_knots)
      hal_design <- make_design_matrix(
        X = as.matrix(data[phase_2_idx, W_name, drop=FALSE]),
        blist = basis_list)
      g1W_fit <- cv.glmnet(x = hal_design, y = A[phase_2_idx],
                           family = "binomial", alpha = 1, foldid = foldid_obs,
                           parallel = parallel)
      g1W <- numeric(nrow(data))
      g1W[Delta == 1] <- as.numeric(predict(g1W_fit, newx = hal_design,
                                            s = "lambda.min", type = "response"))
    }
  }

  # apply truncation
  if (is.null(g_bounds)) {
    g_bounds <- c(5/sqrt(sum(Delta))/log(sum(Delta)), 1-5/sqrt(sum(Delta))/log(sum(Delta)))
  }
  g1W <- .bound(g1W, g_bounds); g1W[Delta == 0] <- 0

  # IPCW-TMLE (TODO: gaussian Y)
  HAW <- ifelse(Delta == 1, A/g1W-(1-A)/(1-g1W), 0)
  H1W <- ifelse(Delta == 1, 1/g1W, 0)
  H0W <- ifelse(Delta == 1, -1/(1-g1W), 0)
  DFullNC <- HAW*(Y-Q$QAW)+Q$Q1W-Q$Q0W; K <- Delta/Pi
  fluc <- as.numeric(
    coef(glm(Y[phase_2_idx] ~ -1+offset(qlogis(Q$QAW[phase_2_idx]))+HAW[phase_2_idx],
             weights = K[phase_2_idx], family = "quasibinomial")))
  if (is.na(fluc)) fluc <- 0
  Q1W_star <- plogis(qlogis(Q$Q1W)+fluc*H1W)
  Q0W_star <- plogis(qlogis(Q$Q0W)+fluc*H0W)
  QAW_star <- Q1W_star*A+Q0W_star*(1-A)
  psi_ipcw_tmle <- weighted.mean(Q1W_star[phase_2_idx]-Q0W_star[phase_2_idx], w = K[phase_2_idx])
  eic_ipcw_tmle <- K*(Q1W_star-Q0W_star-psi_ipcw_tmle+HAW*(Y-QAW_star))
  se_ipcw_tmle <- sqrt(var(eic_ipcw_tmle)/length(eic_ipcw_tmle))
  lower_ipcw_tmle <- psi_ipcw_tmle+qnorm(0.025)*se_ipcw_tmle
  upper_ipcw_tmle <- psi_ipcw_tmle+qnorm(0.975)*se_ipcw_tmle

  # estimate full-data non-centered EIC regression m=E(DFullNC|Delta=1,V)
  data$DFullNC <- DFullNC
  if (DFullNCReg_pool) {
    # DFullNC ~ Delta + V on all observations, evaluate at Delta = 1
    DFullNC_task <- sl3_Task$new(data = data,
                                 covariates = c(Delta_name, V_name),
                                 outcome = "DFullNC", outcome_type = "gaussian")
    DFullNC_pred_task <- sl3_Task$new(data = data_Delta1,
                                      covariates = c(Delta_name, V_name),
                                      outcome = "DFullNC",
                                      outcome_type = "gaussian")
  } else {
    # DFullNC ~ V on Delta = 1 observations
    if (is.character(m_method) && m_method == "sl3") {
      m_method <- get_default_sl3_learners("gaussian")
    }
    if (is.list(m_method)) {
      DFullNC_task <- sl3_Task$new(data = data[phase_2_idx, , drop=FALSE],
                                   covariates = V_name, outcome = "DFullNC",
                                   folds = folds_obs, outcome_type = "gaussian")
      DFullNC_pred_task <- sl3_Task$new(data = data, covariates = V_name,
                                        outcome = "DFullNC", folds = folds,
                                        outcome_type = "gaussian")
      DFullNCReg_lrnr_stack <- Stack$new(m_method)
      DFullNCReg_lrnr <- make_learner(Pipeline, Lrnr_cv$new(DFullNCReg_lrnr_stack),
                                      Lrnr_cv_selector$new(loss_squared_error))
      suppressMessages(DFullNCReg_fit_sl <- DFullNCReg_lrnr$train(DFullNC_task))
      DFullNCReg <- DFullNCReg_fit_sl$predict(DFullNC_pred_task)
    } else if (m_method == "HAL") {
      basis_list <- enumerate_basis(
        x = as.matrix(data[phase_2_idx, V_name, drop=FALSE]),
        max_degree = enumerate_basis_args$max_degree,
        smoothness_orders = enumerate_basis_args$smoothness_orders,
        num_knots = enumerate_basis_args$num_knots)
      hal_design <- make_design_matrix(
        X = as.matrix(data[, V_name, drop=FALSE]), blist = basis_list)
      m_fit <- cv.glmnet(x = hal_design[phase_2_idx, , drop=FALSE],
                         y = DFullNC[phase_2_idx], family = "gaussian",
                         alpha = 1, foldid = foldid_obs, parallel = parallel)
      DFullNCReg <- as.numeric(predict(m_fit, newx = hal_design,
                                       s = "lambda.min", type = "response"))
    }
  }

  # TMLE -----------------------------------------------------------------------
  fn <- function(args) {
    eps <- args[1]
    gamma <- args[2]
    QAW_eps <- plogis(qlogis(Q$QAW)+eps*HAW)
    Q0W_eps <- plogis(qlogis(Q$Q0W)+eps*H0W)
    Q1W_eps <- plogis(qlogis(Q$Q1W)+eps*H1W)
    DFullNC_eps <- HAW*(Y-QAW_eps)+Q1W_eps-Q0W_eps
    DFullNCReg_gamma <- DFullNCReg+gamma*K
    psi <- mean(K*(Q1W_eps-Q0W_eps))
    score_1 <- mean(K*(DFullNC_eps-DFullNCReg_gamma))
    score_2 <- mean(DFullNCReg_gamma)-psi
    return(c(score_1, score_2))
  }

  # solve non-linear system
  sol <- nleqslv(x = c(0, 0), fn = fn, method = "Newton")
  eps <- sol$x[1]; gamma <- sol$x[2]

  # updates
  QAW_star <- plogis(qlogis(Q$QAW)+eps*HAW)
  Q0W_star <- plogis(qlogis(Q$Q0W)+eps*H0W)
  Q1W_star <- plogis(qlogis(Q$Q1W)+eps*H1W)
  DFullNC_star <- HAW*(Y-QAW_star)+Q1W_star-Q0W_star
  DFullNCReg_star <- DFullNCReg+gamma*K

  # point estimate and inference
  psi_tmle <- mean(DFullNCReg_star)
  eic_tmle <- as.numeric(K*(DFullNC_star-DFullNCReg_star)+DFullNCReg_star-psi_tmle)
  se_tmle <- sqrt(var(eic_tmle)/nrow(data))
  lower_tmle <- psi_tmle+qnorm(0.025)*se_tmle
  upper_tmle <- psi_tmle+qnorm(0.975)*se_tmle

  # check PnEIC
  PnEIC_1 <- mean(K*(DFullNC_star-DFullNCReg_star))
  PnEIC_2 <- mean(DFullNCReg_star)-psi_tmle

  # A-IPCW ---------------------------------------------------------------------
  psi_aipcw <- mean(K*(DFullNC-DFullNCReg)+DFullNCReg)
  eic_aipcw <- K*(DFullNC-DFullNCReg)+DFullNCReg-psi_aipcw
  se_aipcw <- sqrt(var(eic_aipcw)/nrow(data))
  lower_aipcw <- psi_aipcw+qnorm(0.025)*se_aipcw
  upper_aipcw <- psi_aipcw+qnorm(0.925)*se_aipcw

  return(list(psi_ipcw_tmle = psi_ipcw_tmle,
              lower_ipcw_tmle = lower_ipcw_tmle,
              upper_ipcw_tmle = upper_ipcw_tmle,
              psi_tmle = psi_tmle,
              lower_tmle = lower_tmle,
              upper_tmle = upper_tmle,
              psi_aipcw = psi_aipcw,
              lower_aipcw = lower_aipcw,
              upper_aipcw = upper_aipcw))
}

# # estimate Q and g
# if (is.null(QAW)) {
#   if (Q_g_method == "MI") {
#     # use multiple imputation to estimate Q and g
#     W.stage2_aug <- as.data.frame(matrix(NA, nrow = length(Y), ncol = ncol(W.stage2)))
#     W.stage2_aug[Delta.W == 1, ] <- W.stage2
#     names(W.stage2_aug) <- names(W.stage2)
#     df_mi <- data.frame(Y = Y,
#                         A = A,
#                         Delta.W = Delta.W,
#                         W,
#                         W.stage2_aug)
#     if (mi_engine == "mice") {
#       init <- mice::mice(df_mi, maxit = 0)
#       predM <- init$predictorMatrix
#       imp  <- mice::mice(df_mi, predictorMatrix = predM, m = Nimp,
#                          maxit = 20, print = FALSE)
#       imp_list <- lapply(seq_len(Nimp), function(x) {
#         mice::complete(imp, x)
#       })
#     } else if (mi_engine == "mixgb") {
#       if (!requireNamespace("mixgb", quietly = TRUE)) {
#         stop("Please install the 'mixgb' package to use mixgb imputation.")
#       }
#       mixgb_args <- c(list(data = df_mi, m = Nimp), mixgb_params)
#       imp_list <- do.call(mixgb::mixgb, mixgb_args)
#     }
#
#     Q1W_all <- Q0W_all <- QAW_all <- g1W_all <- matrix(0, nrow = length(Y), ncol = length(imp_list))
#
#     for (m in 1:length(imp_list)) {
#       comp <- imp_list[[m]]
#       tmle_m <- tmle::tmle(Y = comp$Y,
#                            A = comp$A,
#                            W = cbind(comp[, colnames(W), drop = FALSE],
#                                      comp[, colnames(W.stage2), drop = FALSE]),
#                            family = Q.family,
#                            g.SL.library = argList$g.SL.library,
#                            Q.SL.library = argList$Q.SL.library,
#                            verbose = FALSE)
#       Q1W_m <- tmle_m$Qinit$Q[, "Q1W"]
#       Q0W_m <- tmle_m$Qinit$Q[, "Q0W"]
#       QAW_m <- Q1W_m*A+Q0W_m*(1-A)
#       g1W_m <- tmle_m$g$g1W
#
#       Q1W_all[, m] <- Q1W_m; Q0W_all[, m] <- Q0W_m; QAW_all[, m] <- QAW_m
#       g1W_all[, m] <- g1W_m
#     }
#     Q1W <- rowMeans(Q1W_all); Q0W <- rowMeans(Q0W_all); QAW <- rowMeans(QAW_all)
#     g1W <- rowMeans(g1W_all)
#
#     # compute relevant quantities
#     H1W <- 1/g1W; H0W <- -1/(1-g1W); HAW <- A/g1W-(1-A)/(1-g1W)
#     DFullNC <- HAW*(Y-QAW)+Q1W-Q0W
#     GAW <- H1W-H0W-HAW^2
#     K <- Delta.W/res.twoStage$pi
#   } else if (Q_g_method == "ipcw") {
#     Q1W <- numeric(length(Delta.W))
#     Q0W <- numeric(length(Delta.W))
#     g1W <- numeric(length(Delta.W))
#     Q1W[Delta.W == 1] <- result$tmle$Qinit$Q[, "Q1W"]
#     Q0W[Delta.W == 1] <- result$tmle$Qinit$Q[, "Q0W"]
#     QAW <- Q1W*A+Q0W*(1-A)
#     g1W[Delta.W == 1] <- result$tmle$g$g1W
#     DFullNCReg_pool <- FALSE # when ipcw, can only learn on Delta=1 obs
#
#     # compute relevant quantities
#     H1W <- 1/g1W; H0W <- -1/(1-g1W); HAW <- A/g1W-(1-A)/(1-g1W)
#     H1W[Delta.W == 0] <- 0; H0W[Delta.W == 0] <- 0; HAW[Delta.W == 0] <- 0
#     DFullNC <- HAW*(Y-QAW)+Q1W-Q0W
#     GAW <- H1W-H0W-HAW^2
#     K <- as.numeric(Delta.W/res.twoStage$pi)
#   }

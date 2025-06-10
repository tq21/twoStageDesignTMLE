#' Plug-in imputation-based TMLE
#'
#' - Estimate Qn, gn via inverse weighting, get DFullNC estimate for Delta=1
#' - Estimate E(DFullNC|Delta=1,V) using HAL
#' - epsilon, one-dimensional path through Qn, evaluate P_n Delta/Pi_n DFullNC
#' 1. Regress DFullNC on V, only on Delta=1 observations
#' 2. Generate MI full data, run full data TMLE on each, obtain empirical EIC,
#'    then average to get E(DFullNC|Delta=1,V)
#' 3. Use MI to directly estimate E(DFullNC|Delta=1,V)
library(nleqslv)
library(sl3)
imp_plugin_logistic <- function(Y, A, W, W_all, Delta.W, W.stage2,
                                Q_g_method = "MI",
                                Nimp = 1,
                                QAW=NULL, Q0W=NULL, Q1W=NULL, g1W=NULL,
                                DFullNCReg_pool=TRUE,
                                Z=NULL, Delta = rep(1, length(Y)), pi=NULL, piform=NULL, pi_oracle=NULL,
                                DFullReg_true = NULL,
                                DFullReg_sl_lib = list(Lrnr_glm$new(), Lrnr_glmnet$new(), Lrnr_dbarts$new()),
                                pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet"),
                                #DFullReg.library = c("SL.glm", "SL.glmnet", "SL.gam", "tmle.SL.dbarts2"),
                                #DFullReg.library = c("SL.glm"),
                                V.pi=10, pi.discreteSL = TRUE, condSetNames = c("A","W","Y"), id = NULL,
                                Q.family = "binomial", augmentW = TRUE,
                                augW.SL.library = c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2"),
                                rareOutcome=FALSE,
                                verbose=FALSE, browse=FALSE, ...) {

  if (browse) browser()

  if(is.null(id)){id <- 1:length(Y)}

  if (is.vector(W)){
    W <- as.matrix(W)
    colnames(W) <- "W1"
  }

  if (is.null(colnames(W))){
    colnames(W) <- paste0("W", 1:ncol(W))
  }

  if(is.vector(W.stage2)){
    W.stage2 <- as.matrix(W.stage2)
    colnames(W.stage2) <- "W.stage2"
  }

  # if augmenting W call TMLE to get initial estimate of Q in full sample
  # do this through the tmle function to get cross-validated initial estimates
  # that are not overfit to the data and are appropriately bounded
  if (augmentW){
    W.Q <- evalAugW(Y, A, W, Delta, id, Q.family, augW.SL.library)
  } else {
    W.Q <- NULL
  }

  # Evaluate conditional sampling probabilities
  res.twoStage <- estimatePi(Y=Y, A=A, W=W, condSetNames=condSetNames, W.Q=W.Q,
                             Delta.W = Delta.W, piform = piform, id = id,
                             pi.SL.library = pi.SL.library,  V = V.pi,
                             discreteSL = pi.discreteSL,
                             verbose=verbose, pi = pi)
  if (!is.null(pi_oracle)) {
    res.twoStage$pi <- pi_oracle
  }

  # Bound normalized obsWeights after rescaling to sum to sum(Delta.W)
  obsWeights <- Delta.W/res.twoStage$pi
  obsWeights <- obsWeights / sum(obsWeights) * sum(Delta.W)
  ub <- sqrt(sum(Delta.W)) * log(sum(Delta.W)) / 5
  obsWeights <- .bound(Delta.W/res.twoStage$pi, c(0, ub))

  # Set call to tmle on full data
  argList <- list(...)
  argList$Y <- Y[Delta.W==1]
  argList$A <- A[Delta.W==1]
  if (is.null(W.Q)){
    argList$W <- cbind(W[Delta.W==1,, drop=FALSE], W.stage2)
  } else {
    argList$W <- cbind(W[Delta.W==1,, drop=FALSE], W.Q[Delta.W==1,], W.stage2)
  }
  argList$Z <- Z[Delta.W==1]
  argList$Delta <- Delta[Delta.W==1]
  argList$obsWeights <- obsWeights[Delta.W==1]
  argList$id <- id[Delta.W==1]
  argList$family <- Q.family
  argList$verbose <- verbose
  if (rareOutcome) {
    argList$Q.SL.library <- c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2")
    argList$Q.discreteSL <- TRUE
    argList$V.Q <- 20
  }
  result <- list()
  result$tmle <- try(do.call(tmle::tmle, argList))

  # estimate Q and g
  if (is.null(QAW)) {
    if (Q_g_method == "MI") {
      # use multiple imputation to estimate Q and g
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
                             g.SL.library = argList$g.SL.library,
                             Q.SL.library = argList$Q.SL.library,
                             verbose = FALSE)
        Q1W_m <- tmle_m$Qinit$Q[, "Q1W"]
        Q0W_m <- tmle_m$Qinit$Q[, "Q0W"]
        QAW_m <- Q1W_m*A+Q0W_m*(1-A)
        g1W_m <- tmle_m$g$g1W

        Q1W_all[, m] <- Q1W_m; Q0W_all[, m] <- Q0W_m; QAW_all[, m] <- QAW_m
        g1W_all[, m] <- g1W_m
      }
      Q1W <- rowMeans(Q1W_all); Q0W <- rowMeans(Q0W_all); QAW <- rowMeans(QAW_all)
      g1W <- rowMeans(g1W_all)

      # compute relevant quantities
      H1W <- 1/g1W; H0W <- -1/(1-g1W); HAW <- A/g1W-(1-A)/(1-g1W)
      DFullNC <- HAW*(Y-QAW)+Q1W-Q0W
      GAW <- H1W-H0W-HAW^2
      K <- Delta.W/res.twoStage$pi
    } else if (Q_g_method == "ipcw") {
      Q1W <- numeric(length(Delta.W))
      Q0W <- numeric(length(Delta.W))
      g1W <- numeric(length(Delta.W))
      Q1W[Delta.W == 1] <- result$tmle$Qinit$Q[, "Q1W"]
      Q0W[Delta.W == 1] <- result$tmle$Qinit$Q[, "Q0W"]
      QAW <- Q1W*A+Q0W*(1-A)
      g1W[Delta.W == 1] <- result$tmle$g$g1W
      DFullNCReg_pool <- FALSE # when ipcw, can only learn on Delta=1 obs

      # compute relevant quantities
      H1W <- 1/g1W; H0W <- -1/(1-g1W); HAW <- A/g1W-(1-A)/(1-g1W)
      H1W[Delta.W == 0] <- 0; H0W[Delta.W == 0] <- 0; HAW[Delta.W == 0] <- 0
      DFullNC <- HAW*(Y-QAW)+Q1W-Q0W
      GAW <- H1W-H0W-HAW^2
      K <- as.numeric(Delta.W/res.twoStage$pi)
    }
  } else {
    # compute relevant quantities
    H1W <- 1/g1W; H0W <- -1/(1-g1W); HAW <- A/g1W-(1-A)/(1-g1W)
    DFullNC <- HAW*(Y-QAW)+Q1W-Q0W
    GAW <- H1W-H0W-HAW^2
    K <- Delta.W/res.twoStage$pi
  }

  if (browse) browser()

  # estimate m=E(DFullNC|Delta=1,V)
  if (DFullNCReg_pool) {
    DFullNC_data <- cbind(res.twoStage$d.pi, DFullNC = DFullNC, Delta.W = Delta.W)
    DFullNC_data_pred <- DFullNC_data; DFullNC_data_pred$Delta.W <- 1
    DFullNC_task <- sl3_Task$new(data = DFullNC_data,
                                 covariates = c(names(res.twoStage$d.pi), "Delta.W"),
                                 outcome = "DFullNC", outcome_type = "gaussian")
    DFullNC_data_pred_task <- sl3_Task$new(data = DFullNC_data_pred,
                                           covariates = c(names(res.twoStage$d.pi), "Delta.W"),
                                           outcome = "DFullNC", outcome_type = "gaussian")
  } else {
    DFullNC_data <- cbind(res.twoStage$d.pi, DFullNC = DFullNC)[Delta.W == 1,]
    DFullNC_data_pred <- cbind(res.twoStage$d.pi, DFullNC = DFullNC)
    DFullNC_task <- sl3_Task$new(data = DFullNC_data,
                                 covariates = names(res.twoStage$d.pi),
                                 outcome = "DFullNC", outcome_type = "gaussian")
    DFullNC_data_pred_task <- sl3_Task$new(data = DFullNC_data_pred,
                                           covariates = names(res.twoStage$d.pi),
                                           outcome = "DFullNC", outcome_type = "gaussian")
  }
  lrnr_stack <- Stack$new(DFullReg_sl_lib)
  DFullReg_lrnr <- make_learner(Pipeline,
                                Lrnr_cv$new(lrnr_stack),
                                Lrnr_cv_selector$new(loss_squared_error))
  suppressMessages(DFullNC_fit_sl <- DFullReg_lrnr$train(DFullNC_task))
  DFullNCReg <- DFullNC_fit_sl$predict(DFullNC_data_pred_task)
  DFullNCReg_R2 <- 1 - sum((DFullNC-DFullNCReg)^2)/sum((DFullNC-mean(DFullNC))^2)
  result$DFullNCReg_R2 <- DFullNCReg_R2

  fn <- function(args) {
    eps <- args[1]
    gamma <- args[2]
    QAW_eps <- plogis(qlogis(QAW)+eps*HAW)
    Q0W_eps <- plogis(qlogis(Q0W)+eps*H0W)
    Q1W_eps <- plogis(qlogis(Q1W)+eps*H1W)
    DFullNC_eps <- HAW*(Y-QAW_eps)+Q1W_eps-Q0W_eps
    DFullNCReg_gamma <- DFullNCReg+gamma*K
    psi <- mean(K*(Q1W_eps-Q0W_eps))
    score_1 <- mean(K*(DFullNC_eps-DFullNCReg_gamma))
    score_2 <- mean(DFullNCReg_gamma)-psi
    return(c(score_1, score_2))
  }

  # solve linear system
  sol <- nleqslv(x = c(0, 0), fn = fn, method = "Newton")
  eps <- sol$x[1]; gamma <- sol$x[2]

  # updates
  QAW_star <- plogis(qlogis(QAW)+eps*HAW)
  Q0W_star <- plogis(qlogis(Q0W)+eps*H0W)
  Q1W_star <- plogis(qlogis(Q1W)+eps*H1W)
  DFullNC_star <- HAW*(Y-QAW_star)+Q1W_star-Q0W_star
  DFullNCReg_star <- DFullNCReg+gamma*K

  # point estimate and inference
  psi <- mean(DFullNCReg_star)
  result$eic <- as.numeric(K*(DFullNC_star-DFullNCReg_star)+DFullNCReg_star-psi)
  se <- sqrt(var(result$eic, na.rm = TRUE) / length(Delta.W))
  result$psi <- psi
  result$lower <- psi + qnorm(0.025) * se
  result$upper <- psi + qnorm(0.975) * se

  # check PnEIC
  PnEIC_1 <- mean(K*(DFullNC_star-DFullNCReg_star))
  PnEIC_2 <- mean(DFullNCReg_star)-psi

  # A-IPCW ---------------------------------------------------------------------
  psi_aipcw <- mean(Delta.W/res.twoStage$pi*(DFullNC-DFullNCReg)+DFullNCReg)
  eic_aipcw <- Delta.W/res.twoStage$pi*(DFullNC-DFullNCReg)+DFullNCReg-psi_aipcw
  se_aipcw <- sqrt(var(eic_aipcw, na.rm = TRUE) / length(Delta.W))
  psi_aipcw_lower <- psi_aipcw + qnorm(0.025) * se_aipcw
  psi_aipcw_upper <- psi_aipcw + qnorm(0.925) * se_aipcw
  result$psi_aipcw <- psi_aipcw
  result$lower_aipcw <- psi_aipcw_lower
  result$upper_aipcw <- psi_aipcw_upper

  result$pi <- res.twoStage$pi

  result$twoStage <- res.twoStage
  result$augW <- W.Q
  class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
  return(result)
}

#' A-IPCW-TMLE, plug-in, imputation-based estimator
#'
#' Compare three methods to estimate E(DFullNC|Delta=1,V)
#' 1. Regress DFullNC on V, only on Delta=1 observations
#' 2. Generate MI full data, run full data TMLE on each, obtain empirical EIC,
#'    then average to get E(DFullNC|Delta=1,V)
#' 3. Use MI to directly estimate E(DFullNC|Delta=1,V)
twoStageTMLE_plugin_DFullReg_compare <- function(Y, A, W, Delta.W, W.stage2, Z=NULL,
                                                 Delta = rep(1, length(Y)), pi=NULL, piform=NULL, pi_oracle=NULL,
                                                 DFullReg_true = NULL,
                                                 pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet", "tmle.SL.dbarts.k.5"),
                                                 DFullReg.library = c("SL.glm", "SL.glmnet", "SL.gam", "tmle.SL.dbarts2"),
                                                 V.pi=10, pi.discreteSL = TRUE, condSetNames = c("A","W","Y"), id = NULL,
                                                 Q.family = "gaussian", augmentW = TRUE,
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
  if (rareOutcome){
    argList$Q.SL.library <- c("SL.glm", "SL.glmnet", "tmle.SL.dbarts2")
    argList$Q.discreteSL <- TRUE
    argList$V.Q <- 20
  }
  result <- list()
  result$tmle <- try(do.call(tmle::tmle, argList))

  # estimate non-centered EIC --------------------------------------------------
  # method 1: Regress DFullNC on V, only on Delta=1 observations
  Q1W <- result$tmle$Qinit$Q[, "Q1W"]
  Q0W <- result$tmle$Qinit$Q[, "Q0W"]
  QAW <- Q1W*A[Delta.W == 1]+Q0W*(1-A[Delta.W == 1])
  g1W <- result$tmle$g$g1W
  DFullNC <- (A[Delta.W == 1]/g1W-(1-A[Delta.W == 1])/(1-g1W))*(Y[Delta.W == 1]-QAW)+Q1W-Q0W
  DFullNCReg <- estimateDfullReg(DFull = DFullNC,
                                 Delta = Delta.W,
                                 V = res.twoStage$d.pi,
                                 DFullbounds = c(-Inf, Inf),
                                 DFullform = NULL,
                                 SL.library = DFullReg.library,
                                 verbose = verbose,
                                 discreteSL = TRUE,
                                 Vfold = argList$V.Q)
  DFullNCReg <- DFullNCReg$DFullReg
  DFullNCReg_R2 <- 1 - sum((DFullNC-DFullNCReg[Delta.W == 1])^2)/sum((DFullNC-mean(DFullNC))^2)
  result$DFullNCReg_R2 <- DFullNCReg_R2

  # method 2: Generate MI full data, run full data TMLE on each, obtain empirical EIC,
  #           then average to get E(DFullNC|Delta=1,V)
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
  DFullNC_all <- matrix(0, nrow = length(Y), ncol = Nimp)
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
    DFullNC_all[, m] <- tmle_m$estimates$IC$IC.ATE
  }
  DFullNCReg_MI <- rowMeans(DFullNC_all)
  DFullNCReg_MI_R2 <- 1 - sum((DFullNC-DFullNCReg_MI[Delta.W == 1])^2)/sum((DFullNC-mean(DFullNC))^2)
  result$DFullNCReg_MI_R2 <- DFullNCReg_MI_R2

  # method 3: Use MI to directly estimate E(DFullNC|Delta=1,V)
  DFullNC_aug <- rep(NA, length(Delta.W))
  DFullNC_aug[Delta.W == 1] <- DFullNC
  W.stage2_aug <- as.data.frame(matrix(NA, nrow = length(Y), ncol = ncol(W.stage2)))
  W.stage2_aug[Delta.W == 1, ] <- W.stage2
  names(W.stage2_aug) <- names(W.stage2)
  df_mi <- data.frame(DFullNC = DFullNC_aug, # target of imputation
                      Y = Y,
                      A = A,
                      Delta.W = Delta.W,
                      W,
                      W.stage2_aug)
  init <- mice::mice(df_mi, maxit = 0)
  predM <- init$predictorMatrix
  predM[,] <- 0  # turn off all imputations
  predM["DFullNC", ] <- 1  # use all variables to impute DFullNC
  predM["DFullNC", "DFullNC"] <- 0  # but not itself
  Nimp <- 50
  imp <- mice::mice(df_mi, predictorMatrix = predM, method = "pmm", m = Nimp, maxit = 20, print = FALSE)
  DFullNCReg_all <- sapply(1:Nimp, function(m) mice::complete(imp, m)$DFullNC)
  DFullNCReg_MI_direct <- rowMeans(DFullNCReg_all)
  DFullNCReg_MI_direct_R2 <- 1 - sum((DFullNC-DFullNCReg_MI_direct[Delta.W == 1])^2)/sum((DFullNC-mean(DFullNC))^2)
  result$DFullNCReg_MI_direct_R2 <- DFullNCReg_MI_direct_R2

  # target DFullNCReg ----------------------------------------------------------
  # method 1: Regress DFullNC on V, only on Delta=1 observations
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg[Delta.W == 1]))) # CHECK PT
  H <- as.numeric(Delta.W/res.twoStage$pi)
  epsilon <- coef(glm(DFullNC ~ -1 + offset(DFullNCReg[Delta.W == 1]) + H[Delta.W == 1], family = "gaussian"))
  epsilon[is.na(epsilon)] <- 0
  DFullNCReg[Delta.W == 1] <- DFullNCReg[Delta.W == 1] + epsilon * H[Delta.W == 1]
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg[Delta.W == 1]))) # CHECK PT

  # method 2: Generate MI full data, run full data TMLE on each, obtain empirical EIC,
  #           then average to get E(DFullNC|Delta=1,V)
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg_MI[Delta.W == 1]))) # CHECK PT
  epsilon <- coef(glm(DFullNC ~ -1 + offset(DFullNCReg_MI[Delta.W == 1]) + H[Delta.W == 1], family = "gaussian"))
  epsilon[is.na(epsilon)] <- 0
  DFullNCReg_MI[Delta.W == 1] <- DFullNCReg_MI[Delta.W == 1] + epsilon * H[Delta.W == 1]
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg_MI[Delta.W == 1]))) # CHECK PT

  # method 3: Use MI to directly estimate E(DFullNC|Delta=1,V)
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg_MI_direct[Delta.W == 1]))) # CHECK PT
  epsilon <- coef(glm(DFullNC ~ -1 + offset(DFullNCReg_MI_direct[Delta.W == 1]) + H[Delta.W == 1], family = "gaussian"))
  epsilon[is.na(epsilon)] <- 0
  DFullNCReg_MI_direct[Delta.W == 1] <- DFullNCReg_MI_direct[Delta.W == 1] + epsilon * H[Delta.W == 1]
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg_MI_direct[Delta.W == 1]))) # CHECK PT

  # point estimate and inference -----------------------------------------------
  psi <- mean(DFullNCReg)
  psi_MI <- mean(DFullNCReg_MI)
  psi_MI_direct <- mean(DFullNCReg_MI_direct)
  DFullNC_aug <- numeric(length(DFullNCReg))
  DFullNC_aug[Delta.W == 1] <- DFullNC
  result$eic <- as.numeric(Delta.W/res.twoStage$pi)*(DFullNC_aug-DFullNCReg)+DFullNCReg-psi
  result$eic_MI <- as.numeric(Delta.W/res.twoStage$pi)*(DFullNC_aug-DFullNCReg_MI)+DFullNCReg_MI-psi_MI
  result$eic_MI_direct <- as.numeric(Delta.W/res.twoStage$pi)*(DFullNC_aug-DFullNCReg_MI_direct)+DFullNCReg_MI_direct-psi_MI_direct
  se <- sqrt(var(result$eic, na.rm = TRUE) / length(Delta.W))
  se_MI <- sqrt(var(result$eic_MI, na.rm = TRUE) / length(Delta.W))
  se_MI_direct <- sqrt(var(result$eic_MI_direct, na.rm = TRUE) / length(Delta.W))
  result$psi <- psi
  result$psi_MI <- psi_MI
  result$psi_MI_direct <- psi_MI_direct
  result$lower <- psi + qnorm(0.025) * se
  result$lower_MI <- psi_MI + qnorm(0.025) * se_MI
  result$lower_MI_direct <- psi_MI_direct + qnorm(0.025) * se_MI_direct
  result$upper <- psi + qnorm(0.975) * se
  result$upper_MI <- psi_MI + qnorm(0.975) * se_MI
  result$upper_MI_direct <- psi_MI_direct + qnorm(0.975) * se_MI_direct

  result$twoStage <- res.twoStage
  result$augW <- W.Q
  class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
  return(result)
}

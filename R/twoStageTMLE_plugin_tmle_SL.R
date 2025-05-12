#' Plug-in imputation-based TMLE
#'
#' - Estimate Qn, gn via inverse weighting, get DFullNC estimate for Delta=1
#' - Estimate E(DFullNC|Delta=1,V) using HAL
#' - epsilon, one-dimensional path through Qn, evaluate P_n Delta/Pi_n DFullNC
#' 1. Regress DFullNC on V, only on Delta=1 observations
#' 2. Generate MI full data, run full data TMLE on each, obtain empirical EIC,
#'    then average to get E(DFullNC|Delta=1,V)
#' 3. Use MI to directly estimate E(DFullNC|Delta=1,V)
twoStageTMLE_plugin_tmle_SL <- function(Y, A, W, Delta.W, W.stage2, Z=NULL,
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

  # estimate Q
  # Q_fit <- glm(Y[Delta.W == 1] ~ .,
  #              data = data.frame(A = A[Delta.W==1], cbind(W[Delta.W==1, ,drop=FALSE], W.stage2)),
  #              family = Q.family, weights = 1/res.twoStage$pi[Delta.W == 1])
  # Q1W <- as.numeric(predict(Q_fit, newdata = data.frame(A = 1, cbind(W[Delta.W==1, ,drop=FALSE], W.stage2)), type = "response"))
  # Q0W <- as.numeric(predict(Q_fit, newdata = data.frame(A = 0, cbind(W[Delta.W==1, ,drop=FALSE], W.stage2)), type = "response"))
  # QAW <- Q1W*A[Delta.W == 1]+Q0W*(1-A[Delta.W == 1])
#
  # # estimate g
  # g_fit <- glm(A[Delta.W == 1] ~ ., data = data.frame(cbind(W[Delta.W==1, ,drop=FALSE], W.stage2)),
  #              family = Q.family, weights = 1/res.twoStage$pi[Delta.W == 1])
  # g1W <- as.numeric(predict(g_fit, newdata = data.frame(cbind(W[Delta.W==1, ,drop=FALSE], W.stage2, type = "response"))))
  # gbounds <- c(5/sqrt(length(A))/log(length(A)), 1-5/sqrt(length(A))/log(length(A)))
  # g1W <- .bound(g1W, gbounds)

  # compute clever covariates
  g1W <- result$tmle$g$g1W
  HAW <- A[Delta.W==1]/g1W-(1-A[Delta.W==1])/(1-g1W)
  H1W <- 1/g1W; H0W <- -1/(1-g1W)
  GAW <- H1W-H0W-HAW^2
  K <- Delta.W/res.twoStage$pi
  Q1W <- result$tmle$Qinit$Q[, "Q1W"]
  Q0W <- result$tmle$Qinit$Q[, "Q0W"]
  QAW <- Q1W*A[Delta.W == 1]+Q0W*(1-A[Delta.W == 1])

  # estimate m=E(DFullNC|Delta=1,V)
  DFullNC <- HAW*(Y[Delta.W==1]-QAW)+Q1W-Q0W
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

  # estimate \tilde{m}=E(GAW|Delta=1,V)
  GAWReg <- estimateDfullReg(DFull = GAW,
                             Delta = Delta.W,
                             V = res.twoStage$d.pi,
                             DFullbounds = c(-Inf, Inf),
                             DFullform = NULL,
                             SL.library = DFullReg.library,
                             verbose = verbose,
                             discreteSL = TRUE,
                             Vfold = argList$V.Q)
  GAWReg <- GAWReg$DFullReg
  GAWReg_R2 <- 1 - sum((GAW-GAWReg[Delta.W == 1])^2)/sum((GAW-mean(GAW))^2)
  result$GAWReg_R2 <- GAWReg_R2

  # solve linear system
  gamma <- mean(K[Delta.W==1]*(GAW-GAWReg[Delta.W==1]))/mean(K[Delta.W==1]^2)
  delta <- mean(K[Delta.W==1]*(DFullNC-DFullNCReg[Delta.W==1]))/mean(K[Delta.W==1]^2)
  numer <- mean(K[Delta.W==1]*(Q1W-Q0W)-DFullNCReg[Delta.W==1]-delta*K[Delta.W==1])
  denom <- mean(GAWReg[Delta.W==1]+gamma*K[Delta.W==1]-K[Delta.W==1]*(H1W-H0W))
  epsilon <- numer/denom

  # updates
  QAW_star <- QAW+epsilon*HAW
  Q1W_star <- Q1W+epsilon*H1W
  Q0W_star <- Q0W+epsilon*H0W
  DFullNC_star <- DFullNC+epsilon*GAW
  GAWReg_star <- GAWReg[Delta.W==1]+gamma*K[Delta.W==1]
  DFullNCReg_star <- DFullNCReg[Delta.W==1]+epsilon*GAWReg_star+delta*K[Delta.W==1]

  # check PnEIC
  PnEIC_1 <- mean(K[Delta.W==1]*(DFullNC+epsilon*GAW-DFullNCReg[Delta.W==1]-epsilon*GAWReg[Delta.W==1]-delta*K[Delta.W==1]))
  PnEIC_2 <- mean(DFullNCReg[Delta.W==1]+epsilon*GAWReg[Delta.W==1]+delta*K[Delta.W==1])-mean(K[Delta.W==1]*(Q1W-Q0W+epsilon*(H1W-H0W)))

  # point estimate and inference
  psi <- mean(DFullNCReg_star)
  DFullNC_aug <- numeric(length(Delta.W))
  DFullNC_aug[Delta.W == 1] <- DFullNC_star
  result$eic <- as.numeric(1/res.twoStage$pi[Delta.W==1]*(DFullNC_star-DFullNCReg_star)+DFullNCReg_star-psi)
  se <- sqrt(var(result$eic, na.rm = TRUE) / sum(Delta.W))
  result$psi <- psi
  result$lower <- psi + qnorm(0.025) * se
  result$upper <- psi + qnorm(0.975) * se

  # A-IPCW ---------------------------------------------------------------------
  psi_aipcw <- mean(Delta.W/res.twoStage$pi*(DFullNC_aug-DFullNCReg)+DFullNCReg)
  eic_aipcw <- Delta.W/res.twoStage$pi*(DFullNC_aug-DFullNCReg)+DFullNCReg-psi_aipcw
  se_aipcw <- sqrt(var(eic_aipcw, na.rm = TRUE) / length(Delta.W))
  psi_aipcw_lower <- psi_aipcw + qnorm(0.025) * se_aipcw
  psi_aipcw_upper <- psi_aipcw + qnorm(0.925) * se_aipcw
  result$psi_aipcw <- psi_aipcw
  result$lower_aipcw <- psi_aipcw_lower
  result$upper_aipcw <- psi_aipcw_upper

  result$twoStage <- res.twoStage
  result$augW <- W.Q
  class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
  return(result)
}

#' Plug-in imputation-based TMLE
#'
#' - Estimate Qn, gn via inverse weighting, get DFullNC estimate for Delta=1
#' - Estimate E(DFullNC|Delta=1,V) using HAL
#' - epsilon, one-dimensional path through Qn, evaluate P_n Delta/Pi_n DFullNC
#' 1. Regress DFullNC on V, only on Delta=1 observations
#' 2. Generate MI full data, run full data TMLE on each, obtain empirical EIC,
#'    then average to get E(DFullNC|Delta=1,V)
#' 3. Use MI to directly estimate E(DFullNC|Delta=1,V)
twoStageTMLE_plugin_tmle <- function(Y, A, W, Delta.W, W.stage2, Z=NULL,
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

  # estimate E(DFullNC|Delta=1,V)
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

  # compute Gn(X)=H(1,W)-H(0,W)-H^2(A,W)-PnDelta/Pi(H(1,W)-H(0,W))
  H1W <- A[Delta.W == 1]/g1W
  H0W <- -(1-A[Delta.W == 1])/(1-g1W)
  HAW <- H1W*A[Delta.W == 1]+H0W*(1-A[Delta.W == 1])
  GX <- H1W-H0W-HAW^2-mean(1/res.twoStage$pi[Delta.W == 1]*(H1W-H0W))

  # estimate E(Gn|Delta=1,V)
  GXReg <- estimateDfullReg(DFull = GX,
                            Delta = Delta.W,
                            V = res.twoStage$d.pi,
                            DFullbounds = c(-Inf, Inf),
                            DFullform = NULL,
                            SL.library = DFullReg.library,
                            verbose = verbose,
                            discreteSL = TRUE,
                            Vfold = argList$V.Q)
  GXReg <- GXReg$DFullReg

  # compute fluctuation epsilon
  numer <- mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNCReg[Delta.W == 1]-DFullNC))
  denom <- mean(1/res.twoStage$pi[Delta.W == 1]*(GX-GXReg[Delta.W == 1]))
  eps <- numer/denom

  # compute Qn_eps and point estimate
  QAW_eps <- QAW+eps*HAW
  Q1W_eps <- Q1W+eps*H1W
  Q0W_eps <- Q0W+eps*H0W
  DFullNC_eps <- DFullNC+eps*GX
  DFullNC_eps_aug <- numeric(length(Delta.W)); DFullNC_eps_aug[Delta.W == 1] <- DFullNC_eps
  DFullNCReg_eps <- DFullNCReg+eps*GXReg
  psi <- mean(DFullNCReg_eps[Delta.W == 1])

  # check score conditions
  score_1 <- Delta.W/res.twoStage$pi*(DFullNC_eps_aug-DFullNCReg_eps)
  score_2 <- DFullNCReg_eps-psi
  eic <- score_1+score_2
  print(mean(eic))

  se <- sqrt(var(eic, na.rm = TRUE) / length(Delta.W))
  result$psi <- psi
  result$lower <- psi + qnorm(0.025) * se
  result$upper <- psi + qnorm(0.975) * se

  result$twoStage <- res.twoStage
  result$augW <- W.Q
  class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
  return(result)
}

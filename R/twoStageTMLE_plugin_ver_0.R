#' twoStageTMLE
#'
#' Inverse probability of censoring weighted TMLE for evaluating
#' parameters when the full set of covariates is available on only
#' a subset of observations.
#'
#' @param Y outcome
#' @param A binary treatment indicator
#' @param W covariate matrix observed on everyone
#' @param Delta.W binary indicator of missing second stage covariates
#' @param W.stage2 matrix of second stage covariates observed on subset
#'  of observations
#' @param Z optional mediator of treatment effect for evaluating a
#'  controlled direct effect
#' @param Delta binary indicator of missing value for outcome \code{Y}
#' @param pi optional vector of missingness probabilities for \code{W.stage2}
#' @param piform parametric regression formula for estimating \code{pi} (see Details)
#' @param pi.SL.library super learner library for estimating \code{pi} (see Details)
#' @param V.pi number of cross validation folds for estimating \code{pi}
#'  using super learner
#' @param pi.discreteSL Use discrete super learning when \code{TRUE}, otherwise
#'  ensemble super learning
#' @param condSetNames Variables to include as predictors of missingness
#'  in \code{W.stage2}, any combination of \code{Y, A}, and either \code{W} (for all
#'  covariates in \code{W}), or individual covariate names in \code{W}
#' @param id Identifier of independent units of observation, e.g., clusters
#' @param Q.family Regression family for the outcome
#' @param augmentW When \code{TRUE} include predicted values for the outcome
#'  the set of covariates used to model the propensity score
#' @param augW.SL.library super learner library for preliminary outcome
#' regression model (ignored when \code{augmentW} is \code{FALSE})
#' @param rareOutcome When \code{TRUE} specifies less ambitious SL for Q in call
#' to \code{tmle} (discreteSL, glm, glmnet, bart library, \code{V=20})
#' @param verbose When \code{TRUE} prints informational messages
#' @param ... other parameters passed to the tmle function (not checked)
#'
#' @return object of class 'twoStageTMLE'.
#' \item{tmle}{Treatment effect estimates and summary information}
#' \item{twoStage}{IPCW weight estimation summary, \code{pi} are the probabilities,
#' \code{coef} are SL weights or coefficients from glm fit, \code{type} of
#' estimation procedure, \code{discreteSL} flag indicating whether discrete
#' super learning was used}
#' \item{augW}{Matrix of predicted outcomes based on stage 1 covariates only}
#'
#' @details
#' When using \code{piform} to specify a parametric model for pi that conditions
#' on the outcome use \code{Delta.W} as the dependent variable and \code{Y.orig}
#' on the right hand side of the formula instead of \code{Y}. When writing a
#' user-defined SL wrapper for inclusion in \code{pi.SL.library} use \code{Y}
#' on the left hand side of the formula. If specific covariate names are
#' used on the right hand side use \code{Y.orig} to condition
#' on the outcome.
#'
#' @examples
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rnorm(n)
#' W3 <- rnorm(n)
#' A <- rbinom(n, 1, plogis(-1 + .2*W1 + .3*W2 + .1*W3))
#' Y <- 10 + A + W1 + W2 + A*W1 + W3 + rnorm(n)
#' d <- data.frame(Y, A, W1, W2, W3)

#' # Set 400 with data on W3, more likely if W1 > 1
#' n.sample <- 400
#' p.sample <- 0.5 + .2*(W1 > 1)
#' rows.sample <- sample(1:n, size = n.sample, p = p.sample)
#' Delta.W <- rep(0,n)
#' Delta.W[rows.sample] <- 1
#' W3.stage2 <- cbind(W3 = W3[Delta.W==1])
#' #1. specify parametric models and do not augment W (fast, but not recommended)
#' result1 <- twoStageTMLE(Y=Y, A=A, W=cbind(W1, W2), Delta.W = Delta.W,
#'    W.stage2 = W3.stage2, piform = "Delta.W~ I(W1 > 0) + Y.orig", V.pi= 5,
#'    verbose = TRUE, Qform = "Y~A+W1",gform="A~W1 + W2 +W3", augmentW = FALSE)
#' summary(result1)
#' \donttest{
#' #2. specify a parametric model for conditional missingness probabilities (pi)
#' #   and use default values to estimate marginal effect using \code{tmle}
#' result2 <- twoStageTMLE(Y=Y, A=A, W=cbind(W1, W2), Delta.W = Delta.W,
#'      W.stage2 = cbind(W3)[Delta.W == 1], piform = "Delta.W~ I(W1 > 0)",
#'      V.pi= 5,verbose = TRUE)
#' result2
#' }
#' @seealso
#' * [tmle::tmle()] for details on customizing the estimation procedure
#' * [twoStageTMLEmsm()] for estimating conditional effects
#' * S Rose and MJ van der Laan. A Targeted Maximum Likelihood Estimator for
#' Two-Stage Designs. \emph{Int J Biostat.} 2011 Jan 1; 7(1): 17.
#' \doi{doi:10.2202/1557-4679.1217}
#' @export
twoStageTMLE_plugin_ver_0 <- function(Y, A, W, Delta.W, W.stage2, Z=NULL,
                                      Delta = rep(1, length(Y)), pi=NULL, piform=NULL, pi_oracle=NULL,
                                      DFullReg_true = NULL,
                                      pi.SL.library = c("SL.glm", "SL.gam", "SL.glmnet", "tmle.SL.dbarts.k.5"),
                                      DFullReg.library = c("tmle.SL.dbarts2"),
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

  # target Q -------------------------------------------------------------------
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

  # target DFullNCReg ----------------------------------------------------------
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg[Delta.W == 1]))) # CHECK PT
  H <- as.numeric(Delta.W/res.twoStage$pi)
  epsilon <- coef(glm(DFullNC ~ -1 + offset(DFullNCReg[Delta.W == 1]) + H[Delta.W == 1], family = "gaussian"))
  epsilon[is.na(epsilon)] <- 0
  DFullNCReg[Delta.W == 1] <- DFullNCReg[Delta.W == 1] + epsilon * H[Delta.W == 1]
  #print(mean(1/res.twoStage$pi[Delta.W == 1]*(DFullNC - DFullNCReg[Delta.W == 1]))) # CHECK PT

  # point estimate and inference -----------------------------------------------
  psi <- mean(DFullNCReg)
  DFullNC_aug <- numeric(length(DFullNCReg))
  DFullNC_aug[Delta.W == 1] <- DFullNC
  result$eic <- as.numeric(Delta.W/res.twoStage$pi)*(DFullNC_aug-DFullNCReg)+DFullNCReg-psi
  se <- sqrt(var(result$eic, na.rm = TRUE) / length(Delta.W))
  result$psi <- psi
  result$lower <- psi + qnorm(0.025) * se
  result$upper <- psi + qnorm(0.975) * se

  result$twoStage <- res.twoStage
  result$augW <- W.Q
  class(result) <- "twoStageTMLE"  # for tmleMSM can use same class
  return(result)
}

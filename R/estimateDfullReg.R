#-----------estimateQ----------------
# purpose: estimate Q=E(Y |Z, A,W) data-adaptively,
# unless super learner not available, or user specifies
# initial values or a regression formula
# arguments:
# 	Y - outcome
# 	Z - intermediate variable between A and Y (default= 0 when no int. var.)
#	A - treatment indicator (1=treatment, 0=control)
# 	W - baseline covariates
#	Delta - missingness indicator
#	Q - optional externally estimated values for Q
#	Qbounds - bounds for predicted values
#  	Qform - optional regression formula to use for glm if
#	        non-data adaptive estimation specified
# 	maptoYstar - if TRUE, using logistic fluctuation for bounded, continuous outcomes
# 		estimation inital Q on linear scale, bounded by (0,1),and return on logit scale
#		(will work if family=poisson)
#	SL.library - library of prediction algorithms for Super Learner
#   cvQinit - flag, if TRUE, cross-validate predictions
# discreteSL - if TRUE, use predictions from best learner, otherwise ensemble SL predictions
# 	family - regression family
#	id - subject identifier
# V - number of cross-validation folds
#  obsWeights- typically sampling weights, e.g., 1/p(selected | YAW)
# returns matrix of linear predictors for Q(A,W), Q(0,W), Q(1,W),
#   (for controlled direct effect, 2 additional columns: Q(Z=1,A=0,W), Q(Z=1,A=1,W))
#		family for stage 2 targeting
#		coef, NA, unless Q is estimated using a parametric model
# 		type, estimation method for Q
#----------------------------------------
estimateDfullReg <- function(DFull,
                             Delta,
                             V,
                             DFullbounds,
                             DFullform,
                             SL.library,
                             verbose,
                             discreteSL,
                             Vfold) {
  .expandLib <- function(SL.lib){
    if (is.list(SL.lib)){
      counts <- sapply(SL.lib, length)
      numExtra <-  sum(counts > 2)
      m <- matrix("", nrow = length(SL.lib) + numExtra, ncol = 2)
      rowIndex <- 1
      for (i in 1:length(counts)){
        if(counts[i] ==1 ){
          m[rowIndex, ] <- c(SL.lib[[i]], "All")
        } else{
          m[rowIndex, ] <- SL.lib[[i]][1:2]
        }
        rowIndex <- rowIndex+1
        if (counts[i] > 2){
          for (j in 3:counts[i]){
            m[rowIndex,] <- SL.lib[[i]][c(1,  j)]
            rowIndex <- rowIndex+1
          }
        }
      }
      return(split(m, 1:nrow(m)))
    } else {
      return(SL.lib)
    }
  }
  SL.version <- 2
  m <- NULL
  coef <- NA
  type <- "user-supplied values"
  if(verbose) { cat("\tEstimating initial regression of DFull on V\n")}
  DFullReg <- numeric(length(Delta))
  if(!(is.null(DFullform))){
    if(identical(as.character(as.formula(DFullform)), c("~","DFull", "."))){
      DFullform <- paste("DFull~", paste(colnames(V), collapse="+"))
    }
    DFullform <- as.formula(DFullform)
    m <- suppressWarnings(glm(DFullform, data=data.frame(DFull, V[Delta == 1, , drop=FALSE]), family="gaussian"))
    DFullReg <- as.numeric(predict(m, newdata=data.frame(V)))
    coef <- coef(m)
    type <- "glm, user-supplied model"
    DFullReginit <- list(DFullReg=DFullReg, coef=coef, type=type)
  } else {
    if (verbose) {cat("\t using SuperLearner\n")}
    n <- sum(Delta)
    X <- data.frame(V[Delta == 1, , drop = FALSE])
    newX <- data.frame(V)
    if (packageDescription("SuperLearner")$Version < SL.version){
      arglist <- list(Y=DFull, X=X, newX=newX, SL.library=SL.library,
                      V=Vfold, family="gaussian", save.fit.library=FALSE)
    } else {
      arglist <- list(Y=DFull, X=X, newX=newX, SL.library=SL.library,
                      cvControl=list(V=Vfold), family="gaussian",
                      control = list(saveFitLibrary=FALSE))
    }
    suppressWarnings(m <- try(do.call(SuperLearner, arglist)))
    if(inherits(m, "SuperLearner")) {
      coef <- m$coef
      if (discreteSL){
        keepAlg <- which.min(m$cvRisk)
        SL.coef <- 1
        type <- paste("SuperLearner, discrete, selected",  paste(.expandLib(SL.library)[[keepAlg]], collapse = ", "))
      } else {
        keepAlg <- which(m$coef > 0)
        SL.coef <- m$coef[m$coef > 0]
        type <- "SuperLearner, ensemble"
      }
      if (discreteSL) {
        keepAlg <- which.min(m$cvRisk)
        SL.coef <- 1
        DFullReg <- as.numeric(m$library.predict[,keepAlg])
      } else {
        DFullReg <- as.numeric(as.matrix(m$library.predict[,keepAlg]) %*% SL.coef)
      }
    } else {
      stop("Super Learner failed when estimating DFullReg. Exiting program\n")
    }
  }

  if (any(is.na(DFullReg)) | inherits(m, "try-error")) {
    if (verbose) {cat("\t Warning: \nRunning main terms regression for 'DFull' using glm\n")}
    DFullform <- paste("DFull~", paste(colnames(V), collapse="+"))
    m <- glm(DFullform, data=data.frame(DFull = DFull, V[Delta == 1, , drop = FALSE]), family="gaussian")
    DFullReg <- as.numeric(predict(m, newdata=data.frame(V)))
    coef <- coef(m)
    type="glm, main terms model"
  }
  DFullReg <- .bound(DFullReg, DFullbounds)
  DFullReginit <- list(DFullReg=DFullReg, coef=coef, type=type, SL.library=SL.library)
  return(DFullReginit)
}

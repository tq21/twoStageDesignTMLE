#' Estimate E(DFullNC|Delta=1,V)
#'
#' Method 2:
#' Regress DFullNC on Delta and V, on all observations, evaluate at Delta=1
estimateDFullReg_MI <- function(DFull,
                                  Delta,
                                  V,
                                  DFullbounds,
                                  DFullform,
                                  SL.library,
                                  verbose,
                                  discreteSL,
                                  Vfold) {
  browser()
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
    n <- length(Delta)
    X <- data.frame(V, Delta = Delta)
    newX <- data.frame(V, Delta = 1)
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

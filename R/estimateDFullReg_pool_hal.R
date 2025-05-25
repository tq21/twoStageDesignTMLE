#' Estimate E(DFullNC|Delta=1,V)
#'
#' Method 2:
#' Regress DFullNC on Delta and V, on all observations, evaluate at Delta=1
library(hal9001)
library(glmnet)
estimateDFullReg_pool_hal <- function(DFull,
                                      Delta,
                                      V) {
  basis_list <- enumerate_basis(x = as.matrix(cbind(V, Delta = Delta)),
                                max_degree = 2,
                                smoothness_orders = 1)
  hal_design <- make_design_matrix(X = as.matrix(cbind(V, Delta = Delta)),
                                   blist = basis_list)
  hal_design_pred <- make_design_matrix(X = as.matrix(cbind(V, Delta = 1)),
                                        blist = basis_list)
  fit <- cv.glmnet(x = hal_design,
                   y = DFull,
                   family = "gaussian",
                   alpha = 1)
  pred <- as.numeric(predict(fit, newx = hal_design_pred, s = "lambda.min"))
  return(list(DFullReg = pred))
}

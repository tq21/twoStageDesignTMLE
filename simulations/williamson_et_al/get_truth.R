get_true_ATE <- function(N = 1e6, XScenario = 1, YScenario = 1,
                         lowcor = 0.1, midcor = 0.4, highcor = 0.7, gencor = 0.2) {
  # Step 1: Generate covariates (X is overwritten below)
  covars <- XFunc(N = N, lowcor = lowcor, midcor = midcor,
                  highcor = highcor, gencor = gencor, XScenario = XScenario)

  # Step 2: Generate potential outcomes under treatment (X = 1)
  covars_treat <- covars
  covars_treat$X <- 1
  Y1 <- YFunc(data = covars_treat, YScenario = YScenario)$Y

  # Step 3: Generate potential outcomes under control (X = 0)
  covars_ctrl <- covars
  covars_ctrl$X <- 0
  Y0 <- YFunc(data = covars_ctrl, YScenario = YScenario)$Y

  # Step 4: Return true ATE
  return(mean(Y1) - mean(Y0))
}

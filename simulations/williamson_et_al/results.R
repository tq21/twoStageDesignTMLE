library(dplyr)

truth <- get_true_ATE(N = 1e7,
                      XScenario = 1.6,
                      YScenario = 4.17,
                      lowcor = 0.1,
                      midcor = 0.4,
                      highcor = 0.7,
                      gencor = 0.2)

res_df %>%
  summarize(abs_bias = abs(mean(psi - truth)),
            se = sd(psi),
            mse = mean((psi - truth)^2),
            bias_se = abs_bias / se,
            coverage = mean((lower <= truth) & (truth <= upper)),
            .by = c("est"))

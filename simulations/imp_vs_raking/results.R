library(dplyr)
library(knitr)
library(kableExtra)

res_df <- read.csv("out/mice-tmle_x1_y4.1_m2.2_0611_150836.csv")
res_df <- read.csv("out/rak_x1_y4.1_m2.2_0611_111024.csv")
truth <- get_true_ATE(N = 1e7, XScenario = 1, YScenario = 4.1)
res_df %>%
  summarize(
    abs_bias = abs(mean(psi - truth)),
    se = sd(psi),
    mse = mean((psi - truth)^2),
    bias_se = abs_bias / se,
    cover = mean((lower <= truth) & (truth <= upper)),
    oracle_cover = mean(truth >= psi+qnorm(0.025)*sd(psi) & truth <= psi+qnorm(0.975)*sd(psi)),
    .by = c("n", "est")
  )

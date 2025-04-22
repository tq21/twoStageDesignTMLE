library(dplyr)
library(knitr)
library(kableExtra)
source("sim_data.R")

res_df <- read.csv("out/Y_simple_binomial_Delta_large_0417_131354.csv")
truth <- get_truth(Y_type = "simple_binomial")
res_df %>%
  summarize(
    abs_bias = abs(mean(psi - truth)),
    se = sd(psi),
    mse = mean((psi - truth)^2),
    bias_se = abs_bias / se,
    coverage = mean((lower <= truth) & (truth <= upper)),
    power = mean((lower >= 0) & (0 <= upper)),
    .by = c("n", "est_name")
  )

library(dplyr)
library(knitr)
library(kableExtra)
source("sim_data.R")

res_df <- read.csv("out/miss_type_med_500_0603_212509.csv")
truth <- get_truth()
res_df %>%
  summarize(
    abs_bias = abs(mean(psi - truth)),
    se = sd(psi),
    mse = mean((psi - truth)^2),
    bias_se = abs_bias / se,
    cover = mean((lower <= truth) & (truth <= upper)),
    oracle_cover = mean(truth >= psi+qnorm(0.025)*sd(psi) & truth <= psi+qnorm(0.975)*sd(psi)),
    .by = c("n", "est_name")
  )

library(dplyr)
library(knitr)
library(kableExtra)
library(tidyr)
source("sim_data.R")

res_df <- read.csv("out/Y_complex_miss_strong_0617_100340.csv")
truth <- get_truth(Y_type = "complex")
res_df <- res_df[complete.cases(res_df), ]
res_df %>%
  summarize(
    abs_bias = abs(mean(psi - truth)),
    se = sd(psi),
    mse = mean((psi - truth)^2),
    bias_se = abs_bias / se,
    coverage = mean((lower <= truth) & (truth <= upper)),
    .by = c("n", "est_name")
  )

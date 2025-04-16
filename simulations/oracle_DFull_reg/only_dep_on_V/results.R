library(dplyr)
library(knitr)
library(kableExtra)

res_df <- read.csv("out/0415_141721.csv")
truth <- 0.4
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

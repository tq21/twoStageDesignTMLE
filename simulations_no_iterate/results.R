library(dplyr)
library(knitr)
library(kableExtra)

res_df <- read.csv("out/Y_1_Delta_1.1_0415_233521.csv")
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
  ) %>%
  mutate(across(c(abs_bias, se, mse), ~ formatC(.x, format = "e", digits = 2))) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, caption = "") %>%
  kable_styling(latex_options = c("striped", "hold_position"))

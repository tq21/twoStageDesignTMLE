library(dplyr)

res_df <- read.csv("out/Y_1_Delta_1_0412_213000.csv")
truth <- 0.4
res_df %>%
  summarize(abs_bias = abs(mean(psi - truth)),
            se = sd(psi),
            mse = mean((psi - truth)^2),
            bias_se = abs_bias/se,
            coverage = mean((lower <= truth) & (truth <= upper)),
            power = mean((lower >= 0) & (0 <= upper)),
            .by = c("n", "est_name"))

library(dplyr)
source("src_williamson_et_al/00_utils.R")
source("src_williamson_et_al/01_generate_data.R")
source("src_williamson_et_al/02_methods_design.R")
source("src_williamson_et_al/02_methods_tmle.R")
source("src_williamson_et_al/03_estimate.R")
source("get_truth.R")

set.seed(123)
truth <- get_true_ATE(N = 1e7,
                      XScenario = 1.6,
                      YScenario = 4.17,
                      lowcor = 0,
                      midcor = 0,
                      highcor = 0,
                      gencor = 0)

res_df <- read.csv("out/ver_2_0423_152703.csv")

res_df %>%
  summarize(abs_bias = abs(mean(psi - truth)),
            se = sd(psi),
            mse = mean((psi - truth)^2),
            bias_se = abs_bias / se,
            coverage = mean((lower <= truth) & (truth <= upper)),
            .by = c("n", "est"))

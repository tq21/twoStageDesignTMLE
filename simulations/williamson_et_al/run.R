`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.4",
            .libPaths()))
library(dplyr)
library(purrr)
library(devtools)
library(survey)
library(mice)
library(data.table)
library(marginaleffects)
library(sl3)
load_all()
source("src_williamson_et_al/00_utils.R")
source("src_williamson_et_al/01_generate_data.R")
source("src_williamson_et_al/02_methods_design.R")
source("src_williamson_et_al/02_methods_tmle.R")
source("src_williamson_et_al/03_estimate.R")
set.seed(123)

B <- 500
n_seq <- seq(500, 2000, 500)

run <- function(est_name) {
  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")

      estimators <- "rak"
      if (est_name == "ver_0") {
        estimators <- c("ipcw-tmle-ver-0", estimators)
      } else if (est_name == "ver_1") {
        estimators <- c("ipcw-tmle-ver-1", estimators)
      } else if (est_name == "ver_2") {
        estimators <- c("ipcw-tmle-ver-2", estimators)
      } else if (est_name == "ver_3") {
        estimators <- c("ipcw-tmle-ver-3", estimators)
      } else if (est_name == "plugin_ver_0") {
        estimators <- c("ipcw-tmle-plugin-ver-0", estimators)
      } else if (est_name == "ipcw-tmle-plugin-DFullReg-compare") {
        estimators <- c("ipcw-tmle-plugin-DFullReg-compare", estimators)
      }

      res <- investigate_performance_once(mc_id = 1,
                                          n = .n,
                                          XScenario = 1.6,
                                          YScenario = 4.17,
                                          missScenario = 1,
                                          lowcor = 0, midcor = 0, highcor = 0, gencor = 0,
                                          estimators = estimators,
                                          tmle_args = list(
                                            "g_lib" = c("SL.glm"), "miss_lib" = c("SL.glm"),
                                            "q_lib" = c("SL.glm"), "K" = 5,
                                            "phase1_covars" = c("Y", "X", "Zs", "Zw"),
                                            "phase2_covars" = c("Ws", "Ww")
                                          ),
                                          mi_args = list("n_imp" = 20, maxiter = 25),
                                          raking_args = list("NimpRaking" = 20),
                                          outcome_name = "Y", tx_name = "X",
                                          fam = "binomial", missing_indicator = "is.complete",
                                          cached_datasets = NULL, data_only = FALSE,
                                          plasmode = FALSE,
                                          data_dir = "./", filename_prefix = "data_m1_y1_x1_n10000_id",
                                          filename_suffix = ".rds", read_func = readRDS,
                                          rare_outcome = FALSE)
      res$results <- cbind(res$results, data.frame(n = rep(.n, nrow(res$results))))

      return(res$results)
    })
  })

  return(res_df)
}



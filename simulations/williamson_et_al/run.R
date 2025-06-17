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
library(future)
load_all()
source("src_williamson_et_al/00_utils.R")
source("src_williamson_et_al/01_generate_data.R")
source("src_williamson_et_al/02_methods_design.R")
source("src_williamson_et_al/02_methods_tmle.R")
source("src_williamson_et_al/03_estimate.R")
ncores <- 5
plan(multisession, workers = ncores)

#' function to run simulation
#'
#' @param est_name name of the estimator to run
#' @param xscenario scenario for X
#' @param yscenario scenario for Y
#' @param mscenario scenario for missingness
#' @param B number of Monte-Carlo runs for each sample size
#' @param n_seq sequence of sample sizes to run
#' @param seed seed for random number generation
run <- function(est_name,
                xscenario,
                yscenario,
                mscenario,
                B,
                n_seq,
                seed,
                tmle_args = NULL,
                raking_args = NULL) {

  set.seed(seed)
  print("running estimator: " %+% est_name)
  print("X scenario: " %+% xscenario %+%
        ", Y scenario: " %+% yscenario %+%
        ", missing scenario: " %+% mscenario)

  res_df <- map_dfr(n_seq, function(.n) {
    # set seed
    current_seed <- round(yscenario * 10) +
      round(mscenario * 100) +
      round(xscenario * 100) +
      round(.n * 1000) + round((seed - 1) * 51)
    set.seed(current_seed)

    map_dfr(seq(B), function(.b) {
      print("n: " %+% .n %+% ", b: " %+% .b %+% "...")

      if (!is.null(tmle_args)) {
        tmle_args$phase1_covars <- c("Y", "X", "Zs", "Zw")
        tmle_args$phase2_covars <- c("Ws", "Ww")
      }

      # run given estimator
      res <- tryCatch(
        investigate_performance_once(
          mc_id = .b, n = .n, XScenario = xscenario, fam = "binomial",
          YScenario = yscenario, missScenario = mscenario,
          lowcor = 0.2, midcor = 0.4, highcor = 0.7, gencor = 0.2,
          outcome_name = "Y", tx_name = "X", missing_indicator = "is.complete",
          estimators = est_name, tmle_args = tmle_args, mi_args = NULL,
          raking_args = raking_args, cached_datasets = NULL,
          data_only = FALSE, plasmode = FALSE,
          rare_outcome = FALSE, browse = FALSE
        ), error = function(e) {
          message(conditionMessage(e))
          print(paste0("Error occurred when running Monte-Carlo iteration ", .b))
          traceback()
          output <- list()
          output$results <- data.frame(est = est_name,
                                       psi = NULL,
                                       lower = NULL,
                                       upper = NULL)
          return(output)
        }
      )
      res_df <- cbind(data.frame(n = .n, b = .b), res$results)
      return(res_df)
    })
  })

  return(res_df)
}

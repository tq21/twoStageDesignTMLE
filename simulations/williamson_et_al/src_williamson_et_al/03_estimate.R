# Run the simulation a single time

# Version control log:
# DATE MODIFIED  MODIFIER  PURPOSE
# -------------  --------  -------
# 20230725       bdw       create a file/function that runs the simulation once
# 20230807       bdw       add RRZ, raking, add new options
# 20231115       bdw       harmonize with new functions
# 20231201       eaj       Add RF to list of estimators
# 20240201       eaj       included n_imps to XGB function call, now relies on mice_args value
# 20250502       bdw       add different flavors of TMLE (IPCW-augmented, rare-outcome IPCW)

# run the simulation a single time ---------------------------------------------
# @param mc_id the Monte-Carlo replicate ID
# @param n the sample size
# @param XScenario the covariate scenario
# @param YScenario the outcome regression model
# @param missScenario the missing-data mechanism
# @param lowcor low correlation
# @param midcor middle correlation
# @param highcor high correlation
# @param gencor general correlation
# @param estimators the estimator(s) to run
# @param tmle_args a list of arguments to pass to TMLE (SL libraries and number of cross-validation folds)
# @param mi_args a list of arguments to pass to MI algorithms (number of imputations, number of iterations)
# @param raking_args a list of arguments to pass to raking (number of imputations)
# @param outcome_name the name of the outcome
# @param tx_name the name of the treatment variable
# @param fam the outcome family (for regression models)
# @param missing_indicator the name of the missingness indicator
# @param cached_datasets pre-cached datasets (NULL if we want to generate data)
# @param data_only if TRUE, generates and returns the dataset, performs no other computation
# @param plasmode if TRUE, this is a plasmode simulation
# @param data_dir the directory where datasets should be loaded from (only used in plasmode simulations)
# @param filename_prefix the prefix for filenames (only used in plasmode simulations)
# @param filename_suffix the suffix for filenames (i.e., .Rdata or .rds) (only used in plasmode simulations)
# @param read_func the function to read in datasets (i.e., load or readRDS) (only used in plasmode simulations)
# @param rare_outcome if TRUE (and plasmode = TRUE), simplifies the working outcome regression model by removing an interaction between prior self-harm and PHQ item 9
# @return output from a single simulation
investigate_performance_once <- function(mc_id = 1, n = 100, XScenario = 1, YScenario = 1,
                                         missScenario = 1, lowcor = 0.1, midcor = 0.4, highcor = 0.7,
                                         gencor = 0.2, estimators = c("CC"),
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
                                         rare_outcome = FALSE,
                                         browse = FALSE) {
  if (browse) browser()

  # generate a dataset
  if (is.null(cached_datasets)) {
    df <- gen_data(N = n, XScenario = XScenario, YScenario = YScenario, missScenario = missScenario, lowcor = lowcor, midcor = midcor, highcor = highcor, gencor = gencor)
    df$mc_id <- mc_id
  } else {
    if (is.data.frame(cached_datasets)) {
      df <- cached_datasets[cached_datasets$mc_id == mc_id, -ncol(cached_datasets)]
    } else {
      if (plasmode) {
        if (isTRUE(cached_datasets)) {
          df <- read_func(file = paste0(this_data_dir, filename_prefix, mc_id, filename_suffix))
        } else {
          df <- cached_datasets[[mc_id]]
        }
        if (!any(grepl("AgeAtIndex_cat", names(df)))) {
          df$AgeAtIndex_cat <- factor(
            case_when(
              df$AgeAtIndex <= 24 ~ 0,
              df$AgeAtIndex > 24 & df$AgeAtIndex <= 34 ~ 1,
              df$AgeAtIndex > 34 & df$AgeAtIndex <= 44 ~ 2,
              df$AgeAtIndex > 44 & df$AgeAtIndex <= 54 ~ 3,
              df$AgeAtIndex > 54 & df$AgeAtIndex <= 64 ~ 4,
              df$AgeAtIndex > 64 ~ 5
            ), levels = 0:5, labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+")
          )
        }
        if (rare_outcome) {
          df$Charlson_cat <- ifelse(df$Charlson == 0, 0, 1)
        }
      } else {
        df <- cached_datasets[[mc_id]][, -ncol(cached_datasets[[mc_id]])]
      }
    }
  }
  if (data_only) {
    return(df)
  } else {
    all_outcome_formulas <- get_outcome_formulas(YScenario = YScenario, plasmode = plasmode, rare_outcome = rare_outcome)
    oracle_formula <- all_outcome_formulas$oracle
    population_formula <- all_outcome_formulas$population
    outcome_formula <- all_outcome_formulas$outcome
    outcome_formula_confounded <- all_outcome_formulas$confounded
    outcome_formula_factor <- all_outcome_formulas$outcome_factor
    all_miss_formulas <- get_miss_formula(mScenario = missScenario, xScenario = XScenario, plasmode = plasmode)
    miss_formula <- all_miss_formulas$miss_formula
    miss_formula_factor <- all_miss_formulas$miss_formula_factor
    all_tx_formulas <- get_tx_formula(XScenario = XScenario, plasmode = plasmode, rare_outcome = rare_outcome)
    tx_formula <- all_tx_formulas$tx_formula
    tx_formula_factor <- all_tx_formulas$tx_formula_factor

    if (!plasmode) {
      obs_vars <- c("Y", "X", "Zs", "Zw", "As", "Aw", "Ws_obs", "Ww_obs", "is.complete")
    } else {
      obs_vars <- c("Y", "numEpiType", "sex", "AgeAtIndex", "AgeIndexsq", "AgeAtIndex_cat",
                    "anx_dx_priorYr", "aud_dx_priorYr",
                    "Charlson_cat", "priorSH", "MHIP_prior5yr",
                    "IndexTrtPHQ8_score_cat_obs", "IndexTrtPHQ_item9_score_cat_obs", "is.complete")
    }
    # apply the estimators to it
    all_results <- lapply(as.list(estimators), function(est) {
      run_one_estimator(
        estimator = est, df = df, oracle_formula = oracle_formula,
        population_formula = population_formula,
        outcome_formula = outcome_formula, outcome_formula_factor = outcome_formula_factor,
        outcome_formula_confounded = outcome_formula_confounded,
        tx_formula = tx_formula, tx_formula_factor = tx_formula_factor,
        miss_formula = miss_formula, miss_formula_factor = miss_formula_factor,
        fam = fam, obs_vars = obs_vars, tmle_args = tmle_args, mi_args = mi_args,
        raking_args = raking_args, outcome_name = outcome_name,
        tx_name = tx_name, missing_indicator = missing_indicator, plasmode = plasmode,
        browse = browse
      )
    })
    if ((any(grepl("gr", estimators)) | any(grepl("rak", estimators))) & plasmode) {
      message(paste0("MC id ", mc_id, " complete"))
    }
    # return results
    results_df <- do.call(rbind, lapply(all_results, function(x) x$results))
    coefs <- do.call(rbind, lapply(all_results, function(x) x$coefs))
    results_df$mc_id <- mc_id
    return(list("results" = results_df, "coefs" = coefs, "df" = df))
  }
}

# Run a single estimator
# @param estimator a character string specifying the estimator, e.g., "cc_oracle"
# @param df the dataset
# @param oracle_formula the true outcome regression formula
# @param population_formula the working regression formula (using unobserved variables)
# @param outcome_formula the working outcome regression formula
# @param outcome_formula_factor the working outcome regression formula with explicit factor variables, for TMLE-M
# @param outcome_formula_confounded the confounded outcome regression formula
# @param miss_formula the missing-data formula
# @param miss_formula_factor the missing-data formula with explicit factor variables, for TMLE
# @param tx_formula the treatment assignment formula
# @param tx_formula_factor the treatment assignment formula with explicit factor variables, for TMLE-M
# @param tmle_args a list of arguments for TMLE
# @param mi_args a list of arguments for MI methods
# @param raking_args a list of arguments for raking methods
# @param calOption the calibration option
# @param outcome_name the name of the outcome
# @param tx_name the name of the treatment variable
# @param missing_indicator the name of the missingness indicator
# @param plasmode logical, is this a plasmode simulation?
# @param fam the family, e.g., "binomial"
# @return the results of running the estimator on the dataset
run_one_estimator <- function(estimator, df, oracle_formula,
                              population_formula, outcome_formula, outcome_formula_factor,
                              outcome_formula_confounded, miss_formula, miss_formula_factor,
                              tx_formula, tx_formula_factor, obs_vars, tmle_args, mi_args,
                              raking_args, outcome_name = "Y", tx_name = "X",
                              missing_indicator = "is.complete", plasmode = FALSE,
                              fam = "binomial", browse = FALSE) {
  if (browse) browser()

  # complete-case estimators
  if (grepl("CC", estimator, ignore.case = TRUE) | grepl("population", estimator, ignore.case = TRUE)) {
    if (grepl("oracle", estimator)) {
      this_formula <- oracle_formula
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else if (grepl("population", estimator)) {
      this_formula <- population_formula
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else if (grepl("noW", estimator)) {
      this_formula <- outcome_formula_confounded
      these_data <- df[, !grepl(missing_indicator, names(df))]
    } else {
      this_formula <- outcome_formula
      these_data <- df[complete.cases(df), !grepl(missing_indicator, names(df))]
    }
    output <- run_cc(
      data = these_data, formula = this_formula, procedure = estimator,
      coefficient_of_interest = tx_name, fam = fam
    )
  }
  # raking
  else if (grepl("rak", estimator, ignore.case = TRUE) | grepl("gr", estimator, ignore.case = TRUE)) {
    cal_option <- 1
    if (grepl("X", estimator, ignore.case = TRUE)) {
      cal_option <- 2
    }
    rake_on_y <- grepl("Y", estimator, ignore.case = TRUE)
    start_from_ipw <- grepl("ipw", estimator, ignore.case = TRUE)
    tmp <- run_raking_lr(
      data = df, formula = outcome_formula, miss_formula = miss_formula,
      NimpRaking = raking_args$NimpRaking, calOption = cal_option, fam = fam,
      coefficient_of_interest = tx_name, missing_indicator = missing_indicator,
      start_from_ipw = start_from_ipw, rake_on_y = rake_on_y
    )
    output <- list()
    output$results <- data.frame(est = "raking",
                                 psi = tmp$results[tmp$results$estimand == "RD", "est"],
                                 lower = tmp$results[tmp$results$estimand == "RD", "est"]+qnorm(0.025)*tmp$results[tmp$results$estimand == "RD", "SE"],
                                 upper = tmp$results[tmp$results$estimand == "RD", "est"]+qnorm(0.975)*tmp$results[tmp$results$estimand == "RD", "SE"])
  }
  # Full data Q, g method: MICE; Estimator: Imputation-based plug-in TMLE
  else if (estimator == "mice-tmle") {
    tmp <- imp_plugin_logistic(Y = df$Y,
                               A = df$X,
                               W = df[, setdiff(tmle_args$phase1_covars, c("X", "Y")), drop = FALSE],
                               W.stage2 = df[complete.cases(df), tmle_args$phase2_covars, drop = FALSE],
                               Delta.W = df$is.complete,
                               condSetNames = c("W", "A", "Y"),
                               Q_g_method = tmle_args$Q_g_method,
                               Nimp = tmle_args$Nimp,
                               pi.SL.library = tmle_args$pi.SL.library, V.pi = tmle_args$V.pi,
                               Q.SL.library = tmle_args$Q.SL.library, V.Q = tmle_args$V.Q,
                               g.SL.library = tmle_args$g.SL.library, V.g = tmle_args$V.g,
                               DFullReg_sl_lib = tmle_args$DFullReg_sl_lib,
                               Q.family = "binomial",
                               augmentW = FALSE,
                               verbose = FALSE,
                               browse = FALSE)

    output <- list()
    output$results <- data.frame(est = c("ipcw-tmle", "plugin-tmle"),
                                 psi = c(tmp$tmle$estimates$ATE$psi, tmp$psi),
                                 lower = c(tmp$tmle$estimates$ATE$CI[1], tmp$lower),
                                 upper = c(tmp$tmle$estimates$ATE$CI[2], tmp$upper))
  } else {
    output <- NULL
  }
  return(output)
}

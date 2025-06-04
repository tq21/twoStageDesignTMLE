source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
miss_type <- "med"
res_df <- run(n = 4000,
              B = 1000,
              miss_type = miss_type,
              Q_g_method = "MI",
              pi.SL.library = c("SL.glm", "SL.gam"),
              Q.SL.library = c("SL.glm", "SL.gam"),
              g.SL.library = c("SL.glm", "SL.gam"),
              DFullReg_sl_lib = list(Lrnr_glm$new(),
                                     Lrnr_glmnet$new(),
                                     Lrnr_dbarts$new(),
                                     Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0),
                                     Lrnr_earth$new(degree = 3),
                                     Lrnr_ranger$new()))

write.csv(res_df,
          file = "out/miss_type_" %+% miss_type %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

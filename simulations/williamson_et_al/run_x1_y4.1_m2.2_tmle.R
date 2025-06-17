source("run.R")
res_df <- run(est_name = "mice-tmle",
              xscenario = 1,
              yscenario = 4.1,
              mscenario = 2.2,
              B = 200,
              n_seq = 500,
              seed = 123,
              tmle_args = list(Q_g_method = "MI",
                               Nimp = 1, V.pi = 10, V.Q = 10, V.g = 10,
                               pi.SL.library = c("SL.glm", "SL.gam", "SL.xgboost"),
                               Q.SL.library = c("SL.glm", "SL.gam", "SL.xgboost"),
                               g.SL.library = c("SL.glm", "SL.gam", "SL.xgboost"),
                               DFullReg_sl_lib = list(Lrnr_glm$new(),
                                                      Lrnr_glmnet$new(),
                                                      Lrnr_dbarts$new())))
write.csv(res_df, file = "out/mice-tmle_x1_y4.1_m2.2_" %+% timestamp %+% ".csv", row.names = FALSE)

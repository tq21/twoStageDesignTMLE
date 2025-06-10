source("run.R")

res_df <- run(est_name = "mice-tmle",
              xscenario = 1,
              yscenario = 4.1,
              mscenario = 2.2,
              B = 10,
              n_seq = 2000,
              seed = 123,
              tmle_args = list(Q_g_method = "ipcw",
                               Nimp = 1, V.pi = 10, V.Q = 10, V.g = 10,
                               pi.SL.library = c("SL.glm", "SL.gam", "tmle.SL.dbarts2"),
                               Q.SL.library = c("SL.glm", "SL.gam", "tmle.SL.dbarts2"),
                               g.SL.library = c("SL.glm", "SL.gam", "tmle.SL.dbarts2"),
                               DFullReg_sl_lib = list(Lrnr_glm$new(),
                                                      Lrnr_glmnet$new(),
                                                      Lrnr_dbarts$new())))

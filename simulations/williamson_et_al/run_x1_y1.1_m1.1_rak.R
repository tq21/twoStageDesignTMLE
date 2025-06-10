source("run.R")

res_df <- run(est_name = "rak",
              xscenario = 1,
              yscenario = 4.1,
              mscenario = 2.2,
              B = 10,
              n_seq = 2000,
              seed = 123,
              raking_args = list("NimpRaking" = 20))

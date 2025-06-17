source("run.R")

res_df <- run(est_name = "rak",
              xscenario = 1,
              yscenario = 4.1,
              mscenario = 2.2,
              B = 200,
              n_seq = seq(500, 2000, 500),
              seed = 123,
              raking_args = list("NimpRaking" = 20))
write.csv(res_df, file = "out/rak_x1_y4.1_m2.2_" %+% timestamp %+% ".csv", row.names = FALSE)

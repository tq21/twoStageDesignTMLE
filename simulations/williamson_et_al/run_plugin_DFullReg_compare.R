source("run.R")

est_name <- "ipcw-tmle-plugin-DFullReg-compare"
res_df <- run(est_name = est_name)
write.csv(res_df,
          file = "out/" %+% est_name %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

source("run.R")

est_name <- "ver_3"
res_df <- run(est_name = est_name)
write.csv(res_df,
          file = "out/" %+% est_name %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

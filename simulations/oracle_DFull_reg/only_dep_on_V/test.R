source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
res_df <- run()

write.csv(res_df,
          file = "out/" %+% timestamp %+% ".csv",
          row.names = FALSE)

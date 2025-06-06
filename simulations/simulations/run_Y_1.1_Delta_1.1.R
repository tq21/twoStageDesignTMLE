source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

Y_type <- 1.1
Delta_type <- 1.1
res_df <- run(Y_type = Y_type, Delta_type = Delta_type)

write.csv(res_df,
          file = "out/Y_" %+% Y_type %+% "_Delta_" %+% Delta_type %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

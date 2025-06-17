source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

Y_type <- "complex"
miss_type <- "strong"
res_df <- run(Y_type = Y_type,
              miss_type = miss_type,
              true_Pi = FALSE,
              truth = get_truth(Y_type = Y_type),
              seed = 123)

write.csv(res_df,
          file = "out/Y_" %+% Y_type %+% "_miss_" %+% miss_type %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

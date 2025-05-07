source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

Y_type <- "complex_binomial"
Delta_type <- "med"
Q.family <- "binomial"
res_df <- run(Y_type = Y_type,
              Delta_type = Delta_type,
              true_Pi = FALSE,
              Q.family = Q.family)

write.csv(res_df,
          file = "out/Y_" %+% Y_type %+% "_Delta_" %+% Delta_type %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

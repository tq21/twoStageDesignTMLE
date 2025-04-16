source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

Y_type <- "simple_binomial"
Delta_type <- "med_strong"
Q.family <- "binomial"
res_df <- run(Y_type = Y_type,
              Delta_type = Delta_type,
              true_Pi = TRUE,
              Q.family = "binomial")

write.csv(res_df,
          file = "out/Y_" %+% Y_type %+% "_Delta_" %+% Delta_type %+% "_" %+% timestamp %+% ".csv",
          row.names = FALSE)

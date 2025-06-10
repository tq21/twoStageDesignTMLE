#' Estimator arguments

tmle_args <- list("g_lib" = c("SL.glm"),
                  "miss_lib" = c("SL.glm"),
                  "q_lib" = c("SL.glm"),
                  "K" = 5,
                  "phase1_covars" = c("Y", "X", "Zs", "Zw"),
                  "phase2_covars" = c("Ws", "Ww"))

# Raking estimator arguments
raking_args <- list("NimpRaking" = 20)


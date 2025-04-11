eic <- function(Delta,
                Pi,
                D_full,
                D_full_mean) {
  D_full_aug <- numeric(length(Delta))
  D_full_aug[Delta == 1] <- D_full
  ipcw_wt <- Delta/Pi
  wt_eic <- ipcw_wt*D_full_aug
  proj <- (ipcw_wt-1)*D_full_mean

  return(wt_eic-proj)
}

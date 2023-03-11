##--------------------##
## This file documents the optimization process by Susan Athey's Matrix Completion
## Estimation strategy
##--------------------##


library(MCPanel)

MCoptim = function(dat, T_obs, T0, treat = 1){
  
  # PanelY is the panel data of outcomes arranged time-series by rows and units
  # by columns; MCoptim assume the number of treated units is 1, placed at first 
  # column.
  
  PanelY = as.matrix(unstack(dat, Outcome ~ index))
  mask = matrix(TRUE, nrow(PanelY), ncol(PanelY))
  mask[(T0+1):T_obs,1] = FALSE
  estimated_obj <- mcnnm_cv(PanelY, mask)
  y0hat = estimated_obj$L[,1]
  return(y0hat)
}
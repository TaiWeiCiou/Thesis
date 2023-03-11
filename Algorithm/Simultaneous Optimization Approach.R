library(dplyr)

#Target function 1: Hsiao's Method by SOA
duols <- function(par, Y1, Y0, X1, X0, N_x, N_ctrl, T0){
  beta <- as.matrix(par[1:N_x], N_x, 1)
  a <- as.matrix(par[(N_x+1):(N_ctrl+N_x)], N_ctrl, 1)
  u1 <- Y1 - X1 %*% beta
  u00 <- Y0 - X0 %*% beta
  u0 <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u0[,i] <- u00[((i-1)*T_obs +1):(i*T_obs)]
  }
  target <- sum((Y0 - X0 %*% beta)^2) + sum( (u1 - u0[(1:T0),] %*% a )^2)
}

PanelNormal <- function(pan, T_obs, T0, treat = 1){
  datY = pan
  meanY = mean(as.matrix(datY[1:T0,]))
  meanY_byCol1 = apply(datY, 2, mean)[treat]
  meanY_byCol2 = apply(datY, 2, mean)[-treat]
  meanY_byCol = c(meanY_byCol1,meanY_byCol2)
  meanY_byRow1 = apply(datY[1:T0,], 1, mean)
  meanY_byRow2 = apply(datY[(T0+1):T_obs,-treat],1,mean)
  meanY_byRow = c(meanY_byRow1,meanY_byRow2)
  datY_normal = datY - meanY_byRow - 
    matrix(rep(meanY_byCol,nrow(datY)), nrow(datY), ncol(datY), byrow = T) + 
    meanY
  result = list(NormalPan = datY_normal,
                timefixedeff = meanY_byRow,
                unitfixedeff = meanY_byCol,
                mu = meanY)
  return(result)
}

simuoptim <- function(dat, T_obs, T0, N_ctrl, N_x){
  
  datY = unstack(dat, Outcome ~ index)
  datX1 = unstack(dat, X1~ index)
  datX2 = unstack(dat, X2~ index)
  
  ##normalize data
  ### normalize matrix Y
  YfixedEsti = PanelNormal(datY, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datY_normal = YfixedEsti$NormalPan  
  ### normalize matrix X_1
  X1fixedEsti = PanelNormal(datX1, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datX1_normal = X1fixedEsti$NormalPan
  ### normalize matrix X_2
  X2fixedEsti = PanelNormal(datX2, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datX2_normal = X2fixedEsti$NormalPan
  
  ## vector-ize data for the convenience of calculation
  Y1 <- datY_normal[1:T0,1] %>% as.matrix()
  Y0 <- stack(datY_normal[,-1])[,1] %>% as.matrix()
  X1 <- cbind(datX1_normal[1:T0,1], datX2_normal[1:T0,1]) %>% as.matrix()
  X0 <- cbind(
    stack(datX1_normal[,-1])[,1], stack(datX2_normal[,-1])[,1]
  ) %>% as.matrix()

  
  ## setting parameter
  parr <- rep(0,(N_x+N_ctrl))
  
  ## optimiztion
  optresult <- optim(par = parr,
                     fn = duols,
                     Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
                     N_x = 2, N_ctrl = 10,T0 = T0,
                     method = "SANN")
  
  ## claiming estimated coefficient
  beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
  aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl)], N_ctrl, 1)
  
  ## constructing the counter-factual at treatment period
  X_tr <- cbind(datX1_normal[,1], datX2_normal[,1]) %>% as.matrix()
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  
  ## estimate fixed effect
  ### estimate mu
  meanY = YfixedEsti$mu
  meanX = cbind(X1fixedEsti$mu, X2fixedEsti$mu) %>% as.matrix()
  mu = meanY - meanX %*% beeta
  
  ### estimate alpha
  meanY_byCol = matrix(YfixedEsti$unitfixedeff, (N_treat + N_ctrl),1)
  meanX_byCol = cbind(
    X1fixedEsti$unitfixedeff, 
    X2fixedEsti$unitfixedeff
    ) %>% as.matrix()
  alpha = meanY_byCol - meanX_byCol %*% beeta
  
  ### estimate Xi
  meanY_byRow = matrix(YfixedEsti$timefixedeff, T_obs,1)
  meanX_byRow = cbind(
    X1fixedEsti$timefixedeff, 
    X2fixedEsti$timefixedeff
  ) %>% as.matrix()
  
  Xi = meanY_byRow - meanX_byRow %*% beeta
  
  y0hat <- X_tr %*% beeta + u_co %*% aa + matrix(alpha[1], T_obs,1) + Xi + matrix(mu, T_obs, 1)
  
  
  ## result
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa)
  return(result)
}

## Hsiao's method by SOA
duols.simplex <- function(par, Y1, Y0, X1, X0, N_x, N_ctrl, T0){
  beta <- as.matrix(par[1:N_x], N_x, 1)
  a <- c(1-sum(par[(N_x+1):(N_ctrl+N_x-1)]), par[(N_x+1):(N_ctrl+N_x-1)])
  u1 <- Y1 - X1 %*% beta
  u00 <- Y0 - X0 %*% beta
  u0 <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u0[,i] <- u00[((i-1)*T_obs +1):(i*T_obs)]
  }
  target <- sum((Y0 - X0 %*% beta)^2) + sum( (u1 - u0[(1:T0),] %*% a )^2)
}

simuoptim.simplex <- function(dat, T_obs, T0, N_ctrl, N_x){
  
  datY = unstack(dat, Outcome ~ index)
  datX1 = unstack(dat, X1~ index)
  datX2 = unstack(dat, X2~ index)
  
  ##normalize data
  ### normalize matrix Y
  YfixedEsti = PanelNormal(datY, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datY_normal = YfixedEsti$NormalPan  
  ### normalize matrix X_1
  X1fixedEsti = PanelNormal(datX1, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datX1_normal = X1fixedEsti$NormalPan
  ### normalize matrix X_2
  X2fixedEsti = PanelNormal(datX2, T_obs = T_obs, T0 = T0, treat = 1 ) 
  datX2_normal = X2fixedEsti$NormalPan
  
  ## vector-ize data for the convenience of calculation
  Y1 <- datY_normal[1:T0,1] %>% as.matrix()
  Y0 <- stack(datY_normal[,-1])[,1] %>% as.matrix()
  X1 <- cbind(datX1_normal[1:T0,1], datX2_normal[1:T0,1]) %>% as.matrix()
  X0 <- cbind(
    stack(datX1_normal[,-1])[,1], stack(datX2_normal[,-1])[,1]
  ) %>% as.matrix()
  
  ## setting parameter
  parr <- c(rep(0,N_x), rep(1/N_ctrl, N_ctrl-1))
  
  ## optimiztion
  optresult <- optim(par = parr,
                     fn = duols.simplex,
                     Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
                     N_x = 2, N_ctrl = N_ctrl,T0 = T0,
                     lower = c(rep(-Inf, N_x),rep(0, N_ctrl-1)),
                     upper = c(rep(Inf, N_x),rep(1, N_ctrl-1)),
                     method = "L-BFGS-B")
  
  ## claiming estimated coefficient
  beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
  aa <- as.matrix(c(1-sum(optresult$par[(N_x+1):(N_x+N_ctrl-1)]),
                    optresult$par[(N_x+1):(N_x+N_ctrl-1)]),N_ctrl, 1)
  
  ## constructing the counter-factual at treatment period
  
  X_tr <- cbind(datX1_normal[,1], datX2_normal[,1]) %>% as.matrix()
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  
  ## estimate fixed effect
  ### estimate mu
  meanY = YfixedEsti$mu
  meanX = cbind(X1fixedEsti$mu, X2fixedEsti$mu) %>% as.matrix()
  mu = meanY - meanX %*% beeta
  
  ### estimate alpha
  meanY_byCol = matrix(YfixedEsti$unitfixedeff, (N_treat + N_ctrl),1)
  meanX_byCol = cbind(
    X1fixedEsti$unitfixedeff, 
    X2fixedEsti$unitfixedeff
  ) %>% as.matrix()
  alpha = meanY_byCol - meanX_byCol %*% beeta
  
  ### estimate Xi
  meanY_byRow = matrix(YfixedEsti$timefixedeff, T_obs,1)
  meanX_byRow = cbind(
    X1fixedEsti$timefixedeff, 
    X2fixedEsti$timefixedeff
  ) %>% as.matrix()
  
  Xi = meanY_byRow - meanX_byRow %*% beeta
  
  y0hat <- X_tr %*% beeta + u_co %*% aa + matrix(alpha[1], T_obs,1) + Xi + matrix(mu, T_obs, 1)
  
  
  ## result
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa)
  return(result)
}

# t0 <- Sys.time()
# T_obs = 30; T0 = 20; N_ctrl  =10
# B = 1000
# for (i in 1:B){
#   ## generating a data for testing.
#   dataresult = DGP1(T_obs = 30, T0 = 20, N_ctrl = 10)
#   dat1 =  dataresult$datXu
#   datY0 = dataresult$datY0
#   y0 = datY0[,1]
# 
#   lsresult <- simuoptim(dat = dat1,
#                         T_obs = 30, T0 = 20, N_ctrl = 10, N_x = 2)
#   simplexresult<- simuoptim.simplex(dat = dat1,
#                                     T_obs = 30, T0 = 20, N_ctrl = 10, N_x = 2)
#   
#   ymean = matrix(mean(y0), T_obs, 1)
#   
#   y0hat <- cbind(lsresult$y0hat, simplexresult$y0hat, ymean)
#   
#   mse <- apply(y0hat, 2, function(x) (x - y0)^2)
#   mse_pre <- apply( mse[1:T0,], 2, mean); mse_aft <- apply( mse[((T0+1):T_obs),], 2, mean)
#   mse_ratio <- mse_aft/mse_pre
#   msetable0 <- msetable0 + t(rbind(mse_pre, mse_aft, mse_ratio))
#   msetable <- msetable0/i
# }
# t1 <- Sys.time()
# t1-t0
# 
# View(msetable)


## Hsaio's method with LASSO on both sides
# duols.panalty <- function(par, 
#                           Y1, Y0, X1, X0, N_x, N_ctrl, T0,
#                           lam1, lam2){
#   
#   beta <- as.matrix(par[1:N_x], N_x, 1)
#   a <- as.matrix(par[(N_x+1):(N_ctrl+N_x+1)], (N_ctrl+1), 1)
#   
#   u1 <- Y1 - X1 %*% beta
#   u00 <- Y0 - X0 %*% beta
#   u0 <- matrix(0, T_obs, N_ctrl)
#   for (i in 1:N_ctrl){
#     u0[,i] <- u00[((i-1)*T_obs +1):(i*T_obs)]
#   }
#   u0 <- cbind(rep(1, T_obs), u0)
#   
#   target <- 0.5*sum((Y0 - X0 %*% beta)^2) + 
#     0.5*sum( (u1 - u0[(1:T0),] %*% a )^2) +
#     0.5*lam1 * sum(abs(beta)) + 0.5*lam2* sum(abs(a))
#   
#   return(target)
# }

# lam.max <- function(y,X){
#   max(abs(t(y) %*% X))
# }
# which.min.mat <- function(A){
#   return(which(A == min(A), arr.ind = T))
# }
# nfoldCV <- function(n, n_fold = 5){
#   quo <- n %/% n_fold; r <- n %% n_fold
#   lenlist <- c(rep((quo+1), r), rep(quo, (n_fold-r))) 
#   lenlist2 <- c(0,cumsum(lenlist))
#   sam <- sample(1:n,n); 
#   indicelist <- list()
#   for (i in 1:n_fold){
#     indicelist[[i]] <- sam[(lenlist2[i]+1) : (lenlist2[i+1])]
#   }
#   return(indicelist)
# }
# 
# simuoptim.duolasso <- function(dat, T_obs, T0, N_ctrl, N_x,
#                                D = 10, validnum = 3){
#   
#   datY = unstack(dat, Outcome ~ index)
#   datX1 = unstack(dat, X1~ index)
#   datX2 = unstack(dat, X2~ index)
#   
#   ##normalize data
#   ### normalize matrix Y
#   YfixedEsti = PanelNormal(datY, T_obs = T_obs, T0 = T0, treat = 1 ) 
#   datY_normal = YfixedEsti$NormalPan  
#   ### normalize matrix X_1
#   X1fixedEsti = PanelNormal(datX1, T_obs = T_obs, T0 = T0, treat = 1 ) 
#   datX1_normal = X1fixedEsti$NormalPan
#   ### normalize matrix X_2
#   X2fixedEsti = PanelNormal(datX2, T_obs = T_obs, T0 = T0, treat = 1 ) 
#   datX2_normal = X2fixedEsti$NormalPan
#   
#   ## setting parameter
#   parr <- rep(0,(N_x+N_ctrl))
#   
#   ## vector-ize data for the convenience of calculation
#   Y1 <- datY_normal[1:T0,1] %>% as.matrix()
#   Y0 <- stack(datY_normal[,-1])[,1] %>% as.matrix()
#   X1 <- cbind(datX1_normal[1:T0,1], datX2_normal[1:T0,1]) %>% as.matrix()
#   X0 <- cbind(
#     stack(datX1_normal[,-1])[,1], stack(datX2_normal[,-1])[,1]
#   ) %>% as.matrix()
#   
#   ## generating lambda list
#   lam1.max <- lam.max(y = Y0, X = X0)
#   lam2.max <- lam.max(y = Y1, X = Y0[(1:T0),])
#   lam1list <- exp((1:D)*(log(lam1.max)/D))
#   lam2list <- exp((1:D)*(log(lam2.max)/D)) 
#   
#   #data for CV
#   dat_normal = dat
#   dat_normal['Outcome'] = stack(datY_normal)[,1]
#   dat_normal['X1'] = stack(datX1_normal)[,1]
#   dat_normal['X2'] = stack(datX2_normal)[,1]
#   
#   ## CV
#   
#   msetable <- matrix(0,D,D)
#   
#   t0 <- Sys.time()
#   sam <- sample(2:N, validnum)
#   for (fold in sam){
#     
#     pool <- c(2:N)[-(fold-1)]
#     
#     Y1 <- as.matrix(dat_normal[(dat_normal$index == fold),3][1:T0], T0, 1)
#     Y0 <- as.matrix(dat_normal[(dat_normal$index %in% pool),3], T_obs*(N_ctrl-1), 1)
#     X1 <- as.matrix(
#       dat_normal[(dat_normal$index == fold),c(4,5)][(1:T0),], T0, N_x)
#     X0 <- as.matrix(
#       dat_normal[(dat_normal$index %in% pool),c(4,5)], T_obs*(N_ctrl-1), N_x)
#     
#     for (ii in 1:D){
#       for (jj in 1:D){
#         
#         optresult <- optim(par = rep(0,(N_x+ N_ctrl)),
#                            fn = duols.panalty,
#                            Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
#                            N_x = 2, N_ctrl = N_ctrl-1,T0 = T0,
#                            method = "SANN",
#                            lam1 = lam1list[ii], lam2 = lam2list[jj])
#         
#         beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
#         aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl)], N_ctrl, 1)
#         
#         X_tr <-as.matrix(
#           dat_normal[(dat_normal$index == fold),c(4,5)], T_obs, N_x
#         )
#         
#         u_co0 <- Y0 - X0 %*% beeta
#         u_co <- matrix(0, T_obs, N_ctrl-1)
#         for (i in 1:(N_ctrl-1)){
#           u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
#         }
#         u_co <- cbind(rep(1,T_obs), u_co)
#         
#         y0hat <- X_tr %*% beeta + u_co %*% aa
#         
#         y0 <- dat_normal[(dat_normal$index == fold),3]
#         
#         msetable[ii,jj] <- msetable[ii,jj] + mean(((y0-y0hat)^2)[((T0+1):T_obs)])
#         
#       }
#     }
#     
#   }
#   t1 <- Sys.time()
#   t1-t0
#   
#   i.min <- which.min.mat(msetable)[1];j.min <- which.min.mat(msetable)[2]
#   
#   ## vector-ize data for the convenience of calculation
#   Y1 <- datY_normal[1:T0,1] %>% as.matrix()
#   Y0 <- stack(datY_normal[,-1])[,1] %>% as.matrix()
#   X1 <- cbind(datX1_normal[1:T0,1], datX2_normal[1:T0,1]) %>% as.matrix()
#   X0 <- cbind(
#     stack(datX1_normal[,-1])[,1], stack(datX2_normal[,-1])[,1]
#   ) %>% as.matrix()
#   
#   
#   ## optimiztion
#   optresult <- optim(par = rep(0,(N_x+N_ctrl)),
#                      fn = duols.panalty,
#                      Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
#                      N_x = 2, N_ctrl = N_ctrl,T0 = T0,
#                      method = "SANN",
#                      lam1 = lam1list[i.min],lam2 = lam2list[j.min])
#   
#   ## claiming estimated coefficient
#   beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
#   aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl)], N_ctrl, 1)
#   
#   ## constructing the counter-factual at treatment period
#   
#   X_tr <- cbind(datX1_normal[,1], datX2_normal[,1]) %>% as.matrix()
#   
#   u_co0 <- Y0 - X0 %*% beeta
#   u_co <- matrix(0, T_obs, N_ctrl)
#   for (i in 1:N_ctrl){
#     u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
#   }
#   
#   
#   ## estimate fixed effect
#   ### estimate mu
#   meanY = YfixedEsti$mu
#   meanX = cbind(X1fixedEsti$mu, X2fixedEsti$mu) %>% as.matrix()
#   mu = meanY - meanX %*% beeta
#   
#   ### estimate alpha
#   meanY_byCol = matrix(YfixedEsti$unitfixedeff, (N_treat + N_ctrl),1)
#   meanX_byCol = cbind(
#     X1fixedEsti$unitfixedeff, 
#     X2fixedEsti$unitfixedeff
#   ) %>% as.matrix()
#   alpha = meanY_byCol - meanX_byCol %*% beeta
#   
#   ### estimate Xi
#   meanY_byRow = matrix(YfixedEsti$timefixedeff, T_obs,1)
#   meanX_byRow = cbind(
#     X1fixedEsti$timefixedeff, 
#     X2fixedEsti$timefixedeff
#   ) %>% as.matrix()
#   
#   Xi = meanY_byRow - meanX_byRow %*% beeta
#   
#   y0hat <- X_tr %*% beeta + u_co %*% aa + matrix(alpha[1], T_obs,1) + Xi + matrix(mu, T_obs, 1)
#   
#   
#   ## result
#   msetable.cv <- matrix(0,D*D,3)
#   colnames(msetable.cv) <- c("lam1","lam2","mse")
#   count <- 1
#   for (i in 1:D){
#     for (j in 1:D){
#       msetable.cv[count,] <- c(lam1list[i],lam2list[j],msetable[i,j])
#       count <- count+1
#     }
#   }
#   
#   
#   result <- list( y0hat = y0hat,
#                   beta = beeta,
#                   a = aa,
#                   msetable.cv = msetable.cv)
#   return(result)
# }
# 

## Comparison with simplicity 

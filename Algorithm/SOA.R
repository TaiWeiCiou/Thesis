library(dplyr);library(ggplot2);library(caret);library(glmnet);library(knitr)#kable()
library(gurobi);library(Matrix);library(leaps);library(gsynth);library(Synth)
library(xtable)


## Importing code of dgp
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/DGP.R",
       encoding = "utf-8")

T_obs = 30; T0 = 20; N_ctrl = 10; N_treat = 1; N_x = 2

set.seed(42)
dat1 = DGP1(T_obs = 30, T0 = 20, N_ctrl = 10)$datXu

## core function
regularloss <- function(par, Y, X, p, lam){
  b <- matrix(par, p, 1) 
  loss <- 0.5*sum((Y - X %*% b)^2) + lam * sum(abs(b))
} 

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

# a includes intercept
duols2 <- function(par, Y1, Y0, X1, X0, N_x, N_ctrl, T0){
  beta <- as.matrix(par[1:N_x], N_x, 1)
  a <- as.matrix(par[(N_x+1):(N_ctrl+N_x+1)], (N_ctrl+1), 1)
  u1 <- Y1 - X1 %*% beta
  u00 <- Y0 - X0 %*% beta
  u0 <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u0[,i] <- u00[((i-1)*T_obs +1):(i*T_obs)]
  }
  u0 <- cbind(rep(1,T_obs), u0)
  target <- sum((Y0 - X0 %*% beta)^2) + sum( (u1 - u0[(1:T0),] %*% a )^2)
}


duols.panalty <- function(par, 
                          Y1, Y0, X1, X0, N_x, N_ctrl, T0,
                          lam1, lam2){
  
  beta <- as.matrix(par[1:N_x], N_x, 1)
  a <- as.matrix(par[(N_x+1):(N_ctrl+N_x+1)], (N_ctrl+1), 1)
  
  u1 <- Y1 - X1 %*% beta
  u00 <- Y0 - X0 %*% beta
  u0 <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u0[,i] <- u00[((i-1)*T_obs +1):(i*T_obs)]
  }
  u0 <- cbind(rep(1, T_obs), u0)
  
  target <- 0.5*sum((Y0 - X0 %*% beta)^2) + 
    0.5*sum( (u1 - u0[(1:T0),] %*% a )^2) +
    0.5*lam1 * sum(abs(beta)) + 0.5*lam2* sum(abs(a))
  
  return(target)
}

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


simuoptim <- function(dat, T_obs, T0, N_ctrl, N_x){
  
  ## parsing data
  Y1 <- as.matrix(dat[(dat$index == 1),3][1:T0], T0, 1)
  Y0 <- as.matrix(dat[(dat$index != 1),3], T_obs*N_ctrl, 1)
  X1 <- as.matrix(
    dat[(dat$index == 1),c(4,5)][(1:T0),], T0, N_x)
  X0 <- as.matrix(
    dat[(dat$index != 1),c(4,5)], T_obs*N_ctrl, N_x)
  
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
  
  X_tr <-as.matrix(
    dat[(dat$index == 1),c(4,5)], T_obs, N_x
    )
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  
  y0hat <- X_tr %*% beeta + u_co %*% aa
  
  
  ## result
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa)
  return(result)
}

simuoptim.intercept <- function(dat, T_obs, T0, N_ctrl, N_x){
  
  ## parsing data
  Y1 <- as.matrix(dat[(dat$index == 1),3][1:T0], T0, 1)
  Y0 <- as.matrix(dat[(dat$index != 1),3], T_obs*N_ctrl, 1)
  X1 <- as.matrix(
    dat[(dat$index == 1),c(4,5)][(1:T0),], T0, N_x)
  X0 <- as.matrix(
    dat[(dat$index != 1),c(4,5)], T_obs*N_ctrl, N_x)
  
  ## setting parameter
  parr <- rep(0,(N_x+N_ctrl+1))
  
  ## optimiztion
  optresult <- optim(par = parr,
                     fn = duols2,
                     Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
                     N_x = 2, N_ctrl = N_ctrl,T0 = T0,
                     method = "SANN")
  
  ## claiming estimated coefficient
  beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
  aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl+1)], (N_ctrl+1), 1)
  
  ## constructing the counter-factual at treatment period
  
  X_tr <-as.matrix(
    dat[(dat$index == 1),c(4,5)], T_obs, N_x
  )
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  u_co <- cbind(rep(1,T_obs), u_co)
  
  y0hat <- X_tr %*% beeta + u_co %*% aa
  
  
  ## result
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa)
  return(result)
}

simuoptim.simplex <- function(dat, T_obs, T0, N_ctrl, N_x){
  
  ## parsing data
  Y1 <- as.matrix(dat[(dat$index == 1),3][1:T0], T0, 1)
  Y0 <- as.matrix(dat[(dat$index != 1),3], T_obs*N_ctrl, 1)
  X1 <- as.matrix(
    dat[(dat$index == 1),c(4,5)][(1:T0),], T0, N_x)
  X0 <- as.matrix(
    dat[(dat$index != 1),c(4,5)], T_obs*N_ctrl, N_x)
  
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
  
  X_tr <-as.matrix(
    dat[(dat$index == 1),c(4,5)], T_obs, N_x
  )
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  y0hat <- X_tr %*% beeta + u_co %*% aa
  
  
  ## result
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa)
  return(result)
}


lam.max <- function(y,X){
  max(abs(t(y) %*% X))
}
which.min.mat <- function(A){
  return(which(A == min(A), arr.ind = T))
}
nfoldCV <- function(n, n_fold = 5){
  quo <- n %/% n_fold; r <- n %% n_fold
  lenlist <- c(rep((quo+1), r), rep(quo, (n_fold-r))) 
  lenlist2 <- c(0,cumsum(lenlist))
  sam <- sample(1:n,n); 
  indicelist <- list()
  for (i in 1:n_fold){
    indicelist[[i]] <- sam[(lenlist2[i]+1) : (lenlist2[i+1])]
  }
  return(indicelist)
}

simuoptim.duolasso <- function(dat, T_obs, T0, N_ctrl, N_x,
                               D = 10, validnum = 3){
  
  N <- N_ctrl + 1
  
  ## parsing data
  Y1 <- as.matrix(dat[(dat$index == 1),3][1:T0], T0, 1)
  Y0 <- as.matrix(dat[(dat$index != 1),3], T_obs*N_ctrl, 1)
  X1 <- as.matrix(
    dat[(dat$index == 1),c(4,5)][(1:T0),], T0, N_x)
  X0 <- as.matrix(
    dat[(dat$index != 1),c(4,5)], T_obs*N_ctrl, N_x)
  
  ## setting parameter
  parr <- rep(0,(N_x+N_ctrl))
  
  ## generating lambda list
  lam1.max <- lam.max(y = Y0, X = X0)
  lam2.max <- lam.max(y = Y1, X = Y0[(1:T0),])
  lam1list <- exp((1:D)*(log(lam1.max)/D))
  lam2list <- exp((1:D)*(log(lam2.max)/D)) 

  ## CV
  
  msetable <- matrix(0,D,D)
  
  t0 <- Sys.time()
  sam <- sample(2:N, validnum)
  for (fold in sam){
    
    pool <- c(2:N)[-(fold-1)]
    
    Y1 <- as.matrix(dat[(dat$index == fold),3][1:T0], T0, 1)
    Y0 <- as.matrix(dat[(dat$index %in% pool),3], T_obs*(N_ctrl-1), 1)
    X1 <- as.matrix(
      dat[(dat$index == fold),c(4,5)][(1:T0),], T0, N_x)
    X0 <- as.matrix(
      dat[(dat$index %in% pool),c(4,5)], T_obs*(N_ctrl-1), N_x)
    
    for (ii in 1:D){
      for (jj in 1:D){
        
        optresult <- optim(par = rep(0,(N_x+ N_ctrl)),
                           fn = duols.panalty,
                           Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
                           N_x = 2, N_ctrl = N_ctrl-1,T0 = T0,
                           method = "SANN",
                           lam1 = lam1list[ii], lam2 = lam2list[jj])
        
        beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
        aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl)], N_ctrl, 1)
        
        X_tr <-as.matrix(
          dat[(dat$index == fold),c(4,5)], T_obs, N_x
        )
        
        u_co0 <- Y0 - X0 %*% beeta
        u_co <- matrix(0, T_obs, N_ctrl-1)
        for (i in 1:(N_ctrl-1)){
          u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
        }
        u_co <- cbind(rep(1,T_obs), u_co)
        
        y0hat <- X_tr %*% beeta + u_co %*% aa
        
        y0 <- dat[(dat$index == fold),3]
        
        msetable[ii,jj] <- msetable[ii,jj] + mean(((y0-y0hat)^2)[((T0+1):T_obs)])
        
      }
    }
    
  }
  t1 <- Sys.time()
  t1-t0
  
  i.min <- which.min.mat(msetable)[1];j.min <- which.min.mat(msetable)[2]
  
  ## parsing data
  Y1 <- as.matrix(dat[(dat$index == 1),3][1:T0], T0, 1)
  Y0 <- as.matrix(dat[(dat$index != 1),3], T_obs*N_ctrl, 1)
  X1 <- as.matrix(
    dat[(dat$index == 1),c(4,5)][(1:T0),], T0, N_x)
  X0 <- as.matrix(
    dat[(dat$index != 1),c(4,5)], T_obs*N_ctrl, N_x)
  
  
  ## optimiztion
  optresult <- optim(par = rep(0,(N_x+N_ctrl+1)),
                     fn = duols.panalty,
                     Y1 = Y1, Y0 = Y0,X1 = X1, X0 = X0,
                     N_x = 2, N_ctrl = N_ctrl,T0 = T0,
                     method = "SANN",
                     lam1 = lam1list[i.min],lam2 = lam2list[j.min])
  
  ## claiming estimated coefficient
  beeta <- as.matrix(optresult$par[1:N_x], N_x, 1)
  aa <- as.matrix(optresult$par[(N_x+1):(N_x+N_ctrl+1)], (N_ctrl+1), 1)
  
  ## constructing the counter-factual at treatment period
  
  X_tr <-as.matrix(
    dat[(dat$index == 1),c(4,5)], T_obs, N_x
  )
  
  u_co0 <- Y0 - X0 %*% beeta
  u_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    u_co[,i] <- u_co0[((i-1)*T_obs+1):(i*T_obs)]
  }
  u_co <- cbind(rep(1,T_obs),u_co)
  
  y0hat <- X_tr %*% beeta + u_co %*% aa
  
  
  ## result
  msetable.cv <- matrix(0,D*D,3)
  colnames(msetable.cv) <- c("lam1","lam2","mse")
  count <- 1
  for (i in 1:D){
    for (j in 1:D){
      msetable.cv[count,] <- c(lam1list[i],lam2list[j],msetable[i,j])
      count <- count+1
    }
  }
  
  
  result <- list( y0hat = y0hat,
                  beta = beeta,
                  a = aa,
                  msetable.cv = msetable.cv)
  return(result)
}


mygsynth.duolasso <- function(dat1, T0, T_obs, N_treat = 1, N_ctrl, N_x ){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##--------------------##
  ##demeaning the data 
  ##--------------------##
  datt <- dat1[dat1$index %in% c((N_treat+1):N), c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N_ctrl*T_obs, (N_x+1))) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in (N_treat+1):N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  colnames(datt0) = c("Time","index","Outcome","X1","X2")
  
  dat_tr <- dat1[(dat1$index == 1), c(1,2,3,3+(1:N_x))]
  dat_tr0 <- dat_tr
  obs_imean_tr <- apply(dat_tr[(dat_tr$index == 1 & dat_tr$Time <= T0 ),c(3,3+(1:N_x))], 2, mean)
  for (t in c(1:T0)){
    obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
    obs_dot_tr <- dat_tr[(dat_tr$index == 1 & dat_tr$Time == t),c(3,3+(1:N_x))] - 
      obs_imean_tr - obs_tmean + obs_itmean
    dat_tr0[(dat_tr0$index == 1 & dat_tr0$Time == t),c(3,3+(1:N_x))] <- obs_dot_tr
  }
  colnames(dat_tr0) = c("Time","index","Outcome","X1","X2")

  dat  = rbind(dat_tr0, datt0)


  ##-------------------##
  ## core algo
  ##-------------------##
  
  optresult = simuoptim2(dat = dat, T_obs = T_obs, T0 = T0, N_ctrl = N_ctrl, N_x = N_x)
  
  beta1 = optresult$beta
  
  ##--------------------##
  ## fixed effect
  ##--------------------##
  
  ### mu
  y_dd_bar <- mean(datt[(datt$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(datt[(datt$index %in% c((N_treat+1):N)),3+(1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean)
    xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific effect
  #### for treated unit
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    y_imean0 <- mean(dat1[(dat1$index == i & dat1$Time <= T0),3])
    x_imean0 <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)],2,mean) %>% as.matrix()
    A[i,] <- y_imean0 - t(x_imean0) %*% beta1 - 
      mean(Xi[1:T0])-
      mu 
  }
  
  #### for control units
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - sum(x_imean*beta1) - mu
    A[i,] <- alpha_i
  }
  
  ## estimate the y0
  y0_treat_hat <- matrix(0, T_obs, N_treat)
  for (i in 1:N_treat){
    datX_treat <- dat1[(dat1$index == i), 3+(1:N_x)] %>% as.matrix()
    y0_treat_hat_i <-  optresult$y0hat+
      matrix(mu, T_obs,1)+
      Xi+ 
      matrix(A[i,], T_obs, 1)
    y0_treat_hat[,i] <- y0_treat_hat_i 
  }
      
  resultlist <- list(beta = beta1,
                     Y0_hat = y0_treat_hat,
                     mu = mu,
                     alpha = A,
                     xi = Xi,
                     container = optresult)
  return(resultlist)
}

msetable0 <- data.frame(matrix(0,3,3))
msetable <- data.frame(matrix(0,3,3))
colnames(msetable) <- c("mse_pre","mse_aft","mse_ratio")
row.names(msetable) <- c("duols","simplex","mean")

t0 <- Sys.time()
B = 100
for (i in 1:B){
  ## generating a data for testing.
  dat0 <- DGP1(T_obs = 30, T0 = 20, N_ctrl = 10)
  
  ## setting parameters
  T_obs = 30; T0 = 20; N_ctrl  =10
  
  
  aaa <- simuoptim.intercept(dat = dat0$dat,
                             T_obs = 30, T0 = 20, N_ctrl = 10, N_x = 2)
  aaaa <- simuoptim.simplex(dat = dat0$dat,
                            T_obs = 30, T0 = 20, N_ctrl = 10, N_x = 2)
  
  ymean = rep(mean(dat0$dat[1:T0,3]),T_obs)
  
  
  y0 <- dat0$datY0[,1]
  
  y0hat <- cbind(aaa$y0hat, aaaa$y0hat, ymean)
  
  mse <- apply(y0hat, 2, function(x) (x - y0)^2)
  mse_pre <- apply( mse[1:T0,], 2, mean); mse_aft <- apply( mse[((T0+1):T_obs),], 2, mean)
  mse_ratio <- mse_aft/mse_pre
  msetable0 <- msetable0 + t(rbind(mse_pre, mse_aft, mse_ratio))
  msetable <- msetable0/i
}
t1 <- Sys.time()
t1-t0

View(msetable)


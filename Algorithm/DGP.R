##--------------------##
## Hsiao, Zhou(2018) DGP
##--------------------##
T_obs <- 30; T0 <- 20; T_aft <- T_obs - T0
N_treat <- 1; N_ctrl <- 10; N <- N_treat +N_ctrl;

u_generate <- function(T_obs, N){
  
  v <- matrix(0, T_obs, N)
  for (i in 1:N){
    sig_i <- 0.5*(rchisq(1, df = 1) +1)
    v_i <- rnorm(T_obs, sd = sqrt(sig_i))
    v[,i] <- v_i
  }
  v <- cbind(v[,N],v,v[,i])
  
  u <- matrix(0,T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      u[t,i] <- 2*v[t,i+1] + v[t,i] +v[t,i+2]
    }
  }
  
  return(u)
}

DGP1 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- X1[t,i]*beeta[1,] + X2[t,i]*beeta[2,] + Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
         rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}

#dat <- DGP7(T_obs = 30, T0 = 20, N_ctrl = 10)

DGP2 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2) 
    c2 <- runif(1, min = 1, max = 2)
    c3 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1  +
        c1*Fac[t,1] + c2*Fac[t,2] + c3*Fac[t,3] +
        rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2  +
        c1*Fac[t,1] + c2*Fac[t,2] + c3*Fac[t,3] +
        rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- X1[t,i]*beeta[1,] + X2[t,i]*beeta[2,] + Fac[t,1:2] %*% gama[1:2,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}

DGP2.pure <- function(T_obs, T0, N_ctrl, 
                      N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2) 
    c2 <- runif(1, min = 1, max = 2)
    c3 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1  +
        c1*Fac[t,1] + c2*Fac[t,2] + c3*Fac[t,3] +
        rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2  +
        c1*Fac[t,1] + c2*Fac[t,2] + c3*Fac[t,3] +
        rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <-  Fac[t,1:2] %*% gama[1:2,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}


DGP3 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 2){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    rhoa <- runif(1, min = 0.1, max = 0.9)
    rhob <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1); init2 <- rnorm(1)
    inita <- rnorm(1); initb <- rnorm(1)
    
    for (t in 1:T_obs){
      
      etaa <- rnorm(1)
      x1 <- 1 + rho1*init1  + etaa + rhoa*inita
      X1[t,i] <- x1
      init1 <- x1; inita <- etaa
      
      etab <- rnorm(1)
      x2 <- 1 + rho2*init2  + etab + rhob*initb
      X2[t,i] <- x2
      init2 <- x2; initb <- etab
    }
  }
  
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- X1[t,i]*beeta[1,] + X2[t,i]*beeta[2,] + Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 10))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 10))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","gama1","gama2","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","gama1","gama2","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}

DGP3.pure <- function(T_obs, T0, N_ctrl, 
                      N_treat = 1, N_fac  = 2){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    rhoa <- runif(1, min = 0.1, max = 0.9)
    rhob <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1); init2 <- rnorm(1)
    inita <- rnorm(1); initb <- rnorm(1)
    
    for (t in 1:T_obs){
      
      etaa <- rnorm(1)
      x1 <- 1 + rho1*init1  + etaa + rhoa*inita
      X1[t,i] <- x1
      init1 <- x1; inita <- etaa
      
      etab <- rnorm(1)
      x2 <- 1 + rho2*init2  + etab + rhob*initb
      X2[t,i] <- x2
      init2 <- x2; initb <- etab
    }
  }
  
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 10))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 10))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","gama1","gama2","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","gama1","gama2","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}



DGP4 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  
  beeta <- matrix(c(1,2), 2, 1)
  
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  
  Fac <- matrix(0, T_obs, N_fac)
  init <- rnorm(1)
  for (t in 1:T_obs){
    f2 <- init + rnorm(1)
    f1 <- 0.5*f2 + rnorm(1)
    f3 <- rnorm(1)
    Fac[t,] <- c(f1,f2,f3)
    init <- f2
  }
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- X1[t,i]*beeta[1,] + X2[t,i]*beeta[2,] + Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}


DGP4.pure <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  
  beeta <- matrix(c(1,2), 2, 1)
  
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  
  Fac <- matrix(0, T_obs, N_fac)
  init <- rnorm(1)
  for (t in 1:T_obs){
    f2 <- init + rnorm(1)
    f1 <- 0.5*f2 + rnorm(1)
    f3 <- rnorm(1)
    Fac[t,] <- c(f1,f2,f3)
    init <- f2
  }
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <-  Fac[t,1:2] %*% gama[1:2,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}




DGP5 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 3){
  
  N <- N_treat +N_ctrl
  
  U1 <- u_generate(T_obs = T_obs, N = N_ctrl)
  
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  
  beeta <- matrix(c(1,2), 2, 1)
  
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    for (t in 1:T_obs){
      Y_co[t,i] <- X1[t,i]*beeta[1,] + X2[t,i]*beeta[2,] + Fac[t,] %*% gama[,i] + U1[t,i]
    }
  }
  
  U2 <- u_generate(T_obs = T_obs, N = 1)
  
  beata <- matrix(runif(N_ctrl), N_ctrl,1)
  Y_tr <-  Y_co %*% beata + U2
  
  Y0 <- cbind(Y_tr, Y_co)
  Y <- Y0 + D*delta_bar
  
  
  ## data format parsing
  
  paneldat <- data.frame(matrix(0, T_obs*N, 12))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 12))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","f3","gama1","gama2","gama3","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}


DGP6 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 2){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(rnorm(T_obs*N_fac), T_obs, N_fac)
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 10))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 10))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","gama1","gama2","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","gama1","gama2","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}

DGP7 <- function(T_obs, T0, N_ctrl, 
                 N_treat = 1, N_fac  = 2){
  
  N <- N_treat +N_ctrl
  
  U <- u_generate(T_obs = T_obs, N = N)
  D <- matrix(0, T_obs, N)
  D[((T0+1):T_obs),(1:N_treat)] <-matrix(1,(T_obs-T0),N_treat) 
  delta_bar <- matrix(0, T_obs, N)
  eff <- c(1:10) + rnorm(10)
  delta_bar[((T0+1):T_obs),(1:N_treat)] <- matrix(rep(eff,N_treat), (T_obs - T0), N_treat)
  beeta <- matrix(c(1,2), 2, 1)
  gama <- matrix(rnorm(N_fac*N), N_fac, N )
  Fac <- matrix(0, T_obs, N_fac)
  
  initf1 <- rnorm(1);initf2 <- rnorm(1);
  for (t in 1:T_obs){
    Fac[t,1] <- initf1 +rnorm(1); Fac[t,2] <- initf2 +rnorm(1)
    initf1 <- Fac[t,1]; initf2 <- Fac[t,2]
  }
  
  X1 <- matrix(0, T_obs, N); X2 <- matrix(0, T_obs, N)
  for (i in 1:N){
    
    c1 <- runif(1, min = 1, max = 2); c2 <- runif(1, min = 1, max = 2)
    
    rho1 <- runif(1, min = 0.1, max = 0.9)
    rho2 <- runif(1, min = 0.1, max = 0.9)
    
    init1 <- rnorm(1) 
    init2 <- rnorm(1)
    
    for (t in 1:T_obs){
      
      x1 <- 1 + rho1*init1 + c1*gama[1,i] +c2*Fac[t,1] + rchisq(1, df = 1) - 1
      X1[t,i] <- x1
      init1 <- x1
      
      x2 <- 1 + rho2*init2 + c1*gama[2,i] +c2*Fac[t,2] + rchisq(1, df = 1) - 1
      X2[t,i] <- x2
      init2 <- x2
    }
  }
  
  Y0 <- matrix(0, T_obs, N)
  Y <- matrix(0, T_obs, N)
  for (i in 1:N){
    for (t in 1:T_obs){
      Y0[t,i] <- Fac[t,] %*% gama[,i] + U[t,i]
      Y[t,i] <- Y0[t,i] + D[t,i]*delta_bar[t,i]
    }
  }
  
  paneldat <- data.frame(matrix(0, T_obs*N, 10))
  paneldat0 <- data.frame(matrix(0, T_obs*N, 10))
  count <- 1
  for ( i in 1:N){
    for (t in 1:T_obs){
      paneldat[count,] <- c(t,i,Y[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      paneldat0[count,] <- c(t,i,Y0[t,i],X1[t,i],X2[t,i],Fac[t,],gama[,i],D[t,i])
      count <- count +1
    }
  }
  colnames(paneldat) <- c("Time","index","Outcome","X1","X2","f1","f2","gama1","gama2","D")
  colnames(paneldat0) <- c("Time","index","Y0","X1","X2","f1","f2","gama1","gama2","D")
  
  ## parsing data for adh method
  unitcol <- c()
  for (i in 1:N){
    unitcol <- c(unitcol, rep(paste(i),T_obs))
  }
  unitcol <- as.matrix(unitcol, T_obs*N, 1)
  datADH <- cbind(paneldat[,1:5], unitcol); colnames(datADH)[6] <- "Unit"
  
  ## parsing data for Xu's method
  
  Dd <- c( rep(c(rep(0,T0), rep(1,(T_obs -T0))), N_treat),
           rep(0, T_obs*N_ctrl))
  
  datXu <- cbind(paneldat[,1:5], Dd)
  colnames(datXu)[6] <- "D"
  
  resultlist = list(dat = paneldat,
                    dat0 = paneldat0,
                    datY = Y,
                    datY0 = Y0,
                    datX1 = X1, 
                    datX2 = X2,
                    datADH = datADH,
                    datXu = datXu)
  
  return(resultlist)
}

##--------------------##
## Xu(2017) DGP
##--------------------##

IFEgen <- function(T0,T_obs,N_treat = 1, N_ctrl, sig = 1){
  ## Intercept
  Intercept = 5
  
  ## treatment effect across units by time
  deltabar_t = c(rep(0,T0),c(1:(T_obs - T0)))
  e_it = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  delta_it = matrix(0, T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat+N_ctrl)){
    delta_it[,i] = deltabar_t + e_it[,i]
  }
  
  ## treatment indicators
  D_it_0 = matrix(0,T_obs, N_treat + N_ctrl)
  D_it_1 = matrix(0,T_obs, N_treat + N_ctrl)
  
  for (s in (T0+1):T_obs){
    D_it_1[s,1:N_treat] = 1
  }
  
  ## unobserved factors and factor loading
  sqrt3 = 3^0.5
  w = 0.8
  
  ### factor loading
  lambda_i_1 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3), #treated unit's factor loading
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3) #control unit's factor loading
    )
  lambda_i_2 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3),
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3)
    )
  
  ### unobserved factors
  f_t_1 = rnorm(T_obs)
  f_t_2 = rnorm(T_obs)
  
  ## observed covariates
  beta1 = 1
  beta2 = 3
  
  eta_it_1 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  eta_it_2 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  X_it_1 = matrix(NA, T_obs, N_treat + N_ctrl)
  X_it_2 = matrix(NA, T_obs, N_treat + N_ctrl)
  
  for (i in 1:(N_treat + N_ctrl)){
    for(t in 1:T_obs){
      X_it_1[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_1[t,i]
      X_it_2[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_2[t,i]
    }
  } 
  
  ## unit fixed effect
  alpha_i = rnorm(N_treat+N_ctrl)
  
  ## time fixed effect
  xi_t = rnorm(T_obs)
  
  ## white noise
  epison_it = matrix(rnorm((N_treat+N_ctrl)*T_obs, sd = sig^(1/2)), T_obs, N_treat+N_ctrl)
  
  ## outcome
  Y_it_0 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_0[s,i] = delta_it[s,i]*D_it_0[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] +Intercept
    }
  }
  
  
  Y_it_1 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_1[s,i] = delta_it[s,i]*D_it_1[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] ++Intercept
    }
  }
  
  #ggplot(data = as.data.frame( cbind(Y_it_0[,1],Y_it_1[,1])), aes(x = 1:30)) +
  #  geom_line(aes(y = V1))+
  #  geom_line(aes(y = V2), col = "red")
  
  data_sim1 = as.data.frame(matrix(0,T_obs*(N_treat + N_ctrl),5))
  colnames(data_sim1) = c("Time", "index", "Outcome", "X1", "X2")
  Unitname = c(paste0("V",1:(N_treat+N_ctrl)))
  c = 1
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      data_sim1[c,] = c(s,i, Y_it_1[s,i], X_it_1[s,i], X_it_2[s,i] )
      c = c+1
    }
  }
  data_sim1 = apply(data_sim1, 2, as.numeric) %>% as.data.frame()
  data_sim1$Unit = as.character(data_sim1$index)
  
  D <- c()
  for (i in 1:N_treat){
    D <- c(D, rep(0, T0), rep(1, T_obs-T0))
  }
  for (i in (N_treat+1):(N_treat + N_ctrl)){
    D <- c(D, rep(0,T_obs))
  }
  eff <- c()
  for (i in 1:(N_treat+N_ctrl)){
    eff <- c(eff, rep(0,T0), delta_it[(T0+1):T_obs,i])
  }
  err <- c()
  for (i in 1:(N_treat+N_ctrl)){
    err <- c(err, epison_it[,i])
  }
  
  mu <- rep(5, (N_treat+N_ctrl)*T_obs)
  
  alpha <- c()
  for (i in 1:(N_treat + N_ctrl)){
    alpha <- c(alpha, rep(alpha_i[i], T_obs))
  }
  xi <- rep(xi_t, (N_treat + N_ctrl))
  f1 <- rep(f_t_1, (N_treat + N_ctrl))
  f2 <- rep(f_t_2, (N_treat + N_ctrl))
  
  l1 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l1 <- c(l1, rep(lambda_i_1[i], T_obs))
  }
  
  l2 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l2 <- c(l2, rep(lambda_i_2[i], T_obs))
  }
  
  data_sim2 <- cbind(data_sim1[,-6], eff, err, mu, alpha, xi, f1, l1, f2, l2, D)
  
  
  df = as.data.frame(matrix(NA,T_obs,(N_treat + N_ctrl)))
  
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      df[s,i] = data_sim1[
        which(data_sim1$Time == s & data_sim1$Unit == i),3]
    }
  }
  data = list(
    adh.data = data_sim1,
    su.data = data_sim2,
    y.data = df,
    x1.data = X_it_1,
    x2.data = X_it_2,
    y0.data = Y_it_0
  )
  return(data)
}
AR1gen <- function(T0,T_obs,N_treat = 1, N_ctrl,rho = 0.9, sig = 1){
  ## Intercept
  Intercept = 5
  
  ## treatment effect across units by time
  deltabar_t = c(rep(0,T0),c(1:(T_obs - T0)))
  e_it = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  delta_it = matrix(0, T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat+N_ctrl)){
    delta_it[,i] = deltabar_t + e_it[,i]
  }
  
  ## treatment indicators
  D_it_0 = matrix(0,T_obs, N_treat + N_ctrl)
  D_it_1 = matrix(0,T_obs, N_treat + N_ctrl)
  
  for (s in (T0+1):T_obs){
    D_it_1[s,1:N_treat] = 1
  }
  
  ## unobserved factors and factor loading
  sqrt3 = 3^0.5
  w = 0.8
  
  ### factor loading
  lambda_i_1 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3), #treated unit's factor loading
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3) #control unit's factor loading
    )
  lambda_i_2 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3),
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3)
    )
  
  ### unobserved factors
  f_t_1_0 = rnorm(1)
  f_t_1 = rep(0,T_obs)
  for (s in 1:T_obs) {
    f_t_1[s] = rho*f_t_1_0 + rnorm(1)
    f_t_1_0 = f_t_1[s]
  }
  
  f_t_2 = rnorm(T_obs)
  
  ## observed covariates
  beta1 = 1
  beta2 = 3
  
  eta_it_1 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  eta_it_2 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  X_it_1 = matrix(NA, T_obs, N_treat + N_ctrl)
  X_it_2 = matrix(NA, T_obs, N_treat + N_ctrl)
  
  for (i in 1:(N_treat + N_ctrl)){
    for(t in 1:T_obs){
      X_it_1[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_1[t,i]
      X_it_2[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_2[t,i]
    }
  } 
  
  ## unit fixed effect
  alpha_i = rnorm(N_treat+N_ctrl)
  
  ## time fixed effect
  xi_t = rnorm(T_obs)
  
  ## white noise
  epison_it = matrix(rnorm((N_treat+N_ctrl)*T_obs, sd = sqrt(sig)), T_obs, N_treat+N_ctrl)
  
  ## outcome
  Y_it_0 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_0[s,i] = delta_it[s,i]*D_it_0[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] +Intercept
    }
  }
  
  
  Y_it_1 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_1[s,i] = delta_it[s,i]*D_it_1[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] ++Intercept
    }
  }
  
  #ggplot(data = as.data.frame( cbind(Y_it_0[,1],Y_it_1[,1])), aes(x = 1:30)) +
  #  geom_line(aes(y = V1))+
  #  geom_line(aes(y = V2), col = "red")
  
  data_sim1 = as.data.frame(matrix(0,T_obs*(N_treat + N_ctrl),5))
  colnames(data_sim1) = c("Time", "index", "Outcome", "X1", "X2")
  Unitname = c(paste0("V",1:(N_treat+N_ctrl)))
  c = 1
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      data_sim1[c,] = c(s,i, Y_it_1[s,i], X_it_1[s,i], X_it_2[s,i] )
      c = c+1
    }
  }
  data_sim1 = apply(data_sim1, 2, as.numeric) %>% as.data.frame()
  data_sim1$Unit = as.character(data_sim1$index)
  
  D <- c()
  for (i in 1:N_treat){
    D <- c(D, rep(0, T0), rep(1, T_obs-T0))
  }
  for (i in (N_treat+1):(N_treat + N_ctrl)){
    D <- c(D, rep(0,T_obs))
  }
  eff <- c()
  for (i in 1:(N_treat+N_ctrl)){
    eff <- c(eff, rep(0,T0), delta_it[(T0+1):T_obs,i])
  }
  err <- c()
  for (i in 1:(N_treat+N_ctrl)){
    err <- c(err, epison_it[,i])
  }
  
  mu <- rep(5, (N_treat+N_ctrl)*T_obs)
  
  alpha <- c()
  for (i in 1:(N_treat + N_ctrl)){
    alpha <- c(alpha, rep(alpha_i[i], T_obs))
  }
  xi <- rep(xi_t, (N_treat + N_ctrl))
  f1 <- rep(f_t_1, (N_treat + N_ctrl))
  f2 <- rep(f_t_2, (N_treat + N_ctrl))
  
  l1 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l1 <- c(l1, rep(lambda_i_1[i], T_obs))
  }
  
  l2 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l2 <- c(l2, rep(lambda_i_2[i], T_obs))
  }
  
  data_sim2 <- cbind(data_sim1[,-6], eff, err, mu, alpha, xi, f1, l1, f2, l2, D)
  
  
  df = as.data.frame(matrix(NA,T_obs,(N_treat + N_ctrl)))
  
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      df[s,i] = data_sim1[
        which(data_sim1$Time == s & data_sim1$Unit == i),3]
    }
  }
  data = list(
    adh.data = data_sim1,
    su.data = data_sim2,
    y.data = df,
    x1.data = X_it_1,
    x2.data = X_it_2,
    y0.data = Y_it_0
  )
  return(data)
}
AR2gen <- function(T0,T_obs,N_treat = 1, N_ctrl,rho = 0.9, sig = 1){
  ## Intercept
  Intercept = 5
  
  ## treatment effect across units by time
  deltabar_t = c(rep(0,T0),c(1:(T_obs - T0)))
  e_it = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  delta_it = matrix(0, T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat+N_ctrl)){
    delta_it[,i] = deltabar_t + e_it[,i]
  }
  
  ## treatment indicators
  D_it_0 = matrix(0,T_obs, N_treat + N_ctrl)
  D_it_1 = matrix(0,T_obs, N_treat + N_ctrl)
  
  for (s in (T0+1):T_obs){
    D_it_1[s,1:N_treat] = 1
  }
  
  ## unobserved factors and factor loading
  sqrt3 = 3^0.5
  w = 0.8
  
  ### factor loading
  lambda_i_1 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3), #treated unit's factor loading
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3) #control unit's factor loading
    )
  lambda_i_2 = 
    c(
      runif(N_treat, min = -sqrt3, max = sqrt3),
      runif(N_ctrl, min = sqrt3 -2*w*sqrt3, max = 3*sqrt3 -2*w*sqrt3)
    )
  
  ### unobserved factors
  f_t_1_0 = rnorm(1)
  f_t_1_00 = rnorm(1)
  f_t_1 = rep(0,T_obs)
  for (s in 1:T_obs) {
    f_t_1[s] = rho*f_t_1_00 +(0.1*rho)*f_t_1_0 + rnorm(1)
    f_t_1_0 <- f_t_1_00
    f_t_1_00 <- f_t_1[s]
  }
  
  f_t_2 = rnorm(T_obs)
  
  ## observed covariates
  beta1 = 1
  beta2 = 3
  
  eta_it_1 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  eta_it_2 = matrix(rnorm( (N_treat + N_ctrl)*T_obs ) , T_obs, N_treat + N_ctrl)
  
  X_it_1 = matrix(NA, T_obs, N_treat + N_ctrl)
  X_it_2 = matrix(NA, T_obs, N_treat + N_ctrl)
  
  for (i in 1:(N_treat + N_ctrl)){
    for(t in 1:T_obs){
      X_it_1[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_1[t,i]
      X_it_2[t,i] = 1 + lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] + lambda_i_1[i] + lambda_i_2[i] + f_t_1[s] + f_t_2[s] + eta_it_2[t,i]
    }
  } 
  
  ## unit fixed effect
  alpha_i = rnorm(N_treat+N_ctrl)
  
  ## time fixed effect
  xi_t = rnorm(T_obs)
  
  ## white noise
  epison_it = matrix(rnorm((N_treat+N_ctrl)*T_obs, sd = sqrt(sig)), T_obs, N_treat+N_ctrl)
  
  ## outcome
  Y_it_0 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_0[s,i] = delta_it[s,i]*D_it_0[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] +Intercept
    }
  }
  
  
  Y_it_1 = matrix(NA,T_obs, N_treat + N_ctrl)
  for (i in 1:(N_treat + N_ctrl)){
    for (s in 1:T_obs) {
      Y_it_1[s,i] = delta_it[s,i]*D_it_1[s,i] +
        X_it_1[s,i] * beta1 + X_it_2[s,i] * beta2 + 
        lambda_i_1[i]*f_t_1[s] + lambda_i_2[i]*f_t_2[s] +
        alpha_i[i] + xi_t[s] + epison_it[s,i] ++Intercept
    }
  }
  
  #ggplot(data = as.data.frame( cbind(Y_it_0[,1],Y_it_1[,1])), aes(x = 1:30)) +
  #  geom_line(aes(y = V1))+
  #  geom_line(aes(y = V2), col = "red")
  
  data_sim1 = as.data.frame(matrix(0,T_obs*(N_treat + N_ctrl),5))
  colnames(data_sim1) = c("Time", "index", "Outcome", "X1", "X2")
  Unitname = c(paste0("V",1:(N_treat+N_ctrl)))
  c = 1
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      data_sim1[c,] = c(s,i, Y_it_1[s,i], X_it_1[s,i], X_it_2[s,i] )
      c = c+1
    }
  }
  data_sim1 = apply(data_sim1, 2, as.numeric) %>% as.data.frame()
  data_sim1$Unit = as.character(data_sim1$index)
  
  D <- c()
  for (i in 1:N_treat){
    D <- c(D, rep(0, T0), rep(1, T_obs-T0))
  }
  for (i in (N_treat+1):(N_treat + N_ctrl)){
    D <- c(D, rep(0,T_obs))
  }
  eff <- c()
  for (i in 1:(N_treat+N_ctrl)){
    eff <- c(eff, rep(0,T0), delta_it[(T0+1):T_obs,i])
  }
  err <- c()
  for (i in 1:(N_treat+N_ctrl)){
    err <- c(err, epison_it[,i])
  }
  
  mu <- rep(5, (N_treat+N_ctrl)*T_obs)
  
  alpha <- c()
  for (i in 1:(N_treat + N_ctrl)){
    alpha <- c(alpha, rep(alpha_i[i], T_obs))
  }
  xi <- rep(xi_t, (N_treat + N_ctrl))
  f1 <- rep(f_t_1, (N_treat + N_ctrl))
  f2 <- rep(f_t_2, (N_treat + N_ctrl))
  
  l1 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l1 <- c(l1, rep(lambda_i_1[i], T_obs))
  }
  
  l2 <- c()
  for (i in 1:(N_treat + N_ctrl)){
    l2 <- c(l2, rep(lambda_i_2[i], T_obs))
  }
  
  data_sim2 <- cbind(data_sim1[,-6], eff, err, mu, alpha, xi, f1, l1, f2, l2, D)
  
  
  df = as.data.frame(matrix(NA,T_obs,(N_treat + N_ctrl)))
  
  for (i in 1:(N_treat + N_ctrl)) {
    for (s in 1:T_obs){
      df[s,i] = data_sim1[
        which(data_sim1$Time == s & data_sim1$Unit == i),3]
    }
  }
  data = list(
    adh.data = data_sim1,
    su.data = data_sim2,
    y.data = df,
    x1.data = X_it_1,
    x2.data = X_it_2,
    y0.data = Y_it_0
  )
  return(data)
}



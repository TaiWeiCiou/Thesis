library(dplyr);library(ggplot2);library(caret);library(glmnet);library(knitr)#kable()
library(gurobi);library(Matrix);library(leaps);library(gsynth);library(Synth)
library(xtable)

beeta <- function(par, data){
  mat <- as.matrix(data); beta <- as.matrix(c(1-sum(par),par))
  maty <- mat[,1];matx <- mat[,-1]
  sum((maty - matx %*% beta)^2)
}

coreAlgo <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-14
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        fitB1 <- lm(Outcome ~ .+0, data = datB1)
        datB2 <- cbind(fitB1$residuals, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
      ## estimate the factor loading of the treated unit
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- datt0[(datt0$index == i), 3] %>% as.matrix()
        datX1 <- datt0[(datt0$index == i), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF) %>% as.data.frame()
        fit_i <- lm(V1~.+0, data = dat_treat)
        lambda_i <- coef(fit_i)
        L[,i] <- coef(fit_i)
      }
      
      ## estimate the y0
      y0_dot <- matrix(0, T_obs, N_treat)
      for (i in 1:N_treat){
        datX1 <- datt0[(datt0$index == i),3+(1:N_x)] %>% as.matrix()
        y0_dot_hat_i <- datX1 %*% beta1 + matF %*% L[,i] 
        y0_dot[,i] <- y0_dot_hat_i
      }
      
      y0_treat_hat <- matrix(0,T_obs, N_treat)
      y_itmean <- obs_itmean[1]
      for (i in 1:N_treat){
        y_imean <- mean(datt[(datt$index == i),3])
        for (t in 1:T_obs){
          y_tmean <- mean(datt[(datt$Time == t),3])
          y0_hat <- y0_dot[t,i] + y_imean + y_tmean -y_itmean
          y0_treat_hat[t,i] <- y0_hat
        }
      }
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
      y0_dot <- predict(fit_init, datt0[datt0$index %in% 1:N_treat, 3+(1:N_x)]) %>% as.matrix()
      
      y0_treat_hat <- matrix(0,T_obs, N_treat)
      y_itmean <- obs_itmean[1]
      for (i in 1:N_treat){
        y_imean <- mean(datt[(datt$index == i),3])
        for (t in 1:T_obs){
          y_tmean <- mean(datt[(datt$Time == t),3])
          y0_hat <- y0_dot[t,i] + y_imean + y_tmean -y_itmean
          y0_treat_hat[t,i] <- y0_hat
        }
      }
    }
    
    ## estimating the fixed effect
    
    ### mu
    y_dd_bar <- obs_itmean[1]
    x_dd_bar <- obs_itmean[-1] %>% as.matrix()
    mu <- y_dd_bar - t(x_dd_bar) %*% beta1
    
    ### unit-specific effect
    A <- matrix(0,N,1)
    for(i in 1:N){
      y_imean <- mean(dat1[dat1$index == i,3])
      x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
      alpha_i <- y_imean - sum(x_imean*beta1) - mu
      A[i,] <- alpha_i
    }
    
    ### time-varying effect
    Xi <- matrix(0,T_obs,1)
    for(t in 1:T_obs){
      y_tmean <- mean(dat1[dat1$Time == t,3])
      x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
      xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
      Xi[t,] <- xi_t
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        ## estimate the factor loading of the treated unit w/o s
        L <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- datt0[(datt0$index == i & datt0$Time != s), 3] %>% as.matrix()
          datX1 <- datt0[(datt0$index == i & datt0$Time != s), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[-s,]) %>% as.data.frame()
          fit_i <- lm(V1~.+0, data = dat_treat)
          lambda_i <- coef(fit_i)
          L[,i] <- coef(fit_i)
        }
        
        ## estimate the y0_s_hat
        y0_dot_s_hat <- matrix(0, 1, N_treat)
        for( i in 1:N_treat){
          datX_s <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_dot_s_hat_i <- datX_s %*% beta1 + matF[s,] %*% L[,i]
          y0_dot_s_hat[,i] <- y0_dot_s_hat_i
        }
        
        y0_s_hat <- matrix(0, 1, N_treat)
        y_itmean <- obs_itmean[1]
        y_smean <- mean(datt[(datt$Time == s),3])
        for (i in 1:N_treat){
          y_imean <- mean(datt[(datt$index == i),3])
          y0_s_hat_i <- y0_dot_s_hat + y_imean + y_tmean -y_itmean
          y0_s_hat[,i] <- y0_s_hat_i
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        ## estimate the y0_s_hat
        y0_dot_s_hat <- matrix(0, 1, N_treat)
        for( i in 1:N_treat){
          datX_s <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_dot_s_hat_i <- datX_s %*% beta1 
          y0_dot_s_hat[,i] <- y0_dot_s_hat_i
        }
        
        y0_s_hat <- matrix(0, 1, N_treat)
        y_itmean <- obs_itmean[1]
        y_smean <- mean(datt[(datt$Time == s),3])
        for (i in 1:N_treat){
          y_imean <- mean(datt[(datt$index == i),3])
          y0_s_hat_i <- y0_dot_s_hat + y_imean + y_tmean -y_itmean
          y0_s_hat[,i] <- y0_s_hat_i
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[r_opt+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}
coreAlgo0 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x, R){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-14
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        fitB1 <- lm(Outcome ~ .+0, data = datB1)
        datB2 <- cbind(fitB1$residuals, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
      ## estimate the factor loading of the treated unit
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- datt0[(datt0$index == i), 3] %>% as.matrix()
        datX1 <- datt0[(datt0$index == i), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF) %>% as.data.frame()
        fit_i <- lm(V1~.+0, data = dat_treat)
        lambda_i <- coef(fit_i)
        L[,i] <- coef(fit_i)
      }
      
      ## estimate the y0
      y0_dot <- matrix(0, T_obs, N_treat)
      for (i in 1:N_treat){
        datX1 <- datt0[(datt0$index == i),3+(1:N_x)] %>% as.matrix()
        y0_dot_hat_i <- datX1 %*% beta1 + matF %*% L[,i] 
        y0_dot[,i] <- y0_dot_hat_i
      }
      
      y0_treat_hat <- matrix(0,T_obs, N_treat)
      y_itmean <- obs_itmean[1]
      for (i in 1:N_treat){
        y_imean <- mean(datt[(datt$index == i),3])
        for (t in 1:T_obs){
          y_tmean <- mean(datt[(datt$Time == t),3])
          y0_hat <- y0_dot[t,i] + y_imean + y_tmean -y_itmean
          y0_treat_hat[t,i] <- y0_hat
        }
      }
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
      y0_dot <- predict(fit_init, datt0[datt0$index %in% 1:N_treat, 3+(1:N_x)]) %>% as.matrix()
      
      y0_treat_hat <- matrix(0,T_obs, N_treat)
      y_itmean <- obs_itmean[1]
      for (i in 1:N_treat){
        y_imean <- mean(datt[(datt$index == i),3])
        for (t in 1:T_obs){
          y_tmean <- mean(datt[(datt$Time == t),3])
          y0_hat <- y0_dot[t,i] + y_imean + y_tmean -y_itmean
          y0_treat_hat[t,i] <- y0_hat
        }
      }
    }
    
    ## estimating the fixed effect
    
    ### mu
    y_dd_bar <- obs_itmean[1]
    x_dd_bar <- obs_itmean[-1] %>% as.matrix()
    mu <- y_dd_bar - t(x_dd_bar) %*% beta1
    
    ### unit-specific effect
    A <- matrix(0,N,1)
    for(i in 1:N){
      y_imean <- mean(dat1[dat1$index == i,3])
      x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
      alpha_i <- y_imean - sum(x_imean*beta1) - mu
      A[i,] <- alpha_i
    }
    
    ### time-varying effect
    Xi <- matrix(0,T_obs,1)
    for(t in 1:T_obs){
      y_tmean <- mean(dat1[dat1$Time == t,3])
      x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
      xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
      Xi[t,] <- xi_t
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        ## estimate the factor loading of the treated unit w/o s
        L <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- datt0[(datt0$index == i & datt0$Time != s), 3] %>% as.matrix()
          datX1 <- datt0[(datt0$index == i & datt0$Time != s), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[-s,]) %>% as.data.frame()
          fit_i <- lm(V1~.+0, data = dat_treat)
          lambda_i <- coef(fit_i)
          L[,i] <- coef(fit_i)
        }
        
        ## estimate the y0_s_hat
        y0_dot_s_hat <- matrix(0, 1, N_treat)
        for( i in 1:N_treat){
          datX_s <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_dot_s_hat_i <- datX_s %*% beta1 + matF[s,] %*% L[,i]
          y0_dot_s_hat[,i] <- y0_dot_s_hat_i
        }
        
        y0_s_hat <- matrix(0, 1, N_treat)
        y_itmean <- obs_itmean[1]
        y_smean <- mean(datt[(datt$Time == s),3])
        for (i in 1:N_treat){
          y_imean <- mean(datt[(datt$index == i),3])
          y0_s_hat_i <- y0_dot_s_hat + y_imean + y_tmean -y_itmean
          y0_s_hat[,i] <- y0_s_hat_i
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        ## estimate the y0_s_hat
        y0_dot_s_hat <- matrix(0, 1, N_treat)
        for( i in 1:N_treat){
          datX_s <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_dot_s_hat_i <- datX_s %*% beta1 
          y0_dot_s_hat[,i] <- y0_dot_s_hat_i
        }
        
        y0_s_hat <- matrix(0, 1, N_treat)
        y_itmean <- obs_itmean[1]
        y_smean <- mean(datt[(datt$Time == s),3])
        for (i in 1:N_treat){
          y_imean <- mean(datt[(datt$index == i),3])
          y0_s_hat_i <- y0_dot_s_hat + y_imean + y_tmean -y_itmean
          y0_s_hat[,i] <- y0_s_hat_i
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[R+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}
coreAlgo0_revise <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x, R){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-14
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        fitB1 <- lm(Outcome ~ .+0, data = datB1)
        datB2 <- cbind(fitB1$residuals, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
      ## estimate the factor loading of the treated unit
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- datt0[(datt0$index == i), 3] %>% as.matrix()
        datX1 <- datt0[(datt0$index == i), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF) %>% as.data.frame()
        fit_i <- lm(V1~.+0, data = dat_treat)
        lambda_i <- coef(fit_i)
        L[,i] <- coef(fit_i)
      }
      
      
      ## estimating the fixed effect
      
      ### mu
      y_dd_bar <- obs_itmean[1]
      x_dd_bar <- obs_itmean[-1] %>% as.matrix()
      mu <- y_dd_bar - t(x_dd_bar) %*% beta1
      
      ### unit-specific effect
      A <- matrix(0,N,1)
      for(i in 1:N){
        y_imean <- mean(dat1[dat1$index == i,3])
        x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu
        A[i,] <- alpha_i
      }
      
      ### time-varying effect
      Xi <- matrix(0,T_obs,1)
      for(t in 1:T_obs){
        y_tmean <- mean(dat1[dat1$Time == t,3])
        x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
        xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
        Xi[t,] <- xi_t
      }
      
      y0_treat_hat <- matrix(0, T_obs, N_treat)
      for( i in 1:N_treat){
        datX_treat <- datt0[(datt0$index == i),3+(1:N_x)] %>% as.matrix()
        y0_treat_hat_i <- datX_treat %*% beta1 + 
          matF %*% L[,i] + 
          matrix(mu, T_obs, 1) +  
          matrix(A[i,], T_obs, 1) +  
          matrix(Xi, T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i
      } 
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
      ## estimating the fixed effect
      
      ### mu
      y_dd_bar <- obs_itmean[1]
      x_dd_bar <- obs_itmean[-1] %>% as.matrix()
      mu <- y_dd_bar - t(x_dd_bar) %*% beta1
      
      ### unit-specific effect
      A <- matrix(0,N,1)
      for(i in 1:N){
        y_imean <- mean(dat1[dat1$index == i,3])
        x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu
        A[i,] <- alpha_i
      }
      
      ### time-varying effect
      Xi <- matrix(0,T_obs,1)
      for(t in 1:T_obs){
        y_tmean <- mean(dat1[dat1$Time == t,3])
        x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
        xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
        Xi[t,] <- xi_t
      }
      
      y0_treat_hat <- matrix(0, T_obs, N_treat)
      for( i in 1:N_treat){
        datX_treat <- datt0[(datt0$index == i),3+(1:N_x)] %>% as.matrix()
        y0_treat_hat_i <- datX_treat %*% beta1 + 
          matrix(mu, T_obs, 1) +  
          matrix(A[i,], T_obs, 1) +  
          matrix(Xi, T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i
      } 
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        ## estimate the factor loading of the treated unit w/o s
        L_s <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- datt0[(datt0$index == i & datt0$Time != s), 3] %>% as.matrix()
          datX1 <- datt0[(datt0$index == i & datt0$Time != s), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[-s,]) %>% as.data.frame()
          fit_i <- lm(V1~.+0, data = dat_treat)
          lambda_i <- coef(fit_i)
          L_s[,i] <- coef(fit_i)
        }
        
        y0_treat_hat_s <- matrix(0, T_obs, N_treat)
        for( i in 1:N_treat){
          datXs <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_treat_hat_s_i <- datXs %*% beta1 + 
            matF[s,] %*% L_s[,i] + mu + A[i,] + Xi[s,]
          y0_treat_hat_s[,i] <- y0_treat_hat_s_i
        } 
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_treat_hat_s[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        ## estimate the y0_s_hat
        y0_dot_s_hat <- matrix(0, 1, N_treat)
        for( i in 1:N_treat){
          datX_s <- datt0[(datt0$index == i & datt0$Time == s),3+(1:N_x)] %>% as.matrix()
          y0_dot_s_hat_i <- datX_s %*% beta1 + mu + A[i,] + Xi[s,]
          y0_dot_s_hat[,i] <- y0_dot_s_hat_i
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_dot_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[R+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}

coreAlgo_0502 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x ){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data (only control units)
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
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-14
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        fitB1 <- lm(Outcome ~ .+0, data = datB1)
        datB2 <- cbind(fitB1$residuals, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
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
      
      ## estimate the factor loading of the treated unit
      
      #for(i in 1:N_treat){
      #  y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0) ,3])
      #  x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
      #  alpha_i <- y_imean - sum(x_imean*beta1) - mu
      #  A[i,] <- alpha_i
      #}
      
      A <- matrix(0,N,1)
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
        datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1 + matrix(mu, T0,1) + Xi[1:T0] 
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF[(1:T0),]) %>% as.data.frame()
        fit_i <- lm(V1~., data = dat_treat)
        lambda_i <- coef(fit_i)[-1]
        L[,i] <- lambda_i
        A[i] <- coef(fit_i)[1]
      }
      
      ### unit-specific effect
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
        y0_treat_hat_i <- datX_treat %*% beta1 + matF %*% matrix(L[,i], r,1) +
          matrix(mu, T_obs,1)+
          Xi+ 
          matrix(A[i,], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
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
      #### for treated units
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
        y0_treat_hat_i <- datX_treat %*% beta1  +
          matrix(mu, T_obs,1) +
          Xi +
          matrix(A[i,], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        
        ## estimate the factor loading of the treated unit w/o s
        A_s <- matrix(0, N_treat, 1)
        L_s <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3] %>% as.matrix()
          datX1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1 - matrix(mu, T0-1, 1) - Xi[1:T0][-s]
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[c(1:T0)[-s],]) %>% as.data.frame()
          fit_i <- lm(V1~., data = dat_treat)
          L_s[,i] <- coef(fit_i)[-1]
          A_s[i] <- coef(fit_i)[1]
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 + matF[s,] %*% L_s[,i] +
            mu +
            Xi[s] +
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        
        ## estimate alpha_i
        
        A_s <- matrix(0, N_treat, 1)
        for (i in 1:N_treat){
          y_imean0_s <- mean(dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0),3])
          x_imean0_s <- apply(dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)],2,mean) %>% as.matrix()
          A_s[i] <- y_imean0_s - t(x_imean0_s) %*% beta1 - mean(Xi[1:T0])-mu 
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 +
            mu +
            Xi[s]+ 
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[r_opt+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}
coreAlgo_LS_0502 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ## demeaning the data (only control units)
  datt <- dat1[dat1$Time <= T0 , c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T0, (N_x+1))) %>% as.data.frame()
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T0){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  dat <- datt0
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  beta1 <- beta_init
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    fitA <- lm(V1~.+0,data = as.data.frame(datA))
    
    matF_all <- fitA$fitted.values %>% as.matrix()
    matF_coef <- coef(fitA) %>% as.matrix()
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ### mu
  y_dd_bar <- mean(dat1[(dat1$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(dat1[(dat1$index %in% c((N_treat+1):N)),3 + (1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean) %>% as.matrix()
    xi_t <- y_tmean - t(x_tmean) %*% beta1 - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific fixed effect
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
    datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
    datY1_hat <- datX1 %*% beta1 + matrix(mu, T0, 1) + Xi[1:T0] + matF_all
    e1 <- datY1 - datY1_hat
    A[i] <- mean(e1) 
  }
  
  ### unit-specific effect
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - t(x_imean) %*% beta1 - mu
    A[i,] <- alpha_i
  }
  
  ### estimate y0
  ## estimate the y0
  A_co <- matrix(0, T_obs*N_ctrl,1)
  for (i in 1:N_ctrl) {
    A_co[((i-1)*T_obs +1):(i*T_obs)] <- A[i+1]
  }
  Xi_co <- matrix(Xi, T_obs*N_ctrl, 1)
  mu_co <- matrix(mu, T_obs*N_ctrl, 1)
  daty_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3] %>% as.matrix()
  datx_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3+(1:N_x)]  %>% as.matrix()
  date_co0 <- (daty_co - datx_co %*% beta1 - A_co - Xi_co - mu_co) %>% as.matrix()
  date_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    date_co[,i] <- date_co0[((i-1)*T_obs +1):(i*T_obs)]
  }
  datF_treat <- date_co %*% matF_coef
  datx_treat <- dat1[(dat1$index == 1), 3+(1:N_x)] %>% as.matrix()
  y0_treat_hat <- datx_treat %*% beta1 + datF_treat + A[1] + Xi + matrix(mu, T_obs, 1)  
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_LASSO_0502 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ## demeaning the data (only control units)
  datt <- dat1[dat1$Time <= T0 , c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T0, (N_x+1))) %>% as.data.frame()
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T0){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  dat <- datt0
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  beta1 <- beta_init
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    X_pre <- datA[,-1]; y_pre <- datA[,1]
    
    fit_lasso <- glmnet(X_pre, y_pre, alpha = 1, intercept=FALSE)
    #plot(fit_lasso)
    cv_lasso <- cv.glmnet(X_pre, y_pre, alpha = 1, intercept=FALSE)
    #plot(cv_lasso)
    bestlam <- cv_lasso$lambda.min
    
    matF_all <- predict(fit_lasso, s = bestlam, newx = X_pre)
    matF_coef <-  as.matrix(coef.glmnet(fit_lasso, s = bestlam)[-1])
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  
  ### mu
  y_dd_bar <- mean(dat1[(dat1$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(dat1[(dat1$index %in% c((N_treat+1):N)),3 + (1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean) %>% as.matrix()
    xi_t <- y_tmean - t(x_tmean) %*% beta1 - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific fixed effect
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
    datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
    datY1_hat <- datX1 %*% beta1 + matrix(mu, T0, 1) + Xi[1:T0] + matF_all
    e1 <- datY1 - datY1_hat
    A[i] <- mean(e1) 
  }
  
  ### unit-specific effect
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - t(x_imean) %*% beta1 - mu
    A[i,] <- alpha_i
  }
  
  ### estimate y0
  ## estimate the y0
  A_co <- matrix(0, T_obs*N_ctrl,1)
  for (i in 1:N_ctrl) {
    A_co[((i-1)*T_obs +1):(i*T_obs)] <- A[i+1]
  }
  Xi_co <- matrix(Xi, T_obs*N_ctrl, 1)
  mu_co <- matrix(mu, T_obs*N_ctrl, 1)
  daty_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3] %>% as.matrix()
  datx_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3+(1:N_x)]  %>% as.matrix()
  date_co0 <- (daty_co - datx_co %*% beta1 - A_co - Xi_co - mu_co) %>% as.matrix()
  date_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    date_co[,i] <- date_co0[((i-1)*T_obs +1):(i*T_obs)]
  }
  datF_treat <- date_co %*% matF_coef
  datx_treat <- dat1[(dat1$index == 1), 3+(1:N_x)] %>% as.matrix()
  y0_treat_hat <- datx_treat %*% beta1 + datF_treat + A[1] + Xi + matrix(mu, T_obs, 1)  
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_simplex_0502 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ## demeaning the data (only control units)
  datt <- dat1[dat1$Time <= T0 , c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T0, (N_x+1))) %>% as.data.frame()
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T0){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  dat <- datt0
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  beta1 <- beta_init
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    aveW <- rep(1/N_ctrl,N_ctrl)
    fitA <- optim(par = aveW[-1], 
                  fn = beeta, data = datA,
                  lower = rep(0,N_ctrl-1), upper = rep(1,N_ctrl-1),
                  method = "L-BFGS-B")
    
    matF_coef <-  as.matrix(c(1-sum(fitA$par),fitA$par))
    matF_all <- datA[,-1] %*% matF_coef
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ### mu
  y_dd_bar <- mean(dat1[(dat1$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(dat1[(dat1$index %in% c((N_treat+1):N)),3 + (1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean) %>% as.matrix()
    xi_t <- y_tmean - t(x_tmean) %*% beta1 - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific fixed effect
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
    datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
    datY1_hat <- datX1 %*% beta1 + matrix(mu, T0, 1) + Xi[1:T0] + matF_all
    e1 <- datY1 - datY1_hat
    A[i] <- mean(e1) 
  }
  
  ### unit-specific effect
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - t(x_imean) %*% beta1 - mu
    A[i,] <- alpha_i
  }
  
  ### estimate y0
  ## estimate the y0
  A_co <- matrix(0, T_obs*N_ctrl,1)
  for (i in 1:N_ctrl) {
    A_co[((i-1)*T_obs +1):(i*T_obs)] <- A[i+1]
  }
  Xi_co <- matrix(Xi, T_obs*N_ctrl, 1)
  mu_co <- matrix(mu, T_obs*N_ctrl, 1)
  daty_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3] %>% as.matrix()
  datx_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3+(1:N_x)]  %>% as.matrix()
  date_co0 <- (daty_co - datx_co %*% beta1 - A_co - Xi_co - mu_co) %>% as.matrix()
  date_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    date_co[,i] <- date_co0[((i-1)*T_obs +1):(i*T_obs)]
  }
  datF_treat <- date_co %*% matF_coef
  datx_treat <- dat1[(dat1$index == 1), 3+(1:N_x)] %>% as.matrix()
  y0_treat_hat <- datx_treat %*% beta1 + datF_treat + A[1] + Xi + matrix(mu, T_obs, 1)  
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_MIO_0502 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ## demeaning the data (only control units)
  datt <- dat1[dat1$Time <= T0 , c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T0, (N_x+1))) %>% as.data.frame()
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T0){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  dat <- datt0
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  beta1 <- beta_init
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    fitA <- MIP0_optim_gsynth(df = as.data.frame(datA) , X_pre = datA[,-1])
    
    matF_all <- fitA$Yhat %>% as.matrix()
    matF_coef <- fitA$beta
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }  
  ### mu
  y_dd_bar <- mean(dat1[(dat1$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(dat1[(dat1$index %in% c((N_treat+1):N)),3 + (1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean) %>% as.matrix()
    xi_t <- y_tmean - t(x_tmean) %*% beta1 - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific fixed effect
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
    datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
    datY1_hat <- datX1 %*% beta1 + matrix(mu, T0, 1) + Xi[1:T0] + matF_all
    e1 <- datY1 - datY1_hat
    A[i] <- mean(e1) 
  }
  
  ### unit-specific effect
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - t(x_imean) %*% beta1 - mu
    A[i,] <- alpha_i
  }
  
  ### estimate y0
  ## estimate the y0
  A_co <- matrix(0, T_obs*N_ctrl,1)
  for (i in 1:N_ctrl) {
    A_co[((i-1)*T_obs +1):(i*T_obs)] <- A[i+1]
  }
  Xi_co <- matrix(Xi, T_obs*N_ctrl, 1)
  mu_co <- matrix(mu, T_obs*N_ctrl, 1)
  daty_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3] %>% as.matrix()
  datx_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3+(1:N_x)]  %>% as.matrix()
  date_co0 <- (daty_co - datx_co %*% beta1 - A_co - Xi_co - mu_co) %>% as.matrix()
  date_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    date_co[,i] <- date_co0[((i-1)*T_obs +1):(i*T_obs)]
  }
  datF_treat <- date_co %*% matF_coef
  datx_treat <- dat1[(dat1$index == 1), 3+(1:N_x)] %>% as.matrix()
  y0_treat_hat <- datx_treat %*% beta1 + datF_treat + A[1] + Xi + matrix(mu, T_obs, 1)  
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}

coreAlgo_LS_0507 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ## demeaning the data (only control units)
  datt <- dat1[dat1$Time <= T0 , c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T0, (N_x+1))) %>% as.data.frame()
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T0){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  dat <- datt0
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  beta1 <- beta_init
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    aveW <- rep(1/N_ctrl,N_ctrl)
    fitB1 <- optim(par = aveW[-1], 
                   fn = beeta0, data = datA,
                   method = "L-BFGS-B")
    matF_coef <-  as.matrix(c(0-sum(fitB1$par),fitB1$par))
    matF_all <- datA[,-1] %*% matF_coef
    
    #fitA <- lm(V1~.+0,data = as.data.frame(datA))
    
    #matF_all <- fitA$fitted.values %>% as.matrix()
    #matF_coef <- coef(fitA) %>% as.matrix()
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ### mu
  y_dd_bar <- mean(dat1[(dat1$index %in% c((N_treat+1):N)),3])
  x_dd_bar <- apply(dat1[(dat1$index %in% c((N_treat+1):N)),3 + (1:N_x)],2,mean) %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
    x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean) %>% as.matrix()
    xi_t <- y_tmean - t(x_tmean) %*% beta1 - mu 
    Xi[t,] <- xi_t
  }
  
  ### unit-specific fixed effect
  A <- matrix(0,N,1)
  for (i in 1:N_treat){
    datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
    datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
    datY1_hat <- datX1 %*% beta1 + matrix(mu, T0, 1) + Xi[1:T0] + matF_all
    e1 <- datY1 - datY1_hat
    A[i] <- mean(e1) 
  }
  
  ### unit-specific effect
  for(i in (N_treat+1):N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - t(x_imean) %*% beta1 - mu
    A[i,] <- alpha_i
  }
  
  ### estimate y0
  ## estimate the y0
  A_co <- matrix(0, T_obs*N_ctrl,1)
  for (i in 1:N_ctrl) {
    A_co[((i-1)*T_obs +1):(i*T_obs)] <- A[i+1]
  }
  Xi_co <- matrix(Xi, T_obs*N_ctrl, 1)
  mu_co <- matrix(mu, T_obs*N_ctrl, 1)
  daty_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3] %>% as.matrix()
  datx_co <- dat1[(dat1$index %in% c((N_treat+1):N)),3+(1:N_x)]  %>% as.matrix()
  date_co0 <- (daty_co - datx_co %*% beta1 - A_co - Xi_co - mu_co) %>% as.matrix()
  date_co <- matrix(0, T_obs, N_ctrl)
  for (i in 1:N_ctrl){
    date_co[,i] <- date_co0[((i-1)*T_obs +1):(i*T_obs)]
  }
  datF_treat <- date_co %*% matF_coef
  datx_treat <- dat1[(dat1$index == 1), 3+(1:N_x)] %>% as.matrix()
  y0_treat_hat <- datx_treat %*% beta1 + datF_treat + A[1] + Xi + matrix(mu, T_obs, 1)  
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}


Aalpha0 <- function(par, ssr){
  A <- c(0-sum(par),par)
  sum((ssr-A)^2)
}
Xii0 <- function(par, ssr){
  xi <- c(0-sum(par),par)
  sum((ssr-xi)^2)
}
beeta0 <- function(par, data){
  mat <- as.matrix(data); beta <- as.matrix(c(0-sum(par),par))
  maty <- mat[,1];matx <- mat[,-1]
  sum((maty - matx %*% beta)^2)
}
coreAlgo_0506 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x ){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data (only control units)
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
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-5
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        #fitB1 <- lm(Outcome ~ .+0, data = datB1)
        aveW <- rep(1/r,r)
        fitB1 <- optim(par = aveW[-1], 
                       fn = beeta0, data = datB1,
                       method = "L-BFGS-B")
        matF_coef <-  as.matrix(c(0-sum(fitB1$par),fitB1$par))
        resB <- datY - datF %*% matrix(matF_coef, r,1)
        datB2 <- cbind(resB, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)[1:2]%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
      ### mu
      y_dd_bar <- mean(datt[(datt$index %in% c((N_treat+1):N)),3])
      x_dd_bar <- apply(datt[(datt$index %in% c((N_treat+1):N)),3+(1:N_x)],2,mean) %>% as.matrix()
      mu <- y_dd_bar - t(x_dd_bar) %*% beta1
      
      ### time-varying effect
      rss_Xi <- rep(0,T_obs)
      for(t in 1:T_obs){
        y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
        x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean)
        rss_xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
        rss_Xi[t] <- rss_xi_t
      }
      aveXi <- rep(1/T_obs,T_obs)
      fitXi <- optim(par = aveXi[-1], 
                     fn = Xii0, ssr = rss_Xi,
                     method = "L-BFGS-B")
      Xi <-  as.matrix(c(0-sum(fitXi$par),fitXi$par), T_obs, 1)
      
      
      
      ## estimate the factor loading of the treated unit
      
      #for(i in 1:N_treat){
      #  y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0) ,3])
      #  x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
      #  alpha_i <- y_imean - sum(x_imean*beta1) - mu
      #  A[i,] <- alpha_i
      #}
      
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
        datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1 + matrix(mu, T0,1) + Xi[1:T0] 
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF[(1:T0),]) %>% as.data.frame()
        aveW <- rep(1/r,r)
        fitA <- optim(par = aveW[-1], 
                      fn = beeta0, data = dat_treat,
                      method = "L-BFGS-B")
        matF_coef <-  as.matrix(c(0-sum(fitA$par),fitA$par))
        #fit_i <- lm(V1~.+0, data = dat_treat)
        L[,i] <- matF_coef
        #A[i,] <- coef(fit_i)[1]
      }
      
      
      ssr_A <- rep(0,N)
      for(i in 1:N_treat){
        y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0),3])
        x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu - mean(Xi, T0, 1)
        ssr_A[i] <- alpha_i
      }
      ### unit-specific effect
      for(i in (N_treat+1):N){
        y_imean <- mean(dat1[dat1$index == i,3])
        x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu
        ssr_A[i] <- alpha_i
      }
      
      aveA <- rep(1/N,N)
      fitAa <- optim(par = aveA[-1], 
                     fn = Aalpha0, ssr = ssr_A,
                     method = "L-BFGS-B")
      A <-  as.matrix(c(0-sum(fitAa$par),fitAa$par), N, 1)
      
      ## estimate the y0
      y0_treat_hat <- matrix(0, T_obs, N_treat)
      for (i in 1:N_treat){
        datX_treat <- dat1[(dat1$index == i), 3+(1:N_x)] %>% as.matrix()
        y0_treat_hat_i <- datX_treat %*% beta1 + matF %*% matrix(L[,i], r,1) +
          matrix(mu, T_obs,1)+
          Xi+ 
          matrix(A[i], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
      ### mu
      y_dd_bar <- mean(datt[(datt$index %in% c((N_treat+1):N)),3])
      x_dd_bar <- apply(datt[(datt$index %in% c((N_treat+1):N)),3+(1:N_x)],2,mean) %>% as.matrix()
      mu <- y_dd_bar - t(x_dd_bar) %*% beta1
      
      ### time-varying effect
      rss_Xi <- rep(0,T_obs)
      for(t in 1:T_obs){
        y_tmean <- mean(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)),3])
        x_tmean <- apply(dat1[(dat1$Time == t & dat1$index %in% c((N_treat+1):N)), 3+(1:N_x)], 2, mean)
        rss_xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
        rss_Xi[t] <- rss_xi_t
      }
      aveXi <- rep(1/T_obs,T_obs)
      fitXi <- optim(par = aveXi[-1], 
                     fn = Xii0, ssr = rss_Xi,
                     method = "L-BFGS-B")
      Xi <-  as.matrix(c(0-sum(fitXi$par),fitXi$par), T_obs, 1)
      
      ssr_A <- rep(0,N)
      for(i in 1:N_treat){
        y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0),3])
        x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu - mean(Xi, T0, 1)
        ssr_A[i] <- alpha_i
      }
      ### unit-specific effect
      for(i in (N_treat+1):N){
        y_imean <- mean(dat1[dat1$index == i,3])
        x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
        alpha_i <- y_imean - sum(x_imean*beta1) - mu
        ssr_A[i] <- alpha_i
      }
      
      aveA <- rep(1/N,N)
      fitAa <- optim(par = aveA[-1], 
                     fn = Aalpha0, ssr = ssr_A,
                     method = "L-BFGS-B")
      A <-  as.matrix(c(0-sum(fitAa$par),fitAa$par), N, 1)
      
      ## estimate the y0
      y0_treat_hat <- matrix(0, T_obs, N_treat)
      for (i in 1:N_treat){
        datX_treat <- dat1[(dat1$index == i), 3+(1:N_x)] %>% as.matrix()
        y0_treat_hat_i <- datX_treat %*% beta1  +
          matrix(mu, T_obs,1) +
          Xi +
          matrix(A[i,], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        
        L_s <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3] %>% as.matrix()
          datX1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1 + matrix(mu, T0-1,1) + Xi[1:T0][-s] 
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[c(1:T0)[-s],]) %>% as.data.frame()
          aveW <- rep(1/r,r)
          fitA <- optim(par = aveW[-1], 
                        fn = beeta0, data = dat_treat,
                        method = "L-BFGS-B")
          matF_coef <-  as.matrix(c(0-sum(fitA$par),fitA$par))
          #fit_i <- lm(V1~.+0, data = dat_treat)
          L_s[,i] <- matF_coef
          #A[i,] <- coef(fit_i)[1]
        }
        
        A_s <- matrix(0,N_treat,1)
        for(i in 1:N_treat){
          y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0),3])
          x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
          alpha_i <- y_imean - sum(x_imean*beta1) - mu - mean(Xi, T0, 1)
          A_s[i,] <- alpha_i
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 + matF[s,] %*% L_s[,i] +
            mu +
            Xi[s] +
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        
        ## estimate alpha_i
        
        A_s <- matrix(0, N_treat, 1)
        for (i in 1:N_treat){
          y_imean0_s <- mean(dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0),3])
          x_imean0_s <- apply(dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)],2,mean) %>% as.matrix()
          A_s[i] <- y_imean0_s - t(x_imean0_s) %*% beta1 - mean(Xi[1:T0])-mu 
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 +
            mu +
            Xi[s]+ 
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[r_opt+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}

coreAlgo_0507_2 <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x ){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data (only control units)
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
  
  dat <- datt0[(datt0$index %in% (N_treat +1):N),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## specify factor numbers
  container <- list()
  MSE_r <- c()
  for (r in 0:5){
    ## iteration part
    beta1 <- beta_init
    if (r != 0){
      
      loop <- 0;bdlist <- NULL
      
      betadist <- 1; tor = 1e-5
      while (betadist > tor){
        
        beta0 <- beta1
        
        ### (a) given beta, estimate F
        
        beta0 <- as.matrix(beta0)
        
        datY_hat <- datX %*% beta0
        dat_e <- datY - datY_hat
        
        rssm <- matrix(0,T_obs,T_obs)
        for ( i in 1:N_ctrl){
          e <- dat_e[(T_obs*(i-1)+1):(T_obs*i),]
          rssm_i <- e %*% t(e)
          rssm <- rssm + rssm_i
        }
        rssm <- rssm/(N_ctrl*T_obs)
        
        pca0 <- prcomp(rssm)
        matF <- pca0$x[,1:r] %>% as.matrix()
        
        ### (b) given F, estimate beta
        datF <- NULL
        c <- 0
        while (c < N_ctrl) {
          datF <- rbind(datF, matF)
          c <- c+1
        }
        #datB1 <- cbind(datY, datF) %>% as.data.frame();colnames(datB1)[1] <- "Outcome"
        #fitB1 <- lm(Outcome ~ .+0, data = datB1)
        #aveW <- rep(1/r,r)
        #fitB1 <- optim(par = aveW[-1], 
        #               fn = beeta0, data = datB1,
        #               method = "L-BFGS-B")
        #matF_coef <-  as.matrix(c(0-sum(fitB1$par),fitB1$par))
        
        L <-( t(datF) %*% (datY - datX %*% beta0) )/ T_obs
        
        resB <- datY - datF %*% matrix(L, r,1)
        datB2 <- cbind(resB, datX) %>% as.data.frame(); colnames(datB2)[1] <- "Outcome"
        fitB2 <- lm(Outcome ~ .+0, data = datB2)
        beta1 <- coef(fitB2)[1:2]%>% as.matrix()
        
        ### evaluate the convergence
        betadist <- dist(beta1 - beta0)
        bdlist <- c(bdlist, betadist)
        loop <- loop +1
      }
      
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
      
      ## estimate the factor loading of the treated unit
      
      #for(i in 1:N_treat){
      #  y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0) ,3])
      #  x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
      #  alpha_i <- y_imean - sum(x_imean*beta1) - mu
      #  A[i,] <- alpha_i
      #}
      
      A <- matrix(0,N,1)
      L <- matrix(0, r, N_treat)
      for (i in 1:N_treat){
        datY1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3] %>% as.matrix()
        datX1 <- dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
        datY1_hat <- datX1 %*% beta1 + matrix(mu, T0,1) + Xi[1:T0] 
        e1 <- datY1 - datY1_hat
        dat_treat <- cbind(e1, matF[(1:T0),]) %>% as.data.frame()
        fit_i <- lm(V1~., data = dat_treat)
        lambda_i <- coef(fit_i)[-1]
        L[,i] <- lambda_i
        A[i] <- coef(fit_i)[1]
      }
      
      ### unit-specific effect
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
        y0_treat_hat_i <- datX_treat %*% beta1 + matF %*% matrix(L[,i], r,1) +
          matrix(mu, T_obs,1)+
          Xi+ 
          matrix(A[i,], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
      
    }else{
      beta1 <- as.matrix(beta1)
      matF <- matrix(NA, T_obs, r)
      
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
      #### for treated units
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
        y0_treat_hat_i <- datX_treat %*% beta1  +
          matrix(mu, T_obs,1) +
          Xi +
          matrix(A[i,], T_obs, 1)
        y0_treat_hat[,i] <- y0_treat_hat_i 
      }
    }
    
    ## calculate leave-one-out MSE for CV
    mse_r <- 0
    if (r != 0){
      for (s in 1:T0){
        
        L_s <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3] %>% as.matrix()
          datX1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1 + matrix(mu, T0-1,1) + Xi[1:T0][-s] 
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[c(1:T0)[-s],]) %>% as.data.frame()
          aveW <- rep(1/r,r)
          fitA <- optim(par = aveW[-1], 
                        fn = beeta0, data = dat_treat,
                        method = "L-BFGS-B")
          matF_coef <-  as.matrix(c(0-sum(fitA$par),fitA$par))
          #fit_i <- lm(V1~.+0, data = dat_treat)
          L_s[,i] <- matF_coef
          #A[i,] <- coef(fit_i)[1]
        }
        
        A_s <- matrix(0,N_treat,1)
        for(i in 1:N_treat){
          y_imean <- mean(dat1[(dat1$index == i & dat1$Time <= T0),3])
          x_imean <- apply(dat1[(dat1$index == i & dat1$Time <= T0), 3+(1:N_x)], 2, mean)
          alpha_i <- y_imean - sum(x_imean*beta1) - mu - mean(Xi, T0, 1)
          A_s[i,] <- alpha_i
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 + matF[s,] %*% L_s[,i] +
            mu +
            Xi[s] +
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }else{
      for (s in 1:T0){
        
        ## estimate the factor loading of the treated unit w/o s
        A_s <- matrix(0, N_treat, 1)
        L_s <- matrix(0, r, N_treat)
        for (i in 1:N_treat){
          datY1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3] %>% as.matrix()
          datX1 <- dat1[(dat1$index == i & dat1$Time != s & dat1$Time <= T0), 3+(1:N_x)] %>% as.matrix()
          datY1_hat <- datX1 %*% beta1 - matrix(mu, T0-1, 1) - Xi[1:T0][-s]
          e1 <- datY1 - datY1_hat
          dat_treat <- cbind(e1, matF[c(1:T0)[-s],]) %>% as.data.frame()
          fit_i <- lm(V1~., data = dat_treat)
          L_s[,i] <- coef(fit_i)[-1]
          A_s[i] <- coef(fit_i)[1]
        }
        
        ## estimate y0_{-s}
        for (i in 1:N_treat){
          datxs <- dat1[(dat1$index == i & dat1$Time == s), 3 + (1:N_x)] %>% as.matrix()
          y0_s_hat <- datxs %*% beta1 + matF[s,] %*% L_s[,i] +
            mu +
            Xi[s] +
            A_s[i]
        }
        
        ## calculate MSE w/o s
        mse_s <- 0
        for (i in 1:N_treat){
          mse <- (dat1[(dat1$index == i & dat1$Time == s),3] - y0_s_hat[,i])^2
          mse_s <- mse_s + mse 
        }
        
        ## summing mse
        mse_r <- mse_r + mse_s
      }
    }
    
    list_r <- list(beta = beta1, matF = matF,
                   Y0_hat = y0_treat_hat,
                   mu = mu, alpha = A, xi = Xi)
    
    ## save the result for each r
    MSE_r <- c(MSE_r , mse_r)
    container[[r+1]] <- list_r
  }
  
  r_opt <- which.min(MSE_r)-1
  result <- container[[r_opt+1]]
  
  resultlist <- list(beta = result$beta,
                     matF = result$matF,
                     Y0_hat = result$Y0_hat,
                     mu = result$mu,
                     alpha = result$alpha,
                     xi = result$xi,
                     r_opt = r_opt,
                     container = container)
  return(resultlist)
}

coreAlgo_LS <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$Time <= T0 ),]#[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  #beta1 <- beta_init
  beta1 <- beta_ife
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    fitA <- lm(V1~.,data = as.data.frame(datA))
    
    matF_all <- fitA$fitted.values %>% as.matrix()
    matF_coef <- coef(fitA) %>% as.matrix()
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ## estimate the y0
  datY_all <- datt0[,3] %>% as.matrix()
  datX_all <- datt0[,(3+(1:N_x))] %>% as.matrix()
  dat_resi0 <- datY_all - datX_all %*% beta1
  dat_resi <- matrix(0,T_obs,N)
  for (i in 1:N){
    dat_resi[,i] <- dat_resi0[(T_obs*(i-1)+1):(T_obs*i)]
  }
  dat_resi <- cbind(matrix(1,T_obs,1),dat_resi[,-1])
  y0_dot_hat <- datX_all[(1:T_obs),] %*% beta1 + dat_resi %*% matF_coef
  
  
  y0_treat_hat <- matrix(0,T_obs, N_treat)
  y_itmean <- obs_itmean[1]
  for (i in 1:N_treat){
    y_imean <- mean(datt[(datt$index == i),3])
    for (t in 1:T_obs){
      y_tmean <- mean(datt[(datt$Time == t),3])
      y0_hat <- y0_dot_hat[t,i] + y_imean + y_tmean -y_itmean
      y0_treat_hat[t,i] <- y0_hat
    }
  }
  
  ## estimating the fixed effect
  
  ### mu
  y_dd_bar <- obs_itmean[1]
  x_dd_bar <- obs_itmean[-1] %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### unit-specific effect
  A <- matrix(0,N,1)
  for(i in 1:N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - sum(x_imean*beta1) - mu
    A[i,] <- alpha_i
  }
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[dat1$Time == t,3])
    x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
    xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
    Xi[t,] <- xi_t
  }
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_LASSO <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$Time <= T0 ),]#[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  #beta1 <- beta_init
  #beta1 <- c(0,0)
  beta1 <- beta_ife
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    X_pre <- datA[,-1]; y_pre <- datA[,1]
    
    fit_lasso <- glmnet(X_pre, y_pre, alpha = 1,intercept = FALSE )
    #plot(fit_lasso)
    cv_lasso <- cv.glmnet(X_pre, y_pre, alpha = 1)
    #plot(cv_lasso)
    bestlam <- cv_lasso$lambda.min
    
    matF_all <- predict(fit_lasso, s = bestlam, newx = X_pre)
    matF_coef <-  as.matrix(coef.glmnet(fit_lasso, s = bestlam)[-1])
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ## estimate the y0
  datY_all <- datt0[,3] %>% as.matrix()
  datX_all <- datt0[,(3+(1:N_x))] %>% as.matrix()
  dat_resi0 <- datY_all - datX_all %*% beta1
  dat_resi <- matrix(0,T_obs,N)
  for (i in 1:N){
    dat_resi[,i] <- dat_resi0[(T_obs*(i-1)+1):(T_obs*i)]
  }
  dat_resi <- dat_resi[,-1]
  y0_dot_hat <- datX_all[(1:T_obs),] %*% beta1 + dat_resi %*% matF_coef
  
  
  y0_treat_hat <- matrix(0,T_obs, N_treat)
  y_itmean <- obs_itmean[1]
  for (i in 1:N_treat){
    y_imean <- mean(datt[(datt$index == i),3])
    for (t in 1:T_obs){
      y_tmean <- mean(datt[(datt$Time == t),3])
      y0_hat <- y0_dot_hat[t,i] + y_imean + y_tmean -y_itmean
      y0_treat_hat[t,i] <- y0_hat
    }
  }
  
  ## estimating the fixed effect
  
  ### mu
  y_dd_bar <- obs_itmean[1]
  x_dd_bar <- obs_itmean[-1] %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### unit-specific effect
  A <- matrix(0,N,1)
  for(i in 1:N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - sum(x_imean*beta1) - mu
    A[i,] <- alpha_i
  }
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[dat1$Time == t,3])
    x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
    xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
    Xi[t,] <- xi_t
  }
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_simplex <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$Time <= T0 ),]#[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  #beta1 <- beta_init
  #beta1 <- c(0,0)
  beta1 <- beta_ife
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    aveW <- rep(1/N_ctrl,N_ctrl)
    fitA <- optim(par = aveW[-1], 
                  fn = beeta, data = datA,
                  lower = rep(0,N_ctrl-1), upper = rep(1,N_ctrl-1),
                  method = "L-BFGS-B")
    
    matF_coef <-  as.matrix(c(1-sum(fitA$par),fitA$par))
    matF_all <- datA[,-1] %*% matF_coef
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ## estimate the y0
  datY_all <- datt0[,3] %>% as.matrix()
  datX_all <- datt0[,(3+(1:N_x))] %>% as.matrix()
  dat_resi0 <- datY_all - datX_all %*% beta1
  dat_resi <- matrix(0,T_obs,N)
  for (i in 1:N){
    dat_resi[,i] <- dat_resi0[(T_obs*(i-1)+1):(T_obs*i)]
  }
  dat_resi <- dat_resi[,-1]
  y0_dot_hat <- datX_all[(1:T_obs),] %*% beta1 + dat_resi %*% matF_coef
  
  
  y0_treat_hat <- matrix(0,T_obs, N_treat)
  y_itmean <- obs_itmean[1]
  for (i in 1:N_treat){
    y_imean <- mean(datt[(datt$index == i),3])
    for (t in 1:T_obs){
      y_tmean <- mean(datt[(datt$Time == t),3])
      y0_hat <- y0_dot_hat[t,i] + y_imean + y_tmean -y_itmean
      y0_treat_hat[t,i] <- y0_hat
    }
  }
  
  ## estimating the fixed effect
  
  ### mu
  y_dd_bar <- obs_itmean[1]
  x_dd_bar <- obs_itmean[-1] %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### unit-specific effect
  A <- matrix(0,N,1)
  for(i in 1:N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - sum(x_imean*beta1) - mu
    A[i,] <- alpha_i
  }
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[dat1$Time == t,3])
    x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
    xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
    Xi[t,] <- xi_t
  }
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}
coreAlgo_MIO <- function(dat1, T0, T_obs, N_treat, N_ctrl, N_x){
  
  ## variable defining
  N <- N_treat + N_ctrl
  
  ##demeaning the data
  datt <- dat1[, c(1,2,3,3+(1:N_x))]
  datt0 <- cbind(datt[,c(1,2)], matrix(0, N*T_obs, 3)) %>% as.data.frame()
  
  obs_itmean <- apply(datt[,c(3,3+(1:N_x))], 2, mean)
  for (i in 1:N){
    obs_imean <- apply(datt[(datt$index == i),c(3,3+(1:N_x))], 2, mean)
    for (t in 1:T_obs){
      obs_tmean <- apply(datt[(datt$Time == t),c(3,3+(1:N_x))], 2, mean)
      obs_dot <- datt[(datt$index == i & datt$Time == t),c(3,3+(1:N_x))] - 
        obs_imean - obs_tmean + obs_itmean
      datt0[(datt0$index == i & datt0$Time == t),c(3,3+(1:N_x))] <- obs_dot
    }
  }
  
  dat <- datt0[(datt0$Time <= T0 ),]#[(datt0$index %in% (N_treat +1):N ),]
  
  ## parsing
  datY <- dat[,3] %>% as.matrix()
  datX <- dat[,3+(1:N_x)] %>% as.matrix()
  
  ## init--least square coefficient
  dat_init <- cbind(datY, datX) %>% as.data.frame(); colnames(dat_init)[1] <- "Outcome"
  fit_init <-lm(Outcome ~ .+0, data = dat_init) 
  beta_init <- coef(fit_init)
  
  ## iteration part
  #beta1 <- beta_init
  beta1 <- beta_ife
  
  loop <- 0;bdlist <- NULL
  betadist <- 1; tor = 1e-14
  while (betadist > tor){
    
    beta0 <- beta1
    
    ### (a) given beta, estimate F
    
    beta0 <- as.matrix(beta0)
    
    datY_hat <- datX %*% beta0
    dat_e <- datY - datY_hat
    
    datA <- matrix(0,T0,N)
    for (i in 1:N){
      datA[,i] <- dat_e[(T0*(i-1)+1):(T0*i)]
    }
    
    fitA <- MIP0_optim_gsynth(df = as.data.frame(datA) , X_pre = datA[,-1])
    
    matF_all <- fitA$Yhat %>% as.matrix()
    matF_coef <- fitA$beta
    
    ### (b) given F, estimate beta
    datF <- NULL
    c <- 0
    while (c < N) {
      datF <- rbind(datF, matF_all)
      c <- c+1
    }
    
    datB <- cbind(datY-datF, datX) %>% as.data.frame(); colnames(datB)[1] <- "Outcome"
    fitB <- lm(Outcome ~ .+0, data = datB)
    beta1 <- coef(fitB)%>% as.matrix()
    
    ### evaluate the convergence
    betadist <- dist(beta1 - beta0)
    bdlist <- c(bdlist, betadist)
    loop <- loop +1
  }
  
  ## estimate the y0
  datY_all <- datt0[,3] %>% as.matrix()
  datX_all <- datt0[,(3+(1:N_x))] %>% as.matrix()
  dat_resi0 <- datY_all - datX_all %*% beta1
  dat_resi <- matrix(0,T_obs,N)
  for (i in 1:N){
    dat_resi[,i] <- dat_resi0[(T_obs*(i-1)+1):(T_obs*i)]
  }
  dat_resi <- dat_resi[,-1]
  y0_dot_hat <- datX_all[(1:T_obs),] %*% beta1 + dat_resi %*% matF_coef
  
  
  y0_treat_hat <- matrix(0,T_obs, N_treat)
  y_itmean <- obs_itmean[1]
  for (i in 1:N_treat){
    y_imean <- mean(datt[(datt$index == i),3])
    for (t in 1:T_obs){
      y_tmean <- mean(datt[(datt$Time == t),3])
      y0_hat <- y0_dot_hat[t,i] + y_imean + y_tmean -y_itmean
      y0_treat_hat[t,i] <- y0_hat
    }
  }
  
  ## estimating the fixed effect
  
  ### mu
  y_dd_bar <- obs_itmean[1]
  x_dd_bar <- obs_itmean[-1] %>% as.matrix()
  mu <- y_dd_bar - t(x_dd_bar) %*% beta1
  
  ### unit-specific effect
  A <- matrix(0,N,1)
  for(i in 1:N){
    y_imean <- mean(dat1[dat1$index == i,3])
    x_imean <- apply(dat1[dat1$index == i, 3+(1:N_x)], 2, mean)
    alpha_i <- y_imean - sum(x_imean*beta1) - mu
    A[i,] <- alpha_i
  }
  
  ### time-varying effect
  Xi <- matrix(0,T_obs,1)
  for(t in 1:T_obs){
    y_tmean <- mean(dat1[dat1$Time == t,3])
    x_tmean <- apply(dat1[dat1$Time == t, 3+(1:N_x)], 2, mean)
    xi_t <- y_tmean - sum(x_tmean*beta1) - mu 
    Xi[t,] <- xi_t
  }
  
  resultlist <- list(beta = beta1, 
                     matF = list(matF_all = matF_all, matF_coef = matF_coef), 
                     Y0_hat = y0_treat_hat,
                     mu = mu, alpha = A, xi = Xi)
  
  return(resultlist)
}

MIP0_optim_gsynth <- function(df,X_pre){
  ##data
  y = as.matrix(df[1:T0,1])
  X = X_pre
  
  ## parameters
  cv_fold = 5
  p = ncol(X)
  
  ## splits data
  splits_indice = TSeriesCV_splits(y,X,cv_fold)
  
  ## CV loops
  mse_list = data.frame(matrix(NA,cv_fold,p))
  beta_list = data.frame(matrix(NA,cv_fold*p,2*p))
  #colname(mse_list) = c(rep(paste0("k = _",)))
  c = 1
  for (fold in 1:cv_fold){
    # splits data into training and testing 
    X_train = X[splits_indice$training_indice[[fold]],]
    X_test = X[splits_indice$validation_indice[[fold]],]
    y_train = y[splits_indice$training_indice[[fold]],]
    y_test = y[splits_indice$validation_indice[[fold]],] 
    
    # best subset select: k tuning 
    for (k in 1:p) {
      MIOresult = MIP(y_train, X_train,k)
      beta_list[c,] = MIOresult$x
      
      beta = Matrix(nrow = p, ncol = 1, data = 0, sparse = TRUE)
      beta[1:p] = as.matrix(MIOresult$x[1:p],p,1)
      
      y_predict = X_test %*% beta
      
      mse = sum((y_predict - y_test)^2)
      #pause
      mse_list[fold,k] = mse
      
      c = c+1
    }
  }
  
  mse_k = apply(mse_list, 2 , sum)
  minmsek0 = as.numeric(which(mse_k == min(mse_k))) 
  result = MIP(y, X, k = minmsek0)
  
  beta = as.matrix(result$x[1:p], p,1) 
  
  Yhat_MIO0 = X_pre %*% beta
  
  resultlist = list(Yhat = Yhat_MIO0, beta = beta)
  
  return(resultlist)
}
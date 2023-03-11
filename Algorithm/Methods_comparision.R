library(dplyr);library(ggplot2);library(caret);library(glmnet);library(knitr)#kable()
library(gurobi);library(Matrix);library(leaps);library(gsynth);library(Synth)
library(xtable);library(MCPanel)

## DGP
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/DGP.R",
       encoding = "utf-8")
#source("F:/Oct_CodingSummary/Simulation/DGP.R")
## Methods
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/Methods.R",
       encoding = "utf-8")
#source("F:/Oct_CodingSummary/Simulation/Methods.R")
## Mygsynth and AOA
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/mygsynth.R",
       encoding = "utf-8")
#source("F:/Oct_CodingSummary/Simulation/mygsynth.R")

## Simultaneously Optimization Approach
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/Simultaneous Optimization Approach.R",
       encoding = "utf-8")
#source("F:/Oct_CodingSummary/Simulation/Simultaneous Optimization Approach.R")
##
source("C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/MCoptim.R",
       encoding = "utf-8")
#source("F:/Oct_CodingSummary/Simulation/MCoptim.R")

T0 =30;
T_obs = T0+10;  N_ctrl = 10; N_x = 2

mseratiotable0 <- data.frame(matrix(0,21,3))
mseratiotable <- data.frame(matrix(0,21,3))
B = 1000
loop = 0
count = 1 
set.seed(42)
t0 <- Sys.time()



while(count <= B){
  
  ##Data Generating
  #datatoana <- IFEgen(T0 = 20, T_obs = 30, N_ctrl = 10)
  datatoana <- DGP2.pure(T_obs = T_obs, T0 = T0, N_ctrl =N_ctrl)
  
  
  N =  N_ctrl +1; N_x = 2
  ##Data Parsing
  #data_sim1 <- datatoana$adh.data
  #data_sim2 <- datatoana$su.data
  #df <- datatoana$y.data
  #X_it_1 <- datatoana$x1.data 
  #X_it_2 <- datatoana$x2.data 
  #Y_it_0 <- datatoana$y0.data
  
  data_sim1 <- datatoana$datADH
  data_sim2 <- datatoana$datXu
  df <- datatoana$datY; colnames(df) <- c(paste0("V",1:N))
  X_it_1 <- datatoana$datX1
  X_it_2 <- datatoana$datX2
  Y_it_0 <- datatoana$datY0
  
  
  T_obs <- nrow(df); T0 <- T_obs - 10; N_ctrl <- ncol(df)-1
  
  X_pre <- cbind(rep(1,T0), df[1:T0,-1]) %>% as.matrix()
  #colnames(X_pre)[1] <- "ones"
  y_pre <- as.matrix(df[1:T0, 1])
  
  X_aft <- cbind(rep(1, T_obs-T0 ), df[(T0+1):T_obs, -1]) %>% as.matrix()
  #colnames(X_aft)[1] <- "ones"
  
  dat <- cbind(df[,1], rep(1,T_obs), df[,-1]) %>% as.matrix()
  #colnames(dat)[1:2] <- c("V1", "ones")
  
  ##main algorithm
  
  ### MIO part
  err_messager <- 0
  tryCatch(
    {
      yhat_MIO <- MIP_optimProcess(df)
      yhat_MIO0 <- MIP0_optimProcess(df,X_pre)
    
      w <- MIO_SC(y = df[c(1:T0),1], X =df[c(1:T0),-1])
      yhat_MIOSCM <- df[,-1] %*% matrix(w$x, N_ctrl,1)
    },
    error = function(cnd){err_messager <<- -1; print("skip a loop")}
  )
  
  if (err_messager != -1){
  
    yhat_SCM <- ADH_optimProcess(data_sim1)#SCM
    yhat_SCM0 <- ADH0_optimProcess(data_sim1)#SCM_nocor
    
    yhat_pda <- PDA_optimProcess(df %>% as.data.frame())
    yhat_lasso <- LASSO_optimProcess(X_pre, y_pre, dat)
    yhat_NET <- NET_optimProcess(df)
    
    iferesult <- IFE_optimProcess(data_sim2)
    yhat_ife <- iferesult$yhat
    
    beta_ife <- iferesult$model$beta
    fit_ife <- iferesult$model
    semiresult <- semi_optimProcess(fit_ife,
                                    X_it_1 %>% as.data.frame(), 
                                    X_it_2 %>% as.data.frame())
    yhat_semi1 <- semiresult$yhat1
    yhat_semi3 <- semiresult$yhat3
    
    iferesult_nocor <- IFE_optimProcess_nocor(data_sim2)
    yhat_ife_nocor <- iferesult_nocor$yhat 
    
    ymean <- rep(mean(y_pre),T_obs)
    
    yhat_mygsynth <- coreAlgo_0507_2(data_sim2,T0 = T0, T_obs = T_obs, N_treat =1 , N_ctrl = N_ctrl, N_x = 2)$Y0_hat
    yhat_lsls <- coreAlgo_LS_0502(data_sim2,T0 = T0, T_obs = T_obs, N_treat =1 , N_ctrl = N_ctrl, N_x = 2)$Y0_hat
    yhat_lslasso <-  coreAlgo_LASSO_0502(data_sim2,T0 = T0, T_obs = T_obs, N_treat =1 , N_ctrl = N_ctrl, N_x = 2)$Y0_hat
    yhat_simplex <- coreAlgo_simplex_0502(data_sim2,T0 = T0, T_obs = T_obs, N_treat =1 , N_ctrl = N_ctrl, N_x = 2)$Y0_hat
    yhat_gsynthmio<- coreAlgo_MIO_0502(data_sim2,T0 = T0, T_obs = T_obs, N_treat =1 , N_ctrl = N_ctrl, N_x = 2)$Y0_hat
    
    yhat_hsiaols <- simuoptim(data_sim2, T_obs = T_obs, T0 = T0, N_ctrl = N_ctrl, N_x = N_x)$y0hat
    yhat_hsiaosimplex <- simuoptim.simplex(data_sim2, T_obs = T_obs, T0 = T0, N_ctrl = N_ctrl, N_x = N_x)$y0hat
    
    yhat_MC <- MCoptim(data_sim2, T_obs, T0)
    Y0 = Y_it_0[,1]
    
    df_result = cbind(Y0, ymean, 
                      yhat_SCM, yhat_MIO, yhat_MIOSCM,
                      yhat_pda, yhat_lasso, yhat_NET,yhat_SCM0, yhat_ife_nocor, yhat_MIO0, 
                      yhat_ife, yhat_semi1, yhat_semi3,
                      yhat_mygsynth, yhat_lsls, yhat_lslasso, yhat_simplex, yhat_gsynthmio,
                      yhat_hsiaols,yhat_hsiaosimplex,
                      yhat_MC
                      ) %>% as.data.frame()
    colnames(df_result) <- c("Y0","ymean",
                             "yhat_SCM","yhat_MIO", "yhat_MIOSCM",
                             "yhat_pda","yhat_lasso","yhat_NET","yhat_SCM0", "yhat_ife0","yhat_MIO0",
                             "yhat_ife", "yhat_semi1", "yhat_semi3",
                             "yhat_mygsynth", "yhat_lsls", "yhat_lslasso", "yhat_simplex", "yhat_gsynthmio",
                             "yhat_hsiaols","yhat_hsiaosimplex",
                             "yhat_MC")
    
    df_mse <- apply(df_result, 2, function(x) (df_result$Y0-x)^2)[,-1] %>% as.data.frame()
    colnames(df_mse) <- c("mseMean",
                          "mseSCM", "mseMIO","mseMIOSCM",
                          "msePDA","mseLASSO","mseNET","mseSCM0","mse_ife0","mseMIO0", 
                          "mseIFE", "mse_semi1", "mse_semi3",
                          "mse_mygsynth", "mse_lsls", "mse_lslasso", "mse_simplex", "mse_gsynthmio",
                          "mse_hsiaols","mse_hsiaosimplex",
                          "mse_MC")
    
    
    mse_pre <- apply(df_mse, 2, function(x) mean(x[1:T0]))
    mse_aft <- apply(df_mse, 2, function(x) mean(x[(T0+1):T_obs]))
    mse_ratio <- mse_aft/mse_pre
    
    df_mseratio = cbind(mse_pre, mse_aft, mse_ratio) #%>% as.data.frame()
    colnames(df_mseratio) = c("mse_pre", "mse_aft", "mse_ratio")
    row.names(df_mseratio) = c("Mean",
                               "SCM","BSS-simplex","MIO_SCM",
                               "PDA","LASSO","NET","SCM0","IFE0","MIO0",
                               "IFE","semi1","semi3",
                               "Mygsynth","gsynth_lsls","gsynth_lslasso","gsynth_simplex","gsynth_mio",
                               "Hsiaols", "hsiaosimplex",
                               "MC")
    
    mseratiotable0 <- mseratiotable0 + df_mseratio
    count <- count + 1
    mseratiotable <- mseratiotable0/count
  }

  loop = loop + 1
}

t1 <- Sys.time()
##data tidying
colnames(mseratiotable) = c("mse_pre", "mse_aft", "mse_ratio")
row.names(mseratiotable) = c("Mean",
                             "SCM","BSS-simplex","MIO_SCM",
                             "PDA","LASSO","NET","SCM0","IFE0","MIO0",
                             "IFE","semi1","semi3",
                             "Mygsynth","gsynth_lsls","gsynth_lslasso","gsynth_simplex","gsynth_mio",
                             "Hsiaols","Hsiaosimplex",
                             "MC")

t1-t0

View(mseratiotable)

#write.csv(mseratiotable, "C:/Users/user/Documents/Thesis/Oct_CodingSummary/Simulation/simuresult/msetable_dgppuret30n10.csv")
#write.csv(mseratiotable, "E:/Oct_CodingSummary/Simulation/simuresult/msetable_dgp5t20n10.csv")
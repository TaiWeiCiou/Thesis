library(dplyr);library(ggplot2);library(caret);library(glmnet);library(knitr)#kable()
library(gurobi);library(Matrix);library(leaps);library(gsynth);library(Synth)
library(xtable)


MIP_SC = function(y, X, k){
  # parameters setting
  p = ncol(X)
  
  #model building
  model <- list()
  model$modelnames <- "SC"
  
  #variables type
  model$varnames <- c(
    paste0(rep('b',p), 1:p), # beta as variables
    paste0(rep('z_bar',p), 1:p)  # 1-z as variables
  )
  
  model$vtype <-  c(rep("C", p), rep("B", p)) # beta are cont.v., z_bar are bi.v 
  model$lb <- c(rep(0,p), rep(0,p)) #all beta are non-negative.
  model$ub <- c(rep(Inf,p), rep(1,p))
  
  # objective
  
  ## construct Qc matrix
  Q <- Matrix(nrow = 2*p, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  Q[1:p, 1:p] = 0.5*t(X)%*%X
  Q[(p+1):(2*p), (p+1):(2*p)] = diag(p)
  
  ## setting objective
  model$Qc <- Q #Qc matrix
  model$obj <- -rbind(t(X)%*%y, matrix(rep(0,p),p,1)) #q vectors
  model$objcon <- 0.5*t(y)%*%y %>% as.numeric() # alpha
  
  model$modelsense <- 'min'
  
  # constraints
  ## SOS constraints
  #Make a list with elements like ["SOS1","SOS2",...,"SOSi",..."SOSp"]
  SOS_list = list()
  
  for (i in 1:p) {
    SOS <-  list()           #assign SOSi as a list
    SOS$type <- 1            #SOS-1 constraint
    SOS$index <- c(i,p+i)    # specify bi and z_bari
    SOS$weight <- c(1,2)     #Gurobi's setting: to order the variables.  
    SOS_list[[i]]<- SOS      #Put the constraint setting into a list
  }
  
  model$sos <- SOS_list
  
  ## L1-norm constraint and simplex constraints
  #A <- Matrix(nrow = 1, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  #A[(p+1):(2*p)] <-  c(rep(-1,p))
  #model$A <- A
  #model$rhs <- k-p
  #model$sense <- '<'
  
  A <- Matrix(nrow = 2, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  A[1, 1:p] <- rep(1,p)
  A[2, (p+1):(2*p)] <-  c(rep(-1,p))
  model$A <- A
  model$rhs        <- c(1, k-p)
  model$sense      <- c('=', '<')
  
  
  params = list()
  params$OutputFlag <- 0
  
  result = gurobi(model, params)
  return(result)
}
MIO_SC = function(y, X){
  # parameters setting
  p = ncol(X)
  
  #model building
  model <- list()
  model$modelnames <- "SC"
  
  #variables type
  model$varnames <- c(
    paste0(rep('b',p), 1:p)
  )
  
  model$vtype <-  c(rep("C", p)) # beta are cont.v., z_bar are bi.v 
  model$lb <- c(rep(0,p)) #all beta are non-negative.
  model$ub <- c(rep(1,p))
  
  # objective
  
  ## construct Qc matrix
  Q <- Matrix(nrow = p, ncol = p, data = 0, sparse = TRUE, doDiag = FALSE)
  Q = 0.5*t(X)%*%X
  
  ## setting objective
  model$Qc <- Q #Qc matrix
  model$obj <- -t(X)%*%y #q vectors
  model$objcon <- 0.5*t(y)%*%y %>% as.numeric() # alpha
  
  model$modelsense <- 'min'
  
  # constraints
  
  A <- Matrix(nrow = 1, ncol = p, data = 1, sparse = TRUE, doDiag = FALSE)
  model$A <- A
  model$rhs        <- c(1)
  model$sense      <- c('=')
  
  
  params = list()
  params$OutputFlag <- 0
  
  result = gurobi(model, params)
  return(result)
}

TSeriesCV_splits = function(y,X,nsplits){
  
  n = nrow(X)
  
  q = n %/% (nsplits+1)
  r = n %% (nsplits+1)
  partition_num = c(rep((q+1),r), rep(q,(nsplits+1-r)))
  training_indice = list()
  validation_indice = list()
  for (k in 1:nsplits) {
    training_indice[[k]] = c(1:sum(partition_num[1:k]))
    validation_indice[[k]] = c((sum(partition_num[1:k])+1):sum(partition_num[1:(k+1)]))
  }
  
  indice_list = list(
    "training_indice" = training_indice,
    "validation_indice" = validation_indice
  )
  
  return(indice_list)
}
MIP = function(y, X, k){
  # parameters setting
  p = ncol(X)
  
  #model building
  model <- list()
  model$modelnames <- "regression"
  
  #variables type
  model$varnames <- c(
    paste0(rep('b',p), 1:p), # beta as variables
    paste0(rep('z_bar',p), 1:p)  # 1-z as variables
  )
  
  model$vtype <-  c(rep("C", p), rep("B", p)) # beta are cont.v., z_bar are bi.v 
  model$lb <- c(rep(-Inf,p), rep(0,p)) #all beta are non-negative.
  model$ub <- c(rep(Inf,p), rep(1,p))
  
  # objective
  
  ## construct Qc matrix
  Q <- Matrix(nrow = 2*p, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  Q[1:p, 1:p] = 0.5*t(X)%*%X
  Q[(p+1):(2*p), (p+1):(2*p)] = diag(p)
  
  ## setting objective
  model$Qc <- Q #Qc matrix
  model$obj <- -rbind(t(X)%*%y, matrix(rep(0,p),p,1)) #q vectors
  model$objcon <- 0.5*t(y)%*%y %>% as.numeric() # alpha
  
  model$modelsense <- 'min'
  
  # constraints
  ## SOS constraints
  #Make a list with elements like ["SOS1","SOS2",...,"SOSi",..."SOSp"]
  SOS_list = list()
  
  for (i in 1:p) {
    SOS <-  list()           #assign SOSi as a list
    SOS$type <- 1            #SOS-1 constraint
    SOS$index <- c(i,p+i)    # specify bi and z_bari
    SOS$weight <- c(1,2)     #Gurobi's setting: to order the variables.  
    SOS_list[[i]]<- SOS      #Put the constraint setting into a list
  }
  
  model$sos <- SOS_list
  
  ## L1-norm constraint and simplex constraints
  #A <- Matrix(nrow = 1, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  #A[(p+1):(2*p)] <-  c(rep(-1,p))
  #model$A <- A
  #model$rhs <- k-p
  #model$sense <- '<'
  
  A <- Matrix(nrow = 1, ncol = 2*p, data = 0, sparse = TRUE, doDiag = FALSE)
  A[(p+1):(2*p)] <-  c(rep(-1,p))
  model$A <- A
  model$rhs        <- k-p
  model$sense      <- '<'
  
  
  params = list()
  
  params$OutputFlag <- 0
  
  result = gurobi(model, 
                  params
  )
  return(result)
}
ADH_optimProcess <- function(data_sim1){
  ## pick v by cross-validation
  # data setup for training model
  dataprep.out <-
    dataprep(
      foo = data_sim1 ,
      predictors    = c(4,5),
      dependent     = "Outcome",
      unit.variable = 2,
      time.variable = 1,
      treatment.identifier = "1",
      controls.identifier = unique(data_sim1$index)[-1],
      time.predictors.prior = 1:T0,
      time.optimize.ssr = 1:T0,
      unit.names.variable = 6,
      time.plot = 1:T_obs
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  #
  synth.tables <- synth.tab(
    dataprep.res = dataprep.out,
    synth.res = synth.out
  ) 
  #synth.tables
  
  yhat_SCM <- (dataprep.out$Y0%*%synth.out$solution.w)
  
  return(yhat_SCM)
}
ADH0_optimProcess <- function(data_sim1){
  dat_nocovar <- data_sim1
  dat_nocovar[,c(4,5)] <- matrix(1, T_obs, 2)
  
  dataprep.out <-
    dataprep(
      foo = dat_nocovar ,
      predictors    = c("Outcome"),
      dependent     = "Outcome",
      unit.variable = 2,
      time.variable = 1,
      treatment.identifier = "1",
      controls.identifier = unique(dat_nocovar$index)[-1],
      time.predictors.prior = 1:T0,
      time.optimize.ssr = 1:T0,
      unit.names.variable = 6,
      time.plot = 1:T_obs
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # fit main model with v from training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  
  #
  synth.tables0 <- synth.tab(
    dataprep.res = dataprep.out,
    synth.res = synth.out
  ) 
  #synth.tables
  
  yhat_SCM0 <- (dataprep.out$Y0%*%synth.out$solution.w)
  return(yhat_SCM0)
}
NET_optimProcess <- function(df){
  #####
  # Elastic net
  #unregister_dopar <- function() {
  #  env <- foreach:::.foreachGlobals
  #  rm(list=ls(name=env), pos=env)
  #}
  #unregister_dopar()
  
  cvtype <- trainControl(
    method = "timeslice",
    initialWindow = 5,
    horizon = 3,
    fixedWindow = FALSE,
    search = "random",
    allowParallel = FALSE)
  
  #enetgrid <- expand.grid(alpha = 0.1,lambda = seq(0,5,0.2))
  
  enet.fit <- train(
    V1 ~ ., data = df[1:T0,],
    method = "glmnet",
    trControl = cvtype
    #tuneGrid = enetgrid
  )
  
  fit_elnet_cv = cv.glmnet(X_pre[,-1], y_pre, alpha = enet.fit$bestTune[1], nfolds = 10)
  #fit_elnet_cv$lambda.min
  yhat_NET <- predict(fit_elnet_cv, dat[,-c(1,2)], s = "lambda.min")
  return(yhat_NET)
}
PDA_optimProcess <- function(df){
  # PDA method 
  fit2 <- regsubsets(V1 ~ ., data = df[1:T0,], method = "forward")
  
  summary2 <- summary(fit2)
  selectindex2 <- summary2$which[which.min(summary2$bic),]
  beta2 <- coef(fit2, which.min(summary2$bic))
  x2 <- apply(dat[,-1], 1, function(x) x * selectindex2) 
  x2 <- x2[apply(x2, 1, function(x) all(x != 0)),]
  yhat_pda1.2 <- beta2 %*% x2
  
  yhat_pda <- t(yhat_pda1.2)
  
  return(yhat_pda)
}
LASSO_optimProcess <- function(X_pre, y_pre, dat){
  fit_lasso <- glmnet(X_pre[,-1], y_pre, alpha = 1)
  #plot(fit_lasso)
  cv_lasso <- cv.glmnet(X_pre[,-1], y_pre, alpha = 1)
  #plot(cv_lasso)
  bestlam <- cv_lasso$lambda.min
  yhat_lasso <- predict(fit_lasso, s = bestlam, newx = dat[,-c(1,2)])
  
  return(yhat_lasso)
}
IFE_optimProcess <- function(data_sim2){
  fit_ife <- gsynth(Outcome ~ D + X1 + X2, data = data_sim2,
                    index=c("index","Time"), inference="parametric",
                    se = TRUE, nboots = 1, CV = TRUE,
                    force = "two-way", parallel = FALSE, cores = NULL)
  yhat_ife <- fit_ife$Y.ct
  
  result <- list(model = fit_ife, yhat = yhat_ife)
  
  return(result)
}
IFE_optimProcess_nocor <- function(data_sim2){
  fit_ife <- gsynth(Outcome ~ D, data = data_sim2,
                    index=c("index","Time"), inference="parametric",
                    se = TRUE, nboots = 1, CV = TRUE,
                    force = "two-way", parallel = FALSE, cores = 4)
  yhat_ife <- fit_ife$Y.ct
  
  result <- list(model = fit_ife, yhat = yhat_ife)
  
  return(result)
}
semi_optimProcess <- function(fit_ife, X_it_1, X_it_2){
  b_ife <- fit_ife$beta
  e_ife <- df - (X_it_1*b_ife[1] + X_it_2*b_ife[2])
  fit_semi1 <- lm(V1 ~ ., data = e_ife[(1:T0),])
  yhat_semi1 <- X_it_1[,1]*b_ife[1] + X_it_2[,1]*b_ife[2] + predict(fit_semi1, newdata = e_ife[,-1])
  
  fit_semi3 <- glmnet(e_ife[(1:T0),-1], e_ife[(1:T0),1], alpha = 1 )
  cv.semi <- cv.glmnet( as.matrix(e_ife[(1:T0),-1]), as.matrix(e_ife[(1:T0),1]), alpha = 1 )
  bestlam_semi <- cv.semi$lambda.min
  yhat_semi3 <- X_it_1[,1]*b_ife[1] + X_it_2[,1]*b_ife[2] + predict(fit_semi3,s = bestlam_semi, newx = as.matrix(e_ife[,-1]))
  
  ## optim approach
  #rss_semi <- function(beta, data, N_treat, N_ctrl, T_obs){
  
  #  df = as.data.frame(matrix(NA,T_obs,(N_treat + N_ctrl)))
  #  for (i in 1:(N_treat + N_ctrl)) {
  #    for (s in 1:T_obs){
  #      df[s,i] = data[
  #        which(data$Time == s & data$Unit == i),3]
  #    }
  #  }
  
  #  e_ife <- df - (data$X1*beta[1] + data$X2*beta[2])
  
  #  fit_semi3 <- glmnet(e_ife[(1:T0),-1], e_ife[(1:T0),1], alpha = 1 )
  #  cv.semi <- cv.glmnet( as.matrix(e_ife[(1:T0),-1]), as.matrix(e_ife[(1:T0),1]), alpha = 1 )
  #  bestlam_semi <- cv.semi$lambda.min
  #  yhat_semi3 <- X_it_1[,1]*beta[1] + X_it_2[,1]*beta[2] + predict(fit_semi3,s = bestlam_semi, newx = as.matrix(e_ife[,-1]))
  #  rss = sum((df$V1 - yhat_semi3)^2)
  #  return(rss)
  #}
  
  #initalpha <- c(b_ife[1], b_ife[2])
  #result <- optim(initalpha, rss_semi, 
  #                data = data_sim1, N_treat = 1, N_ctrl = N_ctrl, T_obs = T_obs )
  
  #b_opt <- result$par
  #e_opt <- df - (X_it_1*b_opt[1] + X_it_2*b_opt[2])
  #fit_semi4 <- glmnet(e_opt[(1:T0),-1], e_opt[(1:T0),1], alpha = 1 )
  #cv.semi4 <- cv.glmnet( as.matrix(e_opt[(1:T0),-1]), as.matrix(e_opt[(1:T0),1]), alpha = 1 )
  #bestlam_semi4 <- cv.semi4$lambda.min
  #yhat_semi4 <- X_it_1[,1]*b_opt[1] + X_it_2[,1]*b_opt[2] + predict(fit_semi4,s = bestlam_semi4, newx = as.matrix(e_opt[,-1]))
  
  result <- list(yhat1 = yhat_semi1, yhat3 = yhat_semi3 
                 #yhat4 = yhat_semi4
  )
}
MIP_optimProcess <- function(df){
  ##data
  y = as.matrix(df[1:T0,1])
  X = as.matrix(df[1:T0,-1])
  
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
      MIOresult = MIP_SC(y_train, X_train,k)
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
  minmsek = as.numeric(which(mse_k == min(mse_k))) 
  result = MIP_SC(y, X, k = minmsek)
  
  beta = as.matrix(result$x[1:p], p,1) 
  
  Xa = as.matrix(df[,-1])
  
  Yhat_MIO = Xa %*% beta
  
  return(Yhat_MIO)
}
MIP0_optimProcess <- function(df,X_pre){
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
  
  Xa = as.matrix(dat[,-1])
  
  Yhat_MIO0 = Xa %*% beta
  
  return(Yhat_MIO0)
}
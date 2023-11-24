library('BayesTree')
library('randomForest')
library('gbm')

library(parallel)
library(foreach)
library(doParallel)

# Number of cross-validation folds
NumOfRep <- 20
num_folds <- 5

createFolds <- function(data, k) {
  set.seed(123)  # For reproducibility
  fold_indices <- sample(rep(1:k, length.out = nrow(data)))
  fold_list <- split(1:nrow(data), fold_indices)
  return(fold_list)
}

# Create a data frame with all combinations of hyperparameter values
hyperparametersBART <- expand.grid(
  m=c(50,200),
  nu_q=list(c(3,0.9),c(3,0.99),c(10,0.75)),
  k=c(1,2,3,5)
)

hyperparametersAddiVortes <- expand.grid(
  m=c(50,200),
  nu=c(6),
  q=c(0.85),
  k=c(1,3),
  sd=c(0.8,1.5),
  omega=c(3),
  lambda=c(5,25)
)

hyperparametersBoost <- expand.grid(
  depth=c(1,2,3,4),
  num_trees=c(50,100,200),
  eta=c(0.01, 0.05, 0.10, 0.25)
)



num_cores <- 10 # Specify the number of cores you want to use
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for (l in 1:3){
  set.seed(324)
  
  X<-BenchmarkX[[l]]
  Y<-BenchmarkY[[l]]
  Y<-as.numeric(as.matrix(Y))
  
  n=length(Y)
  TrainSet<-matrix(nrow=trunc(5*n/6),ncol = NumOfRep)
  TestSet<-matrix(nrow=n-trunc(5*n/6),ncol = NumOfRep)
  
  for (h in 1:NumOfRep){
    TrainSet[,h]=sort(sample.int(n,5*n/6));
    TestSetInit<-1:n
    TestSet[,h]=TestSetInit[! TestSetInit %in% TrainSet[,h]];
    
  }  
  print('finished h')
  #Cross validation
  folds <- createFolds(X, k = num_folds)  # 5-fold cross-validation
  
 
  # Initialize a matrix to store the results
  results_matrixAddiVortes <- matrix(ncol = 8, nrow = nrow(hyperparametersAddiVortes))
  rmse_values<-vector(length = num_folds)
  
  # Iterate over hyperparameter combinations
  results_matrixAddiVortes<-foreach(i = 1:nrow(hyperparametersAddiVortes), .combine = rbind) %dopar% {
    library('invgamma')
    library('FNN')
    library('expint')
    
    params <- hyperparametersAddiVortes[i, ]
    rmse_values <- numeric(num_folds)
    
    # Perform cross-validation using a loop
    for (j in 1:num_folds) {
      TrainSetCV <- unlist(folds[-j],use.names = FALSE)
      TestSetCV <- folds[[j]]
      
      m<-params$m
      nu <- params$nu
      q<- params$q
      sd<- params$sd
      omega<- params$omega
      lambda<- params$lambda
      k<- params$k
      
      # Make predictions on the test set
      rmse_values[j] <- AddiVortes_Algorithm(Y[TrainSetCV],as.matrix(X[TrainSetCV,]),m,1200,200,nu,q,k,sd,omega,lambda,Y[TestSetCV],as.matrix(X[TestSetCV,]))$RMSE
    }                   
    # Store the average RMSE from cross-validation
    return(c(m,nu,q,k,sd,omega,lambda,mean(rmse_values)))
    
  }
  
  # Find the best hyperparameter combination based on RMSE
  best_index <- which.min(results_matrixAddiVortes[, 8])
  best_hyperparametersAddiVortes <- results_matrixAddiVortes[best_index, 1:7]
  
  print("Finished AddiVortes Tuning")
  
  # Initialize a matrix to store the results
  results_matrixBART <- matrix(ncol = 5, nrow = nrow(hyperparametersBART))
  
  # Iterate over hyperparameter combinations
  results_matrixBART<-foreach (i = 1:nrow(hyperparametersBART), .combine = rbind) %dopar% {
    library('BayesTree')
    
    params <- hyperparametersBART[i, ]
    rmse_values <- numeric(num_folds)
    
    # Perform cross-validation using a loop
    for (j in 1:num_folds) {
      TrainSetCV <- unlist(folds[-j],use.names = FALSE)
      TestSetCV <- folds[[j]]
      
      m <- params$m
      nu <- params$nu_q[[1]][1]
      q <- params$nu_q[[1]][2]
      k <- params$k
      
      # Train the bart model
      BartOr<-bart(as.matrix(X[TrainSetCV,]),as.numeric(as.matrix(Y[TrainSetCV])),as.matrix(X[TestSetCV,]), ntree = m, sigdf = nu,sigquant = q,k = k)
      
      # Make predictions on the test set
      rmse_values[j] <- (mean((Y[TestSetCV]-BartOr$yhat.test.mean)**2))**0.5
    }
    
    # Store the average RMSE from cross-validation
    results_matrixBART[i,] <- cbind(params$m,params$nu_q[[1]][1],params$nu_q[[1]][2],params$k,mean(rmse_values))

  }
  
  print('Finished BART Tuning')
  
  # Find the best hyperparameter combination based on RMSE
  best_index <- which.min(results_matrixBART[, 5])
  best_hyperparametersBART <- results_matrixBART[best_index, 1:4]
  
  # Initialize a matrix to store the results for each combination
  results_matrixRF <- matrix(ncol = 2, nrow = 4)
  PercentOfVariable = c(0.1,0.25,0.5,1)
  for (i in 1:4) {
    # Perform cross-validation using a loop
    fold_RMSE <- numeric(num_folds)
    for (j in 1:num_folds) {
      TrainSetCV <- unlist(folds[-j],use.names = FALSE)
      TestSetCV <- folds[[j]]
      
      mean.y.hat<-as.vector(randomForest(X[TrainSetCV,],Y[TrainSetCV],X[TestSetCV,],Y[TestSetCV],mtry=round(PercentOfVariable[i]*length(X[1,])))$test$predicted)
      fold_RMSE[j] <- sqrt(mean((Y[TestSetCV]-mean.y.hat)^2))
      
    }
    # Store the average RMSE from cross-validation
    results_matrixRF[i, 1] <- mean(fold_RMSE)
    results_matrixRF[i,2] <- PercentOfVariable[i]
  }
  ordered_indices <- order(results_matrixRF[, 1])
  ordered_results <- results_matrixRF[ordered_indices, ]
  PercentOfVariable<-ordered_results[1,2]
  
  print('Finished RF tuning')
  
  
  CV.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('invgamma')
    library('FNN')
  
    AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),best_hyperparametersAddiVortes[1],1200,200,best_hyperparametersAddiVortes[2],best_hyperparametersAddiVortes[3],best_hyperparametersAddiVortes[4],best_hyperparametersAddiVortes[5],best_hyperparametersAddiVortes[6],best_hyperparametersAddiVortes[7],Y[TestSet[,k]],as.matrix(X[TestSet[,k],]))$RMSE
  }
  
  Default.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('invgamma')
    library('FNN')
    
    AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),200,1200,200,6,0.85,3,0.8,3,25,Y[TestSet[,k]],as.matrix(X[TestSet[,k],]))$RMSE
  }
  
  BART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('BayesTree')
    library('expint')
    
    #bartMAchine
    modelBART<-bart(as.matrix(X[TrainSet[,k],]),as.numeric(as.matrix(Y[TrainSet[,k]])),as.matrix(X[TestSet[,k],]), ntree = best_hyperparametersBART[1], sigdf = best_hyperparametersBART[2],sigquant = best_hyperparametersBART[3],k = best_hyperparametersBART[4])
    BART.RMSE<-(mean((Y[TestSet[,k]]-modelBART$yhat.test.mean)**2))**0.5  
  }
  
  RForest.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('randomForest')
    #randomForests
    
    RandomF<-randomForest(X[TrainSet[,k],],Y[TrainSet[,k]],X[TestSet[,k],],Y[TestSet[,k]],mtry = round(PercentOfVariable*length(X[1,])))
    RForest.RMSE<-(mean((Y[TestSet[,k]]-RandomF$test$predicted)**2))**0.5
  }
    
    
  RMSE<-rbind(CV.AddiVortes.RMSE,Default.AddiVortes.RMSE,BART.RMSE,RForest.RMSE)
  
  print(l)
  
  RRMSE<-apply(RMSE, 2, function(col) col / min(col))
  
  boxplot(RRMSE[1,],RRMSE[2,],RRMSE[3,],RRMSE[4,],horizontal = TRUE, names = unique(c("AddiVortes CV","AddiVortes Default", "BART", "Random Forests")))
  
  if (l==1){
    All.RRMSE<-matrix(nrow=4,ncol=NumOfRep)
    All.RRMSE<-RRMSE
  }
  else{
    All.RRMSE<-cbind(All.RRMSE,RRMSE)
  }
}

stopCluster(cl)

boxplot(All.RRMSE[1,],All.RRMSE[2,],All.RRMSE[3,], All.RRMSE[4,],horizontal = TRUE, names = unique(c("AddiVortes-CV", "AddiVortes-default","BART","Random Forests")))

Smaller.RRMSE.AddiVortesCV<-All.RRMSE[1,All.RRMSE[1,]<1.5]
Smaller.RRMSE.AddiVortesDef<-All.RRMSE[2,All.RRMSE[2,]<1.5]
Smaller.RRMSE.BART<-All.RRMSE[3,All.RRMSE[3,]<1.5]
Smaller.RRMSE.RF<-All.RRMSE[4,All.RRMSE[4,]<1.5]

boxplot(as.numeric(Smaller.RRMSE.AddiVortesCV),as.numeric(Smaller.RRMSE.AddiVortesDef),as.numeric(Smaller.RRMSE.BART), as.numeric(Smaller.RRMSE.RF),horizontal = TRUE, names = unique(c("AddiVortes-CV", "AddiVortes-default","BART","Random Forests")))

length(Smaller.RRMSE.AddiVortesCV)/length(All.RRMSE[1,])
length(Smaller.RRMSE.AddiVortesDef)/length(All.RRMSE[2,])
length(Smaller.RRMSE.BART)/length(All.RRMSE[3,])
length(Smaller.RRMSE.RF)/length(All.RRMSE[4,])


OrderedAddiVortesCV<-sort(as.numeric(Smaller.RRMSE.AddiVortesCV))
OrderedAddiVortesDef<-sort(as.numeric(Smaller.RRMSE.AddiVortesDef))
OrderedRF<-sort(as.numeric(Smaller.RRMSE.RF))
OrderedBART<-sort(as.numeric(Smaller.RRMSE.BART))

OrderedAddiVortesCV[length(Smaller.RRMSE.AddiVortesCV)*0.5]
OrderedAddiVortesCV[length(Smaller.RRMSE.AddiVortesCV)*0.75]

OrderedAddiVortesDef[length(Smaller.RRMSE.AddiVortesDef)*0.5]
OrderedAddiVortesDef[length(Smaller.RRMSE.AddiVortesDef)*0.75]

OrderedBART[length(Smaller.RRMSE.BART)*0.5]
OrderedBART[length(Smaller.RRMSE.BART)*0.75]

OrderedRF[length(Smaller.RRMSE.RF)*0.5]
OrderedRF[length(Smaller.RRMSE.RF)*0.75]

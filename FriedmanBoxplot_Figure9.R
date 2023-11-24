library(parallel)
library(foreach)
library(doParallel)

#100 Friedman datasets
BenchmarkX<-list()
BenchmarkY<-list()
for (i in 1:100){
  
  set.seed(8)
  BenchmarkX[[i]]=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(BenchmarkX[[i]])
  BenchmarkY[[i]]=Ey+sigma*rnorm(n)
}

# Number of cross-validation folds
NumOfRep <- 1
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


hyperparametersBoost <- expand.grid(
  depth=c(1,2,3,4),
  num_trees=c(50,100,200),
  eta=c(0.01, 0.05, 0.10, 0.25)
)

num_cores <- 10  # Specify the number of cores you want to use
cl3 <- makeCluster(num_cores)
registerDoParallel(cl3)


for (l in 1:100){
  set.seed(324)
  
  X<-BenchmarkX[[l]]
  Y<-BenchmarkY[[l]]
  
  Y<-as.numeric(as.matrix(Y))
  
  n=length(Y)
  TrainSet<-matrix(nrow=trunc(n/11),ncol = NumOfRep)
  TestSet<-matrix(nrow=n-trunc(n/11),ncol = NumOfRep)
  
  for (h in 1:NumOfRep){
    TrainSet[,h]=sort(sample.int(n,n/11));
    TestSetInit<-1:n
    TestSet[,h]=TestSetInit[! TestSetInit %in% TrainSet[,h]];
    
  }  
  print('finished h')
  #Cross validation
  folds <- createFolds(X[1:125,], k = num_folds)  # 5-fold cross-validation
  
  
  # Initialize a matrix to store the results
  results_matrixBoost <- matrix(ncol = 4, nrow = nrow(hyperparametersBoost))
  
  # Iterate over hyperparameter combinations
  results_matrixBoost<-foreach(i = 1:nrow(hyperparametersBoost), .combine = rbind) %dopar% {
    library('gbm')
    
    params <- hyperparametersBoost[i, ]
    rmse_values <- numeric(num_folds)
    
    # Perform cross-validation using a loop
    for (j in 1:num_folds) {
      TrainSetCV <- unlist(folds[-j],use.names = FALSE)
      TestSetCV <- folds[[j]]
      
      data_frame_Train<-as.data.frame(cbind(Y[TrainSetCV],X[TrainSetCV,]))
      data_frame_Test<-as.data.frame(cbind(f(X[TestSetCV,]),X[TestSetCV,]))
      
      depth <- params$depth
      num_trees <- params$num_trees
      eta <- params$eta
      
      # Train the bart model
      gradBoost<-gbm(V1~.,data=data_frame_Train,distribution = "gaussian", interaction.depth=depth,n.trees = num_trees, shrinkage = eta)
      
      # Make predictions on the test set
      rmse_values[j] <- (mean((f(X[TestSetCV,])-predict(gradBoost,newdata=data_frame_Test))**2))**0.5
    }
    
    # Store the average RMSE from cross-validation
    results_matrixBoost[i,]<-c(depth,num_trees,eta,mean(rmse_values))
    
  }
  
  # Find the best hyperparameter combination based on RMSE
  best_index <- which.min(results_matrixBoost[, 4])
  best_hyperparametersBoost <- results_matrixBoost[best_index, 1:3]
  
  print('FinishedBoostTuning')
  
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
      BartOr<-bart(as.matrix(X[TrainSetCV,]),as.numeric(as.matrix(Y[TrainSetCV])),as.matrix(X[TestSetCV,]),nskip = 1000,ndpost = 3000, ntree = m, sigdf = nu,sigquant = q,k = k)
      
      # Make predictions on the test set
      rmse_values[j] <- (mean((f(X[TestSetCV,])-BartOr$yhat.test.mean)**2))**0.5
    }
    
    # Store the average RMSE from cross-validation
    results_matrixBART[i,] <- cbind(params$m,params$nu_q[[1]][1],params$nu_q[[1]][2],params$k,mean(rmse_values))
    
  }
  
  print('FinishedBARTTuning')
  
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
      
      mean.y.hat<-as.vector(randomForest(X[TrainSetCV,],Y[TrainSetCV],X[TestSetCV,],f(X[TestSetCV,]),mtry=round(PercentOfVariable[i]*length(X[1,])))$test$predicted)
      fold_RMSE[j] <- sqrt(mean((f(X[TestSetCV,])-mean.y.hat)^2))
      
    }
    # Store the average RMSE from cross-validation
    results_matrixRF[i, 1] <- mean(fold_RMSE)
    results_matrixRF[i,2] <- PercentOfVariable[i]
  }
  ordered_indices <- order(results_matrixRF[, 1])
  ordered_results <- results_matrixRF[ordered_indices, ]
  PercentOfVariable<-ordered_results[1,2]
  print('part 2')
  
  Default.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('invgamma')
    library('FNN')
    
    AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),50,1200,200,6,0.85,3,0.8,3,25,f(X[TestSet[,k],]),as.matrix(X[TestSet[,k],]))$RMSE
  }
  
  BART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('BayesTree')
    library('expint')
    
    #bartMAchine
    modelBART<-bart(as.matrix(X[TrainSet[,k],]),as.numeric(as.matrix(Y[TrainSet[,k]])),as.matrix(X[TestSet[,k],]),nskip = 1000,ndpost = 3000, ntree = best_hyperparametersBART[1], sigdf = best_hyperparametersBART[2],sigquant = best_hyperparametersBART[3],k = best_hyperparametersBART[4])
    BART.RMSE<-(mean((f(X[TestSet[,k],])-modelBART$yhat.test.mean)**2))**0.5  
  }
  
  RForest.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('randomForest')
    #randomForests
    
    RandomF<-randomForest(X[TrainSet[,k],],Y[TrainSet[,k]],X[TestSet[,k],],f(X[TestSet[,k],]),mtry = round(PercentOfVariable*length(X[1,])))
    RForest.RMSE<-(mean((f(X[TestSet[,k],])-RandomF$test$predicted)**2))**0.5
  }
  
  print('hi')
  # Perform cross-validation with the bart function
  Boost.RMSE<-vector(length=NumOfRep)
  Boost.RMSE<-foreach(k = 1:NumOfRep,.combine = cbind) %dopar% {
    library('gbm')
    #Gradient Boosting
    
    data_frame_Train<-as.data.frame(cbind(Y[TrainSet],X[TrainSet,]))
    data_frame_Test<-as.data.frame(cbind(f(X[TestSet,]),X[TestSet,]))
    
    gradientBoost<-gbm(V1~.,data=data_frame_Train,distribution="gaussian",interaction.depth = best_hyperparametersBoost[1],n.trees = best_hyperparametersBoost[2],shrinkage = best_hyperparametersBoost[3])
    Boost.RMSE[k]<-(mean((f(X[TestSet,])-predict(gradientBoost,newdata=data_frame_Test))**2))**0.5
  }
  
  
  RMSE<-rbind(Default.AddiVortes.RMSE,BART.RMSE,RForest.RMSE,Boost.RMSE)
  
  
  print('finished')
  print(l)
  
  boxplot(RMSE[1,],RMSE[2,],RMSE[3,],RMSE[4,],horizontal = TRUE, names = unique(c("AddiVortes Default", "BART", "Random Forests","Gradient Boosting")))
  
  if (l==1){
    All.RMSE<-matrix(nrow=4,ncol=NumOfRep)
    All.RMSE<-RMSE
  }
  else{
    All.RMSE<-cbind(All.RMSE,RMSE)
  }

}

stopCluster(cl3)

All.RMSE<-read.csv(file= "/home/grads/gmwc55/Documents/RRMSEdata/RMSEVAluesP10" )
All.RMSE<-as.matrix(All.RMSE[,2:101],)

for(i in 1:5){
  All.RMSE[i,]<-as.vector(All.RMSE[i,])
}


boxplot(All.RMSE[1,],All.RMSE[2,],All.RMSE[3,], All.RMSE[4,],All.RMSE[5,],horizontal = TRUE, names = unique(c("AddiVortes-CV", "AddiVortes-default","BART","Random Forests", "Boosting")))

boxplot(All.RMSE[1,],All.RMSE[2,],All.RMSE[3,], All.RMSE[4,],horizontal = TRUE, names = unique(c("AddiVortes-Def","BART","Random Forests","Boosting")))

length(All.RRMSE[1,])


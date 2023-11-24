library(parallel)
library(foreach)
library(doParallel)

num_cores <- 11  # Specify the number of cores you want to use
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  xlabels<-c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
  
  par(mfrow = c(2, 5))
  
  # Assuming you have a dataset 'data' and a black-box regression model 'model'
  for (i in 1:10){
    library('invgamma')
    library('FNN')
    library('expint')
    
    # Create a data frame with X1 values and other variables the same as your dataset
    data_range <-rep(list(X),n)
    
    predictions<-vector(length=11)
    
    # Make predictions on data_range using the black-box regression model
    
    predictions<- foreach(j = c(1,10,20,30,40,50,60,70, 80,90,100),.combine = rbind) %dopar% {
      library('invgamma')
      library('FNN')
      library('expint')
      library("plotrix")
      
      h<-sort(X[,i])
      data_range[[j]][,i] <- rep(h[j],n)
      return(
        mean(AddiVortes_Algorithm(Y,X,200,6000,1000,6,0.85,3,0.8,3,25,Y,data_range[[j]])$mean_yhat_Test)
        )
    }
    
    h<-sort(X[,i])
    # Calculate the partial dependence for X1 by averaging predictions
    plot(h[c(1,10,20,30,40,50,60,70, 80,90,100)],predictions,type='b',ylim=c(9,18), xlab =xlabels[i],ylab = "")
    if (i == 1 | i==6) {
      mtext("Partial dependence",  side = 2, line = 2.5, )  # Adjust line parameter as needed
    }
    
  }
  stopCluster(cl)

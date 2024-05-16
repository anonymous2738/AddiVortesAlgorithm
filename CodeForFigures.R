library(parallel)
library(foreach)
library(doParallel)  

figure5<- function(){

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
}

figure6<-function(){
  f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }

  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 510      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)

  n=length(Y)
  TrainSet=sort(sample.int(n,50*n/51));
  TestSet=1:n;
  TestSet=TestSet[! TestSet %in% TrainSet];
  
  #Creating graph
  
  lineCol <- c("#e23224", "orange", "#762d92", "black", "#0928ff")
  lineType <- c(1, 2, 3, 4, 5)  # Different line types for each line
  labels <- c('m=10', 'm=20', 'm=50', 'm=100', 'm=200')  # Labels for each line
  i <- 1
    
    # Create an empty plot with the x-axis range and labels
  plot(1, type = "n", xlim = c(1, 10), ylim = c(0, 100), xlab = "X-axis", ylab = "Percentage of covariate Used")
    
  for (m in c(10, 20, 50, 100, 200)) {
     if (m == 10) {
      plot(AddiVortes_Algorithm(Y[TrainSet], as.matrix(X[TrainSet, ]), m, 6000, 1000, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, type = "b", col = lineCol[1], lty = lineType[1], xlim = c(1, 10),ylim=c(0,0.3), ylab = "Percentage Used",xlab="Covariate",cex.lab=1.5,lwd=2)
    } 
    else {print("hi")
      lines(AddiVortes_Algorithm(Y[TrainSet], as.matrix(X[TrainSet, ]), m, 6000, 1000, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, col = lineCol[i],lty = lineType[i] , type="b",lwd=2)
    }
      
    # Add custom x-axis labels at breaks from 1 to 10
    axis(1, at = 1:10)
      i <- i + 1
    }

    
  # Add a legend to the plot
  legend("topright", legend = labels, col = lineCol, lty = lineType,lwd=2,inset = 0.1)
}

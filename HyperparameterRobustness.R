 hyperparametersAddiVortesRobust <- matrix(
    c(10, 6,3,
      0.75,0.85,0.95,
      4,3,2,
      1.5,0.8,0.8,
      3,3,1,
      25,25,5), ncol=6
  )
  
  num_samples<- 100
  rmse_values<-numeric(3)
  RobustnessValues<-array(dim=c(3,100,8))
  FullRobustnessValues<-array(dim=c(3,100,8))
  
  num_cores <- 8  # Specify the number of cores you want to use
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  for (m in c(1,10,50,100,200,300,400,500)){
    t<-which(m==c(1,10,50,100,200,300,400,500))
    RobustnessValues[,,t]<-foreach (i = 1:num_samples,.combine=cbind) %dopar% {
      library('invgamma')
      library('FNN')
      library('expint')
      
      for (i in 1:nrow(hyperparametersAddiVortesRobust)){
        params<-hyperparametersAddiVortesRobust[i,]
        print(params)
        nu <- params[1]
        q<- params[2]
        k<- params[3]
        sd<- params[4]
        omega<- params[5]
        lambda<- params[6]
        
        rmse_values[i] <- AddiVortes_Algorithm(Y[TrainSet],as.matrix(X[TrainSet,]),m,1200,200,nu,q,k,sd,omega,lambda,f(X[TestSet,]),as.matrix(X[TestSet,]))$RMSE
      }

      print(rmse_values)
      return(rmse_values)
    }
   print(m)
  }
  
  
  stopCluster(cl)
  
  FullMeanValues<-matrix(nrow = 3,ncol = 8)
  FullLowerquantileValues<-matrix(nrow = 3,ncol = 8)
  FullUpperquantileValues<-matrix(nrow = 3,ncol = 8)
  
  for (i in 1:8){
    FullMeanValues[,i]<-rowMeans(RobustnessValues[,,i])
    FullLowerquantileValues[,i]<- apply(RobustnessValues[,,i], 1, quantile, probs = c(0.05, 0.95))[1,]
    FullUpperquantileValues[,i]<- apply(RobustnessValues[,,i], 1, quantile, probs = c(0.05, 0.95))[2,]
  }
  

  plotCI(c(1,10,50,100,200,300,400,500),FullMeanValues[1,],FullUpperquantileValues[1,]-FullMeanValues[1,],FullMeanValues[1,]-FullLowerquantileValues[1,], sfrac=0, scol = 'red',pch = 16,cex = 0.6,col = 'red',xlim = c(1,500),ylim = c(1,5),lwd=2,xlab = "Number of Tessellations, m", ylab = "RMSE" )
  
  points(c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullMeanValues[2,], pch = 16,cex = 0.6,col = 'blue')
  arrows( c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullUpperquantileValues[2,], c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullLowerquantileValues[2,], angle = 90, code = 3, length = 0, col = "blue",lwd= 2)
  
  points(c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullMeanValues[3,], pch = 16,cex = 0.6, col = 'green')
  arrows( c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullUpperquantileValues[3,], c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullLowerquantileValues[3,], angle = 90, code = 3, length = 0, col = "green",lwd = 2)
  
  legend("topright", legend = c("Aggressive","Default","Conservative"), col = c('red','blue','green'), lwd = 2,inset = 0.1,) 

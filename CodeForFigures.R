library("plotrix")

figure2<-function(){

  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  par(mfrow =c(1,3))
  
  AddiVortes_Algorithm_Plot<-function(y,x,m,max_iter,burn_in,nu,q,k,var,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    #Initialize the Prediction Set, Dimension set and Tessellation Set
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,var))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    plotForSigmaSquared<-vector(length = max_iter)
    plotForRMSE<-vector(length = max_iter)
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=sd(y)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(y ~ x)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(y)-length(x[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(y,yScaled,nu,lambda,SumOfAllTess)
      plotForSigmaSquared[i]=(SigmaSquared*(max(y)-min(y))^2)^0.5
      plotForRMSE[i]=(mean((yScaled-SumOfAllTess)^2))^0.5
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,var) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
      print(i)
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test_ColSum=(colSums(TestMatrix)/(length(TestMatrix[,1])))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    LowerConfidenceTRAINValue<-vector(length=length(mean_yhat))
    UpperConfidenceTRAINValue<-vector(length=length(mean_yhat))
    
    for (i in 1:length(mean_yhat)){
      PredictionMatrix[i,]<-sort(PredictionMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
    
    
    
    LowerConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    UpperConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    
    for (i in 1:length(mean_yhat_Test)){
      TestMatrix[i,]<-sort(TestMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
  
    
    
    plotCI(y,mean_yhat, UpperConfidenceTRAINValue-mean_yhat,mean_yhat-LowerConfidenceTRAINValue,sfrac=0, scol = 'grey',ylab = "posterior intervals", xlab = "In-Sample f(x)",  cex.lab = 1.5) 
    abline(0,1)
    
    plotCI(YTest,mean_yhat_Test,UpperConfidenceTESTValue-mean_yhat_Test,mean_yhat_Test-LowerConfidenceTESTValue,sfrac=0, scol = 'grey', ylab = "posterior intervals", xlab = "Out-of-Sample f(x)", cex.lab = 1.5)
    abline(0,1)
    
    #plot(plotForSigmaSquared,col=c(rep('red',burn_in),rep('black',max_iter-burn_in)),xlab="MCMC iteration",ylab="Sigma draw",type="l")
    plot(plotForSigmaSquared, type = "n", col = "black", lwd = 2, xlab="MCMC iteration",ylab="Sigma draw", cex.lab = 1.5)
    
    # Draw the first segment in red
    segments(x0 = 1:(burn_in - 1), y0 = plotForSigmaSquared[1:(burn_in - 1)],
             x1 = 2:burn_in, y1 = plotForSigmaSquared[2:burn_in], col = "red", lwd = 2)
    # Draw the second segment in black
    segments(x0 = (burn_in):(max_iter - 1), y0 = plotForSigmaSquared[(burn_in+1):(max_iter - 1)],
             x1 = (burn_in+1):max_iter, y1 = plotForSigmaSquared[(burn_in):max_iter], col = "black", lwd = 2)
    abline(1,0)
    
    return(
      data.frame(
        RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
      )
    )
  }
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],200,2000,200,6,0.85,3,0.8,3,25,f(X[TestSet,]),X[TestSet,],IntialSigma = "Linear")

}

figure3<-function(){
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  par(mfrow =c(1,2))
  
  AddiVortes_Algorithm<-function(y,x,m,max_iter,nu,q,k,var,Omega,lambda_rate,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    #Initialize:
    #Prediction Set (A list of vectors with the output values for each tessellation), 
    #Dimension set (A list of vectors with the covariates included in the tessellaions);
    #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,var))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    AverageNumberOfCells<-vector(length = max_iter)
    AverageNumberOfDim<-vector(length = max_iter)
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=sd(y)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(y ~ x)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(y)-length(x[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      NumOfCells<-0
      NumOfDim<-0
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(y,yScaled,nu,lambda,SumOfAllTess)
      
      for (j in 1:m){
        
        
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,var) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
        NumOfCells<-NumOfCells+length(Tess[[j]][,1])
        NumOfDim<-NumOfDim+length(Tess[[j]][1,])
      }
      
      AverageNumberOfCells[i]<-NumOfCells/m
      AverageNumberOfDim[i]<-NumOfDim/m
      print(i)
    }
  
    
    plot(AverageNumberOfDim,type="l",xlab="MCMC iteration", ylab="Average Number Of Dimensions per Tessellaion")
    plot(AverageNumberOfCells,type="l",xlab="MCMC iteration", ylab="Average Number Of Centers per Tessellaion")
  
    
  }
  
  AddiVortes_Algorithm(Y[TrainSet],X[TrainSet,],200,200,6,0.85,3,0.8,3,25,IntialSigma = "Linear")
}

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
    # Create a data frame with X1 values and other variables the same as your dataset
    data_range <-rep(list(X),n)
    
    predictions<-vector(length=11)
    
    # Make predictions on data_range using the black-box regression model
    
    predictions<- foreach(j = c(1,10,20,30,40,50,60,70, 80,90,100),.combine = rbind) %dopar% {
      library('invgamma')
      library('FNN')
      library("plotrix")
      
      h<-sort(X[,i])
      data_range[[j]][,i] <- rep(h[j],n)
      return(
        mean(AddiVortes_Algorithm_Plot(Y,X,200,6000,1000,6,0.85,3,0.8,3,25,Y,data_range[[j]])$mean_yhat_Test)
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

AddiVortes_Algorithm_Plot_figure6<-function(y,x,m,max_iter,burn_in,nu,q,k,var,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){

  p=length(x[1,]) #p is the number of covariates used
  
  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }

  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  #Initialize the Prediction Set, Dimension set and Tessellation Set
  
  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,var))))}
  
  SumOfAllTess=rep(mean(yScaled),length(yScaled));
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2;
  LastTessPred=matrix
  acceptedTess=0;
  PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
  TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
  CovariatesUsed<-rep(0,length(x[1,]))
  plotForSigmaSquared<-vector(length = max_iter)
  plotForRMSE<-vector(length = max_iter)
  AverageNumberOfCells<-vector(length = max_iter)
  AverageNumberOfDim<-vector(length = max_iter)
  
  #finding lambda
  if (IntialSigma=="Naive"){
    SigmaSquaredHat=sd(y)
  }
  else{
    MultiLinear<-lm(y ~ x)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(y)-length(x[1,])-1)
  }
  print(SigmaSquaredHat)
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
  
  print(lambda)
  for (i in 1:max_iter){
    NumOfCells<-0
    NumOfDim<-0
    
    #Sample Whole model Sigma squared
    SigmaSquared=SigmaSquaredCalculation(y,yScaled,nu,lambda,SumOfAllTess)
    print((mean((yScaled-SumOfAllTess)^2))^0.5)
    plotForSigmaSquared[i]=(SigmaSquared*(max(y)-min(y))^2)^0.5
    plotForRMSE[i]=(mean((yScaled-SumOfAllTess)^2))^0.5
    
    for (j in 1:m){
      NewTessOutput<-NewTess(xScaled,j,Tess,Dim,var)
      TessStar<-NewTessOutput[[1]]
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred)
      R_ijOld<-ResidualsOutput[[1]]
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]]
      IndexesStar<-ResidualsOutput[[6]]
      Indexes<-ResidualsOutput[[7]]
      
      if (!any(n_ijNew==0)){
        
        LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate);
        
        if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){
          #print(exp(log(runif(n=1, min=0, max=1))))
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
          LastTessPred=Pred[[j]][IndexesStar]
          acceptedTess=acceptedTess+1
        }
        else {
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{
        Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
      #print(j)
      CovariatesUsed[Dim[[j]]]<-CovariatesUsed[Dim[[j]]]+1
      NumOfCells<-NumOfCells+length(Tess[[j]][,1])
      NumOfDim<-NumOfDim+length(Tess[[j]][1,])
    }
    print(i)
    
    AverageNumberOfCells[i]<-NumOfCells/m
    AverageNumberOfDim[i]<-NumOfDim/m
    
    if (i>burn_in){
      PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
      TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
    }
  }
  
  mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)

  CovariataesUsedPercentage<-CovariatesUsed/sum(CovariatesUsed)
  
  return(
    data.frame(
      CovariataesUsedPercentage
    )
  )
}


  
  f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }

  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 510      #number of observations
  set.seed(23)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)

  n=length(Y)
  TrainSet=sort(sample.int(n,50*n/51));
  TestSet=1:n;
  TestSet=TestSet[! TestSet %in% TrainSet];
  par(mfrow=c(1,1))
  #Creating graph
  
  lineCol <- c("#e23224", "orange", "#762d92", "black", "#0928ff")
  lineType <- c(1, 2, 3, 4, 5)  # Different line types for each line
  labels <- c('m=10', 'm=20', 'm=50', 'm=100', 'm=200')  # Labels for each line
  i <- 1
    
    # Create an empty plot with the x-axis range and labels
  plot(1, type = "n", xlim = c(1, 10), ylim = c(0, 100), xlab = "X-axis", ylab = "Percentage of covariate Used")
    
  for (m in c(10, 20, 50, 100, 200)) {
     if (m == 10) {
      plot(AddiVortes_Algorithm_Plot_figure6(Y[TrainSet], as.matrix(X[TrainSet, ]), m, 6000, 1000, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, type = "b", col = lineCol[1], lty = lineType[1], xlim = c(1, 10),ylim=c(0,0.3), ylab = "Percentage Used",xlab="Covariate",cex.lab=1.5,lwd=2)
    } 
    else {print("hi")
      lines(AddiVortes_Algorithm_Plot_figure6(Y[TrainSet], as.matrix(X[TrainSet, ]), m, 6000, 1000, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, col = lineCol[i],lty = lineType[i] , type="b",lwd=2)
    }
      
    # Add custom x-axis labels at breaks from 1 to 10
    axis(1, at = 1:10)
      i <- i + 1
    }

    
  # Add a legend to the plot
  legend("topright", legend = labels, col = lineCol, lty = lineType,lwd=2,inset = 0.1)
}

figure7<-function(){
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
        
        rmse_values[i] <- AddiVortes_Algorithm_Plot(Y[TrainSet],as.matrix(X[TrainSet,]),m,1200,200,nu,q,k,sd,omega,lambda,f(X[TestSet,]),as.matrix(X[TestSet,]))$RMSE
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

  }

figure9<-function(){
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
  }

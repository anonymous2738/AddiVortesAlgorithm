library("plotrix")

AddiVortes_Algorithm<-function(y,x,m,max_iter,burn_in,nu,q,k,var,Omega,lambda_rate,YTest,XTest,o,IntialSigma = "Linear"){
  
  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  #YTest=(YTest-(max(YTest)+min(YTest))/2)/(max(YTest)-min(YTest));
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
  
  print(Tess)
  
  PredictionTest = TestPrediction(XTest,m,Tess,Dim,Pred)*(max(y)-min(y))+(max(y)+min(y))/2;
  AverageRMSETest=sqrt(mean(((YTest-mean_yhat_Test)^2)));

  PredictionTrain = TestPrediction(xScaled,m,Tess,Dim,Pred)*(max(y)-min(y))+(max(y)+min(y))/2;
  AverageRMSETrain=sqrt(mean((y-mean_yhat)^2));

  
  #plots
  
  plot(AverageNumberOfDim,type="l",xlab="MCMC iteration", ylab="Average Number Of Dimensions per Tessellaion")
  plot(AverageNumberOfCells,type="l",xlab="MCMC iteration", ylab="Average Number Of Centers per Tessellaion")
  
  CovariataesUsedPercentage<-CovariatesUsed/sum(CovariatesUsed)
  plot(CovariataesUsedPercentage)
  
  plotCI(y,mean_yhat, UpperConfidenceTRAINValue-mean_yhat,mean_yhat-LowerConfidenceTRAINValue,sfrac=0, scol = 'grey',ylab = "posterior intervals", xlab = "In-Sample f(x)",  cex.lab = 1.5,main = paste("p=", o)) 
  abline(0,1)
  
  plotCI(YTest,mean_yhat_Test,UpperConfidenceTESTValue-mean_yhat_Test,mean_yhat_Test-LowerConfidenceTESTValue,sfrac=0, scol = 'grey', ylab = "posterior intervals", xlab = "Out-of-Sample f(x)", cex.lab = 1.5,main= paste("p=", o))
  abline(0,1)
  
  #plot(plotForSigmaSquared,col=c(rep('red',burn_in),rep('black',max_iter-burn_in)),xlab="MCMC iteration",ylab="Sigma draw",type="l")
  plot(plotForSigmaSquared, type = "n", col = "black", lwd = 2, xlab="MCMC iteration",ylab="Sigma draw", cex.lab = 1.5,  main= paste("p=", o))
  
  # Draw the first segment in red
  segments(x0 = 1:(burn_in - 1), y0 = plotForSigmaSquared[1:(burn_in - 1)],
           x1 = 2:burn_in, y1 = plotForSigmaSquared[2:burn_in], col = "red", lwd = 2)
  # Draw the second segment in black
  segments(x0 = (burn_in):(max_iter - 1), y0 = plotForSigmaSquared[(burn_in+1):(max_iter - 1)],
           x1 = (burn_in+1):max_iter, y1 = plotForSigmaSquared[(burn_in):max_iter], col = "black", lwd = 2)
  abline(1,0)
  
  return(
    data.frame(
      RMSE = sqrt(mean((YTest-mean_yhat_Test)^2)),
      CovariataesUsedPercentage
    )
  )
}

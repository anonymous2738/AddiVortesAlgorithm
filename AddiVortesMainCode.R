AddiVortes_Algorithm<-function(y,x,m,max_iter,burn_in,nu,q,k,var,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){
  
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
  
  #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
  PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
  TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
  
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
  }
  
  #finding the mean of the predition over the iterations and then unscaling the predictions.
  mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  
  return( #Returns the RMSE value for the test samples.
    data.frame(
      RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
    )
  )
}

hello<- function(){
  print('Hello')
  }


library('invgamma')
library('FNN')

SigmaSquaredCalculation<-function(y,yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
  
  n=length(y)
  SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+(max(y)-min(y))^2*sum((yScaled-SumOfAllTess)^2))/2)
  
  return(SigmaSquared/(max(y)-min(y))^2)
}

NewTess<-function(x,j,Tess,Dim,var){ #Propose a new tessellation
  
  p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
  
  DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar=Tess #Similar for the tessellation matrix.
  
  if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
    NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
    DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
    TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,var)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
    Modification="AD"}
  else if (p<0.4){ #Remove a dimension if p is less then 0.4.
    RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
    DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
    TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
    Modification="RD"}
  else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
    TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,var)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
    Modification="AC"}
  else if (p<0.8){ #Add a centre if p is less then 0.8. 
    CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
    TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
    Modification="RC"}
  else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
    TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,var) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
    Modification="Change"}
  else{ #Swop a dimension.
    NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
    DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
    DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
    TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,var) #Add new normally sampled coordinates new dimension added.
    Modification="Swop"}
  
  TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
  
  return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
}

fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
  return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
}

Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
  if (length(Tess[,1])==1){ #only 1 centre
    CellsForGivenTess=rep(1,length(x[,1]))
  }
  else{ #multiple
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
  }
  return(CellsForGivenTess)
}

AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
  
  d=length(Dim[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(Tess[[j]][,1]);

  #The Log Likelihood Ratio in the acceptence ratio
  LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))

  #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(2)}
  }
  else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)

    #Adjustments.
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  }
  else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)

    #Adjustments.
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred}
  else{
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  #Initializing Sizes
  
  R_ijOld=rep(0,length(Pred[[j]]))
  n_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[indexes==i])
    n_ijOld[i]<-sum(indexes==i)
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  n_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i])
    n_ijNew[i]<-sum(IndexesStar==i)}
  
  return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
}

NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
  }
  return(PredSet)
}

TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}


f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
}

sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
n = 200      #number of observations
set.seed(9)
X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(X)
Y=Ey+sigma*rnorm(n)

n=length(Y)
TrainSet=sort(sample.int(n,3*n/6))
TestSet=1:n
TestSet=TestSet[! TestSet %in% TrainSet]

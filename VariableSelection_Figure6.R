#Friedman data for the graph.

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

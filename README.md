AddiVortes
===========

Implementation of the (Bayesian) Additive Voronoi Tessellation (AddiVortes) algorithm in Rstudio.

Copyright (C) 2024
Adam Stone  
Department of Mathematical sciences, Durham University
& 
John Paul Gosling
Department of Mathematical sciences, Durham University
 
Setup Instructions
------------------

To install the AddiVortes algorithm in R, you can run the following code in Rstudio: 

```r

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/AddiVortesMainCode.R")

```


Benchmark Real Datasets
-----------------------------

To import real world datasets in Rstudio one can run the following code:

```r
source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/Datasets.R")

```

For each dataset, this imports the full datasets but also the feature matrix (X_dataset) and output variable (Y_dataset) so that the Addivortes algorithm can be implemented. For example, for the Boston to run the AddiVortes algorithm:

```r
Boston #Full Boston dataset

n=length(Y_Boston)
TrainSet=sort(sample.int(n,5*n/6))
TestSet=1:n
TestSet=TestSet[! TestSet %in% TrainSet]

 AddiVortes_Algorithm(Y_Boston[TrainSet],X_Boston[TrainSet,],200,2000,200,6,0.85,3,0.8,3,25,Y_Boston[TestSet,],X_Boston[TestSet,],IntialSigma = "Linear")

```
Reproducing Figures in paper 
---------------------------


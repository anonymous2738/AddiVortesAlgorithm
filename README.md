AddiVortes
===========

Implementation of the (Bayesian) Additive Voronoi Tessellation (AddiVortes) algorithm in Rstudio.

Copyright (C) 2024

Adam Stone  & John Paul Gosling

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

To import the real world benchmark datasets used in the paper in Rstudio one can run the following code:

```r
source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/Datasets.R")

```

For each dataset, this imports the full datasets but also the feature matrix (X_dataset) and output variable (Y_dataset) so that the Addivortes algorithm can be implemented. For example, to run the AddiVortes algorithm for the Boston dataset, use the following code:

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

To reproduce the figures in the paper, source the following Github page by running:

```r

source_url("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/main/CodeForFigures.R")

```

### Real World datsets Boxplot


**Warning:** This figure use parallel processing using 10 cores at a time, producing this figure is only recommended if you have 12+ cores.

### Friedman Fimulation Figures

To produce figures 3, 4 or 8 you can simply run:

```r
figure3() #Approximate time 2 minutes 30 seconds
figure4() #Approximate time 2 minutes 30 seconds
figure8() #Approximate time 8 minutes
```
The other figures take longer computational time to produce so we give the option to reduce the number of interations or predictions made, compared to ones used in the paper. The results may be slightly worse but still show strength of the AddiVortes model.

**Warning:** These figures use parallel processing using up to 10 cores at a time, producing these figures is only recommended if you have 12+ cores.

```r
figure5(max_iter = 6000 , burn_in = 1000) 
figure6(max_iter = 6000 , burn_in = 1000) 
figure7(max_iter = 6000 , burn_in = 1000, num_samples = 1000) 
figure9() 
```



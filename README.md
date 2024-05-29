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

AddiVortes_algorithm()

```
Reproducing Figures in paper 
---------------------------


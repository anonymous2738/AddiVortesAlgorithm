# Load necessary library if you haven't already
if (!requireNamespace("invgamma", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)

# Import the dataset into R
abalone <- read_csv(https://github.com/anonymous2738/AddiVortesAlgorithm/blob/DataSets/abalone_Data.csv)
X_abalone<-abalone[,2:8]
Y_abalone<-abalone[,9]

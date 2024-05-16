# Load necessary library if you haven't already
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)

# Import the dataset into R
abalone <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/abalone_Data.csv")
X_abalone<-abalone[,2:8]
Y_abalone<-abalone[,9]

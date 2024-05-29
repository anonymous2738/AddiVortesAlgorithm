# Load necessary library if you haven't already
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(readr)

# Import the dataset into R
abalone <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/abalone_Data.csv")
X_abalone<-as.matrix(abalone[,2:8])
Y_abalone<-as.numeric(as.matrix(abalone[,9]))



Baskball <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Baskball_Data.csv")
X_Baskball<-Baskball[,2:5]
Y_Baskball<-Baskball[,6]

Boston <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/BostonHousing_Data.csv")
X_Boston<-Boston[,2:14]
Y_Boston<-Boston[,15]

Enroll <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Enroll_Data.csv")
X_Enroll<-Enroll[,2:7]
Y_Enroll<-Enroll[,8]

Fat <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/FAT_Data.csv")
X_Fat<-Fat[,2:14]
Y_Fat<-Fat[,15]

HatCo <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/HatCo_Data.csv")
X_HatCo<-HatCo[,2:13]
Y_HatCo<-HatCo[,14]

Labour <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Labour_Data.csv")
X_Labour<-Labour[,2:20]
Y_Labour<-Labour[,21]

Mpg <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/MPG_Data.csv")
X_Mpg<-Mpg[,2:8]
Y_Mpg<-Mpg[,9]

Medicare <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Medicare_Data.csv")
X_Medicare<-Medicare[,2:7]
Y_Medicare<-Medicare[,8]

Price <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Price_Data.csv")
X_Price<-Price[,2:16]
Y_Price<-Price[,17]

Rate <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/RATE_Data.csv")
X_Rate<-Rate[,2:10]
Y_Rate<-Rate[,11]

Edu <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/edu_Data.csv")
X_Edu<-Edu[,2:7]
Y_Edu<-Edu[,8]

Ozone <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/ozone_Data.csv")
X_Ozone<-Ozone[,2:8]
Y_Ozone<-Ozone[,9]


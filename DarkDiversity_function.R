install.packages("vegan")
library(vegan)
data( dune)
#import data
as.data.frame(table(fauna2))
indx <- sapply(fauna2, is.factor)
fauna2[indx] <- lapply(fauna2[indx], function(x) as.numeric(as.character(x)))
source( "D:/dark.pred.ca.R")
source("D:/Moumita_Research/MGhosh_2019/Task/R")

# Predict dark diversity for the first 5 cases based on the remaining 15
dpc.1 <- dark.pred.ca( fauna2[6:20,], fauna2[1:5,], "binminpred", 0.7)
dpc.1 <- dark.pred.ca( fauna2[6:20,], fauna2[1:5,], "minpred", 0.8)
k<- dpc.1$pred.mat # predicted abundances in the species pool (rescaled)
dpc.1$pred.mat

# Write the result
write.csv(dpc.1$pred.mat,"D:/mydata_minpred.csv")
pa <- decostand(fauna2, "pa")
boxplot(as.vector(k) ~ unlist(pa), xlab="Presence", ylab="UNO")
 
 
 

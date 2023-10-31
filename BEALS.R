install.packages("vegan")
library(vegan)
\*import data set*\
as.data.frame(table(fauna2))
indx <- sapply(fauna2, is.factor)
fauna2[indx] <- lapply(fauna2[indx], function(x) as.numeric(as.character(x)))

#beals(fauna2)
x <- beals(fauna2)
## Smoothed values against presence or absence of species
pa <- decostand(fauna2, "pa")
boxplot(as.vector(x) ~ unlist(pa), xlab="Presence", ylab="Beals")
 
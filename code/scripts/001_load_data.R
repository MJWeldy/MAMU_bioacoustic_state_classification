#Generate data R script
library(ggplot2)
Y <- rnorm(1000, 10, 2)
X <- rnorm(1000, 5, 1)
df <- data.frame(X = X, Y = Y)
head(df)
return(print(head(df)))
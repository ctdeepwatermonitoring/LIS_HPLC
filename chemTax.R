# Modified from example - package limSolve Chemtax function
# https://cran.r-project.org/web/packages/limSolve/limSolve.pdf
library(BCE)
library(limSolve)

#1.Pigment composition of the various algal groups (based on Li 2004?)

ratio <- read.csv("data/ratios.csv",header = TRUE, row.names=1)
ratio <- ratio *100
ratio <- as.matrix(ratio)

#2.Pigment composition of the samples measured with HPLC. Ratio pigment:Chla calculated
#First try just with the experiments 1 and 2 data

field <- read.csv("data/field.csv", header=TRUE, row.names=1)

i <- 2

# 1. Graphical representation of the chemtax example input data
palette(rainbow(12, s = 0.6, v = 0.75))
mp <- apply(ratio, MARGIN = 2, max)
pstars <- rbind(t(t(ratio)/mp) ,
                sample = field[i,]/max(field[i,]))
stars(pstars, len = 0.9, key.loc = c(7.2, 1.7),scale=FALSE,ncol=4,
      main = "Pigment composition", draw.segments = TRUE,
      flip.labels=FALSE)

# 2. Estimating the algal composition of the field sample
# Nx <-nrow(Chemtax$Ratio) # Example Data
Nx <-nrow(ratio)

# equations that have to be met exactly Ex=f:
# sum of all fraction must be equal to 1.
EE <- rep(1, Nx)
FF <- 1

# inequalities, Gx>=h:
# all fractions must be positive numbers
GG <- diag(nrow = Nx)
HH <- rep(0, Nx)

# equations that must be reproduced as close as possible, Ax ~ b
# = the field data; the input ratio matrix and field data are rescaled
# AA <- Chemtax$Ratio/rowSums(Chemtax$Ratio) #Example Data
# BB <- Chemtax$Field/sum(Chemtax$Field) #Example Data

AA <- ratio/rowSums(ratio)
BB <- as.numeric(field[i,]/sum(field[i,]))


# 1. Solve with lsei method
X <- lsei(t(AA), array(BB), EE, FF, GG, HH)$X
Sample <- data.frame(Algae = rownames(ratio),
                      fraction = X)

# plot results
barplot(X, names = rownames(ratio), col = heat.colors(length(rownames(ratio))),
        cex.names = 0.5, main = "LIS example solved with lsei")


# 2. Bayesian sampling;
# The standard deviation on the field data is assumed to be 0.01
# jump length not too large or NO solutions are found!
xs <- xsample(t(AA), BB, EE, FF, GG, HH, sdB = 0.01, jmp = 0.025)$X
pairs(xs, main= "Chemtax, Bayesian sample")
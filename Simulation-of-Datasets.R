############################################
##   Simulation of the datasets           ##
##       analysed in the simulation study ##
############################################

RNGkind(sample.kind = "Rounding")
set.seed(280617)


##required libraries
library(drc)

#####
##functions needed

##four-parametric log logistic function
my.loglog <- function(x, b, c, d, e){
  return(c + (d-c)/(1+exp(b*(log(x)-log(e)))))
}

#####
##further constant values: the doses where gene expression is measured and the number of samples per dosis
doses <- c(0, 25, 150, 350, 450, 550, 800, 1000)
samples <- 3

##matrix of parameters corresponding to the three scenarios
pars <- matrix(NA, nrow=3, ncol=4)
colnames(pars) <- c("b", "c", "d", "e")
rownames(pars) <- c("Scenario I", "Scenario II", "Scenario III")
pars[1,] <- c(-6, 0, 0.58, 450)
pars[2,] <- c(-3, 0, 4, 900)
pars[3,] <- c(-3, 0, 1.16, 400)


##The standard deviations for the simulations are depending on the range of the curve that the simulation is based on
## they were chosen fitting a linear model that models the standard deviation depending on the absolute value of the ranges.
## "small" and "large" come from these standard deviations by dividing by two or multiplying by two
std.devs <- matrix(NA, nrow=3, ncol=3)
colnames(std.devs) <- c("Small", "Middle", "Large")
rownames(std.devs) <- c("Scenario I", "Scenario II", "Scenario III")
std.devs[1,] <- c(0.095, 0.189, 0.379)
std.devs[2,] <- c(0.131, 0.261, 0.522)
std.devs[3,] <- c(0.107, 0.213, 0.427)


##Calculate the cross product of all combinations. 
##all in all, 9 combinations (3 scenarios with 3 standard deviations each) are considered
situations <- expand.grid(c("Small", "Middle", "Large"), c("Scenario I", "Scenario II", "Scenario III"),
                          stringsAsFactors = FALSE)
colnames(situations) <- c("Var", "Scen")


####
##Simulation of the data:
for(co in 1:9){
  comb <- situations[co,]

  ##the response values calculated for the concentrations considered, based on the choice of true parameter 
    ##in the respective scenario
  resp <-  my.loglog(rep(doses, each=samples), pars[comb$Scen, 1], pars[comb$Scen, 2], pars[comb$Scen, 3], pars[comb$Scen, 4])

  ##prepare the 1000 x 24 matrix, in which the sampled expression values are stored
  mat.expression <- matrix(NA, ncol=length(resp), nrow=1000)
  colnames(mat.expression) <- paste0(rep(doses, each=samples), ".", rep(1:3, 8))

  ##for 1000 simulation runs, add normally distributed noise to the response values calculated above, 
    ##thus yielding the simulated concentration-expression matrix
  for(i in 1:1000){
    mat.expression[i,] <- resp + rnorm(length(resp), mean = 0, sd = std.devs[comb$Scen, comb$Var])
  }

  ##split this before saving to avoid spaces in filename
  scen.save <- strsplit(comb$Scen, " ")[[1]][2]
  ##save this dataset
  save(mat.expression, file = paste0("DataSet.Scenario", scen.save, ".Var", comb$Var, ".RData"))
}










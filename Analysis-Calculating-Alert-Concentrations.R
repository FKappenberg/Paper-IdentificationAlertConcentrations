####################################################
## Calculation of the four alert concentrations   ##
##            ALOEC, LOEC, ALEC, LEC              ##
####################################################

##source file: Functions-Calculating-Alert-Concentrations contains all the functions needed
  ##for calculation of the alert concentrations
source("Functions-Calculating-Alert-Concentrations.R")


##constant values:
doses <- c(0, 25, 150, 350, 450, 550, 800, 1000)
samples <- 3
threshold <- log2(1.5)



##########################
## Variance Small

load("DataSet.ScenarioI.VarSmall.RData")
scen1.small <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioII.VarSmall.RData")
scen2.small <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioIII.VarSmall.RData")
scen3.small <- mat.expression
rm(mat.expression)



##################
## Calculation of all alert concentrations
aloecs.scen1.small <- apply(scen1.small, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen1.small <- apply(scen1.small, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen1.small <- apply(scen1.small, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen1.small <- apply(scen1.small, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen1.small <- t(sapply(models.scen1.small, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


aloecs.scen2.small <- apply(scen2.small, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen2.small <- apply(scen2.small, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen2.small <- apply(scen2.small, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen2.small <- apply(scen2.small, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen2.small <- t(sapply(models.scen2.small, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


aloecs.scen3.small <- apply(scen3.small, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen3.small <- apply(scen3.small, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen3.small <- apply(scen3.small, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen3.small <- apply(scen3.small, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen3.small <- t(sapply(models.scen3.small, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


save(aloecs.scen1.small, aloecs.scen2.small, aloecs.scen3.small,
     loecs.t.scen1.small, loecs.t.scen2.small, loecs.t.scen3.small,
     loecs.dun.scen1.small, loecs.dun.scen2.small, loecs.dun.scen3.small,
     vpll.scen1.small, vpll.scen2.small, vpll.scen3.small,
     file="Benchmark.VarSmall.TwoSided.RData")


##
# Determine the number of negative diagonal entries of the covariance matrix
cov.scen1.small <- sapply(models.scen1.small, function(mod) diag(vcov(mod)))
neg.cov.scen1.small <- union(union(which(cov.scen1.small[1,] < 0), which(cov.scen1.small[2,] < 0)),
                             union(which(cov.scen1.small[3,] < 0), which(cov.scen1.small[4,] < 0)))

cov.scen2.small <- sapply(models.scen2.small, function(mod) diag(vcov(mod)))
neg.cov.scen2.small <- union(union(which(cov.scen2.small[1,] < 0), which(cov.scen2.small[2,] < 0)),
                             union(which(cov.scen2.small[3,] < 0), which(cov.scen2.small[4,] < 0)))

cov.scen3.small <- sapply(models.scen3.small, function(mod) diag(vcov(mod)))
neg.cov.scen3.small <- union(union(which(cov.scen3.small[1,] < 0), which(cov.scen3.small[2,] < 0)),
                             union(which(cov.scen3.small[3,] < 0), which(cov.scen3.small[4,] < 0)))

save(neg.cov.scen1.small, neg.cov.scen2.small, neg.cov.scen3.small,
     file="NegativeCovarianceEntries-VarSmall.RData")


###############################
## Variance Middle

load("DataSet.ScenarioI.VarMiddle.RData")
scen1.middle <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioII.VarMiddle.RData")
scen2.middle <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioIII.VarMiddle.RData")
scen3.middle <- mat.expression
rm(mat.expression)


aloecs.scen1.middle <- apply(scen1.middle, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen1.middle <- apply(scen1.middle, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen1.middle <- apply(scen1.middle, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen1.middle <- apply(scen1.middle, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})


vpll.scen1.middle <- t(sapply(models.scen1.middle, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


aloecs.scen2.middle <- apply(scen2.middle, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen2.middle <- apply(scen2.middle, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen2.middle <- apply(scen2.middle, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen2.middle <- apply(scen2.middle, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen2.middle <- t(sapply(models.scen2.middle, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


aloecs.scen3.middle <- apply(scen3.middle, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen3.middle <- apply(scen3.middle, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen3.middle <- apply(scen3.middle, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen3.middle <- apply(scen3.middle, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen3.middle <- t(sapply(models.scen3.middle, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


save(aloecs.scen1.middle, aloecs.scen2.middle, aloecs.scen3.middle,
     loecs.t.scen1.middle, loecs.t.scen2.middle, loecs.t.scen3.middle,
     loecs.dun.scen1.middle, loecs.dun.scen2.middle, loecs.dun.scen3.middle,
     vpll.scen1.middle, vpll.scen2.middle, vpll.scen3.middle,
     file="Benchmark.VarMiddle.TwoSided.RData")


##
# Determine the number of negative diagonal entries of the covariance matrix
cov.scen1.middle <- sapply(models.scen1.middle, function(mod) diag(vcov(mod)))
neg.cov.scen1.middle <- union(union(which(cov.scen1.middle[1,] < 0), which(cov.scen1.middle[2,] < 0)),
                              union(which(cov.scen1.middle[3,] < 0), which(cov.scen1.middle[4,] < 0)))

cov.scen2.middle <- sapply(models.scen2.middle, function(mod) diag(vcov(mod)))
neg.cov.scen2.middle <- union(union(which(cov.scen2.middle[1,] < 0), which(cov.scen2.middle[2,] < 0)),
                              union(which(cov.scen2.middle[3,] < 0), which(cov.scen2.middle[4,] < 0)))

cov.scen3.middle <- sapply(models.scen3.middle, function(mod) diag(vcov(mod)))
neg.cov.scen3.middle <- union(union(which(cov.scen3.middle[1,] < 0), which(cov.scen3.middle[2,] < 0)),
                              union(which(cov.scen3.middle[3,] < 0), which(cov.scen3.middle[4,] < 0)))


save(neg.cov.scen1.middle, neg.cov.scen2.middle, neg.cov.scen3.middle,
     file="NegativeCovarianceEntries-VarMiddle.RData")


################################
## Variance Large

load("DataSet.ScenarioI.VarLarge.RData")
scen1.large <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioII.VarLarge.RData")
scen2.large <- mat.expression
rm(mat.expression)

load("DataSet.ScenarioIII.VarLarge.RData")
scen3.large <- mat.expression
rm(mat.expression)


aloecs.scen1.large <- apply(scen1.large, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen1.large <- apply(scen1.large, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen1.large <- apply(scen1.large, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen1.large <- apply(scen1.large, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen1.large <- t(sapply(models.scen1.large, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))



aloecs.scen2.large <- apply(scen2.large, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen2.large <- apply(scen2.large, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen2.large <- apply(scen2.large, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen2.large <- apply(scen2.large, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen2.large <- t(sapply(models.scen2.large, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


aloecs.scen3.large <- apply(scen3.large, 1, function(x) aloec.gene(x, doses, alt="two.sided"))
loecs.t.scen3.large <- apply(scen3.large, 1, function(x) loec.t.gene(x, doses, alt="two.sided"))
loecs.dun.scen3.large <- apply(scen3.large, 1, function(x) loec.dun.gene(x, doses, alt="two.sided"))

models.scen3.large <- apply(scen3.large, 1, function(x){
  object <- drm(x ~ rep(doses, each=3), fct=LL2.4())
  return(object)
})

vpll.scen3.large <- t(sapply(models.scen3.large, function(mod){
  erg <- ALEC(mod, alt="two.sided")
  alec <- ifelse(erg$ALEC<=1000, erg$ALEC, NA)
  alec.ci.low <-ifelse(erg$ALEC<=1000, erg$CI.lower, NA)
  alec.ci.up <- ifelse(erg$ALEC<=1000, erg$CI.upper, NA)
  
  erg <- LEC(mod, alt="two.sided")
  lec <- ifelse(erg$Iteration == "Max.Iter", NA, erg$disk.conc)
  
  res <- c(alec, alec.ci.low, alec.ci.up, lec)
  names(res) <- c("ALEC", "ALEC:LowerCI", "ALEC:UpperCI", "LEC")
  return(res)
}))


save(aloecs.scen1.large, aloecs.scen2.large, aloecs.scen3.large,
     loecs.t.scen1.large, loecs.t.scen2.large, loecs.t.scen3.large,
     loecs.dun.scen1.large, loecs.dun.scen2.large, loecs.dun.scen3.large,
     vpll.scen1.large, vpll.scen2.large, vpll.scen3.large,
     file="Benchmark.VarLarge.TwoSided.RData")

##
# Determine the number of negative diagonal entries of the covariance matrix

cov.scen1.large <- sapply(models.scen1.large, function(mod) diag(vcov(mod)))
neg.cov.scen1.large <- union(union(which(cov.scen1.large[1,] < 0), which(cov.scen1.large[2,] < 0)),
                             union(which(cov.scen1.large[3,] < 0), which(cov.scen1.large[4,] < 0)))

cov.scen2.large <- sapply(models.scen2.large, function(mod) diag(vcov(mod)))
neg.cov.scen2.large <- union(union(which(cov.scen2.large[1,] < 0), which(cov.scen2.large[2,] < 0)),
                             union(which(cov.scen2.large[3,] < 0), which(cov.scen2.large[4,] < 0)))

cov.scen3.large <- sapply(models.scen3.large, function(mod) diag(vcov(mod)))
neg.cov.scen3.large <- union(union(which(cov.scen3.large[1,] < 0), which(cov.scen3.large[2,] < 0)),
                             union(which(cov.scen3.large[3,] < 0), which(cov.scen3.large[4,] < 0)))

save(neg.cov.scen1.large, neg.cov.scen2.large, neg.cov.scen3.large,
     file="NegativeCovarianceEntries-VarLarge.RData")

























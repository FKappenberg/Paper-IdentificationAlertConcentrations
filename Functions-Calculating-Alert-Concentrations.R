############################################################################
## All functions needed to calculate the different alert concentrations   ##
##     ALOEC, LOEC, ALEC, LEC                                             ##
############################################################################

##required libraries
library(drc)
library(multcomp) #for glht, specifically the Dunnett procedure


############################################
## Observation-based alert concentrations ##
############################################

###
## t.Test: Calculation of compound-wise fold changes and p-values of a t-test

#Input: 
#gene     ->  Vector with all expression values for a single gene with increasing doses and all replicates
#doses    ->  Vector with all doses in increasing order (single doses)
#nrep     ->  number of replicates per dose (not control), equal number for all doses is assumed
#ncontrol ->  number of replicates for control, may differ from nrep
#mu0      ->  true difference in means that is to be tested. Defaults to log2(1.5)
#alpha    ->  significance level (not really needed since p-values are considered)
#alt      ->  character string specifying the alternative hypothesis (defaults to "greater")

#Output:
#Dataframe that contains two variables and number of doses many observations. Columns of the dataframe:
#MeanDiff -> fold changes between the corresponding dose and the control
#P.Values -> p-value of the two-sample t-tests between the corresponding dose and the control



t.tests.loec <- function(gene, doses, nrep=3, ncontrol=3, mu0=log2(1.5), alpha=0.05, alt="greater"){
  
  ##Calculation of the difference in mean between each dose and the control
  diff <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
    mean(x) - mean(gene[(1:ncontrol)])
  })   
  names(diff) <- paste0(doses[-1], "-Ctrl")
  
  if(alt == "greater"){
    ##Calculation of a p-value for each two-sample t-test for the comparison of each concentration vs control
    pvals <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
      t.test(x, gene[(1:ncontrol)], mu=mu0, alternative="greater", conf.level=1-alpha)$p.value
    })   
  }
  if(alt == "less"){
    ##Calculation of a p-value for each two-sample t-test for the comparison of each concentration vs control
    pvals <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
      t.test(x, gene[(1:ncontrol)], mu=-mu0, alternative="less", conf.level=1-alpha)$p.value
    })   
  }
  if(alt == "two.sided"){
    ##Calculation of a p-value for each two-sample t-test for the comparison of each concentration vs control
    ##Two one-sided p-values are calculated and the final p-values is given as 2*min(p1, p2)
    pvals1 <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
      t.test(x, gene[(1:ncontrol)], mu=mu0, alternative="greater", conf.level=1-alpha)$p.value
    })  
    pvals2 <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
      t.test(x, gene[(1:ncontrol)], mu=-mu0, alternative="less", conf.level=1-alpha)$p.value
    })   
    pvals <- 2*pmin(pvals1, pvals2)
  }
  names(pvals) <- paste0(doses[-1], "-Ctrl")
  
  res <- data.frame(MeanDiff = diff, P.Values=pvals)
  return(res)
}


##Calculation of the LOEC based on the t-test
loec.t.gene <- function(gene, doses, nrep=3, ncontrol=3, mu0=log2(1.5), alpha=0.05, alt="greater"){
  ##apply the t.test function to one specific gene
  erg <- t.tests.loec(gene, doses, nrep=nrep, ncontrol=ncontrol, mu0=mu0, alpha=alpha, alt=alt)
  
  ##only genes with unambiguous profiles are considered (regardless of the type of analysis)
  ##i.e. if any gene exist with significant results in both directions, no LOEC can be calculated
  if((any(erg$P.Values < alpha & erg$MeanDiff>mu0)) & (any(erg$P.Values < alpha & erg$MeanDiff<(-mu0)))){
    loec <- NA
  }
  ##in the case of unambiguous profiles, the loec is the lowest concentration where the p-value is smaller than alpha
  else{
    loec <- ifelse(any(erg$P.Values<alpha), doses[1+min(which(erg$P.Values<alpha))], NA)
  }
  return(loec)
}

##Calculation of the ALOEC
##although no significance is assessed here, calculation of the ALOEC is also based on the t.test function
  ##since the fold changes are returned by this function as well
aloec.gene <- function(gene, doses, nrep=3, ncontrol=3, mu0=log2(1.5), alpha=0.05, alt="greater"){
  ##apply the t.test function to one specific gene
  erg <- t.tests.loec(gene, doses, nrep=nrep, ncontrol=ncontrol, mu0=mu0, alpha=alpha)
  ##only genes with unambiguous profiles are considered (regardless of the type of analysis)
  if((any(erg$MeanDiff > mu0)) & any(erg$MeanDiff < -mu0)){
    aloec <- NA
  }
  ##in the case of an unambiguous profile, the aloec is the lowest concentration where the fold change is smaller than -mu0 / larger than mu0
  else{
    if(alt == "greater") aloec <- ifelse(any(erg$MeanDiff > mu0), doses[1+min(which(erg$MeanDiff > mu0))], NA)
    if(alt == "smaller") aloec <- ifelse(any(erg$MeanDiff < -mu0), doses[1+min(which(erg$MeanDiff < -mu0))], NA)
    if(alt == "two.sided") aloec <- ifelse(any(abs(erg$MeanDiff) > mu0), doses[1+min(which(abs(erg$MeanDiff) > mu0))], NA)
  }
  return(aloec)
}


###
## Dunnett-procedure for the alternative calculation of the LOEC

#Input: 
#gene     -> Vector with all expression values for a single gene with increasing doses and all replicates
#doses    -> Vector with all doses in increasing order
#nrep     -> number of replicates per dose (not control), equal number for all doeses is assumed
#ncontrol -> number of replicates for control, may differ from nrep
#mu0      -> true difference in means that is to be tested
#alpha    -> significance level (not really needed since p-values are considered)
#alt      -> specifying the alternative hypothesis (defaults to "greater")

#Output:
#Dataframe that contains two variables and number of doses many observations. Columns of the dataframe:
#MeanDiff -> fold changes between the corresponding dose and the control
#P.Values -> p-value of the two-sample t-tests between the corresponding dose and the control

dunnett.tests.loec <- function(gene, doses, nrep=3, ncontrol=3, mu0=log2(1.5), alpha=0.05, alt="greater"){
  
  ##Creation of a dataset to which a linear model can be fitted
  df <- data.frame(Gene=gene, Dose=factor(c(rep(doses[1], ncontrol), rep(doses[-1],each=nrep))))
  mod <- lm(Gene ~ Dose, data=df)
  
  ##Application of the "glht" function to the linear model, extraction of the p-values
  if(alt == "greater"){
    pvals <- summary(glht(mod, linfct=mcp(Dose="Dunnett"), rhs=log2(1.5), alternative="greater", alpha=alpha))$test$pvalue
    names(pvals) <- paste0(doses[-1], "-Ctrl")
  }
  if(alt == "less"){
    pvals <- summary(glht(mod, linfct=mcp(Dose="Dunnett"), rhs=-log2(1.5), alternative="less", alpha=alpha))$test$pvalue
    names(pvals) <- paste0(doses[-1], "-Ctrl")
  }
  if(alt == "two.sided"){
    pvals1 <- summary(glht(mod, linfct=mcp(Dose="Dunnett"), rhs=log2(1.5), alternative="greater", alpha=alpha))$test$pvalue
    pvals2 <- summary(glht(mod, linfct=mcp(Dose="Dunnett"), rhs=-log2(1.5), alternative="less", alpha=alpha))$test$pvalue
    pvals <- 2*pmin(pvals1, pvals2)
    names(pvals) <- paste0(doses[-1], "-Ctrl")
  }
  
  ##Calculation of the difference in mean between each dose and the control
  diff <- tapply(gene[-(1:ncontrol)], rep(doses[-1], each=nrep), function(x){
    mean(x) - mean(gene[(1:ncontrol)])
  })   
  names(diff) <- paste0(doses[-1], "-Ctrl")
  
  res <- data.frame(MeanDiff = diff, P.Values=pvals)
  return(res)
}

##Calculation of the LOEC based on the Dunnett procedure
loec.dun.gene <- function(gene, doses, nrep=3, ncontrol=3, mu0=log2(1.5), alpha=0.05, alt="greater"){
  ##apply the dunnett.test function to one specific gene
  erg <- dunnett.tests.loec(gene, doses, nrep=nrep, ncontrol=ncontrol, mu0=mu0, alpha=alpha, alt=alt)
  
  ##only genes with unambiguous profiles are considered (regardless of the type of analysis)
  if((any(erg$P.Values < alpha & erg$MeanDiff>mu0)) & (any(erg$P.Values < alpha & erg$MeanDiff<(-mu0)))){
    loec <- NA
  }
  ##in the case of an unambiguous profile, the loec is the lowest concentration where the p-value is smaller than alpha
  else{
    loec <- ifelse(any(erg$P.Values<alpha), doses[1+min(which(erg$P.Values<alpha))], NA)
  }
  return(loec)
}



######################################
## Model-based alert concentrations ##
######################################


###
## ALEC-function for calculating the alec together with the confidence interval

##Input:
#object -> a fitted drc object (fitted with LL2.4())
#lambda -> the response value for which the corresponding concentration (ALEC) is to be calculated
#alpha  -> The niveau for the confidence interval (on the dosis axis!)
#alt  -> specifying the alternative hypothesis (defaults to "greater")

##Output:
#Data frame with the following information:
#ALEC         -> Estimated concentration
#CI.lower     -> Lower limit of the confidence interval (on the dosis axis!)
#CI.upper     -> Upper limit of the confidence interval (on the dosis axis!)
#sd.ALEC      -> standard deviation of the estimated ALEC, calculated using "our own" delta rule
#log2(lambda) -> The original response value 


ALEC <- function(object, lambda=log2(1.5), alpha=0.05, alt="greater"){
  library("msm")
  library("drc")
  
  ##extract the four parameters from the fitted drc-object
  ##note: since fitting was conducted using "LL2.4()", the parameter e is actually the re-parameterised version log(e)
  b <- coef(object)['b:(Intercept)']
  c <- coef(object)['c:(Intercept)']
  d <- coef(object)['d:(Intercept)']
  e <- coef(object)['e:(Intercept)']
  
  ##Not the intersection of the curve with lambda itself is of interest, but with lower asymptote +/- lambda
  ##f0 is the value of the left asymptote regardless of whether the curve is increasing or decreasing
  f0 <- ifelse(b>0, d, c)
  
  if(alt == "greater"){
    thresh <- f0 + lambda
    ##Calculation of the ALEC via the inverse, note that the parameter e actually is log(e)
    est <- exp(e)*((d-thresh)/(thresh-c))^(1/b)
  } 
  if(alt == "less"){
    thresh <- f0 - lambda
    ##Calculation of the ALEC via the inverse, note that the parameter e actually is log(e)
    est <- exp(e)*((d-thresh)/(thresh-c))^(1/b)
  }
  if(alt == "two.sided"){
    thresh.1 <- f0-lambda
    thresh.2 <- f0+lambda
    est.1 <- exp(e)*((d-thresh.1)/(thresh.1-c))^(1/b)
    est.2 <- exp(e)*((d-thresh.2)/(thresh.2-c))^(1/b)
    ##only one of the calculated thresholds is in the range between the asymptote values, 
      ##thus only for one threshold, a valid estimate is obtained. This is then the ALEC
    est <- ifelse(is.na(est.1), est.2, est.1)
  }

  ##The standard deviation of the estimate is calculated via the deltamethod and then used 
    ##to determine limits of the confidence interval
  seformula <- sprintf("~ log(exp(x4)*((x3-%f)/(%f-x2))^(1/x1))", lambda,lambda)
  se <- deltamethod(as.formula(seformula), coef(object), vcov(object))
  est.lo <- exp(log(est)- qt(1-alpha/2,df.residual(object))*se)
  est.up <- exp(log(est)+ qt(1-alpha/2,df.residual(object))*se)
  est.ci.delta <- as.matrix(cbind(est,est.lo,est.up))
  
  res <- data.frame(ALEC= est, CI.lower = est.lo, CI.upper = est.up, sd.ALEC=se, "log2(Lambda)" =lambda) 
  return(res)
}

###
##Preparations for calculating the LEC:
## the following functions are all needed at some point in the calculation of the LEC

####
##Gradient function - calculates the gradient of a 4pLL model (in which the parameter e is already log-transformed)

##Input:
##object -> a 4pLL model, fitted with LL2.4()
##x.conc -> the concentration for which the gradient is to be calculated

##Output:
##grad, the Gradient (ie the partial derivatives for each of the parameters - e equals log(e)!)

grad.4pLL <- function(object, x.conc){
  
  ##extract the coefficients from the drc object
  b <- coef(object)['b:(Intercept)']
  c <- coef(object)['c:(Intercept)']
  d <- coef(object)['d:(Intercept)']
  e <- coef(object)['e:(Intercept)']
  
  ##if the gradient for concentration 0 is to be calculated, the gradient is a vector with three entries =0 and the entry
  ##corresponding to c or d (depending on the parameter b) =1
  if(x.conc == 0){
    if(b<0){
      ##then c corresponds to the left asymptote
      return(c(0,1,0,0))
    }
    else{
      ##then d corresponds to the left asymptote
      return(c(0,0,1,0))
    }
  }
  else{
    ##these terms can be derived theoretically by applying simple differentiation rules
    g1 <- -((d-c)*(log(x.conc)-e)*exp(b*(log(x.conc)-e)))/((1+exp(b*(log(x.conc)-e)))^2)
    g2 <- 1-(1/(1+exp(b*(log(x.conc)-e))))
    g3 <- (1/(1+exp(b*(log(x.conc)-e))))
    g4 <- (b*(d-c)*exp(b*(log(x.conc)-e)))/((1+exp(b*(log(x.conc)-e)))^2)
    
    grad <- c(g1, g2, g3, g4)
    names(grad) <- c("b", "c", "d", "e")
    return(grad)
  }
}


###
##Theoretical variance - calculates the variance of each single difference based on the delta-method theory

##Input:
##object  -> fitted 4pLL object
##x1, x2  -> two doses, the variance between the difference of their respective response values is to be calculated

##Output:
##The variance of the difference between the response values for the two respective doses. x2 defaults to the control

theor.var <- function(object, x1, x2=0){
  
  ##extract the covariance matrix and calculate the gradients of the object for the two concentrations
  sigma <- vcov(object) 
  g1 <- grad.4pLL(object, x1)
  g2 <- grad.4pLL(object, x2)
  
  var.diff <- g1%*%sigma%*%g1 + g2%*%sigma%*%g2 - 2*g1%*%sigma%*%g2
  return(var.diff)
}


###
##model-test function, calculates p-values for the new 4pLL model based test for a given model and concentration

##Input:
##model   -> fitted 4pLL object (fitted with LL2.4())
##x.conc  -> concentration for which the significant different by a lambda from the lower asymptote is to be calculated
##lambda  -> difference between concentration and asymptote that needs to be exceeded

##Output:
##log2(f(x))      -> response value for concentration x.conc
##CI.L            -> lower limit of the confidence interval
##CI.R            -> upper limit of the confidence interval
##sd              -> standard deviation of (fx - f0)
##test-stat-value -> value of the teststatistic
##P.Value         -> p-value of the test

model.test <- function(model, x.conc, lambda, alpha=0.05, alt="greater"){
  library("msm")
  
  ##extract the coefficients from the model, e is actually log(e)
  bb <- coef(model)['b:(Intercept)']
  cc <- coef(model)['c:(Intercept)']
  dd <- coef(model)['d:(Intercept)']
  ee <- coef(model)['e:(Intercept)']
  
  ##returns the value of a 4pLL function with e being log(e)
  llog <- function(x, b, c, d, e){
    c + (d - c) / (1 + exp(b * (log(x) - e)))
  }
  
  ##response of the concentration in question
  fx <- llog(x = x.conc, bb, cc, dd, ee)
  
  ##value of the left asymptote:
  f0 <- ifelse(bb>0, dd, cc)
  
  ##variance of (fx - f0)
  var.FC <- theor.var(model, x1=x.conc)
  
  if(alt == "greater"){
    est.lo.FC <- fx - qnorm(1-alpha)*sqrt(var.FC)
    est.up.FC <- Inf
    est.ci.delta.FC <- as.matrix(cbind(fx, est.lo.FC, est.up.FC, var.FC))
    
    ##test-statistic and p-value based on the normal distribution
    test.stat <- ((fx-(f0+lambda))/sqrt(var.FC))
    pvalue <- 1-pnorm(test.stat)
  }
  if(alt == "less"){
    est.lo.FC <- -Inf
    est.up.FC <- fx + qnorm(1-alpha)*sqrt(var.FC)
    est.ci.delta.FC <- as.matrix(cbind(fx, est.lo.FC, est.up.FC, var.FC))
    
    ##test-statistic and p-value based on the normal distribution
    test.stat <- ((fx-(f0-lambda))/sqrt(var.FC))
    pvalue <- pnorm(test.stat)
  }
  if(alt == "two.sided"){
    est.lo.FC <- fx - qnorm(1-alpha/2)*sqrt(var.FC)
    est.up.FC <- fx + qnorm(1-alpha/2)*sqrt(var.FC)
    est.ci.delta.FC <- as.matrix(cbind(fx, est.lo.FC, est.up.FC, var.FC))
    
    ##test.statistic and p-value based on the normal distribution
    test.stat.g <- ((fx-(f0+lambda))/sqrt(var.FC))
    test.stat.l <- ((fx-(f0-lambda))/sqrt(var.FC))
    p.g <- 1-pnorm(test.stat.g)
    p.l <- pnorm(test.stat.l)
    pvalue <- 2*min(p.g, p.l)
    test.stat <- ifelse(p.g<p.l, test.stat.g, test.stat.l)
  }
  
  res <- data.frame("log2(f(x))"= fx, CI.lower = est.lo.FC, CI.upper = est.up.FC, "sd(f(x))" = sqrt(var.FC), "test.stat" = test.stat, p = pvalue)
  colnames(res) <- c("logFC", "CI.L", "CI.R", "sd", "test-stat.value", "P.Value")
  return(res)
}



###
##LEC: calculation of the LEC using grid search and the functions implemented above

##Input:
##model   -> a fitted drc object (fitted with LL2.4())
##lambda  -> the response value for which the corresponding concentration (ALEC) is to be calculated
##start   -> lower limit of the first interval used in the bisection method
##finish  -> upper limit of the first interval used in the bisection method
##epsilon -> length of interval in the bisection limit that leads to a stop of the method
##alpha   -> The niveau for the confidence interval (on the dosis axis)
##alt     -> specifying the alternative hypothesis (defaults to "greater")

##Output:
##Dataframe with the following entries:
##erg       -> dataframe with the results from the model.test function
##disk.conc -> LEC
##Iteration -> "Max.Iter" if no LEC could be calculated, "Reached" else

##Calculation of the LEC
LEC <- function(model, lambda=log2(1.5), start=1, finish=1000, epsilon=0.001, alpha=0.05, alt="greater"){

  ##if the result is not significant even for the highest concentration tested, stop
  disk.conc <- finish
  erg <- model.test(model, disk.conc, lambda=lambda, alpha=alpha, alt=alt)
  if(erg$P.Value>alpha | is.na(erg$P.Value)){
    res <- data.frame(erg, disk.conc, "Iteration" ="Max.Iter")
    return(res) 
    break
  }
  
  ##otherwise execute the bisection algorithm as presented in the paper
  lower <- start
  upper <- finish
  repeat{
    disk.conc <- (lower + upper)/2
    erg <- model.test(model, disk.conc, lambda = lambda, alpha=alpha, alt=alt)
    if(is.na(erg$P.Value)){
      res <- data.frame(erg, disk.conc, "Iteration" ="Max.Iter")
      return(res) 
      break
    }
    else{
      if(erg$P.Value<alpha){
        upper <- disk.conc
      }
      else{
        lower <- disk.conc
      }
    }
    if(abs(upper- lower)<epsilon){
      res <- data.frame(erg, disk.conc, "Iteration" ="Reached")
      return(res) 
      break
    }
  }
}





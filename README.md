
# Comparison of observation-based andmodel-based identification of alert concentrationsfrom concentration-expression data

An important goal of concentration-response studies in toxicology is to determine an 'alert' concentration where a critical level of the response variable is exceeded. 
In a classical observation-based approach, only measured concentrations are considered as potential alert concentrations. 
Alternatively, a parametric curve is fitted to the data that describes the relationship between concentration and response.
For a prespecified effect level, both an absolute estimate of the alert concentration and an estimate of the lowest concentration where the effect level is exceeded significantly are of interest. 

In a simulation study for gene expression data, we compared the observation-based and the model-based approach for both absolute and significant exceedance of the prespecified effect level.
Results show that, compared to the observation-based approach, the model-based approach overestimates the true alert concentration less often and more frequently leads to a valid estimate, especially for genes with large variance.

In this GitHub repository, the code needed to reproduce the simulation studies as well as the resulting datasets is stored.

## Simulation-of-Datasets.R

The code to simulate 9 concentration-response datasets.
These are divided into three scenarios introduced in detail in the paper.
For each scenario, three different values of the standard deviation used for sampling response data are used.
Resulting datasets are:

- DataSet.ScenarioI.VarSmall.RData
- DataSet.ScenarioI.VarMiddle.RData
- DataSet.ScenarioI.VarLarge.RData
- DataSet.ScenarioII.VarSmall.RData
- DataSet.ScenarioII.VarMiddle.RData
- DataSet.ScenarioII.VarLarge.RData
- DataSet.ScenarioIII.VarSmall.RData
- DataSet.ScenarioIII.VarMiddle.RData
- DataSet.ScenarioIII.VarLarge.RData

## Functions-Calculating-Alert-Concentrations.R

This document comprises all code needed for calculating the four different alert concentrations introduced in the paper, whereby the LOEC can be calculated based on the $t$-test and based on the Dunnett procedure.
ALOEC and the $t$-test LOEC are calculated with the help of the function "t.tests.loec" , the Dunnett procedure LOEC is calculated based on the function "dunnett.tests.loec".

The ALEC and its confidence interval is directly calculated in the function "ALEC".
For calculation of the LEC, several auxiliary functions are needed. 
Those are based on theoretical results for calculating the gradient of a 4pLL model or the Variance of the difference of two function values. 
Making use of the auxiliary functions "grad.4pLL"", "theor.var", and "model.test", the LEC is calculated via the bisection method in the function "LEC".

## Analysis-Calculating-Alert-Concentrations.R

Here, all alert concentrations are calculated, making use of the code provided in "Functions-Calculating-Alert-Concentrations.R". 
Results are stored according to Variance, i.e. in one of the following files, the respective 5 (ALOEC, 2x LOEC, ALEC, LEC) alert concentrations for all three scenarios are stored.
The resulting RData-Sets are:

- Benchmark.VarSmall.TwoSided.RData
- Benchmark.VarMiddle.TwoSided.RData
- Benchmark.VarLarge.TwoSided.RData

Additionally, for each 4pLL model fitted, it is determined whether the covariance matrix contains negative diagonal entries.
The indices of the simulated genes for which this is the case are stored separately in the three RData-Sets:

- NegativeCovarianceEntries-VarSmall.RData
- NegativeCovarianceEntries-VarMiddle.RData
- NegativeCovarianceEntries-VarLarge.RData

Due to the relatively long runtime of the glht-function used for calculating the Dunnett procedure, this analysis needs a couple of minutes runtime.

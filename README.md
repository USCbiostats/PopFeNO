# PopFeNO
Functions, datas and demo file for package PopFeNO (temp)

This repository stores demo document, R code, simulation results for the package in develepment temporary named PopFeNO.This package is for the following article:

PopFeNO provides organized R functions and walkthrough manual for researchers to choose their desired model to analysis FeNO data. 

### Installation

This package requires installation of JAGS. If not already installed, follow instructions at: https://mcmc-jags.sourceforge.io/
The following R packages were also required for data reformating, model fitting and implement the JAGS in R:

```{r}
require(MASS)
require(lme4)
require(nlme)
require(reshape2)
require(R2jags)
```

### Demostration document walk through

The demo document present how to  to implement six methods for research involving multiple flow exhaled nitric oxide (FeNO) measurements in a study population, where the goal is to relate estimated NO parameters to factors of interest (i.e., covariate(s) X).

In the demo file ("DeMO.pdf"), a simulated FeNO data set was first generated using given parameters (i.e. Population level NO parameters, covariate coefficient, correlation matrix). The NO parameters used here were CANO, log(CawNO), log(DawNO). We used natural log transformation of the later two parameters. 

This is a balanced dataset which means every participants has the same number of maneuvers. Alternative Jags model functions were required for unbalanced data.

The sixed methods presented here are four two-stage approaches (Two stage NLS, HMA, NLME, HB) and two unified approaches (Unified NLME, HB). 

* NLS: Natural log transform-both-sides nonlinear least squares model.
* HMA: Högman and Merilӓinen Algorithm, using an iterative algorithm involving a third-order approximation.
* NLME:Nonlinear mixed effects (NLME) approach, 
* HB:  Hierarchical Bayesian model 

HB methods requires hours to achieve a converged results so we also provide previously fitted results in the repository (use random seed to control the repeatablity).

The summary plot re-transformed the log(CawNO) and log(DawNO) for a more meaningful interpretation.

### Application of this package (in future)

We want to present a direct comparison of the popular methods for FeNO data analysis. Researchers can view the performance and decide which to use based on their own data's characteristivs.

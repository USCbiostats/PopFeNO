# PopFeNO
This repository stores demonstration files, R code, and simulation results related to the following articles:

Weng J, Molshatzki N, Marjoram P, Gauderman WJ, Gilliland FD, Eckel SP. Hierarchical Bayesian models to estimate covariate effects on airway and alveolar nitric oxide. Scientific Reports. 2021; 11(1):17180. PMID:34433846. PMC8387480.

Weng, J., Molshatzki, N., Marjoram, P., Gauderman, W., Gilliland, F. and Eckel, S.P., 2022. Longitudinal Hierarchical Bayesian models of covariate effects on airway and alveolar nitric oxide. medRxiv.

### Installation

This code requires installation of JAGS. If not already installed, follow instructions at: https://mcmc-jags.sourceforge.io/

The following R packages were also required for data reformating, model fitting and implementing JAGS in R:

```{r}
require(MASS)
require(lme4)
require(nlme)
require(reshape2)
require(R2jags)
```

### Demonstration

The demo document shows how to implement six methods for research involving multiple flow exhaled nitric oxide (FeNO) measurements in a study population, where the goal is to relate estimated NO parameters to factors of interest (e.g., environmental exposures, disease status) henceforth called covariate(s) X.

In the demo file ("DeMO.pdf"), a simulated FeNO data set was first generated using given parameters (i.e. Population level NO parameters, covariate coefficient, correlation matrix). The NO parameters used here were CANO, log(CawNO), and log(DawNO). We used natural log transformation of the latter two parameters. 

This is a balanced dataset which means every participant has the same number of maneuvers. Alternative Jags model functions were required for unbalanced data.

The six methods presented here are four two-stage approaches (Two stage NLS, HMA, NLME, HB) and two unified approaches (Unified NLME, HB). 

* NLS: Natural log transform-both-sides nonlinear least squares model.
* HMA: Högman and Merilӓinen Algorithm, using an iterative algorithm involving a third-order approximation.
* NLME:Nonlinear mixed effects (NLME) approach, 
* HB:  Hierarchical Bayesian model 

HB methods require hours to achieve a converged results so we also provide previously fitted results in the repository (use random seed to control the repeatablity).

The summary plot re-transformed the log(CawNO) and log(DawNO) for a more meaningful interpretation.

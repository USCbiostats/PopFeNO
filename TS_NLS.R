library(lme4)
library(nlme)

# nonlinear least square function for FeNO two compartment model
# fit n nls models for n participants
nlsListStart <- function(dat){
  out <- nlsList(logeno ~  log(exp(logCaw) + (Ca- exp(logCaw))*exp(-exp(logDaw)/flow))|id,# grouped by id
                 data=dat, start=list(Ca = 2, logCaw = log(68), logDaw = log(12)), # good starting values can improve convergence
                 control=list(tolerance=0.001))
  out
} 
# obtain the estimated NO parameters for N participants, ordered by id. 
TS_NLS_StageI<-function(dat){
  fitnls <- nlsListStart(dat)# fit NLS function
  estCoef <- coef(fitnls)# obtain estimated NO parameters
  estCoef$id<-as.numeric(rownames(estCoef))# add id column for later grouped regression
  estCoef_o<-estCoef[order(estCoef$id,decreasing=FALSE),]# sort the data by id
  return(estCoef_o)
}


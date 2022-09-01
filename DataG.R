
library(MASS) # needed for mvrnorm()
library(reshape2) # needed for melt()


#----------------- Unbalanced, noised flow, FeNo data -------------------#
# Use given flow or simulated flow, which could be randomly mutated unbalanced and noised



FlowFun<-function(Ndat,P){
  flow=matrix((rep(c(30,30,50,50,100,100,300,300),times=Ndat)+rnorm(Ndat*8,0,1))*rbinom(Ndat*8, 1, 0.9),nrow=Ndat,ncol=8,byrow = T)
  return(flow)
}
Xfn=function(n){rnorm(n,0,1)}

DataGeneratorCS<-function(Ndat      = 10,            # number of study participants
                          popmean   = c(1.5,3.5,2.5),# vector of NO parameter population means, when X=0, in this order: CA, logCaw, logDaw
                          X         = rnorm(Ndat,0,1),# function to generate participant-level covariate X as function of number of study participants
                          beta      = c(0.1,0.1,0.1),
                          flowTarget=c(rep(30,2),rep(50,2),rep(100,2),rep(300,2)),# vector of regression coefficients relating X to each NO parameter, in this order: CA, logCaw, logDaw
                          Flow=FlowFun(10,0.9), # Must be a Ndat*max length of flow matrix, missing values were 0s
                          SD=0.1,              #SD of the error
                          sdalphaCa            = 0.45, # population SD of CA
                          sdalphalogCaw        = 0.65, # population SD of logCaw
                          sdalphalogDaw        = 0.55, # population SD of logDaw
                          coralphalogCawCa     = 0.66, # correlation of NO parameters: logCaw, CA
                          coralphalogCawlogDaw = -0.35,# correlation of NO parameters: logCaw, logDaw
                          coralphalogDawCa     = -0.38 # correlation of NO parameters: logDaw, CA
){
  # generate X
  print(X)

  # create variance/covariance matrix for NO parameters based on SD and correlations
  NOcov <- matrix(c(sdalphaCa^2, coralphalogCawCa*sdalphaCa*sdalphalogCaw, coralphalogDawCa*sdalphaCa*sdalphalogDaw,          
                    coralphalogCawCa*sdalphaCa*sdalphalogCaw, sdalphalogCaw^2, coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw,
                    coralphalogDawCa*sdalphaCa*sdalphalogDaw, coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw, sdalphalogDaw^2),
                  3,3)
  logeno <- matrix(NA,ncol=length(flowTarget),nrow=Ndat) # create empty 2-D array for log(FeNO), each row for one participant
  colnames(logeno) <- flowTarget
  NOparam <- matrix(NA,ncol=3,nrow=Ndat)           # create empty 2-D array for NO parameters, each row for one participant
  colnames(NOparam) <- c("Ca","logCaw","logDaw")
  
  for(i in 1:Ndat){
    repeat{
      param <- mvrnorm(1, X[i] *  beta + popmean, NOcov) # sample NO parameters (Ca, logCaw, logDaw) via multi-normal distribution for participant i
      FeNOi  <- exp(param[2])+(param[1]-exp(param[2]))*exp(-1*exp(param[3])/Flow[i,]) # calculate mean FeNO for participant i using Ca[i], Caw[i], Daw[i]
      logFeNOobsi  <- ifelse(FeNOi>0,log(FeNOi) + rnorm(length(Flow[i,]),0,SD),NA) # calculate log(FeNO) for each measurement using mean FeNO if mean FeNO is non-negative
      if(!any(is.na(logFeNOobsi)) & param[1]>0){break}#satisfy all constrains: Ca[i] is non negative (param[1]) and FeNO[i,j] is non-negative
    }
    logeno[i,]  <- logFeNOobsi # store sampled logFeNO for participant i's multi-flow measurement
    NOparam[i,] <- c(param[1],param[2],param[3]) #store sampled NO parameter for participant i
  }
  # create alternative format dataset for frequentist methods
  logenodata<-melt(data.frame(id=c(1:Ndat),logeno),id="id")
  logenodata<-logenodata[order(logenodata$id,decreasing = FALSE),]
  dat <- data.frame(id=rep(c(1:Ndat),each=length(flowTarget)),eno=NA,logeno=NA,flow=as.vector(t(Flow)))
  dat$logeno  <- logenodata$value
  dat$eno     <-exp(dat$logeno)
  dat$flowTarget <- rep(flowTarget,times=Ndat)
 
  dat<-dat[(dat$flow != 0)>0,]
  # "real" data will have observed flow different from target flow rates
  return(list(dat=dat,datNOparam=NOparam,datX=data.frame(id=1:Ndat,X=X)))
}


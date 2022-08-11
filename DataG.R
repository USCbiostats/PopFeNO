
library(MASS) # needed for mvrnorm()
library(reshape2) # needed for melt()

# Generate Cross-Sectional Data
DataGeneratorCS<-function(Ndat=500,            # number of study participants
                          alpha=c(1.5,3.5,2.5),# vector of NO parameter population means, when X=0, in this order: CA, logCaw, logDaw
                          Xfn=function(n){rnorm(n,0,1)},# function to generate participant-level covariate X as function of number of study participants
                          beta=c(0.1,0.1,0.1), # vector of regression coefficients relating X to each NO parameter, in this order: CA, logCaw, logDaw
                          Flow=c(rep(30,2),rep(50,2),rep(100,2),rep(300,2)), # vector of flow rates in multiple flow FeNO protocol
                          SD=0.1,              #SD of the error
                          sdalphaCa            = 0.45, # population SD of CA
                          sdalphalogCaw        = 0.65, # population SD of logCaw
                          sdalphalogDaw        = 0.55, # population SD of logDaw
                          coralphalogCawCa     = 0.66, # correlation of NO parameters: logCaw, CA
                          coralphalogCawlogDaw = -0.35,# correlation of NO parameters: logCaw, logDaw
                          coralphalogDawCa     = -0.38 # correlation of NO parameters: logDaw, CA
                          ){
  # generate X
  X <- Xfn(Ndat)
  # create variance/covariance matrix for NO parameters based on SD and correlations
  NOcov <- matrix(c(sdalphaCa^2, coralphalogCawCa*sdalphaCa*sdalphalogCaw, coralphalogDawCa*sdalphaCa*sdalphalogDaw,          
                    coralphalogCawCa*sdalphaCa*sdalphalogCaw, sdalphalogCaw^2, coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw,
                    coralphalogDawCa*sdalphaCa*sdalphalogDaw, coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw, sdalphalogDaw^2),
                  3,3)
  logeno <- matrix(NA,ncol=length(Flow),nrow=Ndat) # create empty 2-D array for log(FeNO), each row for one participant
  colnames(logeno) <- Flow
  NOparam <- matrix(NA,ncol=3,nrow=Ndat)           # create empty 2-D array for NO parameters, each row for one participant
  colnames(NOparam) <- c("Ca","logCaw","logDaw")
  
  for(i in 1:Ndat){
    repeat{
      param <- mvrnorm(1, X[i] *  beta + alpha, NOcov) # sample NO parameters (Ca, logCaw, logDaw) via multi-normal distribution for participant i
      FeNOi  <- exp(param[2])+(param[1]-exp(param[2]))*exp(-1*exp(param[3])/Flow) # calculate mean FeNO for participant i using Ca[i], Caw[i], Daw[i]
      logFeNOobsi  <- ifelse(FeNOi>0,log(FeNOi) + rnorm(length(Flow),0,SD),NA) # calculate log(FeNO) for each measurement using mean FeNO if mean FeNO is non-negative
      if(!any(is.na(logFeNOobsi)) & param[1]>0){break}#satisfy all constrains: Ca[i] is non negative (param[1]) and FeNO[i,j] is non-negative
    }
    logeno[i,]  <- logFeNOobsi # store sampled logFeNO for participant i's multi-flow measurement
    NOparam[i,] <- c(param[1],param[2],param[3]) #store sampled NO parameter for participant i
  }
  # final dataset for Bayesian methods in JAGS, required in list format
  datJAGS <- list(NOparam=NOparam,logeno=logeno,X=X)
  
  # create alternative format dataset for frequentist methods
  logenodata<-melt(data.frame(id=c(1:Ndat),logeno),id="id")
  logenodata<-logenodata[order(logenodata$id,decreasing = FALSE),]
  dat <- data.frame(id=rep(c(1:Ndat),each=length(Flow)),eno=NA,logeno=NA,flow=rep(Flow,Ndat))
  dat$logeno  <- logenodata$value
  dat$eno     <-exp(dat$logeno)
  dat$flowTarget <- dat$flow # "real" data will have observed flow different from target flow rates
  return(list(dat=dat,datJAGS=datJAGS,datX=data.frame(id=1:Ndat,X=X)))
}


# Generate Longitudinal Data
DataGeneratorL<-function(beta0,#vector of NO parameter population means, when X=0, in this order: CA, logCaw, logDaw
                         beta1,#association with environmental factor
                         Xfn=function(n,m){matrix(rnorm(n*m,0,1),ncol=n,nrow=m)}, # function to generate participant-visit level covariate X as function of number of study participants and visit times
                         # The covariates are independent and identically distributed normal variables
                         Flow,# vector of flow rates in multiple flow FeNO protocol
                         SD,  #SD of the error
                         sdCa            = 0.50,# personal SD of CA (visit level)
                         sdlogCaw        = 0.40,# personal SD of logCaw
                         sdlogDaw        = 0.33,# personal SD of logDaw
                         corlogCawCa     = 0,# personal correlation of NO parameters: logCaw, CA
                         corlogCawlogDaw = 0,# personal correlation of NO parameters: logCaw, logDaw
                         corlogDawCa     = 0,# personal correlation of NO parameters: logDaw, CA
                         sdalphaCa            = 0.45,# population SD of CA (participant level)
                         sdalphalogCaw        = 0.65,# population SD of logCaw
                         sdalphalogDaw        = 0.55,# population SD of logDaw
                         coralphalogCawCa     = 0.66,# population correlation of NO parameters: logCaw, CA
                         coralphalogCawlogDaw = -0.35,# population correlation of NO parameters: logCaw, logDaw
                         coralphalogDawCa     = -0.38,# population correlation of NO parameters: logDaw, CA
                         Ndat,Visit){# number of study participants and visit times
  # generate X
  X <- Xfn(Ndat,Visit)
  # create visit level variance/covariance matrix for NO parameters based on SD and correlations
  Visit_NOcov<- matrix(c(sdCa^2,corlogCawCa*sdCa*sdlogCaw,corlogDawCa*sdCa*sdlogDaw,
                                  corlogCawCa*sdCa*sdlogCaw,sdlogCaw^2,corlogCawlogDaw*sdlogCaw*sdlogDaw,
                                  corlogDawCa*sdCa*sdlogDaw,corlogCawlogDaw*sdlogCaw*sdlogDaw,sdlogDaw^2),3,3)
  # create participant level variance/covariance matrix for NO parameters based on SD and correlations
  Sub_NOcov<- matrix(c(sdalphaCa^2,coralphalogCawCa*sdalphaCa*sdalphalogCaw,coralphalogDawCa*sdalphaCa*sdalphalogDaw,
                                coralphalogCawCa*sdalphaCa*sdalphalogCaw,sdalphalogCaw^2,coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw,
                                coralphalogDawCa*sdalphaCa*sdalphalogDaw,coralphalogCawlogDaw*sdalphalogCaw*sdalphalogDaw,sdalphalogDaw^2),3,3)
  # create empty 3-D array for log(FeNO)
  logeno <- array(NA,dim=c(length(Flow),Ndat,Visit),dimnames = list(Flow,c(1:Ndat),c(1:Visit)))
  # create empty 3-D array for NO parameters
  NOparam <-array(NA,dim=c(3,Ndat,Visit),dimnames=list(c("Ca","logCaw","logDaw"),c(1:Ndat),c(1:Visit)))
  # create empty 3-D array for average NO parameters for each participant at their visit
  mean<-array(NA,dim=c(3,Ndat,Visit),dimnames=list(c("Ca","logCaw","logDaw"),c(1:Ndat),c(1:Visit)))
  # create array to store the randomness for participants
  alpha<-array(NA,dim=c(3,Ndat),dimnames=list(c("Ca","logCaw","logDaw"),c(1:Ndat)))
 
  for(i in 1:Ndat){
    repeat{
      alpha[,i]<-mvrnorm(1,c(0,0,0),Sub_NOcov) #sample the randomness for participant i
      #create temporary array to store the average NO parameters for participant i
      mean_sub<-array(NA,dim=c(3,Visit),dimnames=list(c("Ca","logCaw","logDaw"),c(1:Visit))) 
      #create temporary arrays to store NO parameters and logFeNO for participant i at visit j
      NO_sub<-array(NA,dim=c(3,Visit),dimnames=list(c("Ca","logCaw","logDaw"),c(1:Visit))) 
      logFeNOobsi <- array(NA,dim=c(length(Flow),Visit),dimnames = list(Flow,c(1:Visit)))
      for(j in 1:Visit ){
        mean_sub[,j]<-beta0+alpha[,i]+X[j,i]*beta1#X:row=visit/column=participant
        NO_sub[,j] <- mvrnorm(1,mean_sub[,j],Visit_NOcov) #add the visit level randomness
        FeNOi  <- exp(NO_sub[2,j])+(NO_sub[1,j]-exp(NO_sub[2,j]))*exp(-1*exp(NO_sub[3,j])/Flow)
        logFeNOobsi[,j]  <- ifelse(FeNOi>0,log(FeNOi) + rnorm(length(Flow),0,SD),NA)
      }
      if((!any(is.na(logFeNOobsi))) & sum(NO_sub[1,]>0)==Visit){# check if all the sampled Ca and FeNO are non-negative for participant i
        break
      }#satisfy all, continue for next participant
    }
    
    logeno[,i,]  <- logFeNOobsi
    NOparam[,i,] <-NO_sub
    mean[,i,]<-mean_sub
    
  }
  # final dataset for Bayesian methods in JAGS, required in list format
  dataJAGS<-list(NOparam=NOparam,logeno=logeno,X=X)
  
  # ---------create alternative format dataset for frequentist methods--------#
  # convert X matrix to long formated, repeated Nflow times for each visit
  convertXflow<-function(Xmatrix,Nflow){# X has Ndat columns and Nvisit rows
    out<-rep(Xmatrix[,1],each=Nflow)
    for(i in 2:Ndat){
      out<-c(out,rep(Xmatrix[,i],each=Nflow))
    }
    return(out)
  }
  X_long                <-convertXflow(X,Nflow)
  
  # reformate the 3-D array to 2-D array
  logenodata   <-melt(data.frame(id=c(1:Ndat),t(logeno[,,1])),id="id")
  for(i in 2:Visit){
    logenodata<-cbind(logenodata,melt(data.frame(id=c(1:Ndat),t(logeno[,,i])),id="id"))
  }
  
  colnames(logenodata)  <-c("id","flow","logeno","visit")
  logenodata$flow       <-rep(c(30,50,100,300),each=Ndat*2,times=3)
  logenodata            <-logenodata[order(logenodata$id,logenodata$visit,decreasing = FALSE),]
  logenodata$eno        <-exp(logenodata$logeno)
  logenodata$idvisit    <-paste(logenodata$id,logenodata$visit,sep="_")
  dat                   <-logenodata
  #-----------------#
  return(list(dat=dat,datJAGS=datJAGS,datX=data.frame(id=dat$id,visit=dat$visit,X=X_long)))
}



#----------------- Unbalanced, noised flow, FeNo data -------------------#
#Use given flow or randomly mutated unbalanced, noised flow

# assume 10% of the flow were 'missing'
# generate random vector as the indicator whether the flow were measured
# Multiply the indicator with the flow matrix


FlowFun<-function(Ndat,P){
  flow=matrix((rep(c(30,30,50,50,100,100,300,300),times=Ndat)+rnorm(Ndat*8,0,1))*rbinom(Ndat*8, 1, 0.9),nrow=Ndat,ncol=8,byrow = T)
  return(flow)
}

DataGeneratorCS_UbN<-function(Ndat=10,            # number of study participants
                          alpha=c(1.5,3.5,2.5),# vector of NO parameter population means, when X=0, in this order: CA, logCaw, logDaw
                          Xfn=function(n){rnorm(n,0,1)},# function to generate participant-level covariate X as function of number of study participants
                          beta=c(0.1,0.1,0.1),
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
  X <- Xfn(Ndat)

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
      param <- mvrnorm(1, X[i] *  beta + alpha, NOcov) # sample NO parameters (Ca, logCaw, logDaw) via multi-normal distribution for participant i
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


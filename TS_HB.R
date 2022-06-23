# For balanced data: n participants each has same number of maneuvers (same multiple flow rate)
# Jage model function
# Ca should be non-negative, only single distribution can add the truncation constrain. 
# Thus, we separate the multi-normal distribution of (Ca,logCaw,logDaw) to:
# conditional normal distribution: Ca|(logCaw, logDaw) and
# bivariate normal distribution: (logCaw, logDaw)
TSHB_JagsModel <- function()  {
  for( i in 1 : Ndat) {# iterate all participants
    for( j in 1 : Nflows) {# iterated all flows (maneuvers)
      logeno[i,j] ~ dnorm(logmu[i,j],tau_c)# add measurement error
      mu[i,j] <- exp(logCaw[i])+(Ca[i]-exp(logCaw[i]))*exp(-1*exp(logDaw[i])/flow[j])# mean estimated FeNO
      logmu[i,j]<-log(max(mu[i,j] ,0.001))#if mean FeNO is non-negative, log transform it, otherwise set the value to 0.001
    }
    # store NO parameters
    Ca[i]       <- param_subject[i, 1]
    logCaw[i]   <- param_subject[i, 2] 
    logDaw[i]   <- param_subject[i, 3] 
    #------------sample bivariate normal (logCaw, logDaw) and zero-truncated normal Ca|(logCaw, logDaw)
    # add randomness of participant i to NO participants
    # use truncated normal distribution for Ca
    param_subject[i,1] ~ dnorm(mu1_22[i],tau1_22);T(0,)
    # calculated the conditional mean for Ca|(logCaw, logDaw)
    mu1_22[i]<-beta0_Ca+part1[1,1]*(param_subject[i,2]-beta0_logCaw)+part1[1,2]*(param_subject[i,3]-beta0_logDaw)
    # sample logCaw and logDaw
    param_subject[i,2:3] ~ dmnorm(c(beta0_logCaw,beta0_logDaw),tau22[1:2,1:2])
    #--------------#
  }
  #------- percision of Ca|(logCaw, logDaw)--------#
  tau1_22<-1/cov1_22
  # variance of Ca Conditioned on (logCaw,logDaw)
  cov1_22[1,1]<-varCa-(covCalogCaw^2*tau22[1,1]+2*covCalogDaw*tau22[2,1]*covCalogCaw+covCalogDaw^2*tau22[2,2])
  # Sigm(Ca, logCaw, logDaw) * tau(logCaw, logDaw)
  part1[1,1]<-covCalogCaw*tau22[1,1]+covCalogDaw*tau22[2,1]
  part1[1,2]<-covCalogCaw*tau22[1,2]+covCalogDaw*tau22[2,2]
  # Sigma (logCaw, logDaw)
  tau22[1:2,1:2]<-inverse(cov22[1:2,1:2])
  #--------------#
  # calculate covariances
  cov22[1,1] <- varlogCaw
  cov22[1,2] <- covlogCawlogDaw
  cov22[2,1] <- covlogCawlogDaw
  cov22[2,2] <- varlogDaw
  
  covCalogCaw       <-corlogCawCa*sqrt(varCa*varlogCaw)
  covCalogDaw       <-corlogDawCa*sqrt(varCa*varlogDaw)
  covlogCawlogDaw   <-corlogCawlogDaw*sqrt(varlogCaw*varlogDaw)
  # calculate variances and SDs
  varCa<-1/tauCa
  varlogCaw<-1/taulogCaw
  varlogDaw<-1/taulogDaw
  sdCa     <- sqrt(varCa)
  sdlogCaw <- sqrt(varlogCaw)
  sdlogDaw <- sqrt(varlogDaw)
  sigma_c<-sqrt(1/tau_c)
  # sample the third correlation, conditioned on the previous two, to assure positive definite
  corlogDawCa   ~ dunif(L_corlogDawCa,U_corlogDawCa)
  L_corlogDawCa <- corlogCawCa*corlogCawlogDaw-sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  U_corlogDawCa <- corlogCawCa*corlogCawlogDaw+sqrt((corlogCawlogDaw^2-1)*(corlogCawCa^2-1))
  
  corlogCawCa        ~ dunif(-1,1)
  corlogCawlogDaw    ~ dunif(-1,1)
  # sample precision
  tauCa ~ dgamma(1.0E-3, 1.0E-3)
  taulogCaw ~ dgamma(1.0E-3, 1.0E-3)
  taulogDaw ~ dgamma(1.0E-3, 1.0E-3)
  tau_c       ~ dgamma(1.0E-3, 1.0E-3)##random error in regression
  # sample population mean
  beta0_Ca     ~ dnorm(beta0_prior[1],1e-6)
  beta0_logCaw ~ dnorm(beta0_prior[2],1e-6)
  beta0_logDaw ~ dnorm(beta0_prior[3],1e-6)
  
} 

# Two stage HB function with default setting values
TSHB_iter<-function(beta0_prior=c(2,4,3),
                    rhat=1.1,addon.iter=6000,Max_update=40,
                    n.final=6000,N.iterT=6000,N.thinM=2,N.burnin=4000,N.chain=3,
                    flow=flow,dat=dat,
                    tracing=c("beta0_Ca","beta0_logCaw","beta0_logDaw",
                               "Ca","logCaw","logDaw",
                               "sdCa","sdlogCaw","sdlogDaw","corlogCawCa","corlogCawlogDaw","corlogDawCa",
                               "sigma_c")
                    ){
  # data needed for JAGS model
  data_TSHB <- list(Ndat=length(dat$X),Nflows=length(flow),flow=flow,logeno=dat$logeno,
                  beta0_prior=beta0_prior)
  # prepare initial values (number of list = number of parallel chains)
  inits_TSHB <- 
    list(
      list(
        tau_c=runif(1,0.05,400), 
        beta0_Ca=runif(1,0,5),
        beta0_logCaw=runif(1,0,5),
        beta0_logDaw=runif(1,0,5)
      ),
      list(
        tau_c=runif(1,0.05,400),
        beta0_Ca=runif(1,0,5),
        beta0_logCaw=runif(1,0,5),
        beta0_logDaw=runif(1,0,5)
      ),
      list(
        tau_c=runif(1,0.05,400),
        beta0_Ca=runif(1,0,5),
        beta0_logCaw=runif(1,0,5),
        beta0_logDaw=runif(1,0,5)
      )
    )
  
  # specify which parameters you want to track
  # Example: Ca for participant i has rowname Ca[i]
  param_TSHB <-c("beta0_Ca","beta0_logCaw","beta0_logDaw",
                        "Ca","logCaw","logDaw",
                        "sdCa","sdlogCaw","sdlogDaw","corlogCawCa","corlogCawlogDaw","corlogDawCa",
                        "sigma_c")
  # simulation start time
  time<-proc.time()[3]
  # Jags function
  fit_TSHB<- jags(data=data_TSHB, parameters.to.save=param_TSHB, inits=inits_TSHB,
                  model.file=TSHB_JagsModel, n.chains=N.chain, n.iter=N.iterT, 
                  n.burnin=N.burnin,n.thin=N.thinM) 
  # model converges when Rhat<1.1
  # record how much iterations and check points are met
  N_updateT=0
  N_final=0
  #---------total iterations: [adaptive phase + short update phase] + [long update phase * N_updateT(+ final phase * (N_final-1))] + [final phase]--------#
  while(N_updateT<Max_update){# we need a maximum iteration to finish the job
    # check point 1
    if(max(fit_TSHB$BUGSoutput$summary[rownames(fit_TSHB$BUGSoutput$summary) %in% tracing,"Rhat" ])<rhat){
      fit_TSHB<-update(fit_TSHB,n.iter=n.final)# if Rhat<1.1, go to the final phase to obtain posterior distributions
      # count how many times we entered the final phase
      N_final=N_final+1
      # check point 2: if met, finish simulation
      if(max(fit_TSHB$BUGSoutput$summary[rownames(fit_TSHB$BUGSoutput$summary) %in% tracing,"Rhat" ])<rhat){
        break
      }
    }# either check point 1 or 2 is not met, repeat update phase
    fit_TSHB<-update(fit_TSHB,n.iter=addon.iter)
    # count how many update phase we entered
    N_updateT<-N_updateT+1
  }
  
  # end time
  timespend<-proc.time()[3]-time
  # trimmed mean
  trimmean<-apply(fit_TSHB$BUGSoutput$sims.matrix,2,function(x) mean(x, 0.2))
  # 99% quantile
  quantile1_99<-t(apply(fit_TSHB$BUGSoutput$sims.matrix,2,function(x) quantile(x, c(0.005,0.995))))
  # output result list
  result<-list(
    timespend=timespend,
    summary=cbind(fit_TSHB$BUGSoutput$summary,"trimmean"=trimmean,quantile1_99),
    Iter_T=N.iterT+N_updateT*addon.iter+n.final*N_final,# total iterations after the main function
    Iter_track<-c(N.iterT=N.iterT,N_updateT=N_updateT,addon.iter=addon.iter,n.final=n.final,N_final=N_final)# total number of update phases and final phases
  )
  return(result)
}





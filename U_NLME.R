library(MASS)
library(lme4)
library(nlme)
library(reshape2)

nonLinInitPlogJ <-function(mCall,LHS,data,...){
  xy <- sortedXyData(mCall[["predictor"]],LHS,data)
  fit <- lm(xy[,"y"]~I(1/xy[,"x"]) + I(1/xy[,"x"]^2))# qualP function
  est <- fit$coef
  CaNO <- est[1]
  JawNO <- est[2] - 2*est[1]*est[3]/est[2]
  DawNO <- -2*est[3]/est[2]
  #log transform Caw, Daw, Jaw, the logged value set to be at least 0.0001
  logDawNO <- log(max(DawNO,0.0001))
  logJawNO <- log(max(JawNO,0.0001))
  logCawNO <- log(max(JawNO/DawNO,0.0001))
  value <- c(CaNO,logCawNO,logDawNO)
  names(value) <- mCall[c("CaNO","logCawNO","logDawNO")]
  value # output NO parameters
}
## quadT self-starter prelim function
nonLinInitTlogJ <-function(mCall,LHS,data,...){
  xy <- sortedXyData(mCall[["predictor"]],LHS,data)
  fit <- lm(I(xy[,"y"]*xy[,"x"])~I(xy[,"x"]) + I(1/xy[,"x"]))#qualP function
  est <- fit$coef
  CaNO <- est[2]
  JawNO <- est[1] - 2*est[2]*est[3]/est[1]
  DawNO <- -2*est[3]/est[1]
  #log transform Caw, Daw, Jaw, the logged value set to be at least 0.0001
  logDawNO <- log(max(DawNO,0.0001))
  logJawNO <- log(max(JawNO,0.0001))
  logCawNO <- log(max(JawNO/DawNO,0.0001))
  value <- c(CaNO,logCawNO,logDawNO)
  names(value) <- mCall[c("CaNO","logCawNO","logDawNO")]
  value # output NO parameters
}

## model without gradient
nonLinLogModellogJ <- function(predictor,CaNO,logCawNO,logDawNO){
  log(exp(logCawNO) + (CaNO - exp(logCawNO))*exp(-exp(logDawNO)/predictor))
}
# model with gradient
nonLinLogModellogJwithGrad <- deriv( ~log(exp(logCawNO) + (CaNO - exp(logCawNO))*exp(-exp(logDawNO)/predictor)),
                                     c("CaNO","logCawNO","DawNO"),
                                     function(predictor,CaNO,logCawNO,logDawNO){})

#self-starting function with gradient
SSnonLinLogPlogJGrad <- selfStart(nonLinLogModellogJwithGrad, nonLinInitPlogJ, c("CaNO","logCawNO","logDawNO"))
#self-starting function without gradient
SSnonLinLogPlogJwoGrad <- selfStart(nonLinLogModellogJ, nonLinInitPlogJ, c("CaNO","logCawNO","logDawNO"))


U_NLME_direct <- function(	dat,       # multiple flow FeNO dataset (id, eno, logeno, flow, flowTarget)
                            X, 	   # optional dataset containing id and person-level X
                            tol,   # tolerance  second fit (usually smaller)
                            outputFit=TRUE # if TRUE, outputs fitted model additionally
){
  
  mydat <- subset(merge(dat,X,by="id"),select=c(logeno,flow,id,X)) 
  dat2 <- groupedData(logeno~ flow | id, data=mydat)
  # nls separately for each group
  fitnlsList <- nlsList(logeno~SSnonLinLogPlogJwoGrad(flow,CaNO,logCawNO,logDawNO),data=dat2,
                        control=list(tolerance = 0.1))
  # starting values for fixed effects
  cf <- coef(fitnlsList)# this coefficients were not ordered by id
  
  startFix <- unlist(lapply(cf, median, na.rm = TRUE))
  idREstart <- matrix(0,nrow=length(unique(mydat$id)),ncol=3)
  colnames(idREstart) <- names(startFix); rownames(idREstart) <- unique(mydat$id)
  idmat <- matrix(as.integer(unlist(strsplit(rownames(coef(fitnlsList)),"/"))),byrow=T,ncol=1)
  for(i in 1:3){
    fillcol <- tapply(cf[,i],as.factor(idmat[,1]),mean,na.rm=T)# re-arrange estimated coefficients by id
    idREstart[which(!is.nan(fillcol)),i] <- fillcol[which(!is.nan(fillcol))]-startFix[i] #get the random effects for each id as coef-median(coef),fill NA with 0
  }
  
  
  
  fitU.nlme <- nlme(logeno~SSnonLinLogPlogJwoGrad(flow,CaNO,logCawNO,logDawNO),
                    data=mydat,
                    fixed = list(CaNO ~ X, logCawNO ~ X, logDawNO ~ X), # modity this part for multiple X, together with the starting values
                    random= pdLogChol(CaNO +logCawNO + logDawNO ~ 1), # specify the correlations for random effect
                    groups=~id,# group by id (participants)
                    # specify starting values of fixed and random effects
                    start=list( fixed=c(startFix[1],0,startFix[2],0,startFix[3],0),# startFix is the start value for intercepts, 0's are slope.
                                random=list(id=idREstart)      
                    ),
                    verbose=F,
                    control=list(tolerance = tol,opt="nlminb",natural=TRUE,niterEM=50)
  )
  return(fitU.nlme)
}

# If the direct approach can't converge, fit TS_NLME first and update with Xs.
U_NLME_update <- function( TSNLME_StageIout,
                           dat,       # multiple flow FeNO dataset (id, eno, logeno, flow, flowTarget)
                           X, 	   # optional dataset containing id and person-level X
                           tol=0.1,# tolerance for second fit (usually smaller)
                           outputFit=TRUE # if TRUE, outputs fitted model additionally
){
  startFix <- fixef(TSNLME_StageIout$fit)
  mydat <- subset(merge(dat,X,by="id"),select=c(logeno,flow,id,X)) 
  fitU.nlme <- update(TSNLME_StageIout$fit,
                      fixed = list(CaNO ~ X, logCawNO ~ X, logDawNO ~ X),
                      start = c(startFix[1],0,startFix[2],0,startFix[3],0),
                      control=list(tolerance=tol)
  )
  return(fitU.nlme)
}

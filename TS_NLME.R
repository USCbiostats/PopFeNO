library(MASS)
library(lme4)
library(nlme)
library(reshape2)

## quadP self-starter prelim function
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



# Stage I: obatin NO parameters via NLME function 
TS_NLME_StageI <- function(	dat,       # multiple flow FeNO dataset (id, eno, logeno, flow, flowTarget)
                            X=NA, 	   # optional dataset containing id and person-level X
                            tol1=0.1,  # tolerance for initial fit (usually larger)
                            tol2=0.01, # tolerance for second fit (usually smaller)
                            outputFit=TRUE # if TRUE, outputs fitted model additionally
){
  if(!is.null(dim(X))){
    dat<-merge(dat,X,by="id")
    mydat <- subset(dat,select=c(logeno,flow,id,X)) 
  }else{
    mydat <- subset(dat,select=c(logeno,flow,id))
  }
   
  
  dat2 <- groupedData(logeno~ flow | id, data=mydat)
  # nls separately for each group
  fitnlsList <- nlsList(logeno~SSnonLinLogPlogJwoGrad(flow,CaNO,logCawNO,logDawNO),data=dat2,
                        control=list(tolerance = 0.1))
  # starting values for fixed effects
  cf <- coef(fitnlsList)# this coefficients were not ordered by id
  
  startfix <- unlist(lapply(cf, median, na.rm = TRUE))
  # ----------starting values for id-level RE-----------# 
  # There are lot of NA in fitnls results
  # We use the medians as starting value for fixed components
  idREstart <- matrix(0,nrow=length(unique(mydat$id)),ncol=3)
  colnames(idREstart) <- names(startfix); rownames(idREstart) <- unique(mydat$id)
  idmat <- matrix(as.integer(unlist(strsplit(rownames(coef(fitnlsList)),"/"))),byrow=T,ncol=1)
  for(i in 1:3){
    fillcol <- tapply(cf[,i],as.factor(idmat[,1]),mean,na.rm=T)# re-arrange estimated coefficients by id
    idREstart[which(!is.nan(fillcol)),i] <- fillcol[which(!is.nan(fillcol))]-startfix[i] #get the random effects for each id as coef-median(coef),fill NA with 0
  }
  #----------------#
  
  # we use the self-starting function without gradient, can change to another if needed
  fit1.nlme <- nlme(logeno~SSnonLinLogPlogJwoGrad(flow,CaNO,logCawNO,logDawNO),
                    data=mydat,
                    fixed=CaNO +logCawNO + logDawNO ~ 1, # specify function for fixed effects
                    random= list(pdLogChol(CaNO +logCawNO + logDawNO ~ 1)), # specify the correlations for random effect
                    groups=~id,# group by id (participants)
                    # specify starting values of fixed and random effects
                    start=list( fixed=c(startfix["CaNO"],startfix["logCawNO"],startfix["logDawNO"]),
                                random=list(id=idREstart)      
                    ),
                    verbose=F,# don't output detailed simulation information
                    control=list(tolerance = tol1,opt="nlminb",natural=TRUE,niterEM=50))# use a larger tolerance value 1 to get an initial model
  # update model with smaller tolerance value 2
  fit2.nlme <- update(fit1.nlme,control=list(tolerance =tol2,maxIter=50),verbose=T)
  fit <- fit2.nlme
  if(outputFit){
    summary(fit2.nlme)
  }
  
  cf <- coef(fit)
  cf <- data.frame(cf,id=as.numeric(rownames(cf)))
  cf_o<-cf[order(cf$id,decreasing=FALSE),]# sort the data by id
  names(cf_o)[which(names(cf_o)=="CaNO")] <- "Ca"
  names(cf_o)[which(names(cf_o)=="logCawNO")] <- "logCaw"
  names(cf_o)[which(names(cf_o)=="logDawNO")] <- "logDaw"
  
  return(list(ests=cf_o,fit=fit2.nlme,dat=mydat))

}


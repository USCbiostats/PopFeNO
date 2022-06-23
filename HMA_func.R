Hogman3rdOrderAlg <- function(fL,fM,fH,eL,eM,eH){  
  # linT approximation on high flows to get S=CaNO and I
  CaNO_est <- S <- (fH*eH - fM*eM)/(fH-fM)
  # solve for I using linT model at medium flow rate (flow=100)
  I <- (fM*eM - S*fM)
  JawNO_est  <- I
  
  # first order approximation starting value
  Dw0_fo <- I/(eL - S)
  
  # iterative algorithm function, 1st order
  italg <- function(SV,nIterAlg=10){
    # object to store Dw values for each iteration
    outIterAlg <- numeric(nIterAlg)
    # set starting value (SV)
    Dw_n <- SV
    for(i in 1:nIterAlg){
      Dw_nplus1 <- outIterAlg[i] <-  -I*(exp(-Dw_n/fL) - exp(-Dw_n/fM))/(eL-eM)
      diffDw <- abs(Dw_n - Dw_nplus1)
      #print(i);print(Dw_nplus1)
      # stop if values of Dw_n and Dw_nplus1 are close enough
      if(diffDw <0.001){break}
      Dw_n <- Dw_nplus1  
    }
    out <- list(DawNO=Dw_nplus1,nIter=i,diffDaw=diffDw)
    out
  }
  DawNO_result_fo <- italg(Dw0_fo)
  DawNO_est_fo <- DawNO_result_fo$DawNO
  
  CawNO_est_fo <- (I + S*DawNO_est_fo)/DawNO_est_fo
  
  Ic <- I/(1-(fM+fH)*DawNO_est_fo/(2*(fM*fH)+(fH^3-fM^3)/(6*fM^2*fH^2*(fH-fM)/DawNO_est_fo^2)))
  Sc <- S-DawNO_est_fo*Ic/2/fM/fH+(fM+fH)*DawNO_est_fo^2*Ic/6/fM^2/fH^2
  
  # iterative algorithm function, 3rd order
  italgto <- function(SV,nIterAlg=10){
    # object to store Dw values for each iteration
    outIterAlg <- numeric(nIterAlg)
    # set starting value (SV)
    Dw_n <- SV
    for(i in 1:nIterAlg){
      Dw_nplus1 <- outIterAlg[i] <-  -Ic*(exp(-Dw_n/fL) - exp(-Dw_n/fM))/(eL-eM)
      diffDw <- abs(Dw_n - Dw_nplus1)
      #print(i);print(Dw_nplus1)
      # stop if values of Dw_n and Dw_nplus1 are close enough
      if(diffDw <0.001){break}
      Dw_n <- Dw_nplus1  
    }
    out <- list(DawNO=Dw_nplus1,nIter=i,diffDaw=diffDw)
    out
  }
  DawNO_result_to <- italgto(DawNO_est_fo)
  DawNO_est_to <- DawNO_result_to$DawNO
  
  DawNO_est_final <- DawNO_est_to  
  CawNO_est_final <- Ic/DawNO_est_final + Sc
  CaNO_est_final <- Sc
  JawNO_est_final <- DawNO_est_final*(CawNO_est_final - CaNO_est_final)
  
  outfinal <- c(  CaNO=unname(CaNO_est_final),
                  DawNO=unname(DawNO_est_final),
                  JawNO=unname(JawNO_est_final),
                  CawNO=unname(CawNO_est_final)
  )
  # data checks (both should be true)
  ## test CaNO estimate is positive
  test1 <- CaNO_est_final > 0
  ## test whether measured data is mathematically consistent with model
  LHS <- (eL-eM)/(eM-eH)
  RHS <- fH/fL*((fM-fL)/(fH-fM))
  test2 <- LHS < RHS
  ## To change the function to only give non-missing values for 
  ## valid datasets that produce positive CaNO estimates, remove the 
  ## comments from the following two lines
  # if(!(test2)) outfinal[1:4] <- NA # leave in negative CaNO estimates
  # if(!(test1 & test2)) outfinal[1:4] <- NA
  c(outfinal,valid=unname(test2))
}

################################################################################
################## APPLICATION 3 ###############################################
################################################################################

## Run RC02_ME2_Fit_Functions_C.R

library(survival)
library(ggplot2)

################################################################################
########################## DATA (DAYS) ##########################################
################################################################################

library(KMsurv)

data(larynx)

str(larynx)

time <- larynx$time
status <- larynx$delta   # 1 = death

x <- time[status == 1]

n <- length(x)

#################################
########## FIT RESULTS ##########
#################################

Results <- FIT.MEG2.RD(x)

##########################
##### Performance measures
##########################

Results$Start.Table

Results$Comparison

##############################
###### Log-likelihood graphics
##############################

print(Results$G.NR)

print(Results$G.NGEM)

##################################
##### Estimates by algorithm #####
##################################

Results$NR.estimates

Results$NGEM.estimates

###################################################
########## COMPARISON USING GOF MEASURES ##########
###################################################

RES.GOF <- Compare.MEG(x)

RES.GOF$Summary

RES.GOF$Estimates

RES.GOF$GOF.Table

RES.GOF$Ranking

RES.GOF$Top.Models

#####################
##### HISTOGRAM #####

print(RES.GOF$Plots$Histogram)

###############################################
##### KM and ESTIMATED SURVIVAL FUNCTIONS #####
###############################################

print(RES.GOF$Plots$Survival)

##############################################################################
## WALD CONFIDENCE INTERVALS
##############################################################################

fit.MEG2a <- fit.MEG2(
  
  x = x,
  
  eta0 = log(2),
  
  nu0 = log(1/mean(x))
  
)

Wald.MEG2 <- data.frame(
  
  Algorithm = c("NR","NR","NGEM","NGEM"),
  
  Parameter = c("alpha","beta","alpha","beta"),
  
  Lower = c(
    
    fit.MEG2a$Wald$NR$CI.alpha[1],
    
    fit.MEG2a$Wald$NR$CI.beta[1],
    
    fit.MEG2a$Wald$NGEM$CI.alpha[1],
    
    fit.MEG2a$Wald$NGEM$CI.beta[1]
    
  ),
  
  Upper = c(
    
    fit.MEG2a$Wald$NR$CI.alpha[2],
    
    fit.MEG2a$Wald$NR$CI.beta[2],
    
    fit.MEG2a$Wald$NGEM$CI.alpha[2],
    
    fit.MEG2a$Wald$NGEM$CI.beta[2]
    
  ),
  
  row.names = NULL
  
)

Wald.MEG2

##############################################################################
## LR CONFIDENCE INTERVALS
##############################################################################

LR.MEG2 <- data.frame(
  
  Algorithm = c("NR","NR","NGEM","NGEM"),
  
  Parameter = c("alpha","beta","alpha","beta"),
  
  Lower = c(
    
    fit.MEG2a$LR$NR$CI.alpha[1],
    
    fit.MEG2a$LR$NR$CI.beta[1],
    
    fit.MEG2a$LR$NGEM$CI.alpha[1],
    
    fit.MEG2a$LR$NGEM$CI.beta[1]
    
  ),
  
  Upper = c(
    
    fit.MEG2a$LR$NR$CI.alpha[2],
    
    fit.MEG2a$LR$NR$CI.beta[2],
    
    fit.MEG2a$LR$NGEM$CI.alpha[2],
    
    fit.MEG2a$LR$NGEM$CI.beta[2]
    
  ),
  
  row.names = NULL
  
)

LR.MEG2


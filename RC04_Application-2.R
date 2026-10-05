################################################################################
################## APPLICATION 2 ###############################################
################################################################################

## Run RC02_ME2_Fit_Functions_C.R

library(survival)
library(ggplot2)

################################################################################
########################## DATA (DAYS) ##########################################
################################################################################

x <- c(0.3, 0.3, 4.0, 5.0, 5.6, 6.2, 6.3, 6.6, 6.8, 7.4, 7.5, 8.4, 8.4, 10.3,11.0, 11.8, 12.2, 12.3, 13.5, 14.4, 14.4, 14.8,
       15.5, 15.7, 16.2, 16.3, 16.5, 16.8, 17.2, 17.3, 17.5, 17.9, 19.8, 20.4, 20.9, 21.0, 21.0, 21.1, 23.0, 23.4, 23.6,
       24.0, 24.0, 27.9, 28.2, 29.1, 30.0, 31.0, 31.0, 32.0, 35.0, 35.0, 37.0, 37.0, 37.0, 38.0, 38.0, 38.0, 39.0, 39.0,
       40.0, 40.0, 40.0, 41.0, 41.0, 41.0, 42.0, 43.0, 43.0, 43.0, 44.0, 45.0, 45.0, 46.0, 46.0, 47.0, 48.0, 49.0, 51.0,
       51.0, 51.0, 52.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 60.0, 60.0, 61.0, 62.0, 65.0, 65.0, 67.0, 67.0, 68.0,
       69.0, 78.0, 80.0,83.0, 88.0, 89.0, 90.0, 93.0, 96.0, 103.0, 105.0, 109.0, 109.0, 111.0, 115.0, 117.0, 125.0,
       126.0, 127.0, 129.0, 129.0, 139.0, 154)

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


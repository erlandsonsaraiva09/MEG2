################################################################################
########## Functions for the MEG2 distribution #################################
################################################################################

################################################################################

##### This file contains functions to:
# - evaluate the density function (dMEG2)
# - evaluate the cumulative distribution function (pMEG2) and 
#   the survival function via lower.tail = FALSE
# - compute quantiles (qMEG2)
# - generate random samples (rMEG2)

# The implemented functions follow the standard R conventions adopted by
# dnorm(), pnorm(), qnorm(), and rnorm().
#
# Additional options:
# - dMEG2(): log-density via log = TRUE
# - pMEG2(): survival probabilities via lower.tail = FALSE
#            and log-probabilities via log.p = TRUE
# - qMEG2(): upper-tail quantiles via lower.tail = FALSE
#            and log-probabilities via log.p = TRUE
# - rMEG2(): random generation by mixture representation (default)
#            or inverse transform method (option = 2)

################################################################################

############################################################
###### Density function of the MEG2 distribution ###########
############################################################

       dMEG2 <- function(x, alpha, beta, log = FALSE){
           if(alpha <= 1)
         stop("alpha must be greater than 1.")
           if(beta <= 0)
         stop("beta must be positive.")
  
# Initialize density vector
 
        dens <- numeric(length(x))
  
# Support: x > 0
  
         ind <- which(x > 0)
  
           if(length(ind) > 0)
             {
          xx <- x[ind]
   dens[ind] <- (beta^(alpha - 1)/(1 + beta^(alpha - 2)))*(1 + beta * xx^(alpha - 1) / gamma(alpha)) *exp(-beta*xx)}
  
# Outside the support
  
dens[x <= 0] <- 0
  
           if(log)
             {
        dens <- log(dens)
             }
  
       return(dens)
             }

############################################################
##### CDF of the MEG2 distribution #########################
############################################################
       
       pMEG2 <- function(x, alpha, beta, lower.tail = TRUE, log.p = FALSE)
             {
           if(alpha <= 1)
         stop("alpha must be greater than 1.")
           if(beta <= 0)
         stop("beta must be positive.")
         
# Initialize CDF vector
  
           F <- numeric(length(x))
         
# Support: x > 0

         ind <- which(x > 0)
           if(length(ind) > 0)
             {
          xx <- x[ind]
      F[ind] <- 1 - (beta^(alpha - 2)*exp(-beta*xx) + pgamma(beta*xx, shape = alpha, lower.tail = FALSE))/(1 + beta^(alpha - 2))
             }
         
# For x <= 0, F(x) = 0
      
   F[x <= 0] <- 0
         
# Upper tail probabilities

           if(!lower.tail)
           F <- 1 - F
         
# Log-probabilities
   
           if(log.p)
           F <- log(pmax(F, .Machine$double.xmin))
       return(F)
             }
       
############################################################
##### Quantile function of the MEG2 distribution ###########
############################################################
       
       qMEG2 <- function(p, alpha, beta, lower.tail = TRUE, log.p = FALSE)
             {
           if(alpha <= 1)
         stop("alpha must be greater than 1.")
           if(beta <= 0)
         stop("beta must be positive.")
         
# Transform probabilities if necessary

           if(log.p)
           p <- exp(p)
         
           if(!lower.tail)
           p <- 1 - p
         
           if(any(p < 0 | p > 1))
         stop("Probabilities must be between 0 and 1.")
         
# Output vector
  
           q <- numeric(length(p))
         
# Boundary cases
  
   q[p == 0] <- 0
   q[p == 1] <- Inf
         
# Internal probabilities
    
         ind <- which(p > 0 & p < 1)
         
for(i in ind){
           u <- p[i]
           
########################################################
# Quantile equation: F(x|alpha,beta) - u = 0
########################################################
           
           h <- function(x){pMEG2(x, alpha = alpha, beta = beta) - u}
           
########################################################
# Search procedure for the interval
########################################################
           
           l <- 1
           
        while(h(l) <= 0 && l < 10000)
           l <- l + 1
           
           if(l == 10000)
         stop("Unable to bracket the root.")
           
########################################################
# Compute the quantile
########################################################
           
        q[i] <- uniroot(h,interval = c(l - 1, l))$root
             }
         
       return(q)
             }
       
############################################################
##### Random generation from the MEG2 distribution #########
# option = 1 : Mixture representation (default) ############
# option = 2 : Inverse transform method ####################
############################################################
       
       rMEG2 <- function(n, alpha, beta, option = 1)
             {
           if(alpha <= 1)
         stop("alpha must be greater than 1.")
           if(beta <= 0)
         stop("beta must be positive.")
         
           if(length(n) != 1 || !is.numeric(n) || n <= 0 || n %% 1 != 0)
         stop("The typed value of 'n' must be a positive integer.")
         
           n <- as.integer(n)
         
           if(!(option %in% c(1, 2)))
         stop("Invalid option. Use option = 1 (mixture representation) or option = 2 (inverse transform).")
         
################
# Mixture weight
################
         
           w <- beta^(alpha - 2)/(1 + beta^(alpha - 2))
         
############################################
# Option 1: Mixture representation (default)
############################################
         
           if(option == 1)
             {
           z <- rbinom(n, size = 1, prob = w)
           x <- ifelse(z == 1, rexp(n, rate = beta), rgamma(n, shape = alpha, rate = beta))
       return(x)
             }
         
####################################
# Option 2: Inverse transform method
####################################
         
           if(option == 2)
             {
           u <- runif(n)
           x <- qMEG2(p = u, alpha = alpha, beta = beta)
       return(x)
             }
         stop("Invalid option. Use option = 1 (mixture representation) or option = 2 (inverse transform).")
             }
       
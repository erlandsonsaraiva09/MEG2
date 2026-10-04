################################################################################
########## Functions for fitting the MEG2 model to a data set ##################
################################################################################

## ## Run RC01_ME2_Functions.R

################################################################################
########## REQUIRED PACKAGES ####################################################
################################################################################
 
               library(ggplot2)
               library(reshape2)
   
############################################################
##### Random generation from the MEG2 distribution #########
# option = 1 : Mixture representation (default) ############
# option = 2 : Inverse transform method ####################
############################################################

                rMEG2 <- function(n, alpha, beta, option = 1)
                      {
                    if(!is.numeric(alpha) || length(alpha)!=1 || !is.finite(alpha) || alpha<=1)
                  stop("alpha must be greater than 1.")

                    if(!is.numeric(beta) || length(beta)!=1 || !is.finite(beta) || beta<=0)
                  stop("beta must be positive.")
  
                    if(length(n) != 1 || !is.numeric(n) || n <= 0 || n %% 1 != 0)
                  stop("The typed value of 'n' must be a positive integer.")
  
                    n <- as.integer(n)
  
                    if(!(option %in% c(1, 2)))
                  stop("Invalid option. Use option = 1 (mixture representation) or option = 2 (inverse transform).")
  
################
# Mixture weight
################
  
                    b <- beta^(alpha-2)
                    w <- b/(1+b)  
                    
############################################
# Option 1: Mixture representation (default)
############################################
  
                    if(option == 1)
                      {
                    z <- rbinom(n, size = 1, prob = w)
                    x <- numeric(n)
               id.exp <- z == 1
           x[id.exp]  <- rexp(sum(id.exp), rate = beta)
           x[!id.exp] <- rgamma(sum(!id.exp), shape = alpha, rate = beta)
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

################################################################################
########## TRANSFORMATIONS (FIXED) ##############################################
################################################################################
############################################################
# Transformation ensuring alpha > 1.
# A small offset is added to prevent numerical instabilities
# in the evaluation of lgamma(alpha - 1) and digamma(alpha - 1)
# when alpha approaches the boundary value 1.
############################################################

       alpha_from_eta <- function(eta)
                      {
                alpha <- 1+exp(eta)
                  pmax(alpha, 1 + 1e-8)
                      }
                
         beta_from_nu <- function(nu)
                      {
                 beta <- exp(nu)
                  pmax(beta,1e-12)
                      }
         
################################################################################
########## PENALIZED LOG-LIKELIHOOD) ###########################################
################################################################################
         
          loglik_pena <- function(alpha,beta,x){
                    n <- length(x)
                    A <- exp(log(beta) + (alpha - 1)*log(x) - lgamma(alpha))              
                   ll <- n*(alpha-1)*log(beta) - n*log(1+beta^(alpha-2)) + sum(log1p(A)) - beta*sum(x) - lgamma(alpha-1)
                return(ll)
                      }
           
################################################################################
########## Score functions ######################################################
################################################################################

          score_alpha <- function(alpha, beta, x){
                   n  <- length(x)
                    A <- exp(log(beta) + (alpha - 1)*log(x) - lgamma(alpha))              
               alpha1 <- pmax(alpha - 1, 1e-8)
                    U <- n*log(beta)/(1 + beta^(alpha - 2)) + sum((A/(1 + A))*(log(x) - digamma(alpha))) - digamma(alpha1)
                return(U)
                      }
            
           score_beta <- function(alpha, beta, x){
                    n <- length(x)
                    A <- exp((alpha - 1)*log(x) - lgamma(alpha))
           
                     n*(alpha - 1)/beta - n*(alpha - 2)*beta^(alpha - 3)/(1 + beta^(alpha - 2)) + sum(A/(1 + beta*A)) - sum(x)
                      }

################################################################################
########## Hessian (diagonal terms) ############################################
################################################################################

##### alpha
    
           hess_alpha <- function(alpha, beta, x){
                   n  <- length(x)
                    A <- exp(log(beta) + (alpha - 1)*log(x) - lgamma(alpha))      
               alpha1 <- pmax(alpha-1,1e-8)
        
                term1 <- -n*beta^alpha*(beta*log(beta)/(beta^2+beta^alpha))^2
                term2 <- sum((A/(1+A))*((log(x)-digamma(alpha))^2/(1+A) - trigamma(alpha)))
                term3 <- -trigamma(alpha1)
                term1 + term2 + term3
                      }
    
##### beta
    
            hess_beta <- function(alpha,beta,x){
                    n <- length(x)
        
                term1 <- -n*(alpha-1)/beta^2
                term2 <- -n*(alpha-2)*(alpha-3)*beta^(alpha-4)/(1+beta^(alpha-2))
                term3 <- n*(alpha-2)^2*beta^(2*alpha-6)/(1+beta^(alpha-2))^2
                    A <- exp((alpha - 1)*log(x) - lgamma(alpha))
                term4 <- -sum((A/(1 + beta*A))^2)        
                term1 + term2 + term3 + term4
                      }
     
################################################################################
########## Coordinate-wise Newton-Raphson algorithm ############################
################################################################################
     
              NR_MEG2 <- function(x,
                                  eta0 = log(2),
                                  nu0 = log(1/mean(x)),
                                  tol = 1e-6,
                                  maxit = 1000,
                                  alpha.fixed = NULL,
                                  beta.fixed  = NULL){
              
###############################################################################
##### Initial values
###############################################################################
              
                  eta <- eta0
                  nu  <- nu0
             ll_trace <- rep(NA, maxit)
            converged <- FALSE
              
###############################################################################
##### Fixed parameters
###############################################################################
              
                    if(!is.null(alpha.fixed))
                  eta <- log(alpha.fixed - 1)
              
                    if(!is.null(beta.fixed))
                   nu <- log(beta.fixed)
              
                   t0 <- proc.time()[3]
              
###############################################################################
##### Newton iterations
###############################################################################
              
     for(k in 1:maxit){
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
          ll_trace[k] <- loglik_pena(alpha, beta, x)
                
###############################################################################
##### Update eta
###############################################################################
                
              U_alpha <- score_alpha(alpha, beta, x)
              H_aa    <- hess_alpha(alpha, beta, x)
                
                U_eta <- exp(eta)*U_alpha
                H_eta <- exp(eta)*U_alpha + exp(2*eta)*H_aa
                
                    if(!is.finite(H_eta)){
                H_eta <- -1e-10
                      }
                  else 
                    if(abs(H_eta) < 1e-10)
                      {
                H_eta <- ifelse(H_eta >= 0, 1e-10, -1e-10)
                      }
                
                    if(is.null(alpha.fixed)){
              eta_new <- eta - U_eta/H_eta
                  
                    if(!is.finite(eta_new))
              eta_new <- eta
                      }
                  else{
              eta_new <- eta
                      }
                
###############################################################################
##### Update nu (Gauss-Seidel)
###############################################################################
                
            alpha_new <- alpha_from_eta(eta_new)
                
               U_beta <- score_beta(alpha_new, beta, x)
               H_bb   <- hess_beta(alpha_new, beta, x)
                
                 U_nu <- exp(nu)*U_beta
                 H_nu <- exp(nu)*U_beta + exp(2*nu)*H_bb
                
                    if(!is.finite(H_nu)){
                 H_nu <- -1e-10
                      }
                  else 
                    if(abs(H_nu) < 1e-10){
                 H_nu <- ifelse(H_nu >= 0, 1e-10, -1e-10)
                      }
           
                    if(is.null(beta.fixed)){
               nu_new <- nu - U_nu/H_nu
                  
                    if(!is.finite(nu_new))
               nu_new <- nu
                      }
                  else{
               nu_new <- nu
                      }
                
###############################################################################
##### Convergence criterion
###############################################################################
                
                 diff <- abs(eta_new - eta) +
                   abs(nu_new - nu)
                
                  eta <- eta_new
                  nu  <- nu_new
                
                    if(is.finite(diff) && diff < tol){
                  
            converged <- TRUE
                break
                      }
                      }
              
###############################################################################
##### Computational time
###############################################################################
              
              time_NR <- proc.time()[3] - t0
              
###############################################################################
##### Final estimates
###############################################################################
              
              eta_hat <- eta
              nu_hat  <- nu
              
            alpha_hat <- alpha_from_eta(eta_hat)
            beta_hat  <- beta_from_nu(nu_hat)
              
###############################################################################
##### Output
###############################################################################
              
                return(
                  list(eta_hat = eta_hat,
                       nu_hat = nu_hat,
                       alpha_hat = alpha_hat,
                       beta_hat  = beta_hat,
                       logLik = loglik_pena(alpha_hat, beta_hat, x),
                       iter = k,
                       ll_trace = ll_trace[1:k],
                       converged = converged,
                       time = time_NR
                      )
                      )
                      }
            
################################################################################            
########## NGEM - eta and nu ###################################################
################################################################################
          
######################################
###### E-step: Posterior probabilities
######################################
            
    tau_update_eta_nu <- function(eta, nu, x){
              
################
##### Parameters
################
              
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
              
####################
##### Mixture weight
####################
              
                    b <- beta^(alpha - 2)
                    w <- b/(1 + b)
              
#########################
##### Component densities
#########################
              
                   f1 <- beta*exp(-beta*x)
                   f2 <- exp(alpha*log(beta) + (alpha - 1)*log(x) - beta*x - lgamma(alpha))
              
#############################
##### Posterior probabilities
#############################
              
                  den <- w*f1 + (1 - w)*f2
                  tau <- (w*f1)/den
              
##########################
##### Numerical protection
##########################
              
 tau[!is.finite(tau)] <- 0.5
                  tau <- pmin(pmax(tau, 1e-12), 1 - 1e-12)
                return(tau)
                      }
            
########################################
##### Penalized complete-data Q-function
########################################
            
          Q_pc_eta_nu <- function(eta, nu, x, tau){
              
################
##### Parameters
################
              
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
              
#################
##### Sample size
#################
              
                    n <- length(x)
              
################
##### Q-function
################
              
                    Q <- (n*alpha - sum(tau))*log(beta) - n*log1p(beta^(alpha - 2)) - sum((1 - tau)*lgamma(alpha)) + (alpha - 1)*sum((1 - tau)*log(x)) - beta*sum(x) - lgamma(alpha - 1)
              
##########################
##### Numerical protection
##########################
              
                    if(!is.finite(Q))
                    Q <- -Inf
                return(Q)
                      }
            
###################################
##### Complete-data score for alpha
###################################
            
       score_alpha_pc <- function(alpha, beta, x, tau){
                    n <- length(x)
               alpha1 <- pmax(alpha - 1, 1e-8)
              
                    U <- n*log(beta) - n*beta^(alpha - 2)*log(beta)/(1 + beta^(alpha - 2)) - sum((1 - tau)*digamma(alpha)) + sum((1 - tau)*log(x)) - digamma(alpha1)
                return(U)
                      }
            
##################################
##### Complete-data score for beta
##################################
            
        score_beta_pc <- function(alpha, beta, x, tau){
                    n <- length(x)
                    U <- (n*alpha - sum(tau))/beta - n*(alpha - 2)*beta^(alpha - 3)/(1 + beta^(alpha - 2)) - sum(x)
                return(U)
                      }
            
#####################################
##### Complete-data Hessian for alpha
#####################################
            
        hess_alpha_pc <- function(alpha, beta, x, tau){
                    n <- length(x)
               alpha1 <- pmax(alpha - 1, 1e-8)
                    J <- -n*log(beta)^2*beta^(alpha - 2)/(1 + beta^(alpha - 2))^2 -  sum((1 - tau)*trigamma(alpha)) - trigamma(alpha1)
                return(J)
                      }
            
####################################
##### Complete-data Hessian for beta
####################################
            
         hess_beta_pc <- function(alpha, beta, x, tau){
                    n <- length(x)
                    J <- -(n*alpha - sum(tau))/beta^2 
                         - n*(alpha - 2)*(alpha - 3)*
                           beta^(alpha - 4)/
                           (1 + beta^(alpha - 2)) +
                           n*(alpha - 2)^2*
                           beta^(2*alpha - 6)/(1 + beta^(alpha - 2))^2
              
                return(J)
                      }
            
#################################
##### Transformed score functions
#################################
            
         score_eta_pc <- function(eta, nu, x, tau){
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
      
                   exp(eta) * score_alpha_pc(alpha, beta, x, tau)
              
                      }
            
###############################################################################
            
          score_nu_pc <- function(eta, nu, x, tau){
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
              
                   exp(nu) * score_beta_pc(alpha, beta, x, tau)
              
                      }
            
###############################################################################
##### Transformed Hessian functions
###############################################################################
            
          hess_eta_pc <- function(eta, nu, x, tau){
              
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
              
              U.alpha <- score_alpha_pc(alpha, beta, x, tau)
              J.alpha <- hess_alpha_pc(alpha, beta, x, tau)
              
                   exp(eta)*U.alpha + exp(2*eta)*J.alpha
                  
                      }
            
###############################################################################
            
           hess_nu_pc <- function(eta, nu, x, tau){
                alpha <- alpha_from_eta(eta)
                beta  <- beta_from_nu(nu)
              
               U.beta <- score_beta_pc(alpha, beta, x, tau)
               J.beta <- hess_beta_pc(alpha, beta, x, tau)
              
                   exp(nu)*U.beta + exp(2*nu)*J.beta
              
                      }
            
################################################################################
##### Newton-GEM algorithm #####################################################
################################################################################
            
          NGEM_eta_nu <- function(x,
                                  eta0 = log(2),
                                  nu0 = log(1/mean(x)),
                                  tol = 1e-6,
                                  max_iter = 1000,
                                  alpha.fixed = NULL,
                                  beta.fixed  = NULL){
              
####################
##### Initial values
####################
              
                  eta <- eta0
                  nu  <- nu0
              
                    if(!is.null(alpha.fixed))
                  eta <- log(alpha.fixed - 1)
              
                    if(!is.null(beta.fixed))
                   nu <- log(beta.fixed)
              
             ll_trace <- rep(NA,max_iter)
              
            converged <- FALSE
              
                   t0 <- proc.time()[3]
              
#####################
##### Main iterations
#####################
              
for(iter in 1:max_iter){
                
#######################################
##### Observed penalized log-likelihood
#######################################
                
        ll_trace[iter] <- loglik_pena(alpha_from_eta(eta), beta_from_nu(nu), x)
                
##################
##### E-step #####
##################
                
                   tau <- tau_update_eta_nu(eta,nu, x)
                
########################
##### M-step : eta #####
########################
                
                     U <- score_eta_pc(eta, nu, x, tau)
                     J <- hess_eta_pc(eta, nu, x, tau)
                
                     if(!is.finite(J))
                     J <- -1e-10
                
                     if(abs(J) < 1e-10)
                     J <- ifelse(J>=0,1e-10,-1e-10)
                
                     d <- -U/J
                
                     if(!is.finite(d))
                     d <- 0
                
                     d <- max(min(d,5),-5)
                
                  step <- 1
                
                 Q.old <- Q_pc_eta_nu(eta, nu, x, tau)
                
                 repeat{
               eta.try <- eta + step*d
                 Q.new <- Q_pc_eta_nu( eta.try, nu, x, tau)
                  
                     if(is.finite(Q.new) && Q.new > Q.old)
                  break
                  
                  step <- step/2
                  
                     if(step < 1e-8)
                  break
                       }
                
                     if(is.null(alpha.fixed))
               eta.new <- eta + step*d
           else
               eta.new <- eta
                
#######################
##### M-step : nu ##### 
#######################
                
                   tau <- tau_update_eta_nu(eta.new, nu, x)
                     U <- score_nu_pc(eta.new, nu, x, tau)
                     J <- hess_nu_pc(eta.new, nu, x, tau)
                
                     if(!is.finite(J))
                     J <- -1e-10
                
                     if(abs(J) < 1e-10)
                     J <- ifelse(J>=0,1e-10,-1e-10)
                
                     d <- -U/J
                
                     if(!is.finite(d))
                     d <- 0
                
                     d <- max(min(d,5),-5)
                
                  step <- 1
                
                 Q.old <- Q_pc_eta_nu(eta.new, nu, x, tau)
                
                 repeat{
                nu.try <- nu + step*d
                 Q.new <- Q_pc_eta_nu(eta.new, nu.try, x, tau)
                  
                     if(is.finite(Q.new) && Q.new > Q.old)
                break
                  step <- step/2
                  
                     if(step < 1e-8)
                break
                       }
                
                     if(is.null(beta.fixed))
                nu.new <- nu + step*d
            else
                nu.new <- nu
                
#######################
##### Convergence #####
#######################
                
                 delta <- abs(eta.new-eta) + abs(nu.new-nu)
                   eta <- eta.new
                   nu  <- nu.new
                
                     if(is.finite(delta) && delta < tol){
             converged <- TRUE
             break
                       }
                       }
              
##############################
##### Computational time #####
##############################
              
             time.NGEM <- proc.time()[3]-t0
              
#####################
##### Estimates #####
#####################
              
               eta_hat <- eta
                nu_hat <- nu
             alpha_hat <- alpha_from_eta(eta_hat)
              beta_hat <- beta_from_nu(nu_hat)
              
##################
##### Output #####
##################
              
                 return(
                   list(eta_hat=eta_hat,
                        nu_hat=nu_hat,
                        alpha_hat=alpha_hat,
                        beta_hat=beta_hat,
                        logLik=loglik_pena(alpha_hat, beta_hat, x),
                        iter=iter,
                        converged=converged,
                        ll_trace=ll_trace[1:iter],
                        time=time.NGEM
                       )
                       )
                       }
            
################################################################################            
########## WALD INTERVALS ######################################################
################################################################################
            
###########################################
##### Penalized log-likelihood in (eta, nu)
###########################################
            
     loglik_pen_eta_nu <- function(eta, nu, x){
                 alpha <- alpha_from_eta(eta)
                  beta <- beta_from_nu(nu)
            loglik_pena(alpha = alpha, beta  = beta, x = x)
                       }
            
#######################################################
##### Numerical Hessian of the penalized log-likelihood
#######################################################
            
        Hessian_eta_nu <- function(fit, x, h = 1e-5){
              
#####################
##### Estimates #####
#####################
              
                   eta <- fit$eta_hat
                   nu  <- fit$nu_hat
              
####################################
##### Penalized log-likelihood #####
####################################
              
                     f <- function(eta, nu){
      loglik_pen_eta_nu(eta = eta, nu  = nu, x   = x)
                       }
              
###############################
##### Adaptive step sizes #####
###############################
              
                    h1 <- h*max(1, abs(eta))
                    h2 <- h*max(1, abs(nu))
              
#################################################
##### Second derivative with respect to eta #####
#################################################
              
                    f0 <- f(eta, nu)
                   H11 <- (f(eta + h1, nu) - 2*f0 +f(eta - h1, nu))/h1^2
              
################################################
##### Second derivative with respect to nu #####
################################################
              
                   H22 <- (f(eta, nu + h2) - 2*f0 + f(eta, nu - h2))/h2^2
              
############################
##### Cross derivative #####
############################
              
                   H12 <- (f(eta + h1, nu + h2) - f(eta + h1, nu - h2) - f(eta - h1, nu + h2) + f(eta - h1, nu - h2))/(4*h1*h2)
              
##########################
##### Hessian matrix #####
##########################
              
                     H <- matrix(c(H11, H12, H12, H22), nrow = 2, byrow = TRUE)
              
##############################
##### Numerical symmetry #####
##############################
              
                     H <- (H + t(H))/2
                 return(H)
                       }
            
#################################################
##### Covariance matrix and standard errors #####
#################################################
            
            SE_eta_nu <- function(fit, x, h = 1e-5){
              
#############################
##### Numerical Hessian #####
#############################
              
                    H <- Hessian_eta_nu(fit = fit, x = x, h = h)
              
#######################################
##### Observed information matrix #####
#######################################
              
                    J <- -H
              
#############################
##### Covariance matrix ##### 
#############################
              
                Sigma <- tryCatch(
                 solve(J),
                error = function(e){
                    if(requireNamespace("MASS", quietly = TRUE)){
            MASS::ginv(J)
                      }else{
                matrix(NA,2,2)
                      }
                      }
                      )
              
###########################
##### Standard errors #####
###########################
              
                   se <- sqrt(pmax(diag(Sigma),0))
            names(se) <- c("eta","nu")
              
##################
##### Output #####
##################
              
                return(
                  list(Hessian = H,
                       Information = J,
                       Covariance = Sigma,
                       eta = se[1],
                       nu = se[2]
                      )
                      )
                      }
            
#####################################
##### Wald confidence intervals #####
#####################################
            
   WaldIntervals.MEG2 <- function(fit, x, level = 0.95, h = 1e-5){
              
###########################
##### Standard errors #####
###########################
              
                   SE <- SE_eta_nu(fit = fit, x = x, h = h)
              
##########################
##### Critical value #####
##########################
              
                    z <- qnorm(1 - (1 - level)/2)
              
#####################
##### Estimates ##### 
#####################
              
              eta.hat <- fit$eta_hat
              nu.hat  <- fit$nu_hat
              
######################################
##### Wald intervals in (eta,nu) #####
######################################
              
               CI.eta <- c(eta.hat - z*SE$eta, eta.hat + z*SE$eta)
                CI.nu <- c(nu.hat - z*SE$nu, nu.hat + z*SE$nu)
              
###############################
##### Back-transformation #####
###############################
              
             CI.alpha <- alpha_from_eta(CI.eta)
              CI.beta <- beta_from_nu(CI.nu)
              
##################
##### Output #####
##################
              
                  out <- list(
                
 #####################
 ##### Estimates #####
######################
                
             Estimate = c(alpha = fit$alpha_hat, beta  = fit$beta_hat),
                
###########################
##### Standard errors #####
###########################
                
                   SE = c(eta = SE$eta, nu  = SE$nu),
                
#############################
##### Covariance matrix #####
#############################
                
           Covariance = SE$Covariance,
                
################################
##### Confidence intervals #####
################################
                
               CI.eta = CI.eta,
                CI.nu = CI.nu,
             CI.alpha = CI.alpha,
              CI.beta = CI.beta,
                
##############################
##### Information matrix #####
##############################
                
          Information = SE$Information,
                
############################
##### Confidence level #####
############################
                
                level = level
                      )
           class(out) <- "WaldIntervals.MEG2"
                return(out)
                      }
            
################################################################################            
########## LR INTERVALS ########################################################
################################################################################
           
############################################
##### Profile penalized log-likelihood #####
############################################
            
 ProfileLogLik_eta_nu <- function(value,
                                  parameter = c("eta","nu"),
                                  fit,
                                  x,
                                  algorithm = c("NGEM","NR"),
                                  tol = 1e-6,
                                  max_iter = 1000){
              
#####################
##### Arguments #####
#####################
              
            parameter <- match.arg(parameter)
            algorithm <- match.arg(algorithm)
              
############################################
##### Maximum penalized log-likelihood #####
############################################
              
           loglik.max <- fit$logLik
              
####################################
##### Constrained optimization #####
####################################
               
                  FIT <- tryCatch({
                    if(parameter == "eta"){
          alpha.fixed <- alpha_from_eta(value)
                     if(algorithm == "NGEM"){
            NGEM_eta_nu(x = x,
                        eta0 = fit$eta_hat,
                        nu0 = fit$nu_hat,
                        alpha.fixed = alpha.fixed,
                        tol = tol,
                        max_iter = max_iter
                       )
                    
                       }else{
                NR_MEG2(x = x,
                        eta0 = fit$eta_hat,
                        nu0 = fit$nu_hat,
                        alpha.fixed = alpha.fixed,
                        tol = tol,
                        maxit = max_iter
                       )
                       }
                       }else{
            beta.fixed <- beta_from_nu(value)
                     if(algorithm == "NGEM"){
            NGEM_eta_nu(x = x,
                        eta0 = fit$eta_hat,
                        nu0 = fit$nu_hat,
                        beta.fixed = beta.fixed,
                        tol = tol,
                        max_iter = max_iter
                       )
                       }else{
                NR_MEG2(x = x,
                        eta0 = fit$eta_hat,
                        nu0 = fit$nu_hat,
                        beta.fixed = beta.fixed,
                        tol = tol,
                        maxit = max_iter
                       )
                       }
                       }
                       },
                 error = function(e) NULL)
              
###############################
##### Failed optimization #####
###############################
              
                     if(is.null(FIT))
                 return(NULL)
              
                     if(!FIT$converged)
                 return(NULL)
              
                     if(!is.finite(FIT$logLik))
                 return(NULL)
              
                     if(!is.finite(FIT$alpha_hat))
                 return(NULL)
              
                     if(!is.finite(FIT$beta_hat))
                 return(NULL)
              
######################################
##### Likelihood-ratio statistic #####
######################################
              
                    LR <- 2*(loglik.max - FIT$logLik)
              
##################
##### Output #####
##################
              
                 return(
                   list(LR = LR,
                        alpha = FIT$alpha_hat,
                        beta = FIT$beta_hat,
                        eta = FIT$eta_hat,
                        nu = FIT$nu_hat
                       )
                       )
                       }
            
#################################################
##### Likelihood-ratio endpoint in (eta,nu) #####
#################################################
            
    LR_endpoint_eta_nu <- function(fit,
                                   x,
                                   parameter = c("eta","nu"),
                                   side = c("lower","upper"),
                                   level = 0.95,
                                   algorithm = c("NGEM","NR"),
                                   tol = 1e-6,
                                   max_iter = 1000,
                                   k = 1,
                                   eps = 1e-4,
                                   max.expand = 20){
              
#####################
##### Arguments #####
#####################
              
             parameter <- match.arg(parameter)
                  side <- match.arg(side)
             algorithm <- match.arg(algorithm)
              
###########################
##### Standard errors #####
###########################
              
                    SE <- SE_eta_nu(fit = fit, x = x)
                     if(parameter == "eta"){
             theta.hat <- fit$eta_hat
              theta.se <- SE$eta
                       }else{
             theta.hat <- fit$nu_hat
              theta.se <- SE$nu
                       }
              
###################################
##### Initial search interval #####
###################################
              
                     if(side == "lower"){
                 lower <- theta.hat - k*theta.se
                 upper <- theta.hat
                       }else{
                 lower <- theta.hat
                 upper <- theta.hat + k*theta.se
                       }
              
##########################
##### Critical value #####
##########################
              
                  crit <- qchisq(level, df = 1)
              
#####################################
##### Likelihood-ratio equation #####
#####################################
              
               LR.fun <- function(theta){
                  RES <- suppressWarnings(
  ProfileLogLik_eta_nu(value = theta,
                       parameter = parameter,
                       fit = fit,
                       x = x,
                       algorithm = algorithm,
                       tol = tol,
                       max_iter = max_iter
                      )
                      )
                
                    if(is.null(RES))
                return(NULL)
                
            RES$Value <- RES$LR - crit
                
                return(RES)
                      }
              
#####################################
##### Evaluate initial interval #####
#####################################
              
            RES.lower <- LR.fun(lower)
            RES.upper <- LR.fun(upper)
              
                    if(is.null(RES.lower) || is.null(RES.upper))
                return(NA_real_)
              
              f.lower <- RES.lower$Value
              f.upper <- RES.upper$Value
              
#######################################################
##### Expand interval until the root is bracketed #####
#######################################################
              
               expand <- 0
              
                 while(is.finite(f.lower) && is.finite(f.upper) && sign(f.lower) == sign(f.upper)){
                
########################################
##### Maximum number of expansions #####
########################################
                
                    if(expand >= max.expand){
               warning("Likelihood-ratio endpoint could not be bracketed.")
                return(NA_real_)
                      }
                
#####################################
##### Increase expansion factor #####
#####################################
                
                    k <- k + 1
                
##################################
##### Expand search interval #####
##################################
                
                    if(side == "lower"){
                lower <- theta.hat - k*theta.se
            RES.lower <- LR.fun(lower)
                  
                    if(is.null(RES.lower))
                return(NA_real_)
                  
              f.lower <- RES.lower$Value
        alpha.current <- RES.lower$alpha
        beta.current  <- RES.lower$beta
                      }else{
                upper <- theta.hat + k*theta.se
            RES.upper <- LR.fun(upper)
                  
                    if(is.null(RES.upper))
                return(NA_real_)
                  
              f.upper <- RES.upper$Value
        alpha.current <- RES.upper$alpha
                  
        beta.current  <- RES.upper$beta
                      }
                
################################
##### Check mixture weight ##### 
################################
                
                    if(!is.finite(alpha.current) || !is.finite(beta.current))
                return(NA_real_)
                
            w.current <- beta.current^(alpha.current - 2) / (1 + beta.current^(alpha.current - 2))
                
                    if(!is.finite(w.current))
                return(NA_real_)
                
                    if(w.current <= eps || w.current >= 1 - eps){
               warning("Search stopped because the mixture weight became numerically degenerate.")
                  
                return(NA_real_)
                      }
                
################################
##### Number of expansions #####
################################
                
               expand <- expand + 1
                      }
              
###################################
##### Root of the LR equation #####
###################################
              
                 ROOT <- tryCatch(
               uniroot(f = function(theta){
                       RES <- ProfileLogLik_eta_nu(
                       value = theta,
                       parameter = parameter,
                       fit = fit,
                       x = x,
                       algorithm = algorithm,
                       tol = tol,
                       max_iter = max_iter
                      )
                    
                    if(is.null(RES))
                return(Inf)
                    
               RES$LR - crit
                    
                      },
                  
                       lower = lower,
                       upper = upper,
                       tol = tol
                      )$root,
                error = function(e) NA_real_
                      )
              
##################
##### Output #####
##################
              
                return(ROOT)
                      }
            
#################################################
##### Likelihood-ratio confidence intervals #####
#################################################
            
     LRIntervals.MEG2 <- function(fit,
                                  x,
                                  level = 0.95,
                                  algorithm = c("NGEM","NR"),
                                  tol = 1e-6,
                                  max_iter = 1000,
                                  eps = 1e-4){
              
#####################
##### Algorithm #####
#####################
              
            algorithm <- match.arg(algorithm)
              
#####################################
##### Confidence limits for eta #####
#####################################
              
            eta.lower <- LR_endpoint_eta_nu(fit = fit,
                                            x = x,
                                            parameter = "eta",
                                            side = "lower",
                                            level = level,
                                            algorithm = algorithm,
                                            tol = tol,
                                            max_iter = max_iter,
                                            eps = eps)
              
            eta.upper <- LR_endpoint_eta_nu(fit = fit,
                                            x = x,
                                            parameter = "eta",
                                            side = "upper",
                                            level = level,
                                            algorithm = algorithm,
                                            tol = tol,
                                            max_iter = max_iter,
                                            eps = eps)
              
####################################
##### Confidence limits for nu #####
####################################
              
             nu.lower <- LR_endpoint_eta_nu(fit = fit,
                                            x = x,
                                            parameter = "nu",
                                            side = "lower",
                                            level = level,
                                            algorithm = algorithm,
                                            tol = tol,
                                            max_iter = max_iter,
                                            eps = eps)
              
             nu.upper <- LR_endpoint_eta_nu(fit = fit,
                                            x = x,
                                            parameter = "nu",
                                            side = "upper",
                                            level = level,
                                            algorithm = algorithm,
                                            tol = tol,
                                            max_iter = max_iter,
                                            eps = eps)
              
################################
##### Confidence intervals #####
################################
              
               CI.eta <- c(Lower = eta.lower, Upper = eta.upper)
                CI.nu <- c(Lower = nu.lower,Upper = nu.upper)
              
###############################
##### Back-transformation #####
###############################
              
             CI.alpha <- alpha_from_eta(CI.eta)
              CI.beta <- beta_from_nu(CI.nu)
              
      names(CI.alpha) <- c("Lower","Upper")
       names(CI.beta) <- c("Lower","Upper")
              
##################
##### Output #####
##################
              
                  OUT <- list(Estimate = c(alpha = fit$alpha_hat, beta = fit$beta_hat),
                              Estimate.eta.nu = c(eta = fit$eta_hat, nu = fit$nu_hat),
                              CI.eta = CI.eta,
                              CI.nu = CI.nu,
                              CI.alpha = CI.alpha,
                              CI.beta = CI.beta,
                              Level = level,
                              Algorithm = algorithm,
                              Tolerance = tol,
                              Max.Iterations = max_iter,
                              Max.Expansions = 20,
                              Weight.Tolerance = eps)
              
           class(OUT) <- "LRIntervals.MEG2"
              
                return(OUT)
              
                      }
            
########################
##### Print method #####
########################
            
print.LRIntervals.MEG2 <- function(x, ...){
              
                    cat("\n")
                    cat("=============================================================\n")
                    cat(" Likelihood-Ratio Confidence Intervals for MEG2\n")
                    cat("=============================================================\n\n")
              
#####################
##### Estimates #####
#####################
              
                    cat("Maximum penalized likelihood estimates\n\n")
                   
                   EST <- data.frame(Parameter = c("alpha","beta"), Estimate  = round(x$Estimate, 6), row.names = NULL)
                  print(EST, row.names = FALSE)
              
###############################################
###### Estimates in the transformed scale #####
###############################################
              
                    cat("\n")
                    cat("Estimates in the transformed scale\n\n")
              
                  EST2 <- data.frame(Parameter = c("eta","nu"), Estimate  = round(x$Estimate.eta.nu, 6), row.names = NULL)
                  print(EST2, row.names = FALSE)
              
##########################################
##### Confidence intervals (eta, nu) #####
##########################################
              
                    cat("\n")
                    cat("Likelihood-ratio confidence intervals (eta, nu)\n\n")
              
                  TAB1 <- data.frame(Parameter = c("eta","nu"),
                                     Lower = round(c(x$CI.eta[1], x$CI.nu[1]), 6),
                                     Upper = round(c(x$CI.eta[2], x$CI.nu[2]), 6),
                                     row.names = NULL)
                  print(TAB1, row.names = FALSE)
              
##############################################
##### Confidence intervals (alpha, beta) #####
##############################################
              
                    cat("\n")
                    cat("Likelihood-ratio confidence intervals (alpha, beta)\n\n")
              
                  TAB2 <- data.frame(Parameter = c("alpha","beta"),
                                     Lower = round(c(x$CI.alpha[1], x$CI.beta[1]), 6),
                                     Upper = round(c(x$CI.alpha[2], x$CI.beta[2]), 6),
                                    row.names = NULL)
                  print(TAB2, row.names = FALSE)
              
##################################
##### Additional information #####
##################################
              
                    cat("\n")
                    cat("Confidence level   :", x$Level, "\n")
                    cat("Algorithm          :", x$Algorithm, "\n")
                    cat("Tolerance          :", x$Tolerance, "\n")
                    cat("Maximum iterations :", x$Max.Iterations, "\n")
                    cat("Maximum expansions :", x$Max.Expansions, "\n")
                    cat("Weight tolerance   :", x$Weight.Tolerance, "\n")
                    cat("\n")
              
              invisible(x)
                       }
            
################################################################################            
########## FIT MEG 2 ###########################################################
################################################################################
            
##############################
##### MEG2 model fitting #####
##############################
            
              fit.MEG2 <- function(x,
                                   eta0 = log(2),
                                   nu0 = log(1/mean(x)),
                                   tol = 1e-6,
                                   max_iter = 1000,
                                   level = 0.95,
                                   benchmark = TRUE,
                                   time.rep = 10){
              
##########################
##### Input checking #####
##########################
              
                     x <- as.numeric(x)
              
                     if(any(!is.finite(x)))
                   stop("The sample contains non-finite observations.")
              
                     if(any(x <= 0))
                   stop("All observations must be positive.")
              
                     if(length(x) < 2)
                   stop("Sample size must be at least two.")
              
                     n <- length(x)
              
############################
##### Point estimation #####
############################
              
                fit.NR <- tryCatch(
                NR_MEG2(x = x,
                        eta0 = eta0,
                        nu0 = nu0,
                        tol = tol,
                        maxit = max_iter
                       ),
                 error = function(e) NULL 
                       )
              
              fit.NGEM <- tryCatch(
            NGEM_eta_nu(x = x,
                        eta0 = eta0,
                        nu0 = nu0,
                        tol = tol,
                        max_iter = max_iter
                       ),
                 error = function(e) NULL
                       )
              
############################
##### Check estimation #####
############################
              
                     if(is.null(fit.NR))
                   stop("Newton-Raphson estimation failed.")
              
                     if(is.null(fit.NGEM))
                   stop("NGEM estimation failed.")
              
###################################
##### Computational benchmark #####
###################################
              
              if(benchmark){
                   time.NR <- numeric(time.rep)
                 time.NGEM <- numeric(time.rep)
for(i in seq_len(time.rep)){
                time.NR[i] <- tryCatch(
                    NR_MEG2(x = x,
                            eta0 = eta0,
                            nu0 = nu0,
                            tol = tol,
                            maxit = max_iter
                           )$time,
                     error = function(e) NA_real_
                           )
                  
              time.NGEM[i] <- tryCatch(
                NGEM_eta_nu(x = x,
                            eta0 = eta0,
                            nu0 = nu0,
                            tol = tol,
                            max_iter = max_iter
                           )$time,
                     error = function(e) NA_real_
                           )
                           }
                
              mean.time.NR <- mean(time.NR, na.rm = TRUE)
                sd.time.NR <- sd(time.NR, na.rm = TRUE)
                
            mean.time.NGEM <- mean(time.NGEM, na.rm = TRUE)
              sd.time.NGEM <- sd(time.NGEM, na.rm = TRUE)
                           }else{
              mean.time.NR <- fit.NR$time
                sd.time.NR <- NA_real_
                
            mean.time.NGEM <- fit.NGEM$time
              sd.time.NGEM <- NA_real_
                           }
              
#############################
##### Convergence flags #####
#############################
              
                     NR.ok <- isTRUE(fit.NR$converged)
                   NGEM.ok <- isTRUE(fit.NGEM$converged)
              
#####################################
##### Wald confidence intervals #####
#####################################
              
                  if(NR.ok){
                   Wald.NR <- WaldIntervals.MEG2(fit = fit.NR, x = x, level = level)
                           }else{
                   Wald.NR <- NULL
                           }
              
                if(NGEM.ok){
                 Wald.NGEM <- WaldIntervals.MEG2(fit = fit.NGEM, x = x, level = level)
                           }else{
                 Wald.NGEM <- NULL
                           }
              
#################################################
##### Likelihood-ratio confidence intervals #####
#################################################.  
              
                  if(NR.ok){
                     LR.NR <- LRIntervals.MEG2(fit = fit.NR, 
                                               x = x, 
                                               level = level, 
                                               algorithm = "NR",
                                               tol = tol,
                                               max_iter = max_iter)
                           }else{
                     LR.NR <- NULL
                           }
              
                if(NGEM.ok){
                   LR.NGEM <- LRIntervals.MEG2(fit = fit.NGEM,
                                               x = x,
                                               level = level,
                                               algorithm = "NGEM",
                                               tol = tol,
                                               max_iter = max_iter)
                
                           }else{
                   LR.NGEM <- NULL
                           }
              
####################################
##### Table 1: Point estimates #####
####################################
              
                    Table1 <- data.frame(Parameter = c("alpha","beta"),
                                         `Newton-Raphson` = round(
                                         c(if(NR.ok) fit.NR$alpha_hat else NA_real_,
                                           if(NR.ok) fit.NR$beta_hat else NA_real_), 4),
                                         NGEM = round(
                  
                                         c(if(NGEM.ok) fit.NGEM$alpha_hat else NA_real_, 
                                           if(NGEM.ok) fit.NGEM$beta_hat else NA_real_), 4),
                 
               check.names = FALSE
                           )
              
#########################################
##### Table 2: Confidence intervals #####
#########################################
              
                    Table2 <- data.frame(Parameter = c("alpha","beta"),
                                         NR_Wald = c(
                                         if(NR.ok)
                                         sprintf("(%.4f, %.4f)",
                                         Wald.NR$CI.alpha[1],
                                         Wald.NR$CI.alpha[2])
                                         else "(NA, NA)",
                  
                                         if(NR.ok)
                                         sprintf("(%.4f, %.4f)",
                                         Wald.NR$CI.beta[1],
                                         Wald.NR$CI.beta[2])
                                         else "(NA, NA)"),
                
                                         NR_LR = c(
                                         if(NR.ok)
                                         sprintf("(%.4f, %.4f)",
                                         LR.NR$CI.alpha[1],
                                         LR.NR$CI.alpha[2])
                                         else "(NA, NA)",
                                         if(NR.ok)
                                         sprintf("(%.4f, %.4f)",
                                         LR.NR$CI.beta[1],
                                         LR.NR$CI.beta[2])
                                         else "(NA, NA)"),
                
                                         NGEM_Wald = c(
                                         if(NGEM.ok)
                                         sprintf("(%.4f, %.4f)",
                                         Wald.NGEM$CI.alpha[1],
                                         Wald.NGEM$CI.alpha[2])
                                         else "(NA, NA)",
                                         if(NGEM.ok)
                                         sprintf("(%.4f, %.4f)",
                                         Wald.NGEM$CI.beta[1],
                                         Wald.NGEM$CI.beta[2])
                                         else "(NA, NA)"),
                
                                         NGEM_LR = c(
                                         if(NGEM.ok)
                                         sprintf("(%.4f, %.4f)",
                                         LR.NGEM$CI.alpha[1],
                                         LR.NGEM$CI.alpha[2])
                                         else "(NA, NA)",
                                         if(NGEM.ok)
                                         sprintf("(%.4f, %.4f)",
                                         LR.NGEM$CI.beta[1],
                                         LR.NGEM$CI.beta[2])
                                         else "(NA, NA)"),
                
               check.names = FALSE
                           )
              
##########################################
##### Table 3: Algorithm performance #####
##########################################
              
                    Table3 <- data.frame(Algorithm = c("Newton-Raphson", "NGEM"),
                                         Converged = c(NR.ok, NGEM.ok),
                                         Iterations = c(
                                         if(NR.ok) fit.NR$iter else NA_integer_,
                                         if(NGEM.ok) fit.NGEM$iter else NA_integer_),
                
                                         Mean_Time = round(c(mean.time.NR, mean.time.NGEM), 4),
                                         SD_Time = round(
                                         c(sd.time.NR, sd.time.NGEM), 4),
                
                                         Penalized_LogLik = round(c(
                                         if(NR.ok) fit.NR$logLik else NA_real_,
                                         if(NGEM.ok) fit.NGEM$logLik else NA_real_), 4),
                
               check.names = FALSE
                           )
              
#########################
##### Output object #####
#########################
              
                       out <- list(
                
###############################
##### General information #####
###############################
                
               Sample.Size = n,
                      Call = match.call(),
                
##########################
##### Summary tables #####
##########################
                
           Point.Estimates = Table1,
      Confidence.Intervals = Table2,
     Algorithm.Performance = Table3,
                
############################
##### Point estimation #####
############################
                
                        NR = fit.NR,
                      NGEM = fit.NGEM,
                
#####################################
##### Wald confidence intervals #####
#####################################
                
                      Wald = list(NR = Wald.NR, NGEM = Wald.NGEM),
                
#################################################
##### Likelihood-ratio confidence intervals #####
#################################################
                
                        LR = list(NR = LR.NR, NGEM = LR.NGEM))
              
#################
##### Class #####
#################
              
               class(out) <- "fit.MEG2"
                    return(out)
                          }
            
            

################################################################################
################################################################################
################################################################################
############################ FIT MEG2 FOR REAL DATA ############################
################################################################################
################################################################################
################################################################################
            
              FIT.MEG2.RD <- function(x, runs.time = 50){
              
#######################
##### CHECK INPUT #####
#######################
              
                        if(missing(x))
                      stop("'x' must be supplied.")
              
                        if(!is.numeric(x))
                      stop("'x' must be numeric.")
              
                        if(any(is.na(x)))
                      stop("'x' contains missing values.")
              
                        if(any(x <= 0))
                      stop("All observations must be positive.")
               
                        x <- as.numeric(x)
              
                        n <- length(x)
              
####################################
##### STANDARD INITIALIZATIONS #####
####################################
              
                   starts <- list(list(eta0 = log(0.2), nu0  = log(0.2)),
                                  list(eta0 = log(2), nu0  = log(1/mean(x))),
                                  list(eta0 = log(7), nu0  = log(3)))
              
#################################
##### STARTING VALUES TABLE #####
#################################
              
              Start.Table <- data.frame(Start = 1:3,
                                        eta0 = sapply(starts,function(z) z$eta0),
                                        nu0 = sapply(starts,function(z) z$nu0),
                                        alpha0 = round(sapply(starts,
                                                       function(z)
                                                       alpha_from_eta(z$eta0)),
                                                       4),
                
                                        beta0 = round(sapply(
                                                      starts,
                                                      function(z)
                                                      beta_from_nu(z$nu0)),
                                                      4))
              
######################################
##### AVERAGE COMPUTATIONAL TIME #####
######################################
              
              AverageTime <- function(FUN,eta0, nu0, runs = runs.time){
                    TIMES <- numeric(runs)
   for(i in seq_len(runs)){
                       t0 <- proc.time()[3]
                 invisible(
                       FUN(x = x,eta0 = eta0, nu0 = nu0)
                          )
                 TIMES[i] <- proc.time()[3]-t0
                          }
                         c(Mean = mean(TIMES), SD = sd(TIMES))
                          }
              
################################
##### GENERIC FIT FUNCTION #####
################################
              
              RunAlgorithm <- function(FUN, Algorithm){
                      Fits <- vector("list", length(starts))
                   Summary <- data.frame()
for(i in seq_along(starts)){
                         s <- starts[[i]]
                  
###############
##### Fit #####
###############
                  
                       fit <- FUN(x = x, eta0 = s$eta0, nu0 = s$nu0)
                 Fits[[i]] <- fit
                  
######################################
##### Average computational time #####
######################################
                  
                        tm <- AverageTime(FUN = FUN, eta0 = s$eta0, nu0 = s$nu0)
                  
#########################
##### Summary table #####
#########################
                  
                   Summary <- rbind(Summary,
                                    data.frame(Algorithm = Algorithm,
                                    Start = i,
                                    alpha0 = round(alpha_from_eta(s$eta0),4),
                                    beta0 = round(beta_from_nu(s$nu0),4),
                                    Converged = fit$converged,
                                    Iterations = fit$iter,
                                    Mean_Time_ms = round(1000*tm["Mean"],4),
                                   SD_Time_ms = round(1000*tm["SD"],4),
                                   Time_per_Iteration_ms =
                                   if(fit$converged)
                                   round(1000*tm["Mean"]/fit$iter,6)
                                   else
                                   NA,
                                   Penalized_LogLik = round(fit$logLik,6)))
                           }
                
                   rownames(Summary) <- NULL
                       list(Fits = Fits, Summary = Summary)
                           }
              
##############################
##### RUN NEWTON-RAPHSON #####
##############################
              
                    NR.out <- RunAlgorithm(FUN = NR_MEG2, Algorithm = "Newton-Raphson")
              
####################
##### RUN NGEM #####
####################
              
                  NGEM.out <- RunAlgorithm(FUN = NGEM_eta_nu, Algorithm = "NGEM")
              
###########################
##### BENCHMARK TABLE #####
###########################
              
           Benchmark.Table <- rbind(NR.out$Summary, NGEM.out$Summary)
                        
############################
##### COMPARISON TABLE #####
############################
              
                  NR.Table <- NR.out$Summary
                NGEM.Table <- NGEM.out$Summary
          Comparison.Table <- data.frame(Start = 1:length(starts),
                                         alpha0 = NR.Table$alpha0,
                                         beta0 = NR.Table$beta0,
                                         NR.Conv = NR.Table$Converged,
                                         NGEM.Conv = NGEM.Table$Converged,
                                         NR.Iter = NR.Table$Iterations,
                                         NGEM.Iter = NGEM.Table$Iterations,
                                         NR.Time = NR.Table$Mean_Time_ms,
                                         NGEM.Time = NGEM.Table$Mean_Time_ms,
                                         NR.LogLik = NR.Table$Penalized_LogLik,
                                         NGEM.LogLik = NGEM.Table$Penalized_LogLik)
              
###############################
##### PARAMETER ESTIMATES #####
###############################
              
              NR.estimates <- do.call(rbind,
                               lapply(seq_along(starts), function(i){
                       fit <- NR.out$Fits[[i]]
                 data.frame(Start = i,
                            alpha0 = round(alpha_from_eta(starts[[i]]$eta0),4),
                            beta0 = round(beta_from_nu(starts[[i]]$nu0),4),
                            Converged = fit$converged,
                            alpha_hat =
                            if(fit$converged)
                            round(fit$alpha_hat,6)
                        else
                            NA,
                            beta_hat =
                            if(fit$converged)
                            round(fit$beta_hat,6)
                        else
                            NA,
                            LogLik =
                            if(fit$converged)
                            round(fit$logLik,6)
                        else
                           NA)
                           }))
              
                   rownames(NR.estimates) <- NULL
              
###############################
########## NGEM ESTIMATES #####
###############################
              
            NGEM.estimates <- do.call(rbind,
                               lapply(seq_along(starts), function(i){
                       fit <- NGEM.out$Fits[[i]]
                 data.frame(Start = i,
                            alpha0 = round(alpha_from_eta(starts[[i]]$eta0),4),
                            beta0 = round(beta_from_nu(starts[[i]]$nu0),4),
                            Converged = fit$converged,
                            alpha_hat =
                            if(fit$converged)
                            round(fit$alpha_hat,6)
                        else
                            NA,
                            beta_hat =
                            if(fit$converged)
                            round(fit$beta_hat,6)
                        else
                            NA,
                            LogLik =
                            if(fit$converged)
                            round(fit$logLik,6)
                        else
                            NA
                           )}))
                
                   rownames(NGEM.estimates) <- NULL
              
#########################################
##### PREPARE LOG-LIKELIHOOD TRACES #####
#########################################
              
              PrepareTrace <- function(Fits){
                     valid <- which(sapply(
                              Fits,
                   function(f)
                     isTRUE(f$converged) &&
                   !is.null(f$ll_trace) &&
                     length(f$ll_trace) > 0))
                     Trace <- do.call(rbind,
                     lapply(valid,function(i){
                       fit <- Fits[[i]]
                 data.frame(Iteration = seq_along(fit$ll_trace),
                            LogLik = fit$ll_trace,
                            Start = factor(i))}))
                     Final <- do.call(rbind,
                              lapply(valid,function(i){
                       fit <- Fits[[i]]
                 data.frame(Iteration = length(fit$ll_trace),
                            LogLik = tail(fit$ll_trace,1),
                            Start = factor(i))}))
                       list(valid = valid,
                            trace = Trace,
                            final = Final)}
              
########################
##### BUILD TRACES #####
########################
              
                  NR.Trace <- PrepareTrace(NR.out$Fits)
                NGEM.Trace <- PrepareTrace(NGEM.out$Fits)
              
#########################
##### LEGEND LABELS #####
#########################
              
              LegendLabels <- function(valid){
                     lapply(valid,
                function(i){
                     bquote(alpha[0]==.(round(alpha_from_eta(starts[[i]]$eta0),3))
                            ~ ","
                            ~
                            beta[0]==.(round(beta_from_nu(starts[[i]]$nu0),3)))
                           }
                           )
                           }
              
                 NR.labels <- LegendLabels(NR.Trace$valid)
               NGEM.labels <- LegendLabels(
                 NGEM.Trace$valid
                           )
              
###############################
##### NEWTON-RAPHSON PLOT #####
###############################
              
                     G.NR <- ggplot(NR.Trace$trace,
                       aes(x = Iteration,
                           y = LogLik,
                           colour = Start,
                           group = Start)) +
                 geom_line(linewidth = 1.2) +
                geom_point(size = 2.2) +
                geom_point(data = NR.Trace$final,
                           shape = 21,
                           fill = "white",
                           colour = "black",
                           stroke = 0.9,
                           size = 3.2) +
       scale_colour_manual(values = c("black","red2","chartreuse3"), labels = NR.labels, name = NULL) +
                      labs(title = "Evolution of the penalized log-likelihood across iterations (Newton-Raphson)",
                           x = "Iteration",
                           y = "Penalized log-likelihood") +
             theme_minimal(base_size = 14) +
                     theme(legend.position = "top",
                           legend.direction = "horizontal",
                           legend.text = element_text(size = 10),
                           plot.title = element_text(
                           hjust = .5,
                           face = "bold",
                           size = 12),      
               axis.title = element_text(face = "bold", size = 13),
                                         axis.text = element_text(size = 11),
                                         panel.grid.major = element_line(
                                         colour = "grey85",
                                         linewidth = .5),
         panel.grid.minor = element_line(colour = "grey92",
                                          linewidth = .3))
              
#####################
##### NGEM PLOT #####
#####################
              
          G.NGEM <- ggplot(NGEM.Trace$trace,
                       aes(x = Iteration,
                           y = LogLik,
                           colour = Start,
                           group = Start)) +
                 geom_line(linewidth = 1.2) +
                geom_point(size = 2.2) +
                geom_point(data = NGEM.Trace$final,
                           shape = 21,
                           fill = "white",
                           colour = "black",
                           stroke = 0.9,
                           size = 3.2) +
      scale_colour_manual(values = c("black","red2","chartreuse3"),
                          labels = NGEM.labels,
                          name = NULL) +
                     labs(title = "Evolution of the penalized log-likelihood across iterations (NGEM)", 
                          x = "Iteration",
                          y = "Penalized log-likelihood") +
                          theme_minimal(base_size = 14) +
                    theme(legend.position = "top", 
                          legend.direction = "horizontal",
                          legend.text = element_text(size = 10),
                          plot.title = element_text(
                          hjust = .5,
                          face = "bold",
                          size = 12),
              axis.title = element_text(face = "bold", size = 13),
               axis.text = element_text(size = 11),
        panel.grid.major = element_line(colour = "grey85", linewidth = .5),
        panel.grid.minor = element_line(colour = "grey92", linewidth = .3))
              
###############################
##### BUILD OUTPUT OBJECT #####
###############################
              
                 Results <- list(
                
####################################
##### Standard initializations #####
####################################
                
             Start.Table = Start.Table,
                
#############################
##### Benchmark results #####
#############################
                
               Benchmark = Benchmark.Table,
                
############################
##### Comparison table #####
############################
                
              Comparison = Comparison.Table,
                
###############################
##### Parameter estimates #####
###############################
                
            NR.estimates = NR.estimates,
          NGEM.estimates = NGEM.estimates,
                
#############################
##### Convergence plots #####
#############################
                
                    G.NR = G.NR,
                  G.NGEM = G.NGEM)

                   return(Results)
                         }

################################################################################
################################################################################
################################################################################
############################ COMPARIOSN ########################################
####### MEG2, Exponential, Gamma, Weibull, Lindley, XGamma, Log-Normal #########
################################################################################
################################################################################
################################################################################
            
              Compare.MEG <- function(x, Top = 3){
              
#######################
##### CHECK INPUT #####
#######################
              
                        if(missing(x))
                      stop("'x' must be supplied.")
              
                        if(!is.numeric(x))
                      stop("'x' must be numeric.")
              
                        if(any(is.na(x)))
                      stop("'x' contains missing values.")
              
                        if(any(x <= 0))
                      stop("All observations must be positive.")
              
                        if(length(x) < 5)
                      stop("At least five observations are required.")
              
                        x <- as.numeric(x)
                        n <- length(x)
               
                      Top <- max(1L, as.integer(Top))
                      Top <- min(Top, 7L)
              
###############################
##### GRID USED FOR PLOTS #####
###############################
              
                       xx <- seq(from = 0, to = max(x), length.out = 500)
              
################################
##### GENERIC GOF FUNCTION #####
################################
              
             GOF.measures <- function(logLik, npar, cdf.fun){
                       KS <- ks.test(x, cdf.fun)
                       W2 <- goftest::cvm.test(x, null = cdf.fun)
                data.frame(LogLik = logLik,
                           Deviance = -2*logLik,
                           AIC = -2*logLik + 2*npar,
                           BIC = -2*logLik + log(n)*npar,
                           KS = unname(KS$statistic),
                           KS.pvalue = KS$p.value,
                           W2 = unname(W2$statistic),
                           W2.pvalue = W2$p.value,
                           row.names = NULL)}
              
#######################################################
##### GENERIC MLE FOR ONE-PARAMETER DISTRIBUTIONS #####
#######################################################
                     
         Fit.OneParameter <- function(x, density, start = 1){
                negloglik <- function(theta){
                        if(theta <= 0)
                    return(Inf)
                        f <- density(x, theta)
                        if(any(!is.finite(f)))
                    return(Inf)
                  
                        if(any(f <= 0))
                    return(Inf)
                  
                      -sum(log(f))
                          }
                
                      fit <- optim(par = start, fn = negloglik, method = "L-BFGS-B", lower = 1e-8)
                
                      list(theta = fit$par,
                           logLik = -fit$value,
                           convergence = fit$convergence)
                          }
                     
#########################
##### COLOR PALETTE #####
#########################
              
             Model.Colors <- c(Exponential = "black", 
                               Gamma = "red3",
                               Weibull = "blue3",
                               Lognormal = "forestgreen",
                               Lindley = "brown4",
                               XGamma = "purple",
                               MEG2 = "orange3")
              
###########################
##### FIT EXPONENTIAL #####
###########################
              
                  fit.exp <- fitdistrplus::fitdist(x,"exp")
              
#####################
##### FIT GAMMA #####
#####################
              
                fit.gamma <- fitdistrplus::fitdist(x, "gamma")
              
#######################
##### FIT WEIBULL #####
#######################
              
              fit.weibull <- fitdistrplus::fitdist(x, "weibull")
              
#########################
##### FIT LOGNORMAL #####
#########################
              
                fit.lnorm <- fitdistrplus::fitdist(x, "lnorm")
              
#######################
##### FIT LINDLEY #####
#######################
              
              fit.lindley <- Fit.Lindley(x)
              
######################
##### FIT XGAMMA #####
######################
              
               fit.xgamma <- Fit.XGamma(x)
              
####################
##### FIT MEG2 #####
####################
              
                fit.MEG2a <- fit.MEG2(x = x, eta0 = log(2), nu0 = log(1/mean(x)))
              
##################################
##### EXTRACT NGEM ESTIMATES #####
##################################
              
                alpha.hat <- fit.MEG2a$NGEM$alpha_hat
                 beta.hat <- fit.MEG2a$NGEM$beta_hat
              
#############################################
##### CUMULATIVE DISTRIBUTION FUNCTIONS #####
#############################################
              
                  CDF.exp <- function(z){pexp(z, rate = fit.exp$estimate)}
                CDF.gamma <- function(z){pgamma(z, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])}
              CDF.weibull <- function(z){pweibull(z, shape = fit.weibull$estimate["shape"], scale = fit.weibull$estimate["scale"])}
                CDF.lnorm <- function(z){plnorm(z, meanlog = fit.lnorm$estimate["meanlog"], sdlog = fit.lnorm$estimate["sdlog"])}
              CDF.lindley <- function(z){pLindley(z,theta = fit.lindley$theta)}
               CDF.xgamma <- function(z){pXGamma(z, theta = fit.xgamma$theta)}
                 CDF.MEG2 <- function(z){pMEG2(z, alpha = alpha.hat, beta = beta.hat)}
              
#####################
##### DENSITIES #####
#####################
              
                  Density <- data.frame(x = xx, 
                                        Exponential = dexp(xx, rate = fit.exp$estimate),
                                        Gamma = dgamma(xx, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"]),
                                        Weibull = dweibull(xx, shape = fit.weibull$estimate["shape"], scale = fit.weibull$estimate["scale"]),
                                        Lognormal = dlnorm(xx, meanlog = fit.lnorm$estimate["meanlog"], sdlog = fit.lnorm$estimate["sdlog"]),
                                        Lindley = dLindley(xx, theta = fit.lindley$theta),
                                        XGamma = dXGamma(xx, theta = fit.xgamma$theta),
                                        MEG2 = dMEG2(xx, alpha = alpha.hat, beta = beta.hat))
              
##############################
##### SURVIVAL FUNCTIONS #####
##############################
              
                 Survival <- data.frame(time = xx,
                                        Exponential = pexp(xx, rate = fit.exp$estimate, lower.tail = FALSE),
                                        Gamma = pgamma(xx, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"], lower.tail = FALSE),
                                        Weibull = pweibull(xx, shape = fit.weibull$estimate["shape"], scale = fit.weibull$estimate["scale"], lower.tail = FALSE),
                                        Lognormal = plnorm(xx, meanlog = fit.lnorm$estimate["meanlog"], sdlog = fit.lnorm$estimate["sdlog"], lower.tail = FALSE),
                                        Lindley = pLindley(xx, theta = fit.lindley$theta, lower.tail = FALSE ),
                                        XGamma = pXGamma(xx, theta = fit.xgamma$theta, lower.tail = FALSE),
                                        MEG2 = 1 - pMEG2(xx, alpha = alpha.hat, beta = beta.hat))
              
####################################
##### GOODNESS-OF-FIT MEASURES #####
####################################
              
                  GOF.exp <- GOF.measures(logLik = fit.exp$loglik, npar = 1, cdf.fun = CDF.exp)
                GOF.gamma <- GOF.measures(logLik = fit.gamma$loglik, npar = 2, cdf.fun = CDF.gamma)
              GOF.weibull <- GOF.measures(logLik = fit.weibull$loglik, npar = 2, cdf.fun = CDF.weibull)
                GOF.lnorm <- GOF.measures(logLik = fit.lnorm$loglik, npar = 2, cdf.fun = CDF.lnorm)
              GOF.lindley <- GOF.measures(logLik = fit.lindley$logLik, npar = 1, cdf.fun = CDF.lindley)
               GOF.xgamma <- GOF.measures(logLik = fit.xgamma$logLik, npar = 1, cdf.fun = CDF.xgamma)
                 GOF.MEG2 <- GOF.measures(logLik = fit.MEG2a$NGEM$logLik, npar = 2, cdf.fun = CDF.MEG2)
              
#####################
##### GOF TABLE #####
#####################
              
                GOF.Table <- rbind(Exponential = GOF.exp,
                                   Gamma = GOF.gamma,
                                   Weibull = GOF.weibull,
                                   Lognormal = GOF.lnorm,
                                   Lindley = GOF.lindley,
                                   XGamma = GOF.xgamma,
                                   MEG2 = GOF.MEG2)
              
                GOF.Table <- round(GOF.Table, 4)
              
#########################
##### MODEL RANKING #####
#########################
              
                  Ranking <- data.frame(Model = rownames(GOF.Table), Rank.AIC = rank(GOF.Table$AIC, ties.method = "min"),
                                        Rank.BIC = rank(GOF.Table$BIC, ties.method = "min"),
                                        Rank.KS = rank(GOF.Table$KS, ties.method = "min"),
                                        Rank.W2 = rank(GOF.Table$W2, ties.method = "min"),
                                        row.names = NULL)
              
            Ranking$Total <- Ranking$Rank.AIC +
                            Ranking$Rank.BIC +
                            Ranking$Rank.KS +
                            Ranking$Rank.W2
               
                  Ranking <- Ranking[order(Ranking$Total),]
                  rownames(Ranking) <- NULL
             Ranking$Rank <- seq_len(nrow(Ranking))
              
######################
##### TOP MODELS #####
######################
              
               Top.Models <- Ranking[1:Top, ]
              Best.Models <- Top.Models$Model
              
####################################
##### DATA USED IN THE FIGURES #####
####################################
              
              Density.Top <- Density[,c("x", Best.Models)]
             Survival.Top <- Survival[, c("time", Best.Models)]
              
########################################
##### HISTOGRAM + FITTED DENSITIES #####
########################################
              
             Density.Long <- tidyr::pivot_longer(Density.Top, cols = -x, names_to = "Model", values_to = "Density")
               
                Histogram <- ggplot() +
                             geom_histogram(
                             data = data.frame(x = x),
                             aes(x = x, y = after_stat(density)),
                             bins = 30,
                             fill = "grey90",
                             colour = "black",
                             linewidth = 0.4) +
                             geom_line(
                             data = Density.Long,
                             aes(x = x, y = Density, colour = Model),
                             linewidth = 1.3) +
                             scale_colour_manual(values = Model.Colors[Best.Models]) +
                             labs(title = paste0("Histogram with fitted densities (Top ", Top, " models)"),
                             x = "Time",
                             y = "Density",
                             colour = "Model") +
                             theme_bw(base_size = 14) +
                             theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                             axis.title = element_text(face = "bold", size = 13),
                             axis.text = element_text(size = 11),
                             legend.position = "top",
                             legend.title = element_blank(),
                             legend.text = element_text(
                             size = 11),
                             panel.grid.major = element_line(colour = "grey85"),
                             panel.grid.minor = element_line(colour = "grey92"))
              
##################################
##### KAPLAN-MEIER ESTIMATOR #####
##################################
               
                       KM <- survival::survfit(survival::Surv(x, rep(1, n)) ~ 1)
              
                  KM.Data <- data.frame(time = KM$time, Survival = KM$surv)
              
            Survival.Long <- tidyr::pivot_longer(Survival.Top, cols = -time, names_to = "Model", values_to = "Survival")
              
           Survival.Graph <- ggplot() +
                             geom_step(data = KM.Data, aes(time, Survival),
                             colour = "black",
                             linewidth = 1.4,
                             linetype = "solid") +
                             geom_line(data = Survival.Long, aes(time, Survival, colour = Model),
                             linewidth = 1.3) +
                             scale_colour_manual(values = Model.Colors[Best.Models]) +
                             labs(title = paste0("Empirical and fitted survival functions (Top ", Top, " models)"),
                             x = "Time",
                             y = "Survival probability",
                             colour = "Model") +
                             theme_bw(base_size = 14) +
                             theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                             axis.title = element_text(face = "bold", size = 13),
                             axis.text = element_text(size = 11),
                             legend.position = "top",
                             legend.title = element_blank(),
                             legend.text = element_text(size = 11),
                             panel.grid.major = element_line(colour = "grey85"),
                             panel.grid.minor = element_line(colour = "grey92"))
              
###############################
##### PARAMETER ESTIMATES #####
###############################
              
                Estimates <- rbind(
                             data.frame(Model = "Exponential",
                             Parameter = "rate",
                             Estimate = fit.exp$estimate),
                             
                             data.frame(Model = "Gamma", 
                                        Parameter = c("shape","rate"), 
                                        Estimate = c(fit.gamma$estimate["shape"], 
                                                     fit.gamma$estimate["rate"])),
                
                             data.frame(Model = "Weibull",
                                        Parameter = c("shape","scale"),
                                        Estimate = c(fit.weibull$estimate["shape"], fit.weibull$estimate["scale"])),
                
                             data.frame(Model = "Lognormal",
                                        Parameter = c("meanlog","sdlog"),
                                        Estimate = c(fit.lnorm$estimate["meanlog"], fit.lnorm$estimate["sdlog"])),
                
                             data.frame(Model = "Lindley", 
                                        Parameter = "theta", 
                                        Estimate = fit.lindley$theta),
                
                             data.frame(Model = "XGamma",
                                        Parameter = "theta",
                                        Estimate = fit.xgamma$theta),
                
                             data.frame(Model = "MEG2",
                                        Parameter = c("alpha","beta"),
                                        Estimate = c(fit.MEG2a$NGEM$alpha_hat, fit.MEG2a$NGEM$beta_hat)))
              
       Estimates$Estimate <- round(Estimates$Estimate,6)
              
###########################################
##### SUMMARY OF THE MODEL COMPARISON #####
###########################################
              
                  Summary <- list(Sample.Size = n,
                             Number.of.Models = nrow(GOF.Table),
                             Top = Top,
                             Best.AIC = rownames(GOF.Table)[which.min(GOF.Table$AIC)],
                             Best.BIC = rownames(GOF.Table)[which.min(GOF.Table$BIC)],
                             Best.KS = rownames(GOF.Table)[which.min(GOF.Table$KS)],
                             Best.W2 = rownames(GOF.Table)[which.min(GOF.Table$W2)],
                             Overall.Best = Top.Models$Model[1])
              
#########################
##### FITTED MODELS #####
#########################
              
                     Fits <- list(Exponential = fit.exp,
                                  Gamma = fit.gamma,
                                  Weibull = fit.weibull,
                                  Lognormal = fit.lnorm,
                                  Lindley = fit.lindley,
                                  XGamma = fit.xgamma,
                                  MEG2 = fit.MEG2a)
              
#################
##### PLOTS #####
#################
              
                    Plots <- list(Histogram = Histogram, Survival = Survival.Graph)
              
##########################
##### RESULTS OBJECT #####
##########################
              
                  Results <- list(Summary = Summary,
                                  Estimates = Estimates,
                                  GOF.Table = GOF.Table,
                                  Ranking = Ranking,
                                  Top.Models = Top.Models,
                                  Fits = Fits,
                                  Plots = Plots)
              
#########################
##### PRINT SUMMARY ##### 
#########################
               
                       cat("\n")
                       cat("=============================================================\n")
                       cat("                 Compare.MEG\n")
                       cat("=============================================================\n\n")
              
                       cat("Sample size          :", Summary$SampleSize, "\n")
                       cat("Models compared      :", Summary$Number.of.Models, "\n")
                       cat("Top models requested :", Summary$Top, "\n\n")
              
                       cat("Best according to AIC :", Summary$Best.AIC, "\n")
                       cat("Best according to BIC :", Summary$Best.BIC, "\n")
                       cat("Best according to KS  :", Summary$Best.KS, "\n")
                       cat("Best according to W²  :", Summary$Best.W2, "\n\n")
              
                       cat("Overall best model :", Summary$Overall.Best, "\n\n")
              
                       cat("Top models:\n")
              
                     print(Top.Models, row.names = FALSE)
              
                       cat("\n")
                
                       cat("Available outputs:\n\n")
              
                       cat("Results$Estimates\n")
                       cat("Results$GOF.Table\n")
                       cat("Results$Ranking\n")
                       cat("Results$Top.Models\n")
                       cat("Results$Fits\n")
                       cat("Results$Plots$Histogram\n")
                       cat("Results$Plots$Survival\n\n")
              
                       cat("=============================================================\n")
               
                    return(Results)
                          }
            
################################
##### LINDLEY DISTRIBUTION #####
################################
            
                 dLindley <- function(x, theta, log = FALSE){
                        f <- (theta^2/(1 + theta)) * (1 + x) * exp(-theta*x)
                        if(log)
                    return(log(f))
                        f }
            
################################################################################
            
                 pLindley <- function(q, theta, lower.tail = TRUE){
                        F <- 1 - (1 + theta*q/(1 + theta)) * exp(-theta*q)
                        if(lower.tail)
                    return(F)
                    1 - F }
            
################################################################################
            
                 qLindley <- function(p, theta){
                    sapply(p,
                           function(pp)
                   uniroot(
                           function(x)
                           pLindley(x,theta)-pp,
                           interval=c(0,1000)
                          )$root)
                          }
            
################################################################################
            
                 rLindley <- function(n, theta){qLindley(runif(n), theta)}
            
#######################
##### FIT LINDLEY #####
#######################
            
              Fit.Lindley <- function(x){
                        n <- length(x)
                   loglik <- function(theta){
                        if(theta <= 0)
                    return(Inf)
                       ll <- 2*n*log(theta) - n*log(1 + theta) + sum(log1p(x)) - theta*sum(x)
                    return(-ll)
                          }
                   theta0 <- max(1e-6, 1/mean(x))
                      fit <- optim(par = theta0, fn = loglik, method = "L-BFGS-B", lower = 1e-8)
                      list(theta = fit$par, logLik = -fit$value, convergence = fit$convergence)}
            
###############################
##### XGAMMA DISTRIBUTION #####
###############################
            
                  dXGamma <- function(x, theta, log = FALSE){
                        f <- theta^2/(1 + theta) * (1 + theta*x^2/2) * exp(-theta*x)
                        if(log)
                    return(log(f))
                       f  }
            
################################################################################
            
                  pXGamma <- function(q, theta, lower.tail = TRUE){
                        F <- 1 - (1 + theta*q + theta^2*q^2/(2*(1 + theta))) * exp(-theta*q)
                        if(lower.tail)
                    return(F)
                    1 - F }
            
################################################################################
            
                  qXGamma <- function(p, theta){
                             sapply(p,
                             function(pp)
                             uniroot(
                             function(x)
                             pXGamma(x,theta)-pp,
                             interval=c(0,1000))$root)}
            
################################################################################
            
                  rXGamma <- function(n, theta){qXGamma(runif(n), theta)}
            
######################
##### FIT XGAMMA #####
######################
            
               Fit.XGamma <- function(x){
                        n <- length(x)
                   loglik <- function(theta){
                        if(theta <= 0)
                    return(Inf)
                       ll <- 2*n*log(theta) - n*log(1 + theta) + sum(log(1 + theta*x^2/2)) -  theta*sum(x)
                        if(!is.finite(ll))
                    return(Inf)
                    return(-ll)
                          }
              
                   theta0 <- max(1e-6, 1/mean(x))
                      fit <- optim(par = theta0, fn = loglik, method = "L-BFGS-B", lower = 1e-8)
                      list(theta = fit$par, logLik = -fit$value, convergence = fit$convergence)
                          }
            
            
            
            

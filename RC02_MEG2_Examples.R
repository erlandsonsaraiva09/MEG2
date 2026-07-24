############################################################
##### Example: MEG2 distribution functions #################
############################################################

# The following example illustrates the use of the MEG2
# distribution functions for alpha = 3 and beta = 0.5.
# Specifically, it demonstrates the evaluation of the
# density, distribution, quantile, and random generation
# functions, as well as a graphical comparison between
# simulated observations and the theoretical density.

############################################################

##### Parameter values

  alpha <- 4
  beta  <- 0.88

######################
##### Density function
######################

      x <- c(1, 2, 5, 10)

   dMEG2(x, alpha = alpha, beta = beta)

##### Log-density

   dMEG2(x, alpha = alpha, beta = beta, log = TRUE)

############################################################
# Cumulative distribution function
############################################################

   pMEG2(x, alpha = alpha, beta = beta)

##### Survival probabilities

   pMEG2(x, alpha = alpha, beta = beta, lower.tail = FALSE)

############################################################
# Quantile function
############################################################

      q <- c(0.025, 0.50, 0.975)

   qMEG2(q, alpha = alpha, beta = beta)

############################################################
# Random sample generation
############################################################

set.seed(19)

##### Default method: mixture representation

     x1 <- rMEG2(n = 10, alpha = alpha, beta = beta)
     x1

###### Inverse transform method

     x2 <- rMEG2(n = 10, alpha = alpha, beta = beta, option = 2)
     x2

##########################################################
##### Verification: empirical histogram and fitted density
##########################################################

      x <- rMEG2(n = 10000, alpha = alpha, beta = beta)

    hist(x, probability = TRUE, breaks = 30, main = "Histogram and fitted MEG2 density", xlab = "x")
   curve(dMEG2(x, alpha = alpha, beta = beta), from = 0, to = max(x), add = TRUE, lwd = 2)

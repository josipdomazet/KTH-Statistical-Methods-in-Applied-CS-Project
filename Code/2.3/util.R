# TASK 2.3
loglikelihood <- function(real.y, sigma, beta, N, T, resampling.function){
  all.results <- c()
  for(i in 1:10){
    all.results <- c(all.results, bootstrap.filter(real.y, sigma, beta, N, T, resampling.function)$logl)
  }
  return(mean(all.results))
}

pmh <- function(y.values, sigma2, beta2, Nparticles, T, Niteration, step.size, resampling.function){
  
  x.hat <- matrix(0, nrow = T, ncol = Niteration)
  x.hat.proposed <- matrix(0, nrow = T, ncol = Niteration)
  
  theta <- matrix(0, nrow = 2, ncol = Niteration)
  theta.proposed <- matrix(0, nrow = 2, ncol = Niteration)
  
  # log.likelihood is NOT stored, we only keep it for the current time t
  log.likelihood <- matrix(0, nrow = 1, ncol = Niteration)
  log.likelihood.proposed <- matrix(0, nrow = 1, ncol = Niteration)
  # same for acceptance, it is only kept for time = t
  proposed.accepted <- matrix(0, nrow = 1, ncol = Niteration)
  
  # initial values
  theta[1, 1] <- sigma2
  theta[2, 1] <- beta2
  
  results0 <- bootstrap.filter(y.values, sigma=sqrt(theta[1, 1]), beta=sqrt(theta[2, 1]), Nparticles, T, resampling.function)
  log.likelihood[1, 1] <- results0$logl 

  
  for(k in 2:Niteration){
    
  
    theta.proposed[1,k]<- abs(rnorm(1, theta[1,k-1], step.size[1])) # squared Normal distribution
    theta.proposed[2,k]<- abs(rnorm(1, theta[2,k-1], step.size[2]))

    

    results <- bootstrap.filter(y.values, sigma=sqrt(theta.proposed[1, k]), beta=sqrt(theta.proposed[2, k]), Nparticles, T, resampling.function)
    log.likelihood.proposed[1, k] <- results$logl
  
    
    prior.sigma <- invgamma::dinvgamma(theta.proposed[1, k], shape = 0.01, rate = 0.01, log = T)
    diff.sigma <- prior.sigma - invgamma::dinvgamma(theta[1, k-1], shape = 0.01, rate = 0.01, log = T)

    prior.beta <- invgamma::dinvgamma(theta.proposed[2, k], shape = 0.01, rate = 0.01, log = T)
    diff.beta <- prior.beta - invgamma::dinvgamma(theta[2, k-1], shape = 0.01, rate= 0.01, log = T)
    prior.diff.sum <- diff.sigma + diff.beta
    
    likelihood.dif <- log.likelihood.proposed[1, k] - log.likelihood[1, k-1]
    accept.probability <- exp(prior.diff.sum + likelihood.dif)
  
    
    uniform <- runif(1)
    if(uniform < accept.probability){ # accept new parameters
      theta[, k] <- theta.proposed[, k] 
      log.likelihood[1, k] <- log.likelihood.proposed[1, k]
      #  x.hat[, k] <- x.hat.proposed[, k]
      proposed.accepted[1, k] <- 1
    } else {
      theta[, k] <- theta[, k-1] # note! it is not PROPOSED theta
      log.likelihood[1, k] <- log.likelihood[1, k-1]
      # x.hat[, k] <- x.hat[, k-1]
      proposed.accepted[1, k] <- 0
    }
  }
  # theta is matrix 
  return(list(theta = theta, proposed.accepted = proposed.accepted))
  
}
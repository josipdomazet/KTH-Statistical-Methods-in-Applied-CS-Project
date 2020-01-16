# TASK 2.3 Question 9
loglikelihood <- function(real.y, sigma, beta, N, T, resampling.function){
  all.results <- c()
  for(i in 1:10){
    all.results <- c(all.results, bootstrap.filter(real.y, sigma, beta, N, T, resampling.function)$logl)
  }
  return(mean(all.results))
}


# TASK 2.3 Question 10
pmh <- function(y.values, sigma, beta, Nparticles, TIME, Niteration, step.size, resampling.function){
  
  theta <- matrix(0, nrow = 2, ncol = Niteration)
  theta.proposed <- matrix(0, nrow = 2, ncol = Niteration)

  log.likelihood <- matrix(0, nrow = 1, ncol = Niteration)
  log.likelihood.proposed <- matrix(0, nrow = 1, ncol = Niteration)
  # same for acceptance, it is only kept for time = t
  proposed.accepted <- matrix(0, nrow = 1, ncol = Niteration)
  
  # initial values
  theta[1, 1] <- sigma
  theta[2, 1] <- beta
  
  results0 <- bootstrap.filter(y.values, sigma=sqrt(theta[1, 1]), beta=sqrt(theta[2, 1]), Nparticles, TIME, resampling.function)
  log.likelihood[1, 1] <- results0$logl 

  
  for(k in 2:Niteration){

    theta.proposed[1,k]<- max(0.0001, rnorm(1, theta[1,k-1], step.size[1]))
    theta.proposed[2,k]<- max(0.01, rnorm(1, theta[2,k-1], step.size[2]))
    
    results <- bootstrap.filter(y.values, sigma=sqrt(theta.proposed[1, k]),
                                beta=sqrt(theta.proposed[2, k]), Nparticles, TIME, resampling.function)
    log.likelihood.proposed[1, k] <- results$logl
    
    
    prior.sigma <- nimble::dinvgamma(theta.proposed[1, k], shape = 0.01, scale = 0.01, log = TRUE)
    diff.sigma <- prior.sigma - nimble::dinvgamma(theta[1, k-1], shape = 0.01, scale = 0.01, log = TRUE)

    prior.beta <- nimble::dinvgamma(theta.proposed[2, k], shape = 0.01, scale = 0.01, log = TRUE)
    diff.beta <- prior.beta - nimble::dinvgamma(theta[2, k-1], shape = 0.01, scale = 0.01, log = TRUE)
    prior.diff.sum <- diff.sigma + diff.beta
    
    likelihood.dif <- log.likelihood.proposed[1, k] - log.likelihood[1, k-1]
    accept.probability <- exp(prior.diff.sum + likelihood.dif)

    uniform <- runif(1)
    if(uniform <= accept.probability){ # accept new parameters
      theta[1:2, k] <- theta.proposed[1:2, k] 
      log.likelihood[1, k] <- log.likelihood.proposed[1, k]
      proposed.accepted[1, k] <- 1
    } else {
      theta[1:2, k] <- theta[1:2, k-1] # note! it is not PROPOSED theta
      log.likelihood[1, k] <- log.likelihood[1, k-1]
      proposed.accepted[1, k] <- 0
    }
    
    if(k %% 50 == 0){
      
     cat(sprintf("Current posterior mean for sigma and beta: %f %f \n",mean(theta[1, 1:k]),mean(theta[2, 1:k])))
      
    }
  }
  # theta is matrix 
  return(list(theta = theta, proposed.accepted = proposed.accepted))
  
}
#############################################################################
# Bootstrap particle filter
#############################################################################

bootstrap.filter <- function(real.y, sigma, beta, N, TIME, resampling.function){
  particles <- matrix(0, nrow = TIME, ncol = N)
  x.means <- matrix(0, nrow= TIME, ncol = 1)
  ancestor.indices <- matrix(1, nrow = TIME, ncol = N)
  weights <- matrix(1, nrow = TIME, ncol = N)
  weights.norm <- matrix(1, nrow = TIME, ncol = N)

  particles[1, ] <- rnorm(N, 0, sigma)
  ancestor.indices[1, ] <- 1:N

  weights[1, ] <- dnorm(real.y[1], mean = 0, sd = sqrt(beta^2 * exp(particles[1, ])), log = TIME)
  max.weight <- max(weights[1, ]) # max of log values
  weights[1, ] <- exp(weights[1, ] - max.weight)
  sum_log_weights <- sum(weights[1, ])
  weights.norm[1, ] <- weights[1, ] / sum_log_weights

  log.likelihood <- 0
  log.likelihood <- log.likelihood + max.weight + log(sum_log_weights) - log(N)

  x.means[1, 1] <- sum(weights.norm[1, ] * particles[1, ])

  for(t in 2:TIME){
    new.ancestors <- resampling.function(weights.norm[t-1, ])
    #particles[1:(t-1), ] <- particles[1:(t-1), new.ancestors]
    particles[t, ] <- rnorm(N, mean = 1.0*particles[t-1, new.ancestors], sigma) # using the most important particles

    weights[t, ] <- dnorm(real.y[t], 0, sqrt(beta**2 * exp(particles[t, ])), log = TRUE)
    max.weight <- max(weights[t, ]) # max of log values
    weights[t, ] <- exp(weights[t, ] - max.weight)
    sum_log_weights <- sum(weights[t, ])

    weights.norm[t, ] <- weights[t, ] / sum_log_weights

    x.means[t, 1]  <- sum(weights.norm[t, ] * particles[t, ])

    log.likelihood <- log.likelihood + max.weight +  log(sum_log_weights) - log(N)
  }

  ancestor.index  <- sample(1:N, size=1, prob = weights.norm[TIME, ])
  x.hat.filtered <- particles[, ancestor.index]

  return(list(x.means= x.means, x.hat.filtered = x.hat.filtered, weights.norm = weights.norm, logl = log.likelihood))
}

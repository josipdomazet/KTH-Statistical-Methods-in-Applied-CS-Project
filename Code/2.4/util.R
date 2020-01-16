# TASK 2.4
################################################################################
# Function to perform conditional SMC
################################################################################
csmc <- function(x.values, y.values, sigma, beta, num.particles, TI){

  # create data structures
  particles <- matrix(0, nrow = TI, ncol = num.particles)
  weights <- matrix(0, nrow = TI, ncol = num.particles)
  weights.norm <- matrix(0, nrow = TI , ncol = num.particles)
  ancestor.indices <- matrix(0, nrow = TI, ncol = num.particles)

  particles[1, ] <- rnorm(num.particles, 0, sigma)
  particles[1, num.particles] <- x.values[1] # set the retained particle
  ancestor.indices[1, ] <- 1:num.particles

  weights[1, ] <- dnorm(y.values[1], 0, sd=sqrt(beta^2 * exp(particles[1, ])), log = T)
  max.weight <- max(weights[1, ])
  weights[1, ] <- exp(weights[1, ] - max.weight)
  weights.norm[1, ] <- weights[1, ] / sum(weights[1, ])
  log.likelihood <- 0
  log.likelihood <- log.likelihood + max.weight + log(sum(weights[1, ])) - log(num.particles)

  for(t in 2:TI){

    new.ancestors <- multinomial.resampling(weights.norm[t-1, ])
    new.ancestors[num.particles] <- num.particles

    ancestor.indices[t,] <- ancestor.indices[t-1, new.ancestors]

    particles[t, ] <- rnorm(num.particles, particles[t-1, ancestor.indices[t,]], sigma)
    particles[t, num.particles] <-  x.values[t] # replace last one
    weights[t, ] <- dnorm(y.values[t], 0, sd = sqrt(beta^2 * exp(particles[t, ])), log = T)

    max.weight <- max(weights[t, ])
    weights[t, ] <- exp(weights[t, ] - max.weight)
    weights.norm[t, ] <- weights[t, ] / sum(weights[t, ])

    log.likelihood <- log.likelihood + max.weight + log(sum(weights[t, ])) - log(num.particles)
  }

    index <- sample(1:num.particles, size = 1, prob = weights.norm[TI, ])
    x.star <- rep(0, TI)
    x.star[TI] <- particles[TI, index]

    for(t in TI:2){ #  trace the ancestral path of particle index
      index <- ancestor.indices[t, index]
      x.star[t-1] <-  particles[t-1, index]

    }
    list(x.star = x.star, ancestors = ancestor.indices, logl = log.likelihood)

}


################################################################################
# Function to perform Gibbs sampling
################################################################################
particle.gibbs <- function(y.values, sigma2, beta2, TI, num.particles, num.iteration){
  current.sigma2 <- sigma2
  current.beta2 <- beta2

  # set x values arbitrarily
  bf <- bootstrap.filter(y.values, sqrt(sigma2), sqrt(beta2), num.particles, TI, multinomial.resampling)
  current.x <- bf$x.means


  csmc.result <- csmc(current.x, y.values, sqrt(current.sigma2), sqrt(current.beta2),
                      num.particles, TI)
  current.x <- csmc.result$x.star


  sigma2seq <- tibble(sigma2 = numeric())
  beta2seq <- tibble(beta2 = numeric())

  sigma2seq  <- add_row(sigma2seq, sigma2 = current.sigma2)
  beta2seq<- add_row(beta2seq, beta2 = current.beta2)

  log.likelihood.iters <- rep(NA, num.iteration)
  log.likelihood.iters[1] <- bf$logl
  for(k in 2:num.iteration){

    current.sigma2 <- nimble::rinvgamma(1, shape = 0.01 + TI/2,
                                scale = 0.01 + 0.5 * sum( ( current.x[2:TI] - current.x[1:(TI-1)] )^2))

    sigma2seq <- add_row(sigma2seq, sigma2 = current.sigma2)

    csmc.result <- csmc(current.x, y.values, sqrt(current.sigma2), sqrt(current.beta2),
                                     num.particles, TI)

    current.x <- csmc.result$x.star
    current.beta2 <- nimble::rinvgamma(1, shape = 0.01 + TI/2,
                                       scale = 0.01 + 0.5 * sum(exp(-current.x)*(y.values^2)))
    beta2seq <- add_row(beta2seq, beta2 = current.beta2)
    log.likelihood.iters[k] <- csmc.result$logl
    if(k %% 150 == 0){
      cat(sprintf("Log-likelihood = %.4f \n", csmc.result$logl))
    }
  }

  list(sigma.density = sigma2seq$sigma2, beta.density = beta2seq$beta2, logl = log.likelihood.iters)

}

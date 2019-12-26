##############################################################################
# Sequential importance sampling filter for stochastic volatility model
# N = number of particles
# T = number of time steps
##############################################################################
sis <- function(real.y, N, T){
  x.pf.sis <- matrix(0, nrow = T, ncol=N) 
  w.pf.sis <- matrix(0, nrow = T, ncol=N)
  w.norm.sis <- matrix(0, nrow = T, ncol=N)
  
  x.pf.sis[1, ] <- rnorm(N, mean = 0, sd = 0.16) # Initial, t = 1
  
  w.pf.sis[1, ] <- dnorm(real.y[1], mean = 0, sd = sqrt(0.64**2 * exp(x.pf.sis[1, ]))) # 
  w.norm.sis[1, ] <- w.pf.sis[1, ] / sum(w.pf.sis[1, ])
  
  x.means.sis <- rep(NA, T) 
  x.means.sis[1] <- sum(w.norm.sis[1, ] * x.pf.sis[1,])
  
  for(t in 2:T){
    x.pf.sis[t, ] <- rnorm(N, mean = 1.00 * x.pf.sis[t-1, ], sd = 0.16)
    w.pf.sis[t, ] <- dnorm(real.y[t], mean = 0, sd = sqrt(0.64^2 * exp(x.pf.sis[t, ]))) * w.pf.sis[t-1, ]
    
    w.norm.sis[t, ] <- w.pf.sis[t, ] / sum(w.pf.sis[t, ])
    x.means.sis[t] <- sum(w.norm.sis[t, ] * x.pf.sis[t, ])
  }
  return(list(x.means = x.means.sis, wnorm = w.norm.sis))
}
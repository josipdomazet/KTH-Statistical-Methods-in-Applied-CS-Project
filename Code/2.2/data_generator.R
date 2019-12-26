#############################################################################
# Function used to generate stochastic-volatility data 
############################################################################## 
generate_data <- function(fi, sigma, beta, total_len){
  library(tibble)
  
  results <- tibble(x = numeric(), y = numeric(), t = numeric())
  x <- rnorm(1, 0, sigma) # x1
  y <- rnorm(1, 0, sqrt(beta^2 * exp(x)))
  results <- tibble::add_row(results, x = x, y = y, t = 1)
  for(i in 2:total_len){
    x <- rnorm(1, fi*x, sigma)
    y <- rnorm(1, 0, sqrt(beta^2 * exp(x)))
    results <- tibble::add_row(results, x = x, y = y, t = i)
  }
  return(results)
}

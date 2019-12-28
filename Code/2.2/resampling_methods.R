#############################################################################
# Function used to perform multinomial resampling
# on the given probabilities in weights.vector
#############################################################################
multinomial.resampling <- function(weights.vector, N = NA){
  # generate N values between 0 and 1
  N <- ifelse(is.na(N), length(weights.vector))
  N <- length(weights.vector)
  u <- rep(NA, N)
  u.tilda <- runif(N)^(1/(1:N))
  u[N] <- u.tilda[N] # set the value for the last element
  k <- N-1
  
  while(k > 0){
    u[k] <- u[k+1] * u.tilda[k]
    k <- k-1
  }
  
  out <- rep(NA, N)
  total <- 0
  i <- 1
  j <- 1
  
  while(j <= N && i <= N){
    total <- total + weights.vector[i]
    while(j <= N && total > u[j]){
      out[j] <- i
      j <- j + 1
    }
    i <- i + 1
  }
  out
} 


#############################################################################
# Function used to perform stratified resampling
# on the given probabilities in weights.vector
#############################################################################
stratified.resampling <- function(weights.vector, N = NA){
  N <- ifelse(is.na(N), length(weights.vector))
  output <- rep(NA, N)
  
  total <- weights.vector[1]
  j <- 1
  
  for(i in 1:N){
    u <- (runif(1) + i-1) / N
    while(total < u){
      j <- j + 1
      total <- total + weights.vector[j]
    }
    output[i] <- j
  }
  
  output
}
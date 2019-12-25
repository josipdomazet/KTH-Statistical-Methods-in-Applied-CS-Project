
# function used to extract magic words from sequences
extract.magic.words <- function(data, current.state, N, W){
  magic.words <- matrix(nrow = N, ncol = W)
  end.positions <- current.state + W - 1
  for(i in 1:N){
  magic.words[i, ] <- as.matrix(data[i, current.state[i]:end.positions[i]])[1, ]
  }
  magic.words
}

extract.backgrounds <- function(data, current.state, N, W){
  backgrounds <- matrix(nrow = N, ncol = M - W)
  end.positions <- current.state + W - 1
  for(i in 1:N){
    backgrounds[i, ] <- as.matrix(data[i, -c(current.state[i]:end.positions[i])])[1, ]
  }
  backgrounds
}


gibbs.sampler <- function(my.data, alphabet, num.iterations, N, M, W, K, alpha.bg, alpha.mw){

  # the states are vectors of start positions
  initial.positions <- sample(1:(M - W + 1), N, replace = T)
  samples <- matrix(nrow = 0, ncol = N)
  samples <- rbind(samples, initial.positions)
  
  # do this only ONCE!
  magic.words <- extract.magic.words(my.data, initial.positions, N, W)
  backgrounds <- extract.backgrounds(my.data, initial.positions, N, W)
  
  
  for(it in 1:num.iterations){
    current.state <- samples[it, ]
    end.positions <- current.state + W - 1
    
    new.starts <- matrix(nrow = 1, ncol = N)
    
    # using current stat
    # precompute counts used for posterior
    
    ## counts for background: all characters except the magic words
    background.count.matrix <- matrix(nrow = N, ncol = K)
    
    for(n in 1:N){
      temp <- table(backgrounds[n, ])[alphabet] # use all, then modify in the loop
      temp[is.na(temp)] <- 0
      background.count.matrix[n, ] <- temp
  
    }
    
    for(n in 1:N){
      
      posterior.tilda <- posterior(my.data[n, ], current.state, n, magic.words, backgrounds,
                                   background.count.matrix)
      posterior.norm <- exp.normalize(posterior.tilda)
      new.start.for.n <- sample(1:(M - W + 1), size = 1, prob=posterior.norm)
      current.state[n] <- new.start.for.n
      magic.words[n, ] <- my.data[n, new.start.for.n:(new.start.for.n+W-1)]
      backgrounds[n, ] <- my.data[n, -c(new.start.for.n:(new.start.for.n+W-1))]
      new.starts[1, n] <- new.start.for.n 
    }
    samples <- rbind(samples, new.starts[1, ])
    
  }
  
  # samples represents STATES for EACH iteration! 
  # it's actually called chain because it represents the chain of STATES 
  # state = vector of starting positions for N sequences
  samples
}

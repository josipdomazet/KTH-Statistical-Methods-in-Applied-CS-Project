##############################################################################
# Function used to generate date according to instructions given in task 2.1
##############################################################################
generate.data <- function(N, alphabet, M, W, alpha.bg, alpha.mw){

  theta.bg <- rdirichlet(1, alpha.bg)
  theta.mw <- rdirichlet(W, alpha.mw) # for each magic word position

  data <- matrix(nrow = N, ncol = M)
  all.starts <- rep(NA, N)
  for(n in 1:N){

    generated.seq <- rep(NA, M)
    mw.start <- sample(1:(M - W + 1), 1)
    all.starts[n] <- mw.start
    mw.end <- mw.start + W - 1

    for(i in 1:M){

      if(i %in% mw.start:mw.end){ # it's a magic word part
        generated.seq[i] <- sample(alphabet, 1, prob = theta.mw[i - mw.start + 1, ])
      } else { # if it's background
        generated.seq[i] <- sample(alphabet, 1, prob = theta.bg)
      }
    }
    data[n, ] <- generated.seq
  }

  list(data = data, starts = all.starts)
}


##############################################################################
# Function used to extract magic words from sequences based on
# the current state - the current state represents starting positions
# of magic word for each of the sequences
##############################################################################
extract.magic.words <- function(data, current.state, N, W){
  magic.words <- matrix(nrow = N, ncol = W)
  end.positions <- current.state + W - 1
  for(i in 1:N){
    magic.words[i, ] <- as.matrix(data[i, current.state[i]:end.positions[i]])[1, ]
  }
  magic.words
}


##############################################################################
# Function used to extract background from sequences based on
# the current state - the current state represents starting positions
# of magic word for each of the sequences
##############################################################################
extract.backgrounds <- function(data, current.state, N, W){
  backgrounds <- matrix(nrow = N, ncol = M - W)
  end.positions <- current.state + W - 1
  for(i in 1:N){
    backgrounds[i, ] <- as.matrix(data[i, -c(current.state[i]:end.positions[i])])[1, ]
  }
  backgrounds
}


##############################################################################
# Function used to perform Gibbs sampling.
# It returns the chain, i.e. the sequence of states from all iterations
##############################################################################
gibbs.sampler <- function(my.data, alphabet, num.iterations, N, M, W, K, alpha.bg, alpha.mw){

  # the states are vectors of start positions
  initial.positions <- sample(1:(M - W + 1), N, replace = T)
  chain <- matrix(nrow = 0, ncol = N)
  chain <- rbind(chain, initial.positions)

  # do this only ONCE!
  magic.words <- extract.magic.words(my.data, initial.positions, N, W)
  backgrounds <- extract.backgrounds(my.data, initial.positions, N, W)


  for(it in 1:num.iterations){
    current.state <- chain[it, ]
    end.positions <- current.state + W - 1

    next.state <- rep(NA, N)

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
      current.state[n] <- new.start.for.n # it's importat for the current iteration
      magic.words[n, ] <- my.data[n, new.start.for.n:(new.start.for.n+W-1)]
      backgrounds[n, ] <- my.data[n, -c(new.start.for.n:(new.start.for.n+W-1))]
      next.state[n] <- new.start.for.n
    }
    chain <- rbind(chain, next.state)

  }

  # chain represents sequence STATES for EACH iteration!
  # state = vector of starting positions for N sequences
  chain
}


##############################################################################
# Function used to evaluate marginal likelihood of the background
##############################################################################
bg.prob <- function(counts){
  part2 <- sum(lgamma(counts + alpha.bg) - loggamma.bg)
  log.p <- part1.bg.prob + part2
  log.p
}


##############################################################################
# Function used to evaluate marginal likelihood of the magic words
##############################################################################
mw.prob <- function(counts.position){
  part2 <- sum(lgamma(counts.position + alpha.mw) - loggamma.mw)
  log.p <- part1.mw.prob + part2
  log.p
}


##############################################################################
# Normalize vector of values to [0, 1]
##############################################################################
exp.normalize <- function(prob){
  prob <- exp(prob - max(prob))
  return(prob/sum(prob))
}


##############################################################################
# Function used to evaluate posterior, i.e. full conditional
# probability p(r_n | R_{-n}, D)
##############################################################################
posterior <- function(data, current.state, seq.id, magic.words, backgrounds, matrix.bg){
  prob.distribution <- rep(0, M - W + 1)

  # precompute for the magic word
  precomputed.mw <- matrix(nrow = W, ncol = K) # without seq
  mw.prob.no.first.last <- 0

  for(j in 1:W){
    temp <- table(magic.words[-seq.id, j])[alphabet] # index all seq. except seq.id
    temp[is.na(temp)] <- 0
    precomputed.mw[j, ] <- temp
    # for seq.id
    temp <- table(magic.words[seq.id, j])[alphabet] # separately seq.id
    temp[is.na(temp)] <- 0
    mw.prob.no.first.last <- mw.prob.no.first.last + mw.prob(temp + precomputed.mw[j, ])
  }

  precomputed.bg.sum <- matrix.bg[-seq.id, ] %>% apply(2, sum)
  magic.words[seq.id, ] <- data[1:W]

  oldTotal <- table(data)[alphabet]
  oldTotal[is.na(oldTotal)] <- 0

  for(new.position in 1:(M - W + 1)){


    end.position <- new.position + W - 1

    # calculate old first position
    temp.first <- table(magic.words[seq.id, 1])[alphabet]
    temp.first[is.na(temp.first)] <- 0
    # calculate old last position
    temp.last <- table(magic.words[seq.id, W])[alphabet]
    temp.last[is.na(temp.last)] <- 0
    # subtract probabilitiy values for old positions from ALL other sequences
    mw.prob.no.first.last <- mw.prob.no.first.last - mw.prob(temp.first + precomputed.mw[1, ]) -
      mw.prob(temp.last + precomputed.mw[W, ])
    # now you have posterior part without first and last position

    # permanently
    magic.words[seq.id, ] <- data[new.position:end.position] # modify one magic word
   # backgrounds[seq.id, ] <- data[-c(new.position:end.position)] # modify one background
    # permanently
    magic.count <- table(magic.words[seq.id, ])[alphabet]
    magic.count[is.na(magic.count)] <- 0

    counts.bg <- oldTotal - magic.count

    bg.prob.value <- bg.prob(counts.bg + precomputed.bg.sum)
    # try to replace old code
    # first and last position is changed
    temp.first <- table(magic.words[seq.id, 1])[alphabet]
    temp.first[is.na(temp.first)] <- 0

    temp.last <- table(magic.words[seq.id, W])[alphabet]
    temp.last[is.na(temp.last)] <- 0

    mw.prob.no.first.last <- mw.prob.no.first.last + mw.prob(temp.first + precomputed.mw[1, ]) +
      mw.prob(temp.last + precomputed.mw[W, ])

    prob.distribution[new.position] <- bg.prob.value + mw.prob.no.first.last
  }
  prob.distribution
}

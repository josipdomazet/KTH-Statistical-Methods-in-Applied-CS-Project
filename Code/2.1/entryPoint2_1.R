library(MCMCpack) # dirichet distribution
library(ggplot2) # vizualizations
library(magrittr) # %>%

# import functions
source("util.R")

# Question 1
# Specify arbitrary parameters
N <-5
alphabet <-  c("a", "b", "c", "d")
M <- 30
W<- 10
K <- 4
alpha.bg <-  c(1, 1, 1, 1)
alpha.mw<- c(0.8 ,0.8, 0.8, 0.8)


gen.data.list <- generate.data(N, alphabet = c("a", "b", "c", "d"), M, W, alpha.bg,
               alpha.mw)

# Extract starting positions
starts <- gen.data.list$starts

B <- N * (M - W)
part1.bg.prob <- lgamma(sum(alpha.bg)) - lgamma(B + sum(alpha.bg))  # checked: it works
part1.mw.prob <- lgamma(sum(alpha.mw)) - lgamma(N*W + sum(alpha.mw))

loggamma.mw <- lgamma(alpha.mw)
loggamma.bg <- lgamma(alpha.bg)

num.chains <- 10
all.chains <- matrix(nrow = 0, ncol = N)
#results.gibbs <- gibbs.sampler(gen.data.list$data, alphabet,  150, N, M, W, length(alphabet), alpha.bg, alpha.mw)
# 1 gibbs.sampler returns 1 chain which was initialized randomly
# but I need to have multiple chains because the initial state is important
for(c in 1:num.chains){
  all.chains <- rbind(all.chains, gibbs.sampler(gen.data.list$data, alphabet, 500, N, M, W, length(alphabet), alpha.bg, alpha.mw))
}

tibble(counts = as.matrix(all.chains[, 1])) %>% ggplot(aes(x = counts)) + geom_bar() +
  scale_x_continuous(labels = 1:(M - W + 1), breaks = 1:(M - W + 1))

k <- 8

tibble(counts = as.matrix(all.chains[((k-1)*(501)+1):(k*501)])) %>% ggplot(aes(x = counts)) + geom_bar() +
  scale_x_continuous(labels = 1:(M - W + 1), breaks = 1:(M - W + 1))







### Real task
my.data <- fread("data_q2/input.csv", sep = " ", stringsAsFactors = T, header = F, nrows = 25) %>% as.matrix()
K <- 4 # alphabet
W <- 10 # length of a magic word
N <- 25 # number of rows in data == number of sequences
M <- 1000 # length of one sequence
B <- N * (M - W)
alphabet <- c("a", "c", "g", "t")

alpha.bg <- read.table("data_q2/alphaBg.csv", stringsAsFactors = T, header = F, col.names = "alpha") 
alpha.bg <- alpha.bg$alpha # turn into vector

alpha.mw <- read.table("data_q2/alphaMw.csv", stringsAsFactors = T, header = F, col.names = "alpha")
alpha.mw <- alpha.mw$alpha

part1.bg.prob <- lgamma(sum(alpha.bg)) - lgamma(B + sum(alpha.bg))
part1.mw.prob <- lgamma(sum(alpha.mw)) - lgamma(N * W + sum(alpha.mw))

loggamma.mw <- lgamma(alpha.mw)
loggamma.bg <- lgamma(alpha.bg)


system.time({
results.gibbs <- gibbs.sampler(my.data, alphabet, 1, N, M, W, length(alphabet), alpha.bg, alpha.mw)
})
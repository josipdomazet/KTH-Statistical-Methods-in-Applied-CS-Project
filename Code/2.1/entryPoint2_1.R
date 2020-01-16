library(MCMCpack) # dirichet distribution
library(ggplot2) # vizualizations
library(magrittr) # %>%
library(tibble)
library(dplyr)
library(data.table)
library(tidyr)
# import functions
source("util.R")

# Question 1
# Specify arbitrary parameters
N <- 10
alphabet <-  c("a", "b", "c", "d")
M <- 15
W<- 4
K <- 4
alpha.bg <-  c(1, 1, 1, 1)
alpha.mw<- c(0.8 ,0.8, 0.8, 0.8)

gen.data.list <- generate.data(N, alphabet = c("a", "b", "c", "d"), M, W, alpha.bg,
               alpha.mw)

# Extract starting positions
starts <- gen.data.list$starts

B <- N * (M - W)
part1.bg.prob <- lgamma(sum(alpha.bg)) - lgamma(B + sum(alpha.bg))
part1.mw.prob <- lgamma(sum(alpha.mw)) - lgamma(N*W + sum(alpha.mw))

loggamma.mw <- lgamma(alpha.mw)
loggamma.bg <- lgamma(alpha.bg)

num.chains <- 10
all.chains <- matrix(nrow = 0, ncol = N)

# 1 gibbs.sampler returns 1 chain which was initialized randomly
# but we need to have multiple chains because the initial state is randomised
for(c in 1:num.chains){
  all.chains <- rbind(all.chains, gibbs.sampler(gen.data.list$data, alphabet, 100, N, M, W, length(alphabet), alpha.bg, alpha.mw))
}

all.chains.for.plot <- gather(as_tibble(all.chains), key = "seqNum", value = "startPos", 1:25)
all.chains.for.plot$seqNum <- readr::parse_number(all.chains.for.plot$seqNum)
all.chains.for.plot$seqNum <- all.chains.for.plot$seqNum %>%  factor(levels = 1:25)


all.chains.for.plot$chain <- rep(rep(1:10, each = 21),25) # T
all.chains.for.plot$iteration <- rep(rep(1:21, 10), 25)
all.chains.for.plot <- rename(all.chains.for.plot, seq = seqNum)
all.chains.for.plot <- rename(all.chains.for.plot, rn = startPos)



### Real task
K <- 4 # alphabet
W <- 10 # length of a magic word
N <- 25 # number of rows in data == number of sequences
M <- 1000 # length of one sequence
B <- N * (M - W)
my.data <- fread("data_q2/input.csv", sep = " ", stringsAsFactors = T, header = F, nrows = N) %>% as.matrix()
alphabet <- c("a", "c", "g", "t")

alpha.bg <- read.table("data_q2/alphaBg.csv", stringsAsFactors = T, header = F, col.names = "alpha")
alpha.bg <- alpha.bg$alpha # turn into vector

alpha.mw <- read.table("data_q2/alphaMw.csv", stringsAsFactors = T, header = F, col.names = "alpha")
alpha.mw <- alpha.mw$alpha

part1.bg.prob <- lgamma(sum(alpha.bg)) - lgamma(B + sum(alpha.bg))
part1.mw.prob <- lgamma(sum(alpha.mw)) - lgamma(N * W + sum(alpha.mw))

loggamma.mw <- lgamma(alpha.mw)
loggamma.bg <- lgamma(alpha.bg)

num.chains <- 10
all.chains <- matrix(nrow = 0, ncol = N)
#results.gibbs <- gibbs.sampler(gen.data.list$data, alphabet,  150, N, M, W, length(alphabet), alpha.bg, alpha.mw)
# 1 gibbs.sampler returns 1 chain which was initialized randomly
# but I need to have multiple chains because the initial state is important
#system.time({gibbs.sampler(my.data, alphabet, 1, N, M, W, length(alphabet), alpha.bg, alpha.mw)})
system.time({
for(c in 1:num.chains){
  all.chains <- rbind(all.chains, gibbs.sampler(my.data, alphabet, 20, N, M, W, length(alphabet), alpha.bg, alpha.mw))
}
})

for.work <- as_tibble(all.chains)

# convergence
for.work$chain <- rep(1:10, each = 21)
for.work.gathered <- gather(for.work, key  = "seq", value = "ID", 1:N)

for.work.gathered <- rename(for.work.gathered, rn = ID)
#you can index for.work.gathered 1, 1:2, 1:3, 1:4, ..., 1:nrow(for.work.gathered) for
#plotting convergence

plotPSRD <- function(partly){

Ntotal <- max(partly$iteration)

avg.for.each.chain.and.seq <- partly %>% group_by(seq, chain) %>% summarise(
  avg.position = mean(rn), s2 = 1/(Ntotal - 1) *sum((rn - mean(rn))**2)
) %>% rename(x = seq)
avg.for.each.chain.and.seq # mean(theta)_m

W <- avg.for.each.chain.and.seq %>% group_by(x) %>% summarise(
  W.value = mean(s2)
)

avg.for.each.seq <- avg.for.each.chain.and.seq %>% group_by(x) %>% summarize(
  avg.dot = mean(avg.position)
)
avg.for.each.seq # mean(theta)_dot
avg.for.each.seq <- avg.for.each.seq %>% rename( naziv = x)

B <- avg.for.each.chain.and.seq %>% group_by(x) %>% summarise( B.value = Ntotal/(10 - 1) * sum(
  (avg.position - (avg.for.each.seq %>% filter(naziv == x) %>% pull(avg.dot)))**2
)
)

both <- inner_join(B, W) %>% mutate(varHat = (Ntotal-1)/Ntotal * W.value + 1/Ntotal * B.value) %>%
  mutate(Rhat = sqrt(varHat/W.value))

for (i in 1:nrow(both)){
  finalPlot <- add_row(finalPlot, iteration = Ntotal, seqq = both[i, "x"], psrf = both[i, "Rhat"])
}
finalPlot
}

for.work.gathered$iteration <- rep(rep(1:21, 10), 25)

loop100 <- plotPSRD(all.chains.for.plot %>% filter(iteration <= 2))

for( i in 3:101){
  loop100 <- rbind(loop100, plotPSRD(all.chains.for.plot %>% filter(iteration <= i)))

}
loop100 <- mutate(loop100, psrf = as.double(psrf)) %>% rename(seq.num = seqq)
loop100$seq.num <- unlist(loop100$seq.num)
#loop100$seq.num <- readr::parse_number(loop100$seq.num)

#loop100$seq.num <- loop100$seq.num %>% factor(levels  = 1:20)
loop100  %>% ggplot(aes(iteration, psrf, group = seq.num, color = seq.num)) +
  geom_line()
ggsave("q2-psrfn25.png")

source("util.R")
source("../2.2/bootstrap_filter.R")
source("../2.2/resampling_methods.R")
source("../2.2/data_generator.R")

library(tibble)
library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(invgamma)

search.grid <- tibble(
   sigma = numeric(),
   beta = numeric(),
   logl = numeric()
 )

 for(s in seq(0.001,2,by=0.08)){
   for(b in seq(0.001,2,by=0.08)){
   search.grid <- add_row(search.grid, sigma = s, beta = b)
   }
 }

 my.data <- generate_data(1.00, 0.16, 0.64, 100)

 for(row in 1:nrow(search.grid)){
   search.grid[row,"logl"] <- loglikelihood(real.y = my.data$y,
                                            sigma = search.grid[row,"sigma"] %>% pull,
                                            beta = search.grid[row,"beta"] %>% pull,
                                            200, 100, multinomial.resampling)
 }

 # check for highest values
 search.grid %>% arrange(-logl) %>% top_n(20)



 mlogl <- search.grid %>%  mutate(logl = ifelse(logl > -250, logl, NA)) %>% filter(!is.na(logl)) %>%
   summarize(medianlogl = mean(logl)) %>% pull
 search.grid %>% mutate(logl = ifelse(logl > -100, logl, NA)) %>%
   ggplot(aes(x = sigma, y = beta, fill = logl)) +
   geom_tile() +scale_fill_gradient2(midpoint = mlogl, low = "blue", mid = "white",
                                     high = "red", space = "Lab") +
   xlab(expression(sigma)) +
   ylab(expression(beta))  +
   theme(axis.title.y = element_text(angle = 0))

 # Question 9

 set.seed(100)
 runif(1)
 runif(1)
 my.data <- generate_data(1.00, 0.16, 0.64, 100)

 varyT <- tibble(
   Tvalue = numeric(),
   logl = numeric()
 )
 # Fix N = 100
 for(t in c(25, 50, 60, 70, 80, 100)){
   for(i in 1:100){ # to get enough samples for each T
     varyT <- add_row(varyT, Tvalue = t, logl = bootstrap.filter(my.data$y, 0.16,0.64, 200, t, multinomial.resampling)$logl)
   }
 }

 ggplot(varyT, aes(x = Tvalue, group = Tvalue, y = logl)) +
   geom_boxplot() + ylab("log-likelihood") + xlab("t")
 ggsave("boxplot_vary_t.png")

 varyN <- tibble(
   Nvalue = numeric(),
   logl = numeric()
 )
 # Fix N = 100
 for(n in c(25, 50, 100, 200, 300, 500)){
    for(i in 1:100){ # to get enough samples for each T
     varyN <- add_row(varyN, Nvalue = n, logl = bootstrap.filter(my.data$y, 0.16,0.64, n, 100, multinomial.resampling)$logl)
 }
 }

 ggplot(varyN, aes(x = Nvalue, group = Nvalue, y = logl)) +
   geom_boxplot() + ylab("log-likelihood") + xlab("N")
 ggsave("boxplot_vary_n.png")

########################################################################
# PARTICLE METROPOLIS HASTINGS
########################################################################
# 2000, 0.01, 0.2 c(0.001, 0.01) works fine for BETA
my.data <- generate_data(1.00, 0.16, 0.64, 100)
Niteration <- 10000
burnIn <- 5000
pmh.results <- pmh(my.data$y, 0.1, 0.1, 200, 100, Niteration, c(0.01, 0.2), multinomial.resampling)

PMH.posterior.sigma2 <- tibble(sigma2 = pmh.results$theta[1,burnIn:Niteration])
PMH.posterior.beta2 <- tibble(beta2 = pmh.results$theta[2,burnIn:Niteration])

acceptance.ratio <- sum(pmh.results$proposed.accepted)/length(pmh.results$proposed.accepted)
# visualize distributions

# SIGMA^2
ggplot(PMH.posterior.sigma2, aes(x = sigma2)) +
  geom_histogram(bins = 30, aes(y = ..density..))  +
  geom_density(alpha = 0.2, fill = "lightblue") +
  geom_vline(xintercept = mean(PMH.posterior.sigma2$sigma2), size = 1.2, linetype = 3, color = "blue") +
  xlab(expression(sigma^{2})) +
  ylab("") +
  scale_y_discrete(breaks = NULL) +
  annotate('text', x = 0.06, y = 30,
           label = paste("bar(sigma^2)==~",round(mean(PMH.posterior.sigma2$sigma2), 4)),
           parse = TRUE,size=10)

ggsave("./sigma2.png")

# BETA^2
ggplot(PMH.posterior.beta2, aes(x = beta2)) +
  geom_histogram(bins = 30, aes(y = ..density..))  +
  geom_density(alpha = 0.2, fill = "lightblue") +
  geom_vline(xintercept = mean(PMH.posterior.beta2$beta2), size = 1.2, linetype = 3, color = "blue") +
  xlab(expression(beta^{2})) +
  ylab("") +
  annotate('text', x = 1, y =2,
           label = paste("bar(beta^2)==~",round(mean(PMH.posterior.beta2$beta2), 4)),
           parse = TRUE,size=10)  +
  scale_y_continuous(breaks = NULL)
ggsave("./beta2.png")

# line plot
tibble(sigma2 = pmh.results$theta[2,]) %>% mutate(iterations = 1:Niteration) %>% ggplot(aes(x = iterations,
                                                                                            y = sigma2)) +
  geom_line()
                                                               

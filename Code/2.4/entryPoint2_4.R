source("../2.2/data_generator.R")
source("../2.2/resampling_methods.R")
source("util.R")
source("../2.2/bootstrap_filter.R")

library(tidyverse)

set.seed(12793129)
my.data <- generate_data(1.00, 0.16, 0.64, 100)

set.seed(5055)
iter <- 20000
results.gibbs <- particle.gibbs(my.data$y, 0.1, 0.1, 100, 100, iter)
burn.in <- 5000
plot(results.gibbs$sigma.density)
plot(results.gibbs$beta.density)



PG.posterior.sigma2 <- tibble(sigma2 = results.gibbs$sigma.density[burn.in:20000])


ggplot(PG.posterior.sigma2, aes(x = sigma2)) +
  geom_histogram(bins = 30, aes(y = ..density..))  +
  geom_density(alpha = 0.2, fill = "lightblue") +
  geom_vline(xintercept = mean(PG.posterior.sigma2$sigma2), size = 1.2, linetype = 3, color = "blue") +
  xlab(expression(sigma^{2})) +
  ylab("") +
  scale_y_discrete(breaks = NULL) +
  #ggtitle( expression(paste("Distribution of ", sigma^{2}))) +
  annotate('text', x = 0.06, y = 30,
           label = paste("bar(sigma^2)==~",round(mean(PG.posterior.sigma2$sigma2), 4)),
           parse = TRUE,size=10)
ggsave("./pg-sigma2.png")

PG.posterior.beta2 <- tibble(beta2 = results.gibbs$beta.density[burn.in:20000])
ggplot(PG.posterior.beta2, aes(x = beta2)) +
  geom_histogram(bins = 30, aes(y = ..density..))  +
  geom_density(alpha = 0.2, fill = "lightblue") +
  geom_vline(xintercept = mean(PG.posterior.beta2$beta2), size = 1.2, linetype = 3, color = "blue") +
  xlab(expression(beta^{2})) +
  ylab("") +
  #ggtitle( expression(paste("Distribution of ", sigma^{2}))) +
  annotate('text', x = 0.85, y = 2,
           label = paste("bar(beta^2)==~",round(mean(PG.posterior.beta2$beta2), 4)),
           parse = TRUE,size=10) +
  scale_y_continuous(breaks = NULL)
ggsave("./pg-beta2.png")

both <- tibble(x = PG.posterior.sigma2$sigma2, y = PG.posterior.beta2$beta2, t = burn.in:20000)

ggplot(both, aes(x = x, y = y, color = t)) +
  geom_path() +
  xlab(expression(sigma^{2})) +
  ylab(expression(beta^{2})) +
  geom_vline(xintercept = 0.16**2, color = "red", linetype = 2, size = 1) +
  geom_hline(yintercept = 0.64**2, color = "red", linetype = 2, size = 1) +
  scale_color_continuous(name = "iteration") +
  theme(axis.title.y = element_text(angle = 0))
ggsave("./pg-both.png")

tibble(logl = results.gibbs$logl, k = 1:iter) %>% ggplot(aes(k, logl)) +
  geom_line()

source("data_generator.R")
source("SIS.R")

library(tibble)
library(ggplot2)
library(magrittr)

set.seed(2723)
my.data <- generate_data(1.00, 0.16, 0.64, 100)
sis.filter <- sis(my.data$y, 200, 100)
# estimated value of X_100 = -1.472395
# real value = -1.503388

# histogram of weights
tibble(
  value = as.numeric(sis.filter$wnorm[100, ])) %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100)
  ggsave("sis_histogram2.png")

# visualize results
tibble(t = 1:100, value = sis.filter$x.means) %>% 
  ggplot(aes(x = t, y = value)) +
  geom_line(aes(color = "a")) +
  geom_line(data = tibble(t = 1:100, x = my.data$x), aes(x = t, y = x, color = "b")) +
  ggtitle("Predicted vs. real values") +
  scale_color_manual(name = NULL, breaks = c("a", "b"), 
                     labels = c("SIS", "Real data"), values = c("red", "blue"))


mse.tibble.sis <- tibble(N = as.numeric(), mse = as.numeric())
for(N in seq(10,1500, 5)){
  sis.filter <- sis(my.data$y, N, 100)
  mse.tibble.sis <- add_row(mse.tibble.sis, N = N, mse = (sis.filter$x.means[100] - my.data$x[100])^2)
}

ggplot(mse.tibble.sis, aes(x = N, y = mse)) +
  geom_line() + ylab("MSE")

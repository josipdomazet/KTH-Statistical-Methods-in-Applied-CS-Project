library(tibble)
library(ggplot2)
library(magrittr)
library(RColorBrewer)

source("bootstrap_filter.R")
source("resampling_methods.R")
source("data_generator.R")

set.seed(2723)
my.data <- generate_data(1.00, 0.16, 0.64, 100)
bpf.multinomial <- bootstrap.filter(my.data$y, 0.16, 0.64, 200, 100, multinomial.resampling)

# visualize
tibble(t = 1:100, value = bpf.multinomial$x.means) %>% 
  ggplot(aes(x = t, y = value)) +
  geom_line(aes(color = "a")) +
  geom_line(data = tibble(t = 1:100, x = my.data$x), aes(x = t, y = x, color = "b"),
            alpha = 0.4) +
  ggtitle("Predicted vs real values") +
  scale_color_manual(name = NULL, breaks = c("a", "b"), 
                     labels = c("BPF", "real data"), values = c("red", "blue"))


# point estimate
bpf.multinomial$x.means[100]
my.data$x[100]

# MSE

mse.tibble.bpf.multinomial <- tibble(N = as.numeric(), mse = as.numeric())
for(N in seq(10,1500, 5)){
  bpf.multinomial <- bootstrap.filter(my.data$y, 0.16, 0.64, N, 100, multinomial.resampling)
  mse.tibble.bpf.multinomial <- add_row(mse.tibble.bpf.multinomial, N = N, mse = (bpf.multinomial$x.means[100] - my.data$x[100])^2)
}

ggplot(mse.tibble.bpf.multinomial, aes(x = N, y = mse)) +
  geom_line() + ylab("MSE")

# visualise weights
tibble(value = as.numeric(bpf.multinomial$weights.norm[100, ])) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 30)
ggsave("bpf_multinomial_histogram.png")





#############
bpf.stratified <- bootstrap.filter(my.data$y, 0.16, 0.64, 200, 100, stratified.resampling)
tibble(t = 1:100, value = bpf.stratified$x.means) %>% 
  ggplot(aes(x = t, y = value)) +
  geom_line(aes(color = "a")) +
  geom_line(data = tibble(t = 1:100, x = my.data$x), aes(x = t, y = x, color = "b"),
            alpha = 0.4) +
  ggtitle("Predicted vs real values") +
  scale_color_manual(name = NULL, breaks = c("a", "b"), 
                     labels = c("BPF", "real data"), values = c("red", "blue"))

# point estimate
bpf.stratified$x.means[100]

#MSE
mse.tibble.bpf.stratified <- tibble(N = as.numeric(), mse = as.numeric())
for(N in seq(10,1500, 5)){
  bpf.stratified <- bootstrap.filter(my.data$y, 0.16, 0.64, N, 100, stratified.resampling)
  mse.tibble.bpf.stratified <- add_row(mse.tibble.bpf.stratified, N = N, mse = (bpf.stratified$x.means[100] - my.data$x[100])^2)
}

ggplot(mse.tibble.bpf.stratified, aes(x = N, y = mse)) +
  geom_line() + ylab("MSE")

# weights histogram
tibble(value = as.numeric(bpf.stratified$weights.norm[100, ])) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 5) 
  ggsave("bpf_stratified_histogram.png")


# Compare all three methods
sis.filter <- sis(my.data$y, 200, 100)
tibble(value = as.numeric(sis.filter$wnorm[100, ])) %>% ggplot(aes(x = value)) +
  geom_histogram(bins = 30) 
ggsave("sis_histogram.png")

ggplot(data = tibble(t = 1:100, x= bpf.multinomial$x.means), aes(x = t, y = x, color = "c")) +
  geom_line(aes(color = "b")) +
  geom_line(data = tibble(t = 1:100, x = my.data$x), aes(x = t, y = x, color = "a") , linetype = "dashed") +
  geom_line(data = tibble(t = 1:100, x = sis.filter$x.means), aes(t, x, color = "b")) +
  geom_line(data = tibble(t = 1:100, x = bpf.stratified$x.means), aes(t, x)) +
  geom_line(aes(color = "d")) +
  scale_color_manual(name = NULL,
    breaks = c("a", "b", "c", "d"),
                     labels = c("Real data", "SIS", "BPF multinomial", "BPF stratified"),
                     values = c("#666666", "#D95F02", "#66A61E", "#7570B3")) +
  theme(legend.text = element_text(size = 18))
ggsave("comparison_filters.png", units = c("cm"), width = 30)


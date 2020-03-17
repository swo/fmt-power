#!/usr/bin/env Rscript --vanilla

library(tidyverse)

logit <- function(p) log(p / (1 - p))
pdf <- function(p, sigma) dnorm(logit(p), 0, sigma) / (p * (1 - p))

# Check the integral
integral <- integrate(partial(pdf, sigma = 0.25), 0, 1)
if (abs(integral$value - 1) > integral$abs.error) stop("bad pdf")

n <- 1e4
eps <- 0.0001

dat <- crossing(
  sigma = c(0.5, 1.0, 2.0, 2.5),
  p = seq(eps, 1 - eps, length.out = n)
) %>%
  mutate(y = pdf(p, sigma))

plot <- dat %>%
  ggplot(aes(p, y, color = factor(sigma))) +
  geom_line() +
  cowplot::theme_half_open()

ggsave("fig/sigma-spread.pdf")

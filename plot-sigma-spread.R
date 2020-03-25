#!/usr/bin/env Rscript

library(tidyverse)

logit <- function(p) log(p / (1 - p))
pdf <- function(p, sigma) dnorm(logit(p), 0, sigma) / (p * (1 - p))

# Check the integral
integral <- integrate(partial(pdf, sigma = 0.25), 0, 1)
if (abs(integral$value - 1) > integral$abs.error) stop("bad pdf")

n <- 1e4
eps <- 0.0001

dat <- crossing(
  sigma = c(0.3, 1.0, round(sqrt(2), 2), 2.0),
  p = seq(eps, 1 - eps, length.out = n)
) %>%
  mutate(y = pdf(p, sigma))

plot <- dat %>%
  mutate_at("sigma", factor) %>%
  ggplot(aes(p, y, color = sigma)) +
  geom_line() +
  labs(color = expression(sigma[LO])) +
  cowplot::theme_half_open()

ggsave("fig/sigma-spread.pdf")

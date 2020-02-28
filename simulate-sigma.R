#!/usr/bin/env Rscript --vanilla

# Simulate power under the good/bad model. Note that you need a larger Δp than
# suggested from the Gaussian LO model, since that model has bigger outliers:
# you can count on getting a donor with really crazy p.

library(tidyverse)
source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))

simulate1 <- function(n_donors, patients_per_donor, sigma) {
  p <- invlogit(rnorm(n_donors, mean = 0, sd = sigma))
  success <- rbinom(n_donors, patients_per_donor, p)
  fail <- patients_per_donor - success
  mat <- matrix(c(success, fail), ncol = 2)
  test <- safe_fisher_test(mat)
  test$p.value
}

n_iter <- 1e3

results <- crossing(
  sigma = c(0, 0.25, 0.5, 0.75, 1.0, 2.0),
  n_donors = c(2, 4, 6, 12),
  n_patients = c(12, 24, 48, 96, 192, 384),
) %>%
  mutate(
    patients_per_donor = n_patients / n_donors,
    results = pmap(
      list(n_donors, patients_per_donor, sigma),
      ~ simulate(simulate1, n_iter, ..1, ..2, ..3)
    ),
    n_sig = map_int(results, ~ sum(. < 0.05)),
    power = n_sig / n_iter
  )

# Plot power results --------------------------------------------------

display_results <- results %>%
  mutate_at(vars(starts_with("n_")), factor)

plot <- display_results %>%
  mutate_at(vars(starts_with("n_")), factor) %>%
  ggplot(aes(n_donors, n_patients)) +
  facet_wrap(~ sigma) +
  geom_tile(aes(fill = power)) +
  geom_text(aes(label = round(power * 100))) +
  geom_tile(
    data = filter(display_results, power >= 0.80),
    fill = NA, color = "black"
  ) +
  scale_fill_gradient(
    name = "Power (%)",
    labels = function(x) round(x * 100),
    low = "white", high = "red"
  ) +
  cowplot::theme_half_open()

ggsave("power-sigma.pdf")

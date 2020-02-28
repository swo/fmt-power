#!/usr/bin/env Rscript --vanilla

# Simulate power under the good/bad model. Note that you need a larger Δp than
# suggested from the Gaussian LO model, since that model has bigger outliers:
# you can count on getting a donor with really crazy p.

library(tidyverse)

safe_fisher_test <- function(x) {
  try(return(fisher.test(x)), silent = TRUE)
  chisq.test(x)
}

simulate1 <- function(n_donors, patients_per_donor, phi, p0, p1) {
  quality <- rbinom(n_donors, 1, phi)
  p <- p0 + quality * (p1 - p0)
  success <- rbinom(n_donors, patients_per_donor, p)
  fail <- patients_per_donor - success
  mat <- matrix(c(success, fail), ncol = 2)
  test <- safe_fisher_test(mat)
  test$p.value
}

simulate <- function(n_donors, patients_per_donor, phi, p0, p1, n_iter) {
  control <- as.list(environment()) %>%
    `[[<-`("function", first(sys.call()))

  hash <- digest::digest(control)
  cache_fn <- str_glue("cache/{hash}.rds")
  if (file.exists(cache_fn)) return(readRDS(cache_fn))

  results <- map(1:n_iter, ~ simulate1(n_donors, patients_per_donor, phi, p0, p1))
  saveRDS(results, cache_fn)

  results
}

n_iter <- 1e3

results <- crossing(
  phi = c(0.50),
  # delta_p = c(0.0, 0.25, 0.5, 0.75, 1.0),
  delta_p = seq(0, 1, length.out = 6),
  n_donors = c(2, 4, 6, 12),
  n_patients = c(12, 24, 48, 96, 192),
) %>%
  mutate(
    p0 = 0.5 - delta_p / 2,
    p1 = 0.5 + delta_p / 2,
    patients_per_donor = n_patients / n_donors,
    results = pmap(
      list(n_donors, patients_per_donor, phi, p0, p1),
      partial(simulate, n_iter = n_iter)
    ),
    n_sig = map_int(results, ~ sum(. < 0.05)),
    power = n_sig / n_iter
  )

display_results <- results %>%
  mutate_at(vars(starts_with("n_")), factor)

plot <- display_results %>%
  mutate_at(vars(starts_with("n_")), factor) %>%
  ggplot(aes(n_donors, n_patients)) +
  facet_wrap(~ delta_p) +
  geom_tile(aes(fill = power)) +
  geom_text(aes(label = scales::percent(power))) +
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

ggsave("tmp.pdf")

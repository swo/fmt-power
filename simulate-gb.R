#!/usr/bin/env Rscript

source("utils.R")

simulate_base <- function(n_donors, patients_per_donor, phi, p0, p1) {
  quality <- rbinom(n_donors, 1, phi)
  p <- p0 + quality * (p1 - p0)
  success <- rbinom(n_donors, patients_per_donor, p)
  chisq_p(cbind(success, patients_per_donor - success))
}

simulate_f <- function(n_donors, patients_per_donor, p1) {
  p0 <- 1 - p1
  simulate_base(n_donors, patients_per_donor, 0.5, p0, p1)
}

results <- results_f(simulate_f, 0.5, 1, "cache/gb")

results %>%
  select(n_donors, n_patients, effect_size, x, n, estimate, lci, uci) %>%
  write_tsv("results/gb.tsv")

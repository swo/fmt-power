#!/usr/bin/env Rscript --vanilla

source("utils.R")

simulate_base <- function(n_donors, patients_per_donor, phi, p0, p1) {
  quality <- rbinom(n_donors, 1, phi)
  p <- p0 + quality * (p1 - p0)
  success <- rbinom(n_donors, patients_per_donor, p)
  chisq_p(cbind(success, patients_per_donor - success))
}

simulate_f <- function(n_donors, patients_per_donor, delta_p) {
  p0 <- 0.5 - delta_p / 2
  p1 <- 0.5 + delta_p / 2
  simulate_base(n_donors, patients_per_donor, 0.5, p0, p1)
}

f <- effect_size_factory("goodbad", simulate_f, c(0.0, 1.0))

results <- results_f(f)
results

plot <- plot_f(results)
ggsave("fig/gb.pdf")

#!/usr/bin/env Rscript --vanilla

source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))

simulate_f <- function(n_donors, patients_per_donor, sigma) {
  donor_p <- invlogit(rnorm(n_donors, mean = 0, sd = sigma))
  success <- rbinom(n_donors, patients_per_donor, donor_p)
  chisq_p(cbind(success, patients_per_donor - success))
}

f <- effect_size_factory("sigma", simulate_f, c(1e-4, 1e2))

results <- results_f(f)
results

plot <- plot_f(results)
ggsave("fig/sigma.pdf")

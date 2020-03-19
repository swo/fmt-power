#!/usr/bin/env Rscript --vanilla

source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))

simulate_f <- function(n_donors, patients_per_donor, sigma) {
  donor_p <- invlogit(rnorm(n_donors, mean = 0, sd = sigma))
  success <- rbinom(n_donors, patients_per_donor, donor_p)
  chisq_p(cbind(success, patients_per_donor - success))
}

results <- results_f(simulate_f, 1e-4, 3, "cache/sigma")

results %>%
  select(n_donors, n_patients, effect_size, x, n, estimate, lci, uci) %>%
  write_tsv("results/sigma.tsv")

plot <- plot_f(results) +
  scale_x_continuous(
    name = expression(paste("Effect size (", sigma[LO], ")"))
  )

ggsave("fig/sigma.pdf")

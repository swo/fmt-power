#!/usr/bin/env Rscript

source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))
logit <- function(p) log(p / (1 - p))
pdf <- function(p, sigma_lo) dnorm(logit(p), 0, sigma_lo) / (p * (1 - p))

# Given σ_LO, find σ_p
slo_to_sp <- function(sigma_lo) {
  f <- function(x) (x - 0.5) ** 2 * pdf(x, sigma_lo)
  sqrt(integrate(f, 0, 1)$value)
}

# Reverse
# (The integral fails for σ_LO > 12 or so)
sp_to_slo <- function(sigma_p, interval = c(1e-6, 12)) {
  f <- function(x) slo_to_sp(x) - sigma_p
  uniroot(f, interval)$root
}

simulate_f <- function(n_donors, patients_per_donor, sigma_p) {
  sigma_lo <- sp_to_slo(sigma_p)
  donor_p <- invlogit(rnorm(n_donors, mean = 0, sd = sigma_lo))
  success <- rbinom(n_donors, patients_per_donor, donor_p)
  chisq_p(cbind(success, patients_per_donor - success))
}

# this allows for sourcing the script without running this simulations
if (sys.nframe() == 0) {
  results <- results_f(simulate_f, 0, 0.4, "cache/sigma")

  results %>%
    select(n_donors, n_patients, effect_size, x, n, estimate, lci, uci) %>%
    write_tsv("results/sigma.tsv")
}

library(tidyverse)

# Global simulation parameters ----------------------------------------

n_iter <- 1e3
alpha <- 0.05
target_power <- 0.8
uniroot_tol <- 1e-2

# Simulation functions ------------------------------------------------

# APIs for the key functions:
# simulate_f :: n_donors, patients_per_donor, effect_size -> p_value
# power_f :: simulate_f, n_iter -> power
# effect_size_f :: power_f, target_power -> min_effect_size

fisher_p <- function(x) {
  try(return(fisher.test(x)$p.value), silent = TRUE)
  chisq.test(x)$p.value
}

chisq_p <- function(x) {
  p <- chisq.test(x)$p.value
  ifelse(is.na(p), 1.0, p)
}

power_f <- function(simulate_f, n_donors, patients_per_donor, effect_size) {
  p_values <- map(1:n_iter, ~ simulate_f(n_donors, patients_per_donor, effect_size))
  mean(p_values <= alpha)
}

effect_size_f <- function(simulate_f, interval, n_donors, patients_per_donor) {
  obj <- function(x) {
    power_f(simulate_f, n_donors, patients_per_donor, x) - target_power
  }

  try(return(uniroot(obj, interval, tol = uniroot_tol)$root), silent = TRUE)
  NA
}

effect_size_factory <- function(simulation_name, simulate_f, interval) {
  function(n_donors, patients_per_donor) {
    control <- list(simulation = simulation_name, n_donors = n_donors, patients_per_donor = patients_per_donor, interval = interval)
    cache_fn <- str_glue("cache/{digest::digest(control)}.rds")

    if (file.exists(cache_fn)) {
      readRDS(cache_fn)
    } else {
      print(str_glue("Calling simulation {simulation_name} with {n_donors} donor and {patients_per_donor} PPD"))
      result <- effect_size_f(simulate_f, interval, n_donors, patients_per_donor)
      saveRDS(result, cache_fn)
      result
    }
  }
}

results_f <- function(f) {
  base <- crossing(
    n_donors = c(2, 4, 6, 12),
    n_patients = c(12, 24, 48, 96, 192),
  ) %>%
    mutate(patients_per_donor = n_patients / n_donors)

  top <- base %>%
    filter(n_patients > n_donors) %>%
    mutate(effect_size = map2_dbl(n_donors, patients_per_donor, f))

  bottom <- base %>%
    filter(n_patients == n_donors) %>%
    mutate(effect_size = NA)

  bind_rows(top, bottom)
}

plot_f <- function(results) {
  results %>%
    mutate_at(c("n_donors", "n_patients"), factor) %>%
    ggplot(aes(n_donors, n_patients)) +
    geom_tile(aes(fill = effect_size)) +
    geom_text(aes(label = round(effect_size, 2))) +
    scale_fill_gradient(
      name = "Effect size",
      labels = function(x) round(x, 2),
      low = "white", high = "red"
    ) +
    cowplot::theme_half_open()
}

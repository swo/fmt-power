library(tidyverse)

# Probabilistic bisection algorithm -----------------------------------
# https://people.orie.cornell.edu/shane/theses/ThesisRolfWaeber.pdf

prob_bisect <- function(f, interval, pc = 0.75, n_bin = 1e3, min_conf = 0.99, trace = FALSE) {
  # check for a priori failure conditions
  message <- NULL
  if (f(interval[2]) <= f(interval[1])) message <- "Requires an increasing function"
  if (f(interval[2]) <= 0) message <- "Upper interval not high enough"
  if (f(interval[1]) > 0) message <- "Lower interval not low enough"
  if (!is.null(message)) return(list(x = NA, conf = NA, iter = NA, trace = NULL, message = message))

  qc <- 1 - pc

  # set up a series of breakpoints at the endpoints and between the bins
  break_x <- seq(interval[1], interval[2], length.out = n_bin + 1)
  # set up the bins between the breakpoints
  bin_x <- zoo::rollmean(break_x, 2)

  # start with a flat prior
  pdf <- rep(1 / n_bin, n_bin)
  i <- 1

  if (trace) {
    trace <- tibble(i = i, x = bin_x, y = pdf)
  } else {
    trace <- NULL
  }

  while (max(pdf) < min_conf) {
    # get the median point of the PDF
    cdf <- c(0, cumsum(pdf))
    stopifnot(all.equal(last(cdf), 1))
    median_i <- abs(cdf - 0.5) %>%
      {  match(min(.), .) }

    xn <- break_x[median_i]
    yn <- f(xn)

    if (yn >= 0) {
      pdf[bin_x >= xn] <- qc * pdf[bin_x >= xn]
      pdf[bin_x < xn] <-  pc * pdf[bin_x < xn]
    } else {
      pdf[bin_x >= xn] <- pc * pdf[bin_x >= xn]
      pdf[bin_x < xn] <-  qc * pdf[bin_x < xn]
    }

    # renormalize
    pdf <- pdf / sum(pdf)

    i <- i + 1
    if (!is.null(trace)) trace <- bind_rows(trace, tibble(i = i, x = bin_x, y = pdf))
  }

  list(x = xn, conf = max(pdf), iter = i, trace = trace, message = "converged")
}


# Safe test functions -------------------------------------------------

fisher_p <- function(x) {
  try(return(fisher.test(x)$p.value), silent = TRUE)
  chisq.test(x)$p.value
}

chisq_p <- function(x) {
  p <- chisq.test(x)$p.value
  ifelse(is.na(p), 1.0, p)
}

# Simulation functions ------------------------------------------------

# APIs for the key functions:
# simulate_f :: n_donors, patients_per_donor, effect_size -> p_value
# power_f :: simulate_f, n_iter -> power
# effect_size_f :: power_f, target_power -> min_effect_size

obj_factory <- function(simulate_f, n_donors, n_patients, n_iter = 1e2, target_power = 0.8, alpha = 0.05) {
  function(effect_size) {
    p_values <- map_dbl(1:n_iter, ~ simulate_f(n_donors, n_patients, effect_size))
    mean(p_values <= alpha) - target_power
  }
}

results_f <- function(simulate_f, interval) {
  base <- crossing(
    n_donors = c(2, 4, 6, 12),
    n_patients = c(12, 24, 48, 96, 192),
  ) %>%
    mutate(patients_per_donor = n_patients / n_donors)

  # check that patients per donor are integers
  stopifnot(all.equal(base$patients_per_donor, as.integer(base$patients_per_donor)))

  top <- base %>%
    filter(n_patients > n_donors) %>%
    mutate(
      obj = map2(n_donors, n_patients, ~ obj_factory(simulate_f, .x, .y)),
      bisect = map(obj, ~ prob_bisect(., interval)),
      effect_size = map_dbl(bisect, ~ .$x)
    )

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

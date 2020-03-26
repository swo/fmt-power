#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

find_root <- function(x, y, target = 0.8) {
  # stop if there are no values to fit
  if (any(is.na(x) | is.na(y))) return(NA)
  f <- approxfun(x, y)
  # stop if the function doesn't cross the target value
  if (f(max(x)) < target) return(NA)
  # make a new function with a root at the target value
  g <- function(xx) f(xx) - target
  uniroot(g, range(x))$root
}

results <- tibble(fn = args) %>%
  mutate(data = map(fn, read_tsv)) %>%
  unnest(cols = c(data)) %>%
  mutate(
    simulation = if_else(
      is.na(simulation),
      str_split_fixed(basename(fn), "\\.", 2)[, 1],
      simulation
    )
  ) %>%
  select_at(c("simulation", "n_donors", "n_patients", "effect_size", "estimate")) %>%
  chop(c(effect_size, estimate)) %>%
  mutate(min_effect_size = map2_dbl(effect_size, estimate, find_root)) %>%
  select(simulation, n_donors, n_patients, min_effect_size)

write_tsv(results, "results/min-effect-sizes.tsv")

#!/usr/bin/env Rscript --vanilla

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

find_root <- function(data, target = 0.8) {
  model <- loess(value ~ effect_size, data = data)
  f <- function(x) predict(model, newdata = list(effect_size = x)) - target
  if (f(max(model$x)) < 0) return(NA)
  uniroot(f, range(model$x))$root
}

results <- tibble(fn = args) %>%
  mutate(
    data = map(fn, read_tsv),
    simulation = str_split_fixed(basename(fn), "\\.", 2)[, 1]
  ) %>%
  unnest(cols = c(data)) %>%
  select_at(c("simulation", "n_donors", "n_patients", "effect_size", "estimate", "lci", "uci")) %>%
  pivot_longer(c("estimate", "lci", "uci")) %>%
  nest(data = c(effect_size, value)) %>%
  mutate(
    min_effect_size = round(map_dbl(data, find_root), 2)
  ) %>%
  select(simulation, n_donors, n_patients, name, min_effect_size) %>%
  pivot_wider(values_from = min_effect_size)

write_tsv(results, "results/min-effect-sizes.tsv")

#!/usr/bin/env Rscript --vanilla

# Imagine an ANOVA: each patient exhibits some biomarker value that is due to
# the donor and then to some random noise. What is the minimum ratio of the
# across-donor variance to the within-donor, across-patient variance?  E.g., 3
# means across-donor must exceed the across-patient variance by 3x, while 0.1
# means across-patient variance can be 10x the across-donor variance.

library(tidyverse)

results <- crossing(
  n_donors = c(2, 4, 6, 12),
  n_patients = c(12, 24, 48, 96, 192, 384),
) %>%
  filter(!(n_donors == 12 & n_patients == 12)) %>%
  mutate(
    patients_per_donor = n_patients / n_donors,
    test = pmap(
      list(groups = n_donors, n = patients_per_donor, within.var = 1, power = 0.8),
      power.anova.test
    ),
    effect_size = map_dbl(test, ~ .$between.var)
  )

plot <- results %>%
  mutate_at(vars(starts_with("n_")), factor) %>%
  ggplot(aes(n_donors, n_patients)) +
  geom_tile(aes(fill = effect_size)) +
  geom_text(aes(label = round(effect_size, 2))) +
  scale_fill_gradient(
    name = "Effect size",
    labels = function(x) round(x, 2),
    low = "white", high = "red"
  ) +
  cowplot::theme_half_open()

ggsave("power-anova.pdf")

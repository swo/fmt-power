#!/usr/bin/env Rscript

source("utils.R")

# Effect size is σ_D / (σ_D + σ_P). Convert that to σ_D / σ_P as the "between
# variance" in the ANOVA power function.

results <- results_base %>%
  crossing(effect_size = seq(0, 1 - 1e-12, length.out = global_n_grid)) %>%
  mutate(
    between_var = 1 / ((1 / effect_size) - 1),
    estimate = pmap_dbl(
      list(n_donors, patients_per_donor, between_var),
      ~ power.anova.test(groups = ..1, n = ..2, within.var = 1, between.var = ..3)$power
    ),
    lci = NA,
    uci = NA
  )

results %>%
  select(n_donors, n_patients, effect_size, estimate, lci, uci) %>%
  write_tsv("results/anova.tsv")

plot <- results %>%
  mutate_at(c("n_patients"), ~ fct_rev(factor(.))) %>%
  ggplot(aes(effect_size, estimate)) +
  facet_grid(n_patients ~ n_donors) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  geom_hline(yintercept = c(0, 1)) +
  geom_line(size = 0.7) +
  scale_x_continuous(
    name = expression(paste("Effect size ", sigma[D] / (sigma[D] + sigma[P]))),
    expand = c(0, 0),
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2")
  ) +
  scale_y_continuous(
    name = "Statistical power",
    labels = scales::percent,
    breaks = c(0, 0.5, 0.8, 1)
  ) +
  cowplot::theme_half_open() +
  theme(panel.spacing = unit(1, "lines"))

ggsave("fig/anova.pdf")

#!/usr/bin/env Rscript

source("utils.R")

results <- results_base %>%
  crossing(effect_size = seq(0, 2, length.out = 10)) %>%
  mutate(
    estimate = pmap_dbl(
      list(n_donors, patients_per_donor, effect_size),
      ~ power.anova.test(groups = ..1, n = ..2, within.var = 1, between.var = ..3)$power
    ),
    lci = estimate,
    uci = estimate
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
    name = expression(paste("Effect size (", sigma[donor] / sigma[patient], ")")),
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

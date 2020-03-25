#!/usr/bin/env Rscript --vanilla

library(tidyverse)

# Load data
results <- read_tsv("results/16s.tsv")

# Plot
plot <- results %>%
  ggplot(aes(effect_size, estimate)) +
  facet_grid(simulation ~ n_patients) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  geom_hline(yintercept = c(0, 1)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "gray") +
  geom_line(size = 0.7) +
  scale_x_continuous(
    name = expression(paste("Effect size (", Delta * p, ", %)")),
    labels = function(x) x * 100,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Statistical power",
    labels = scales::percent,
    breaks = c(0, 0.5, 0.8, 1)
  ) +
  cowplot::theme_half_open() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

ggsave("fig/16s.pdf")

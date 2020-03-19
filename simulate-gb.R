#!/usr/bin/env Rscript --vanilla

source("utils.R")

simulate_base <- function(n_donors, patients_per_donor, phi, p0, p1) {
  quality <- rbinom(n_donors, 1, phi)
  p <- p0 + quality * (p1 - p0)
  success <- rbinom(n_donors, patients_per_donor, p)
  chisq_p(cbind(success, patients_per_donor - success))
}

simulate_f <- function(n_donors, patients_per_donor, delta_p) {
  p0 <- 0.5 - delta_p / 2
  p1 <- 0.5 + delta_p / 2
  simulate_base(n_donors, patients_per_donor, 0.5, p0, p1)
}

results <- results_f(simulate_f, 0, 1, "cache/gb")

plot <- plot_f(results) +
  scale_x_continuous(
    name = expression(paste("Effect size (", Delta * p, ", %)")),
    labels = function(x) x * 100,
    expand = c(0, 0)
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("fig/gb.pdf")

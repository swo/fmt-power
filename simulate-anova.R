#!/usr/bin/env Rscript --vanilla

source("utils.R")

f <- function(n_donors, patients_per_donor) {
  power.anova.test(groups = n_donors, n = patients_per_donor, within.var = 1, power = 0.8)$between.var
}

results <- results_f(f)
plot <- plot_f(results)
ggsave("fig/anova.pdf")

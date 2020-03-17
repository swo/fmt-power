#!/usr/bin/env Rscript --vanilla

# Imagine an ANOVA: each patient exhibits some biomarker value that is due to
# the donor and then to some random noise. What is the minimum ratio of the
# across-donor variance to the within-donor, across-patient variance?  E.g., 3
# means across-donor must exceed the across-patient variance by 3x, while 0.1
# means across-patient variance can be 10x the across-donor variance.

library(tidyverse)
source("utils.R")

anova_effect_size_f <- function(n_donors, patients_per_donor) {
  power.anova.test(groups = n_donors, n = patients_per_donor, within.var = 1, power = 0.8)$between.var
}

results <- results_f(anova_effect_size_f)
plot <- plot_f(results)
ggsave("fig/anova.pdf")

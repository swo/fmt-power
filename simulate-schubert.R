#!/usr/bin/env Rscript --vanilla

source("utils.R")
library(vegan)

# Load and process ----------------------------------------------------

otu <- read_tsv("data/raw/cdi_schubert_results/cdi_schubert.otu_table.100.denovo")
raw_metadata <- read_tsv("data/raw/cdi_schubert_results/cdi_schubert.metadata.txt")

metadata <- raw_metadata %>%
  filter(DiseaseState %in% c("H", "nonCDI")) %>%
  mutate(quality = recode(DiseaseState, H = 1, nonCDI = 0)) %>%
  select(sample_id, quality) %>%
  arrange(sample_id)

dissim_fn <- "cache/schubert/dissim.rds"
if (file.exists(dissim_fn)) {
  dissim <- readRDS(dissim_fn)
} else {
  otu_matrix <- otu %>%
    select_at(metadata$sample_id) %>%
    as.matrix() %>%
    { .[rowSums(.) > 0, ] } %>%
    t()

  dissim <- vegdist(otu_matrix)
  saveRDS(dissim, dissim_fn)
}

# Simulate ------------------------------------------------------------

# Which indices in the original data are controls/cases?
q0idx <- which(metadata$quality == 0)
q1idx <- which(metadata$quality == 1)

subset_dissim <- function(dissim, idx) {
  dissim %>%
    as.matrix() %>%
    `[`(idx, idx) %>%
    as.dist()
}

simulate_f <- function(n, delta_p, phi = 0.5) {
  # draw donor qualities
  n_good <- rbinom(1, n, phi)
  n_bad <- n - n_good
  quality <- c(rep(1, n_good), rep(0, n_bad))

  # subset dissim matrix
  good_idx <- sample(q1idx, n_good)
  bad_idx <- sample(q0idx, n_bad)
  Y <- subset_dissim(dissim, c(good_idx, bad_idx))

  # simulate patient outcomes
  donor_p <- 0.5 - delta_p / 2 + delta_p * quality
  X <- rbinom(n, 1, donor_p)

  # if outcomes are all 0 or all 1, pval = 1.0
  if (length(unique(X)) == 1) return(1.0)

  # ignore complete enumeration warnings
  suppressWarnings(adonis(Y ~ X)$aov.tab$`Pr(>F)`[1])
}

# multi draw
simulate_trials <- function(n_patients, effect_size) {
  map_dbl(1:global_n_trials, ~ simulate_f(n_patients, effect_size))
}

# memoize
f <- memoise(simulate_trials, cache = cache_filesystem("cache/schubert"))

results <- crossing(
  n_patients = c(5, 10, 25, 50, 75),
  effect_size = seq(0, 1, length.out = 10)
) %>%
  mutate(
    p_values = pmap(list(n_patients, effect_size), f),
    x = map_dbl(p_values, ~ sum(. <= 0.05)),
    n = map_dbl(p_values, length),
    test = map2(x, n, binom.test),
    estimate = x / n,
    lci = map_dbl(test, ~ .$conf.int[1]),
    uci = map_dbl(test, ~ .$conf.int[2])
  )

results %>%
  select(n_patients, effect_size, x, n, estimate, lci, uci) %>%
  write_tsv("results/schubert.tsv")

plot <- results %>%
  # mutate_at(c("n_patients"), ~ fct_rev(factor(.))) %>%
  ggplot(aes(effect_size, estimate)) +
  facet_grid(. ~ n_patients) +
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

ggsave("fig/schubert.pdf")

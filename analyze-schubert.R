#!/usr/bin/env Rscript --vanilla

# Look at the diarrhea case/control 16S data. Assume "good" donors look like
# controls and "bad" donors look like cases. Given some N of good and N of bad
# donors, and some Δp between their efficacies, what's the statistical power of
# a typical PERMANOVA?

library(tidyverse)
library(vegan)
source("utils.R")

# Load and process ----------------------------------------------------

otu <- read_tsv("data/raw/cdi_schubert_results/cdi_schubert.otu_table.100.denovo")
raw_metadata <- read_tsv("data/raw/cdi_schubert_results/cdi_schubert.metadata.txt")

metadata <- raw_metadata %>%
  filter(DiseaseState %in% c("H", "nonCDI")) %>%
  mutate(quality = recode(DiseaseState, H = 1, nonCDI = 0)) %>%
  select(sample_id, quality) %>%
  arrange(sample_id)

dissim_fn <- "dissim.rds"
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

  adonis(Y ~ X)$aov.tab$`Pr(>F)`[1]
}

root_factory <- function(n, n_iter = 1e2) {
  function(x) {
    print(str_glue("n={n} x={x}"))
    p_values <- map_dbl(1:n_iter, ~ simulate_f(n, x))
    mean(p_values <= 0.05) - 0.8
  }
}

results <- tibble(n = c(5, 10, 25, 50, 75)) %>%
  mutate(
    root_f = map(n, root_factory),
    bisect = map(root_f, ~ prob_bisect(., c(0, 1), n_bin = 1e2, min_conf = 0.8)),
    effect_size = map_dbl(bisect, ~ .$x)
  )

results %>%
  select(n, effect_size) %>%
  write_tsv("cache/schubert.tsv")

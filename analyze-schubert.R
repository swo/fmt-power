#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(vegan)

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

# Simulate outcomes ---------------------------------------------------

# p0 <- 0.01
# delta_p <- 0.98

# p <- p0 + metadata$quality * delta_p
# outcome <- rbinom(nrow(metadata), 1, p)

# Simulate ------------------------------------------------------------

q0idx <- which(metadata$quality == 0)
q1idx <- which(metadata$quality == 1)

subset_dissim <- function(dissim, idx) {
  dissim %>%
    as.matrix() %>%
    `[`(idx, idx) %>%
    as.dist()
}

simulate <- function(n, p0, p1) {
  X <- c(rbinom(n, 1, p0), rbinom(n, 1, p1))

  idx <- c(sample(q0idx, n), sample(q1idx, n))
  Y <- subset_dissim(dissim, idx)

  adonis(Y ~ X)
}

results_fn <- "results.rds"

if (file.exists(results_fn)) {
  results <- readRDS(results_fn)
} else {
  results <- crossing(
    n = c(5, 10, 25, 50),
    delta_p = c(0.01, 0.25, 1.0),
    iter = 1:20
  ) %>%
    mutate(
      test = map2(n, delta_p, ~ simulate(.x, 0.5 - .y / 2, 0.5 + .y / 2)),
      p = map_dbl(test, ~ .$aov.tab$`Pr(>F)`[1])
    )

  saveRDS(results, results_fn)
}

powers <- results %>%
  group_by(n, delta_p) %>%
  summarize(
    n_total = n(),
    n_sig = sum(p <= 0.05)
  ) %>%
  mutate(
    test = map2(n_sig, n_total, binom.test),
    estimate = map_dbl(test, ~ .$estimate),
    cil = map_dbl(test, ~ .$conf.int[1]),
    ciu = map_dbl(test, ~ .$conf.int[2])
  )

# plot <- powers %>%
#   ggplot(aes(n, y = NULL, color = factor(delta_p))) +
#   geom_crossbar(aes(y = estimate, ymin = cil, ymax = ciu)) +
#   geom_point(data = results, aes(y = p), position = "jitter")

# plot <- results %>%
#   ggplot(aes(factor(n), p, color = factor(delta_p))) +
#   geom_violin() +
#   geom_point()

plot <- results %>%
  ggplot(aes(p, fill = factor(n))) +
  facet_grid(delta_p ~ .) +
  geom_histogram(position = "dodge", bins = 20)

ggsave("tmp.pdf")

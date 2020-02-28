#!/usr/bin/env Rscript --vanilla

# Look at the diarrhea case/control 16S data. Assume "good" donors look like
# controls and "bad" donors look like cases. Given some N of good and N of bad
# donors, and some Δp between their efficacies, what's the statistical power of
# a typical PERMANOVA?

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

simulate1 <- function(n, p0, p1) {
  X <- c(rbinom(n, 1, p0), rbinom(n, 1, p1))

  idx <- c(sample(q0idx, n), sample(q1idx, n))
  Y <- subset_dissim(dissim, idx)

  adonis(Y ~ X)
}

simulate <- function(n, p0, p1, n_iter = 1) {
  control <- list(n = n, p0 = p0, p1 = p1, n_iter = n_iter)

  hash <- digest::digest(control)
  cache_fn <- str_glue("cache/{hash}.rds")
  if (file.exists(cache_fn)) return(readRDS(cache_fn)$results)

  results <- map(1:n_iter, ~ simulate1(n, p0, p1))
  saveRDS(list(control = control, results = results), cache_fn)

  results
}

n_iter <- 50
results <- crossing(
  n = c(5, 10, 25, 50, 75),
  delta_p = c(0.01, 0.25, 0.35, 0.50, 1.0),
) %>%
  mutate(test = pmap(list(n, delta_p, n_iter), ~ simulate(..1, 0.5 - ..2 / 2, 0.5 + ..2 / 2, ..3))) %>%
  unnest() %>%
  mutate(p = map_dbl(test, ~ .$aov.tab$`Pr(>F)`[1]))

# histogram of p's against N and Δp
plot <- results %>%
  ggplot(aes(p)) +
  facet_grid(delta_p ~ n, scales = "free_y") +
  geom_histogram(position = "dodge", bins = 20)

# plot of powers
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

plot <- powers %>%
  ggplot(aes(factor(n), estimate, color = factor(delta_p))) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_hline(yintercept = 0.80, linetype = 2) +
  geom_pointrange(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.25))

ggsave("tmp.pdf")

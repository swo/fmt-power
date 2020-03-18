#!/usr/bin/env Rscript --vanilla

source("utils.R")

invlogit <- function(lo) 1 / (1 + exp(-lo))

simulate_f <- function(n_donors, patients_per_donor, sigma) {
  donor_p <- invlogit(rnorm(n_donors, mean = 0, sd = sigma))
  success <- rbinom(n_donors, patients_per_donor, donor_p)
  chisq_p(cbind(success, patients_per_donor - success))
}

interval <- c(1e-4, 3)

results <- results_f(simulate_f, interval)
results

plot <- plot_f(results)
ggsave("fig/sigma.pdf")

# Run spot checks
data_f <- function(f) {
  tibble(x = seq(interval[1], interval[2], length.out = 10)) %>%
    mutate(y = map_dbl(x, f))
}

check_plot <- results %>%
  filter(patients_per_donor > 1) %>%
  mutate(data = map(obj, data_f)) %>%
  select(n_donors, n_patients, data) %>%
  unnest(cols = c(data)) %>%
  ggplot(aes(x, y, color = factor(n_patients))) +
  facet_wrap(~ n_donors) +
  geom_line()

ggsave("fig/sigma-spot-check.pdf")

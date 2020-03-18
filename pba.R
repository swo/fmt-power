#!/usr/bin/env Rscript --vanilla

# https://people.orie.cornell.edu/shane/theses/ThesisRolfWaeber.pdf

library(tidyverse)

stochastic_bisection <- function(f, interval, pc = 0.75, n_bin = 1e2, min_conf = 0.99, trace = FALSE) {
  if (f(interval[2]) <= f(interval[1])) stop("Requires an increasing function")

  qc <- 1 - pc

  # set up a series of breakpoints at the endpoints and between the bins
  break_x <- seq(interval[1], interval[2], length.out = n_bin + 1)
  # set up the bins between the breakpoints
  bin_x <- zoo::rollmean(break_x, 2)

  # start with a flat prior
  pdf <- rep(1 / n_bin, n_bin)
  i <- 1

  if (trace) {
    trace <- tibble(i = i, x = bin_x, y = pdf)
  } else {
    trace <- NULL
  }

  while (max(pdf) < min_conf) {
    # get the median point of the PDF
    cdf <- c(0, cumsum(pdf))
    stopifnot(all.equal(last(cdf), 1))
    median_i <- abs(cdf - 0.5) %>%
      {  match(min(.), .) }

    xn <- break_x[median_i]
    yn <- f(xn)

    if (yn >= 0) {
      pdf[bin_x >= xn] <- qc * pdf[bin_x >= xn]
      pdf[bin_x < xn] <-  pc * pdf[bin_x < xn]
    } else {
      pdf[bin_x >= xn] <- pc * pdf[bin_x >= xn]
      pdf[bin_x < xn] <-  qc * pdf[bin_x < xn]
    }

    # renormalize
    pdf <- pdf / sum(pdf)

    i <- i + 1
    if (trace) trace <- bind_rows(trace, tibble(i = i, x = bin_x, y = pdf))
  }

  list(x = xn, conf = max(pdf), n_iter = i, trace = trace)
}

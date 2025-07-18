library(tidyverse)
library(latex2exp)
library(entropart)

setwd("/mnt/p/DataSciencePoster")

# Fringe trimming ----
# Define parameter grid
N_vals <- seq(100, 500000, by = 100) # Total reads
p_vals <- c(0.001, 0.005, 0.01, 0.02) # Minimum detectable frequencies
alpha_vals <- c(0.1, 0.05, 0.01) # Confidence levels

# Generate a grid of all combinations
grid <- expand.grid(
  N = N_vals,
  p = p_vals,
  alpha = alpha_vals
)

# Calculate the minimum count for each combo
grid <- grid %>%
  mutate(
    n_min = qbinom(1 - alpha, size = N, prob = p),
    label = paste0("p = ", p)
  )

# Plotting
ggplot(grid, aes(x = N, y = n_min, color = factor(p))) +
  geom_line(size = 1.5) +
  facet_wrap(~alpha, labeller = label_bquote(cols = alpha == .(alpha))) +
  scale_y_continuous("Minimum count to retain haplotype") +
  scale_x_continuous("Total number of reads in sample", n.breaks = 4) +
  scale_color_brewer(palette = "Dark2", name = "Frequency threshold (p)") +
  theme_minimal(base_size = 24) +
  ggtitle("Count Thresholds for Haplotype Retention") +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )
ggsave("MinCount.png", dpi = 600, width = 11, height = 8.5, units = "in")

# SNV filter ----
# Minimum minor base count
# Parameters
N_vals <- seq(10, 1000, by = 10)
b_hat_vals <- c(0.001, 0.005, 0.01)
alpha <- c(0.1, 0.05, 0.01)

# Generate data frame
df <- expand.grid(N = N_vals, b_hat = b_hat_vals, alpha = alpha) %>%
  mutate(
    k_min = qbinom(1 - alpha, size = N, prob = b_hat) + 1
  )

# Plot
ggplot(df, aes(x = N, y = k_min, color = factor(b_hat))) +
  geom_line(size = 1.5) +
  scale_color_brewer(
    palette = "Dark2",
    name = TeX("Allele Frequency $(\\hat{b})$")
  ) +
  labs(
    title = "Minimum Minor Base Count for SNV Detection",
    x = "Coverage (Total reads at position)",
    y = "Min variant base count (to call SNV)"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(~alpha, labeller = label_bquote(cols = alpha == .(alpha)))
ggsave("MinBaseCount.png", dpi = 600, width = 11, height = 8.5, units = "in")

# Both major and minor
# Parameters
alpha <- 0.05
error_rate <- 0.05
b_hat <- 0.01

# Range of depths to evaluate
n_vals <- 10:1000

# Compute thresholds
thresholds <- tibble(
  depth = n_vals,
  dominant = qbinom(alpha, size = depth, prob = 1 - error_rate),
  minor = sapply(depth, function(n) {
    k <- 0
    while (1 - pbinom(k - 1, n, b_hat) >= alpha && k <= n) {
      k <- k + 1
    }
    k
  })
)

# Pivot to long format for ggplot
thresholds_long <- thresholds %>%
  tidyr::pivot_longer(
    cols = c(dominant, minor),
    names_to = "ThresholdType", values_to = "MinBaseCount"
  )

# Plot
ggplot(
  thresholds_long,
  aes(x = depth, y = MinBaseCount, color = ThresholdType)
) +
  geom_line(size = 1.5) +
  xlim(10, 1000) +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  scale_color_manual(
    values = c("dominant" = "#1B9E77", "minor" = "#D95F02"),
    labels = c("Dominant Base Threshold", "Minor Base Threshold")
  ) +
  labs(
    title = "Binomial Thresholds for SNV Filtering",
    x = "Total Read Depth (n)",
    y = "Minimum Base Count Required",
    color = "Test Type"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    legend.position = "top"
  )
ggsave("SNVFilter.png", dpi = 600, width = 11, height = 8.5, units = "in")

# Hill number ----
# Simulate haplotype counts and SNV counts
set.seed(42)
n_haps <- 20000

hap_df <- tibble(
  haplotype = sprintf("hap%03d", 1:n_haps),
  count = c(
    sample(1000:20000, 2, replace = TRUE), # high count (real)
    sample(50:999, 20, replace = TRUE), # moderate (ambiguous)
    sample(2:70, 100, replace = TRUE), # low (likely noise)
    rep(1, n_haps - 122) # singletons (to be filtered)
  ),
  snvs = c(
    sample(3:10, 2, replace = TRUE), # real
    sample(2:8, 20, replace = TRUE), # moderate
    sample(0:5, 100, replace = TRUE), # noise
    sample(0:2, n_haps - 122, replace = TRUE) # singletons
  )
)

# Create grid of min_count x min_snv thresholds
grid <- expand_grid(
  min_count = seq(1, 1000, by = 25),
  min_snv = 0:10
)

# Compute Hill numbers across the grid
diversity_grid <- grid %>%
  rowwise() %>%
  mutate(
    Hill0 = {
      filtered <- hap_df %>%
        filter(count >= min_count, snvs >= min_snv)
      if (nrow(filtered) == 0) {
        NA_real_
      } else {
        nrow(filtered)
      }
    },
    Hill1 = {
      filtered <- hap_df %>%
        filter(count >= min_count, snvs >= min_snv)
      if (nrow(filtered) < 2) {
        NA_real_
      } else {
        freqs <- filtered$count / sum(filtered$count)
        exp(-sum(freqs * log(freqs)))
      }
    },
    Hill2 = {
      filtered <- hap_df %>%
        filter(count >= min_count, snvs >= min_snv)
      if (nrow(filtered) < 2) {
        NA_real_
      } else {
        freqs <- filtered$count / sum(filtered$count)
        1 / sum(freqs^2)
      }
    }
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("Hill"),
    names_to = "HillOrder", values_to = "value"
  ) %>%
  mutate(HillOrder = factor(substr(HillOrder, 5, 5),
    levels = c("0", "1", "2"), labels = c("q = 0", "q = 1", "q = 2")
  ))

# Plot contour plots for Hill0, Hill1, and Hill2
ggplot(diversity_grid, aes(x = min_count, y = min_snv, z = value)) +
  geom_contour_filled(
    bins = 18,
    breaks = c(
      0, 1, 2, 5, 10, 15, 20, 50, 75, 100,
      200, 500, 1000, 1500, 2000, 5000, 10000, 20000
    )
  ) +
  facet_wrap(~HillOrder, scales = "free", labeller = label_value) +
  scale_fill_viridis_d(name = "Retained\nHaplotypes") +
  labs(
    title = "Diversity (Hill Numbers) Across Filtering Thresholds",
    x = "Minimum Count Threshold",
    y = "Minimum SNV Threshold"
  ) +
  ylim(c(0, 10)) +
  xlim(c(0, 1000)) +
  theme_minimal(base_size = 24)
ggsave("HillNums.png", dpi = 600, width = 22, height = 8.5, units = "in")

# Real data ----
haps <- read_delim("haplotype_counts.txt", delim = "\t") %>%
  filter(count > 1) |>
  rename(Haplotype = barcode)
table(haps$count)
summary(haps$count)

fringe_trim <- function(counts, p_thresh = 0.001, conf = 0.95) {
  n <- sum(counts)
  n_min <- qbinom(conf, size = n, prob = p_thresh)
  counts[which(counts >= n_min)]
}

hap_counts <- fringe_trim(haps$count)
haps_trimmed <- haps |> filter(count %in% hap_counts)
table(haps$count)
summary(haps$count)
ggplot(haps_trimmed, aes(x = log(count))) +
  geom_histogram() +
  theme_minimal()

# Function to detect SNVs at each aligned position
detect_snvs <- function(seqs, error_rate = 0.05, alpha = 0.05) {
  # seqs: character vector of equal-length aligned haplotype strings
  mat <- do.call(rbind, strsplit(seqs, ""))
  len <- ncol(mat)
  # count bases at each position
  counts_list <- lapply(seq_len(len), function(j) table(mat[, j]))

  # total coverage per site Ni
  n <- sapply(counts_list, sum)

  # first test: nmax_i = max count at site i
  nmax <- sapply(counts_list, max)
  p1 <- pbinom(nmax, size = n, prob = 1 - error_rate)
  snv1 <- which(p1 < alpha)

  # estimate b from normal sites C: second-dominant base freq
  c_inds <- setdiff(seq_len(len), snv1)
  b_vals <- sapply(c_inds, function(j) {
    cs <- sort(counts_list[[j]], decreasing = TRUE)
    if (length(cs) >= 2) cs[2] / sum(cs) else 0
  })
  b <- sum(b_vals * n[c_inds]) / sum(n[c_inds])

  # second test: second-dominant count n2_i, right-tail
  n2 <- sapply(counts_list, function(tbl) {
    cs <- sort(tbl, decreasing = TRUE)
    if (length(cs) >= 2) cs[2] else 0
  })

  p2 <- 1 - pbinom(n2 - 1, size = n, prob = b)
  snv2 <- which(p2 < alpha)

  # final SNV calls = union of sites significant in both tests
  # intersection of snv1 and snv2
  intersect(snv1, snv2)
}

filter_haplotypes_by_snvs <- function(df, snv_sites, reference = NULL) {
  # df: data.frame with columns Haplotype and count
  seqs <- df$Haplotype
  seq_mat <- do.call(rbind, strsplit(seqs, ""))
  if (is.null(reference)) {
    reference <- seqs[which.max(df$count)]
  }
  ref_vec <- strsplit(reference, "")[[1]]

  # Check for differences from reference at any SNV site
  has_real_snv <- apply(seq_mat[, snv_sites, drop = FALSE], 1, function(row) {
    any(row != ref_vec[snv_sites])
  })
  has_real_snv[which.max(df$count)] <- TRUE

  df[has_real_snv, ]
}

snvs <- detect_snvs(haps_trimmed$Haplotype)
haps_snv <- filter_haplotypes_by_snvs(haps_trimmed, snvs)

summary(haps_snv$count)

# build_toy_data.R --------------------------------------------------
# Generates the bundled `toy_*` datasets used by examples and the
# user-cohort vignette. 50 genes/proteins x 24 samples, 4 stages
# (Normal / Early / Mid / Late) x 2 cohorts. The first 10 (RNA) /
# first 8 (protein) features are stage-progressive (`Early_Burst_Up`
# shape); the rest are noise.
#
# Reproducible: `set.seed(2026)` is the only RNG control used here.
# Re-run via `source("data-raw/build_toy_data.R")`.

set.seed(2026L)
n_g <- 50L
n_s <- 24L

# ---- coldata ------------------------------------------------------
toy_coldata <- data.frame(
  sample = paste0("S", sprintf("%02d", seq_len(n_s))),
  stage  = rep(c("Normal", "Early", "Mid", "Late"), each = n_s / 4L),
  cohort = rep(c("CohortA", "CohortB"), times = n_s / 2L),
  row.names = NULL,
  stringsAsFactors = FALSE)

# ---- counts (RNA) -------------------------------------------------
toy_counts <- matrix(
  rnbinom(n_g * n_s, mu = 50, size = 4),
  nrow = n_g, ncol = n_s)
rownames(toy_counts) <- paste0("TGENE", sprintf("%02d", seq_len(n_g)))
colnames(toy_counts) <- toy_coldata$sample

# Inject 10 Early_Burst_Up genes (lift = +50 reads from Stage I)
up_lift <- c(Normal = 0, Early = 60, Mid = 80, Late = 90)
for (i in seq_len(10L)) {
  for (j in seq_len(n_s)) {
    toy_counts[i, j] <- toy_counts[i, j] +
      stats::rpois(1L, up_lift[toy_coldata$stage[j]])
  }
}
storage.mode(toy_counts) <- "integer"

# ---- protein intensity (log2) -------------------------------------
n_p <- 50L
toy_protein <- matrix(
  rnorm(n_p * n_s, mean = 12, sd = 1.2),
  nrow = n_p, ncol = n_s)
rownames(toy_protein) <- paste0("TPROT", sprintf("%02d", seq_len(n_p)))
colnames(toy_protein) <- toy_coldata$sample

# Inject 8 stage-progressive proteins (Early_Burst_Up shape, log2 lift)
prot_lift <- c(Normal = 0, Early = 0.6, Mid = 0.9, Late = 1.0)
for (i in seq_len(8L)) {
  for (j in seq_len(n_s)) {
    toy_protein[i, j] <- toy_protein[i, j] +
      prot_lift[toy_coldata$stage[j]]
  }
}

# ---- save ---------------------------------------------------------
usethis::use_data(toy_counts,  overwrite = TRUE, compress = "xz")
usethis::use_data(toy_coldata, overwrite = TRUE, compress = "xz")
usethis::use_data(toy_protein, overwrite = TRUE, compress = "xz")

cat(sprintf(
  "toy_counts: %d x %d (integer)\ntoy_protein: %d x %d (log2)\ntoy_coldata: %d rows\n",
  nrow(toy_counts), ncol(toy_counts),
  nrow(toy_protein), ncol(toy_protein),
  nrow(toy_coldata)))

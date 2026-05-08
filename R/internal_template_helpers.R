# Private helpers shared by plot_template_atlas() and
# plot_gene_template(). Not exported.
#
# Cache for per-layer 12-template best-match assignments. Computing
# argmax rho across the 10,113-gene atlas is ~1 s; we memoize on
# (layer, atlas-pointer) so repeat calls during a single
# plot_template_atlas() walk pay this once.
.template_assign_cache <- new.env(parent = emptyenv())

# Resolve the 12-template assignment for every gene/protein in a
# layer, by argmax of the per-template Pearson rho computed on z-scored
# 4-point trajectories.
.assign_templates_12 <- function(layer, reference = NULL) {
  layer <- match.arg(layer, c("rna", "protein"))
  ref <- .template_layer_data(layer, reference)
  cache_key <- paste0(layer, "_",
                       digest_pointer(ref))
  if (exists(cache_key, envir = .template_assign_cache,
              inherits = FALSE)) {
    return(get(cache_key, envir = .template_assign_cache))
  }

  e <- new.env()
  utils::data("default_templates", package = "pdactrace", envir = e)
  tpl_mat <- do.call(rbind, e$default_templates)

  beta_cols <- if (layer == "rna") {
    c("rna_beta_N", "rna_beta_E", "rna_beta_M", "rna_beta_L")
  } else {
    c("prot_beta_N", "prot_beta_E", "prot_beta_M", "prot_beta_L")
  }
  M <- as.matrix(ref[, beta_cols, with = FALSE])
  storage.mode(M) <- "numeric"

  rmean <- rowMeans(M, na.rm = FALSE)
  rsd   <- apply(M, 1, stats::sd)
  Z <- (M - rmean) / rsd
  ok <- is.finite(rsd) & rsd > 0 & complete.cases(Z)

  rho_mat <- matrix(NA_real_, nrow = nrow(M), ncol = nrow(tpl_mat),
                     dimnames = list(NULL, rownames(tpl_mat)))
  if (any(ok)) {
    Tn <- t(scale(t(tpl_mat), center = TRUE, scale = TRUE))
    Zn <- Z[ok, , drop = FALSE]
    Zn <- t(scale(t(Zn), center = TRUE, scale = TRUE))
    rho_mat[ok, ] <- (Zn %*% t(Tn)) / 3
  }

  best_idx <- max.col(rho_mat, ties.method = "first")
  best_idx[!ok] <- NA_integer_
  best_template <- colnames(rho_mat)[best_idx]
  best_rho <- rho_mat[cbind(seq_len(nrow(rho_mat)), best_idx)]

  out <- data.table::data.table(
    gene_symbol     = ref$gene_symbol,
    template_argmax = best_template,
    rho_argmax      = best_rho,
    z_N = Z[, 1], z_E = Z[, 2], z_M = Z[, 3], z_L = Z[, 4])
  data.table::setkey(out, gene_symbol)

  assign(cache_key, out, envir = .template_assign_cache)
  out
}

# Pull the right beta-table for a layer.
.template_layer_data <- function(layer, reference = NULL) {
  if (layer == "rna") {
    .get_reference(reference)
  } else {
    e <- new.env()
    utils::data("pdactrace_protein_betas", package = "pdactrace",
                 envir = e)
    e$pdactrace_protein_betas
  }
}

# Cheap content fingerprint of an atlas, so the assignment cache keys
# off "is this the same data?" without holding the whole frame.
digest_pointer <- function(x) {
  paste0("ptr-", format(utils::object.size(x)),
          "-", nrow(x), "-", ncol(x))
}

# Long-format aggregate for one template: returns cohort genes' z-scored
# trajectories ready to feed geom_line / geom_ribbon.
.template_aggregate <- function(layer, template, reference = NULL) {
  layer <- match.arg(layer, c("rna", "protein"))
  asg <- .assign_templates_12(layer, reference = reference)
  cohort <- asg[!is.na(template_argmax) & template_argmax == template]
  if (nrow(cohort) == 0L) {
    return(data.table::data.table(
      gene_symbol = character(0),
      stage = factor(character(0),
                      levels = c("Normal","Early","Mid","Late")),
      z = numeric(0), template = character(0)))
  }
  long <- data.table::rbindlist(list(
    cohort[, .(gene_symbol, stage = "Normal", z = z_N)],
    cohort[, .(gene_symbol, stage = "Early",  z = z_E)],
    cohort[, .(gene_symbol, stage = "Mid",    z = z_M)],
    cohort[, .(gene_symbol, stage = "Late",   z = z_L)]))
  long[, stage := factor(stage,
                          levels = c("Normal","Early","Mid","Late"))]
  long[, template := template]
  long[]
}

# The single ggplot used by both public functions.
.plot_template_panel <- function(agg,
                                  highlight_gene  = NULL,
                                  highlight_color = "#C62828",
                                  template_label  = NULL,
                                  subtitle        = NULL) {
  n_cohort <- length(unique(agg$gene_symbol))
  show_individual <- n_cohort >= 2L
  show_ribbon     <- n_cohort >= 5L
  n_label <- if (n_cohort == 1L) "n = 1"
             else if (!show_ribbon) sprintf("n = %d (no ribbon)", n_cohort)
             else sprintf("n = %d", n_cohort)

  summ <- agg[, .(z_mean = mean(z, na.rm = TRUE),
                  z_sd   = stats::sd(z, na.rm = TRUE)),
              by = stage]
  summ[is.na(z_sd), z_sd := 0]
  summ[, z_lo := z_mean - z_sd]
  summ[, z_hi := z_mean + z_sd]

  p <- ggplot2::ggplot()

  if (show_individual) {
    p <- p + ggplot2::geom_line(
      data = agg,
      mapping = ggplot2::aes(x = stage, y = z, group = gene_symbol),
      color = "grey55", alpha = 0.05, linewidth = 0.18)
  }

  if (show_ribbon) {
    p <- p + ggplot2::geom_ribbon(
      data = summ,
      mapping = ggplot2::aes(x = stage, ymin = z_lo, ymax = z_hi,
                              group = 1),
      fill = "#1565C0", alpha = 0.18, color = NA)
  }

  p <- p + ggplot2::geom_line(
    data = summ,
    mapping = ggplot2::aes(x = stage, y = z_mean, group = 1),
    color = "#0D47A1", linewidth = 0.55) +
    ggplot2::geom_point(
      data = summ,
      mapping = ggplot2::aes(x = stage, y = z_mean),
      color = "#0D47A1", size = 0.9)

  if (!is.null(highlight_gene)) {
    hl <- agg[gene_symbol == highlight_gene]
    if (nrow(hl) == 4L) {
      p <- p + ggplot2::geom_line(
        data = hl,
        mapping = ggplot2::aes(x = stage, y = z, group = gene_symbol),
        color = highlight_color, linewidth = 0.9) +
        ggplot2::geom_point(
          data = hl,
          mapping = ggplot2::aes(x = stage, y = z),
          color = highlight_color, size = 1.4)
    }
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(-2.2, 2.2),
                                  breaks = c(-2, -1, 0, 1, 2)) +
    ggplot2::labs(x = NULL, y = "z-score",
                  title = template_label,
                  subtitle = subtitle) +
    ggplot2::annotate("text", x = 0.55, y = 2.05,
                       label = n_label, hjust = 0,
                       size = 1.7, color = "grey45") +
    pdactrace_axes_theme()

  p
}

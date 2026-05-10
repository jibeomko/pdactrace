#' Project a user-supplied serum proteomics cohort (binary contrast)
#'
#' The serum side of the framework only requires **binary group
#' contrasts** -- PDAC vs healthy control (HC) and optionally
#' Pancreatitis vs HC -- not stage-aware data. This makes user
#' projection asymmetric to the tissue side: where
#' [project_user_cohort()] needs Normal / Early / Mid / Late stage
#' labels for its 12-template trajectory matching, this helper
#' accepts the kind of serum DIA dataset most public PDAC serum
#' studies actually publish (PDAC + HC +/- pancreatitis, no clinical
#' stage labels needed).
#'
#' @section What this function does:
#' \enumerate{
#'   \item Splits samples into PDAC / HC / (optional) Pancreatitis
#'     groups by `coldata[[group_col]]`.
#'   \item Computes `serum_log2fc_PDAC_vs_HC = mean(PDAC) - mean(HC)`
#'     and (if `pan_label` is set)
#'     `serum_log2fc_Pan_vs_HC = mean(Pan) - mean(HC)` per gene.
#'   \item Runs a per-gene two-group test (default `limma_eBayes`;
#'     fallback to `wilcox` or `t_test`) to attach BH-adjusted
#'     p-values for each contrast.
#'   \item Optionally joins the bundled
#'     [pdactrace_reference]'s `rna_pattern` direction to derive a
#'     per-gene `translation_class` (A / B / C):
#'     \itemize{
#'       \item **Class A** -- tissue and serum agree on direction
#'         ('same dir.', concordant translation).
#'       \item **Class B** -- tissue and serum oppose
#'         ('opposite', inverse translation; the rare,
#'         pancreatitis-discriminating case).
#'       \item **Class C** -- decoupled (tissue direction unclear,
#'         or serum log2FC very small).
#'     }
#'     Set `link_to_atlas = FALSE` to skip the atlas join and emit
#'     `translation_class = NA` (deferring class assignment to
#'     downstream code).
#' }
#'
#' Output is shaped to plug straight into
#' [assemble_user_evidence()] via the `serum_summary` argument.
#'
#' @section Input format:
#' \describe{
#'   \item{`intensity`}{A numeric matrix or `data.frame` of
#'     **log2 protein abundance**, with **genes (HGNC symbols) in
#'     rows** and **samples in columns**. Pre-log2 transformation
#'     is required; raw MS intensities will produce nonsensical
#'     log2FC values.}
#'   \item{`coldata`}{A `data.frame` with at least one column
#'     (`group_col`) carrying group labels matching `pdac_label`,
#'     `hc_label`, and (optionally) `pan_label`. `nrow(coldata)`
#'     must equal `ncol(intensity)`. Sample order must match
#'     column order of `intensity`.}
#' }
#' Minimum 3 samples per group is recommended (`limma` and
#' `wilcox` need at least that for a stable variance estimate /
#' rank distribution). Fewer samples emit a `warning`.
#'
#' @param intensity Log2 protein abundance matrix (genes x samples).
#' @param coldata `data.frame` of sample metadata.
#' @param group_col Name of the column in `coldata` carrying group
#'   labels. Default `"group"`.
#' @param pdac_label Group label for PDAC samples. Default `"PDAC"`.
#' @param hc_label Group label for healthy-control samples.
#'   Default `"HC"`.
#' @param pan_label Group label for pancreatitis samples (optional).
#'   Default `NULL` (PDAC vs HC only; no Class A / B assignment
#'   without pancreatitis context, but Class C and serum_detected
#'   are still set).
#' @param test One of `"limma_eBayes"` (default; requires the
#'   `limma` package), `"wilcox"`, or `"t_test"`.
#' @param padj_cutoff Significance threshold for `serum_detected`.
#'   Default `0.05`. Genes with `padj_PDAC_vs_HC <= padj_cutoff`
#'   are flagged as `serum_detected = TRUE`.
#' @param link_to_atlas Logical. If `TRUE` (default), join the
#'   bundled `pdactrace_reference`'s `rna_pattern` direction to
#'   assign `translation_class`. If `FALSE`, leave `translation_class`
#'   as `NA` (defer to downstream code).
#' @param reference Optional `data.table` to inject in place of the
#'   bundled atlas (used by tests).
#' @return A `data.table` with one row per gene in `intensity` and
#'   columns:
#'   \itemize{
#'     \item `gene_symbol`
#'     \item `serum_log2fc_PDAC_vs_HC`, `serum_padj_PDAC_vs_HC`
#'     \item `serum_log2fc_Pan_vs_HC`, `serum_padj_Pan_vs_HC`
#'       (NA if `pan_label = NULL`)
#'     \item `translation_class` (`"A"` / `"B"` / `"C"` / `NA`)
#'     \item `serum_detected` (logical)
#'   }
#'   Suitable for direct use as `assemble_user_evidence(
#'   serum_summary = ...)`.
#' @examples
#' \dontrun{
#'   # Synthetic example: 50 genes x 24 samples (8 PDAC, 8 HC,
#'   # 8 Pancreatitis) with PDAC-up + Pan-flat signal.
#'   set.seed(1)
#'   g <- 50; n <- 24
#'   intensity <- matrix(rnorm(g * n, mean = 5),
#'                       nrow = g, ncol = n,
#'                       dimnames = list(paste0("GENE", seq_len(g)),
#'                                        paste0("S", seq_len(n))))
#'   intensity[1:10, 1:8] <- intensity[1:10, 1:8] + 2  # PDAC-up
#'   coldata <- data.frame(
#'     sample = colnames(intensity),
#'     group  = rep(c("PDAC", "HC", "Pancreatitis"), each = 8))
#'
#'   serum_summary <- project_user_serum_cohort(
#'     intensity, coldata,
#'     group_col  = "group",
#'     pdac_label = "PDAC",
#'     hc_label   = "HC",
#'     pan_label  = "Pancreatitis")
#'   head(serum_summary)
#'
#'   # Feed into the downstream evidence assembler:
#'   ev <- assemble_user_evidence(serum_summary = serum_summary)
#' }
#' @seealso [project_user_cohort()] for the **tissue** side
#'   (stage-aware), [assemble_user_evidence()] for combining
#'   tissue + serum into the audit-rule input,
#'   [plot_serum_direction()] for visualisation.
#' @export
project_user_serum_cohort <- function(intensity,
                                       coldata,
                                       group_col   = "group",
                                       pdac_label  = "PDAC",
                                       hc_label    = "HC",
                                       pan_label   = NULL,
                                       test = c("limma_eBayes",
                                                 "wilcox", "t_test"),
                                       padj_cutoff   = 0.05,
                                       link_to_atlas = TRUE,
                                       reference     = NULL) {
  test <- match.arg(test)

  # ---- input validation -----------------------------------------------
  if (!is.matrix(intensity) && !is.data.frame(intensity)) {
    stop("`intensity` must be a numeric matrix or data.frame ",
         "(genes by samples).", call. = FALSE)
  }
  X <- as.matrix(intensity)
  if (!is.numeric(X)) {
    stop("`intensity` must contain numeric (log2 abundance) ",
         "values.", call. = FALSE)
  }
  if (is.null(rownames(X))) {
    stop("`intensity` must have row names = HGNC gene symbols.",
         call. = FALSE)
  }
  if (!is.data.frame(coldata)) {
    stop("`coldata` must be a data.frame.", call. = FALSE)
  }
  if (nrow(coldata) != ncol(X)) {
    stop(sprintf(
      "nrow(coldata) (%d) must equal ncol(intensity) (%d).",
      nrow(coldata), ncol(X)), call. = FALSE)
  }
  if (!group_col %in% colnames(coldata)) {
    stop(sprintf(
      "coldata has no column '%s'. Available: %s",
      group_col, paste(colnames(coldata), collapse = ", ")),
      call. = FALSE)
  }

  grp <- as.character(coldata[[group_col]])
  pdac_idx <- which(grp == pdac_label)
  hc_idx   <- which(grp == hc_label)
  pan_idx  <- if (is.null(pan_label)) integer(0L)
              else which(grp == pan_label)

  # Friendly check: required groups present?
  unmatched <- setdiff(unique(grp),
                       c(pdac_label, hc_label, pan_label))
  if (length(pdac_idx) == 0L || length(hc_idx) == 0L) {
    stop(sprintf(
      "`%s` column must contain group labels matching ",
      group_col),
      sprintf("`pdac_label = \"%s\"` and `hc_label = \"%s\"`. ",
              pdac_label, hc_label),
      sprintf("Found: %s. ",
              paste(sort(unique(grp)), collapse = ", ")),
      "Pass the labels you actually use, e.g. ",
      "`pdac_label = \"Cancer\"`, `hc_label = \"Control\"`.",
      call. = FALSE)
  }
  if (length(pdac_idx) < 3L || length(hc_idx) < 3L) {
    warning(sprintf(
      "Small group sizes: PDAC = %d, HC = %d. ",
      length(pdac_idx), length(hc_idx)),
      "limma / wilcox tests are unstable with < 3 per group.",
      call. = FALSE)
  }
  if (!is.null(pan_label) && length(pan_idx) < 3L) {
    warning(sprintf(
      "Pancreatitis group size = %d (< 3). Skipping Pan vs HC ",
      length(pan_idx)),
      "contrast.", call. = FALSE)
    pan_idx <- integer(0L)
    pan_label <- NULL
  }

  # ---- log2FC contrasts ----------------------------------------------
  pdac_mean <- rowMeans(X[, pdac_idx, drop = FALSE], na.rm = TRUE)
  hc_mean   <- rowMeans(X[, hc_idx,   drop = FALSE], na.rm = TRUE)
  log2fc_pdac <- pdac_mean - hc_mean
  if (length(pan_idx) > 0L) {
    pan_mean <- rowMeans(X[, pan_idx, drop = FALSE], na.rm = TRUE)
    log2fc_pan <- pan_mean - hc_mean
  } else {
    log2fc_pan <- rep(NA_real_, nrow(X))
  }

  # ---- per-gene tests -------------------------------------------------
  pvec_pdac <- .ps_pvec(X, pdac_idx, hc_idx, test)
  pvec_pan  <- if (length(pan_idx) > 0L)
    .ps_pvec(X, pan_idx, hc_idx, test) else rep(NA_real_, nrow(X))
  padj_pdac <- stats::p.adjust(pvec_pdac, method = "BH")
  padj_pan  <- if (length(pan_idx) > 0L)
    stats::p.adjust(pvec_pan, method = "BH") else
    rep(NA_real_, nrow(X))

  # ---- translation_class via atlas RNA direction ----------------------
  trans_class <- if (isTRUE(link_to_atlas)) {
    .ps_translation_class(rownames(X), log2fc_pdac, reference)
  } else {
    rep(NA_character_, nrow(X))
  }

  # ---- assemble -------------------------------------------------------
  out <- data.table::data.table(
    gene_symbol              = rownames(X),
    serum_log2fc_PDAC_vs_HC  = log2fc_pdac,
    serum_padj_PDAC_vs_HC    = padj_pdac,
    serum_log2fc_Pan_vs_HC   = log2fc_pan,
    serum_padj_Pan_vs_HC     = padj_pan,
    translation_class        = trans_class,
    serum_detected           = !is.na(padj_pdac) &
                                 padj_pdac <= padj_cutoff)
  data.table::setattr(out, "test",         test)
  data.table::setattr(out, "n_pdac",       length(pdac_idx))
  data.table::setattr(out, "n_hc",         length(hc_idx))
  data.table::setattr(out, "n_pan",        length(pan_idx))
  data.table::setattr(out, "padj_cutoff",  padj_cutoff)
  if (length(unmatched) > 0L) {
    message(sprintf(
      "Note: %d sample(s) with unmatched group labels (",
      sum(grp %in% unmatched)),
      paste(unmatched, collapse = ", "),
      ") were ignored.")
  }
  out[]
}

# ---- internal helpers --------------------------------------------------

.ps_pvec <- function(X, group_a_idx, group_b_idx, test) {
  if (test == "limma_eBayes") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      warning("`limma` not installed; falling back to wilcox.",
              call. = FALSE)
      return(.ps_pvec_wilcox(X, group_a_idx, group_b_idx))
    }
    return(.ps_pvec_limma(X, group_a_idx, group_b_idx))
  }
  if (test == "t_test")  return(.ps_pvec_ttest(X, group_a_idx,
                                                  group_b_idx))
  .ps_pvec_wilcox(X, group_a_idx, group_b_idx)
}

.ps_pvec_limma <- function(X, idx_a, idx_b) {
  # Two-group limma with explicit contrast a - b.
  Y <- X[, c(idx_a, idx_b), drop = FALSE]
  grp <- factor(c(rep("a", length(idx_a)),
                   rep("b", length(idx_b))),
                 levels = c("b", "a"))   # b is reference -> coef "a"
  design <- stats::model.matrix(~ grp)
  fit  <- limma::lmFit(Y, design)
  fit2 <- limma::eBayes(fit)
  tt   <- limma::topTable(fit2, coef = "grpa", number = Inf,
                          sort.by = "none")
  tt[match(rownames(X), rownames(tt)), "P.Value"]
}

.ps_pvec_wilcox <- function(X, idx_a, idx_b) {
  vapply(seq_len(nrow(X)), function(i) {
    a <- X[i, idx_a]; b <- X[i, idx_b]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 2L || length(b) < 2L) return(NA_real_)
    suppressWarnings(stats::wilcox.test(a, b)$p.value)
  }, numeric(1L))
}

.ps_pvec_ttest <- function(X, idx_a, idx_b) {
  vapply(seq_len(nrow(X)), function(i) {
    a <- X[i, idx_a]; b <- X[i, idx_b]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 2L || length(b) < 2L ||
        stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
    suppressWarnings(stats::t.test(a, b)$p.value)
  }, numeric(1L))
}

# Translation-class assignment by joining bundled atlas RNA direction.
# A = tissue & serum agree; B = oppose; C = decoupled (tissue dir
# unclear, or serum log2FC near zero).
.ps_translation_class <- function(genes, serum_log2fc,
                                    reference = NULL) {
  ref <- .get_reference(reference)
  hits <- match(genes, ref$gene_symbol)
  rna_dir <- rep(NA_character_, length(genes))
  matched <- !is.na(hits)
  if (any(matched)) {
    pat <- ref$rna_pattern[hits[matched]]
    rna_dir[matched] <- vapply(pat, function(p) {
      if (is.na(p)) return(NA_character_)
      if (grepl("Up|Burst|Peak|Plateau_Up", p)) return("up")
      if (grepl("Down|Loss|Trough|Plateau_Down", p)) return("down")
      NA_character_
    }, character(1L))
  }
  ser_dir <- ifelse(is.na(serum_log2fc) | abs(serum_log2fc) < 0.1,
                     NA_character_,
                     ifelse(serum_log2fc > 0, "up", "down"))
  cls <- rep(NA_character_, length(genes))
  cls[!is.na(rna_dir) & !is.na(ser_dir) & rna_dir == ser_dir] <- "A"
  cls[!is.na(rna_dir) & !is.na(ser_dir) & rna_dir != ser_dir] <- "B"
  # everything that has serum data but tissue direction unknown / borderline
  cls[!is.na(serum_log2fc) & is.na(cls)] <- "C"
  cls
}

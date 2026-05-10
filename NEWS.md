# pdactrace 0.99.5

**Evidence Math layer.** Adds a third interpretability layer between
the existing plain-English provenance text (`format_provenance()`) and
the audit-score decomposition (`explain_score()`): the *mathematical
evidence values* that actually fed those decisions, exposed per axis.
Per-axis-first by design — there is no new black-box composite score.

## New

- `evidence_math(gene_symbol)` — returns a structured list with one
  element per axis (trajectory_fit, effect_magnitude,
  cohort_consistency, rna_protein_coupling, serum_bridge,
  cell_specificity, filter_survival, clinical_role). All values are
  read directly from the bundled atlas or derived by a single named
  arithmetic operation. New derived quantities: `delta_rho =
  rho_best - rho_runner_up` (template specificity margin), `||beta||_2`
  (Euclidean norm of the per-stage beta vector over E,M,L) and
  `cosine(beta_RNA, beta_prot)` (RNA-protein direction concordance).
- `explain_gene(gene_symbol, view = c("evidence", "math", "both"))` —
  text formatter that prints either the plain-English provenance, the
  Evidence Math sections, or both. Mirrors `explain_score()` in
  pattern (verbose-print + invisible structured return).
- `compare_genes(gene_symbols, axes, wide)` — multi-gene tidy-table
  pivot of `evidence_math()` output. Long form by default
  (`gene, axis, metric, value`); `wide = TRUE` for one row per gene
  with `axis.metric` columns suitable for manuscript tables.

## Documentation

- README.md: new **Evidence Math** subsection demonstrating
  `evidence_math("LTBP1")`, `explain_gene("LTBP1", view = "math")`,
  and `compare_genes(...)` with truncated example output.
- README.md: Reproducibility section rewritten to position the
  bundled atlas as the **default offline path** (Layer 1) and the
  user-cohort projection (Layer 2) as the second core scope. The
  optional processed-input rebuild (Layer 3) and raw-data
  reanalysis (Layer 4) are clearly framed as *advanced
  reproducibility / data provenance* outside the package's primary
  user path. `download_phase_csvs()` is documented as not evaluated
  during package checks.
- README.md: note that the appropriate Bioconductor mechanism for
  hosting the larger raw quantification matrices is a companion
  `pdactraceData` ExperimentHub package, planned for a future
  release.

# pdactrace 0.99.4

**Reproducibility vignette + on-demand phase33/34 fetcher.** Closes
the v0.99.3 self-contained-rebuild story: the bundled `inst/extdata`
plus the new `download_phase_csvs()` helper plus a 4-layer
`vignette("05_reproducibility")` give end users a complete documented
path from the prebuilt atlas all the way back to the public
phase-CSV inputs, without cloning the manuscript-monorepo.

## New

- `download_phase_csvs(target, ref, cache, verbose)` — fetches the
  two large upstream phase CSVs (`phase33_deseq2_coef_12template.csv`,
  `phase34_protein_pooled_12template.csv`) on demand from
  `raw.githubusercontent.com/jibeomko/PDAC_biomarker` (Zenodo
  10.5281/zenodo.20067849) and caches them via
  [`BiocFileCache::BiocFileCache()`]. Returns a named character
  vector of cached file paths. `BiocFileCache` is a soft
  (`Suggests`) dependency picked up at runtime via
  `requireNamespace()` — no impact on packages that only need the
  query / score / report API.
- `vignettes/05_reproducibility.Rmd` — 4-layer reproducibility
  guide:
  1. Layer 1 — offline, bundled (use `data/*.rda` directly; no
     network).
  2. Layer 2 — user cohort (project a count / intensity matrix
     through `project_user_cohort()`).
  3. Layer 3 — re-derive the bundled atlas from the package's own
     `inst/extdata` plus the two upstream CSVs cached via
     `download_phase_csvs()`.
  4. Layer 4 — full FASTQ → counts → fits pipeline (out of scope;
     pointer to the manuscript-monorepo).

## Changed

- `DESCRIPTION` — `Suggests:` gains `BiocFileCache` for the new
  helper. Version bumped 0.99.3 → 0.99.4. Imports / Depends are
  unchanged.

# pdactrace 0.99.3

**Self-contained reproducibility — small phase CSVs bundled in
`inst/extdata`.** The package now ships the six small downstream
phase tables (phase2c / phase42 / phase77 / phase29 / phase60 /
phase80; ~52 KB compressed total) so that the full
`build_reference.R` re-derivation chain runs without the companion
manuscript-monorepo. Two large upstream CSVs (phase33 RNA fit,
phase34 protein fit) remain external because they push the tarball
over Bioconductor's 5 MB limit; the build scripts now fall through
a clean lookup chain (`inst/extdata/` → `$PDAC_BASE_DIR/...` →
informative error pointing at the manuscript Zenodo archive
[10.5281/zenodo.20067849](https://doi.org/10.5281/zenodo.20067849)).

## New

- `inst/extdata/phase{2c,29,42,60,77,80}_*.csv.xz` — bundled,
  xz-compressed copies of the six small phase outputs that
  `build_reference.R` joins onto the RNA/protein base. Total ~52 KB
  added to the tarball; no impact on R CMD check or BiocCheck.
- `data-raw/bundle_phase_csvs.R` — bundler that re-creates the
  `inst/extdata/*.csv.xz` set from the manuscript-monorepo when
  `PDAC_BASE_DIR` is set. Strips phase33 to the 16 columns that
  `build_reference.R` actually uses; everything else is a verbatim
  trim + xz pass.

## Changed

- `data-raw/build_reference.R` — adds an internal `read_phase()`
  helper that prefers the bundled `inst/extdata` copy, falls back to
  `$PDAC_BASE_DIR/.../*.csv` (and `.csv.xz`), and errors with a
  clear "set PDAC_BASE_DIR or download from Zenodo
  10.5281/zenodo.20067849" message if neither exists. Six of the
  seven `fread()` calls in the script now route through this helper.
  The two large CSVs (phase33, phase34) still come via
  `read_phase()` but resolve via the fallback chain because the
  tarball can't accommodate them.
- `data-raw/build_protein_betas.R` — same lookup chain via
  `.find_phase_csv()` + xz-aware fread. Re-running the script now
  works against the bundled CSV (set
  `data/pdactrace_protein_betas.rda` rebuild path entirely
  in-package once the phase34 input is available).

## Documentation

- README "Reproducibility" section update is queued — the current
  v0.99.0..v0.99.2 wording still claims phase scripts "live in a
  separate companion repository". A subsequent commit will rewrite
  this section to acknowledge the bundled `inst/extdata` chain and
  the Zenodo download path for the two remaining external CSVs.

# pdactrace 0.99.2

**Human-readable provenance + larger Arial-bold typography on all
plots.** Two reviewer-facing polish passes: the `phaseXX` tags that
the bundled atlas carries in `provenance` are now relabelled to
plain-English evidence sources on the user-facing print path
(phase IDs are preserved verbatim as a `Technical:` footer for
reproducibility), and the NCS-grade ggplot theme that backs every
`plot_*()` function gains larger, sharper Arial-bold typography for
clearer single-column embedding.

## New

- `format_provenance(provenance, style)` (exported) — maps the
  internal phase tags onto plain-English evidence-source labels.
  Three styles: `"compact"` (one-line plus-separated summary),
  `"verbose"` (multi-line bulleted list with explanations),
  `"raw"` (the original phase tags, for the technical footer).
  Mapping: `phase33` → "RNA trajectory", `phase34` → "Tissue
  protein trajectory", `stouffer_consistency` → "Multi-cohort RNA
  consistency", `phase2c` → "scRNA cell origin", `phase42` →
  "Serum / pancreatitis comparison", `phase60` → "7-step serum
  filter audit", `phase77` → "Strict RNA-protein-serum bridge",
  `phase29` → "Resectable-stage marker screen", `phase80` →
  "Predeclared panel member". Unknown tags fall through unchanged
  with an `(unmapped)` suffix in verbose mode.

## Changed

- `query_gene()` print method — replaces the `Provenance:` line
  with two lines: `Evidence:` (plus-separated human labels) and
  `Technical:` (the raw phase IDs).
- `query_gene_detailed()` print method — surfaces the verbose
  bulleted list so that for an interactive query, the reviewer
  sees a per-source one-line explanation immediately, with the
  raw IDs folded below.
- `summarize_gene_evidence()` — same two-line restructuring of the
  trailing provenance block; the rest of the summary text is
  unchanged.
- `R/theme.R` — `pdactrace_axes_theme()` and
  `pdactrace_panel_theme()` switch font family to **Arial** with a
  sans fallback, bump axis titles 6.5 → 8.5 pt (bold), axis text
  5.8 → 7.5 pt (now also bold), plot title 8 → 10 pt, subtitle 5.5
  → 7.5 pt, legend text 5.4 → 7.0 pt, strip text 5.6 → 7.5 pt,
  tighten the axis line / tick contrast (grey15 instead of
  grey30). Plots render the same shape as before but at clearly
  legible sizes for both single-column and double-column embedding.
  No API change — every existing `plot_*()` function picks up the
  new typography automatically.

# pdactrace 0.99.1

**Per-template trajectory PDFs (RNA + protein).** Two new exported
functions render the bundled atlas's 12-template trajectory
catalogue as compact, NCS-grade PDF panels — one panel per template
showing the cohort of genes / proteins that match it best
(z-scored Pearson-rho argmax across all 12 templates) as thin
translucent lines, a mean line, and a ±1 SD ribbon. Reproduces the
historical `fig2C_RNA_*.pdf` figures inside the package and extends
them to the four templates introduced in v0.4.0
(`Late_Burst_Up`, `Late_Loss_Down`, `Monotonic_Up`,
`Monotonic_Down`).

## New

- `plot_template_atlas(layer, templates, output_dir, ...)` —
  atlas-wide reference, returns 12 ggplot panels (or a user-chosen
  subset). When `output_dir` is set, writes one cairo-PDF per
  template via `pdactrace_save()` and reports paths in
  `attr(out, "files")`.
- `plot_gene_template(gene_symbol, layer, output_file, ...)` —
  per-gene overlay: identifies the gene's matched template via the
  same 12-template Pearson-rho argmax that backs
  `classify_trajectory()` / `score_trajectory()`, then draws that
  template's cohort with the gene's own z-scored trajectory
  highlighted on top. Genes whose best-match is a non-Early
  template (i.e., `pdactrace_reference$rna_pattern` is `NA` because
  only Early × 4 is surfaced) get a visible
  *"non-Early best-match"* sub-title.
- `pdactrace_protein_betas` — bundled `data.table` (5,917 measurable
  proteins × 8 columns) carrying per-stage tissue-protein effect
  sizes (`prot_beta_N/E/M/L`), the full 12-template best-match
  label, matched ρ, and BH-adjusted F-test padj. Generated by
  `data-raw/build_protein_betas.R` from
  `phase34_protein_pooled_12template.csv`. Drives the protein layer
  of the two new plot functions; the main atlas
  (`pdactrace_reference`) only carries categorical
  `prot_pattern` / `prot_tier`, which is why this separate object
  is needed for trajectory visualisation.

## Edge-case handling

`.plot_template_panel()` (the private layout helper used by both
public functions) degrades gracefully on tiny cohorts: cohorts of
1 drop the SD ribbon, cohorts of 2–4 drop the ribbon but keep
individual lines, cohorts of ≥ 5 show the full thin-lines + mean +
ribbon. This keeps every one of the 12 panels renderable —
relevant for `template_protein_Late_Loss_Down` (n = 1 in the
current bundled data).

## Tests

- `tests/testthat/test-plot-template-atlas.R` — 6 blocks:
  length-12 / named return shape (RNA + protein), every RNA cohort
  non-empty, `templates` subsetting, invalid-name error, `output_dir`
  writes one PDF per template.
- `tests/testthat/test-plot-gene-template.R` — 5 blocks: ggplot
  return for an Early-surface gene, non-Early fallback subtitle for
  GAPDH (housekeeping), PDF write-out, protein-layer round-trip,
  missing-gene error.
- Total testthat coverage: 262 → **317 PASS** on Linux.

# pdactrace 0.99.0

**First Bioconductor submission.** Per Bioconductor convention, new
packages submit at 0.99.0; subsequent versions in the review queue go
0.99.1, 0.99.2, ... and on acceptance into the release branch the
version becomes 1.0.0. The audit framework, the bundled atlas, the
12-template trajectory catalog, and the external anchor evaluation
results are all carried forward unchanged from the v0.4.1
development line. v0.99.0 is the first version made visible to
Bioconductor and consolidates everything below in a single release.

## Companion repositories and DOIs

* **Software (this package):**
  [github.com/jibeomko/pdactrace](https://github.com/jibeomko/pdactrace),
  Zenodo concept DOI
  [10.5281/zenodo.20069896](https://doi.org/10.5281/zenodo.20069896).
* **Manuscript reproducibility monorepo:**
  [github.com/jibeomko/PDAC_biomarker](https://github.com/jibeomko/PDAC_biomarker),
  Zenodo
  [10.5281/zenodo.20067849](https://doi.org/10.5281/zenodo.20067849).

## New — single-gene + panel evidence reports

- `report_gene(gene_symbol, ...)` — accepts either a single HGNC
  symbol or a multi-gene character vector. Length 1 produces a
  self-contained one-page HTML evidence report (audit-class tag,
  audit-component table, six-axis evidence radar versus the
  `high_confidence_mean` reference, stage-trajectory plot,
  per-cohort breakdown, per-stage detail table, filter trace,
  atlas-version provenance footer). Length 2+ produces a panel
  template with a `compare_candidates()` table at the top, a
  multi-gene evidence radar, and one `explain_score()`-derived
  rationale card per gene. Both templates use a minimal inline
  CSS (no Bootstrap) and weigh ~230–410 KB.
  Templates: `inst/rmd/gene_report.Rmd`, `inst/rmd/panel_report.Rmd`.
  `rmarkdown` and `knitr` are picked up at runtime from `Suggests`
  via `requireNamespace()` — no hard dependency added.

## New — audit-framework explainability

- `explain_score(gene_symbol, verbose = TRUE)` — decomposes one
  gene's `audit_score` into the three weighted evidence axes
  (`evidence_strength` 40%, `biological_coherence` 35%,
  `translational_relevance` 25%) and the two reliability gates
  (`leakage_gate`, `heterogeneity_gate`). Returns a list with the
  structured breakdown plus a one-paragraph plain-English rationale
  derived from which gate state actually fired (housekeeping flag,
  plasma-high-abundance flag, max meta I² bracket). E.g.,
  *"LTBP1 lands in `supported_uncertain` because the
  heterogeneity_gate = 0.7 — max meta I² = 75% in [70%, 90%)"*.
- `compare_candidates(gene_symbols)` — one row per input gene with
  audit class + score, RNA / protein trajectory pattern, translation
  class, dominant scRNA cell origin, serum detectability, max meta
  I², plus a `redundancy_with` column flagging genes that share both
  `rna_pattern` and `cell_origin_top` in the input set. Sorted by
  `audit_score` descending; missing genes returned as NA rows.
- `project_user_cohort(rna, coldata, stage_col, cohort_col, ...)` —
  end-to-end wrapper around `fit_stage_de()` →
  `classify_trajectory()` → `assemble_user_evidence()` →
  `compute_audit_score()`. Accepts either a count matrix +
  `coldata`, or a `SummarizedExperiment` whose `colData()` carries
  the stage / cohort columns. Optional `protein` argument (matrix or
  SE) layers in tissue-protein evidence on the same design.
  Returns a `pdactrace_user_projection` list with each intermediate
  object (`rna_fit`, `rna_pattern`, `prot_fit`, `prot_pattern`,
  `evidence`, `audit`) plus a one-row `summary` of audit-class
  counts. Designed to make "apply pdactrace to my own staged
  cohort" a single line.

## New — single-patient trajectory alignment

- `align_patient_profile(rna_logfc, prot_logfc = NULL, ...)` —
  sample-level analog of `score_trajectory()`. Computes
  audit-weighted Pearson rho between the patient's per-gene log2FC
  vector and each of the three atlas stage axes (`rna_beta_E`,
  `rna_beta_M`, `rna_beta_L`), with Fisher-z 95% CI, plus a
  secondary audit-weighted vote share that surfaces
  disagreement-with-rho as a diagnostic. Returns a
  `pdactrace_patient_alignment` list with `$rna` (4-row data.table),
  `$prot` (categorical concordance, optional), `$summary`
  (one-line), and `$attrs` (provenance). Misalignment with all three
  stage axes is a legitimate output and does not imply Normal.
  Discipline-aligned: deliberately not a supervised classifier — an
  RF / GLM trained on TCGA stage labels would break the package's
  "no supervised fitting under limited ground truth" core principle.

## New — Bioconductor-native input and atlas view

- `fit_stage_de(object, ...)` and `fit_stage_de_protein(object, ...)`
  are S4 generics with two methods each. The `ANY` method preserves
  the original matrix / data.frame signature
  (`fit_stage_de(counts, stage, cohort)`); a new
  `SummarizedExperiment` method accepts
  `fit_stage_de(se, stage_col, cohort_col, assay_name = "counts")`
  (and `assay_name = "intensity"` for the protein side). Both paths
  share the same DESeq2 LRT (RNA) / limma (protein) kernel and
  produce numerically identical output (verified by
  `tests/testthat/test-se-input.R`).
- `as_summarized_experiment(reference = NULL)` — converts the bundled
  10,113-gene atlas into a `SummarizedExperiment` with two assays
  (`rna_beta`, `rna_lfcSE`; both 10,113 × 4 stages), a 4-row
  `colData` (`stage`, `reference_level`), a ~109-column `rowData`
  (audit components, RNA / protein trajectory pattern, translation
  class, scRNA cell origin, serum direction, meta-analysis I², ...),
  and atlas-provenance `metadata`. The bundled object remains a
  `data.table` for fast query-based use cases; this constructor is
  the Bioc-native view on demand.

## New — bundled toy cohort and provenance helpers

- `toy_counts` (50 × 24 integer matrix; 10 stage-progressive
  rows), `toy_coldata` (24 rows: stage, cohort, sample), and
  `toy_protein` (50 × 24 log2 matrix; 8 stage-progressive rows).
  Generated by `data-raw/build_toy_data.R` with `set.seed(2026)`.
  The user-cohort vignette uses these in place of inline synthetic
  generation so examples run end-to-end with no external download.
- `pdactrace_data_sources` — bundled `data.table`, one row per
  public dataset (GEO / PRIDE / PDC / ICGC / GTEx / MassIVE / SRA)
  contributing to the v0.99 atlas (26 rows: 6 RNA + 3 Protein +
  10 scRNA + 3 Serum + 2 Pancreatitis + 2 Validation). Generated by
  `data-raw/build_data_sources.R` from the manifest TSVs.
- `list_data_sources(layer = NULL)` — accessor for the data-source
  table, with optional layer filter.
- `atlas_provenance()` — one-call provenance dossier returning
  package version, atlas snapshot, both Zenodo DOIs, both repo URLs,
  cohort count, and per-layer dataset count.

## New — short function aliases

Five exported aliases alongside their fully-spelled originals to
reduce typing friction. Both forms refer to the same function and
share one Rd page; existing scripts using the long names continue
to work unchanged.

| Long name (kept) | Short alias (new) |
|---|---|
| `evaluate_anchor_enrichment()` | `anchor_enrichment()` |
| `early_pattern_names()` | `early_patterns()` |
| `mid_pattern_names_excluded()` | `mid_patterns()` |
| `align_patient_profile()` | `align_patient()` |
| `classify_protein_trajectory()` | `classify_prot_trajectory()` |

## Changed

- `plot_gene_hexagon()` migrated from bare `theme_void()` to
  `pdactrace_panel_theme()` so all 10 plot functions now share the
  same NCS-grade typography (single-column 88 mm,
  `NCS_W_15COL = 4.72 in`, `NCS_W_DOUBLE = 7.08 in`).
  `geom_text` sizes reduced to 1.7–2.2 mm (≈5–6.5 pt) so the
  hexagon embeds cleanly at single-column width.
- `pdactrace_save()` replaces `cat()` with `message()` so package
  diagnostics route through the standard message stream
  (BiocCheck-compliant).
- The internal `.audit_mc_table()` Monte-Carlo helper replaces the
  manual `set.seed()` + `on.exit()` save/restore block with a
  single `withr::local_seed(seed)` call. The seed parameter is
  scoped to the helper and the global RNG state is restored on
  function exit — same behaviour, BiocCheck-compliant.
- README adds a horizontal 6-panel workflow figure at the top
  (`man/figures/pdactrace_overview.jpg`, 1800 × 1344, 280 KB,
  hand-drawn, license-clean): staged omics evidence → 12-template
  trajectory matching → Early-onset atlas surface → multi-layer
  evidence integration → 3-axis + 2-gate audit scoring → user
  outputs (`query_gene` · `explain_score` · `compare_candidates` ·
  `trace_filters` · `project_user_cohort`). Closing tagline:
  "Transparent prioritization, not a supervised diagnostic
  classifier."

## DESCRIPTION

- `Version: 0.99.0`. `Authors@R` is canonical; the redundant
  `Author:` / `Maintainer:` lines are omitted.
- `Depends: R (>= 4.5.0)` (Bioconductor 3.22's minimum supported R).
- `Imports:` data.table, ggplot2, patchwork, glue, methods, stats,
  utils, withr, SummarizedExperiment.
- `Suggests:` DESeq2, limma, S4Vectors, testthat (>= 3.0.0), knitr,
  rmarkdown, pkgdown.
- `LazyData: true`, `LazyDataCompression: xz`, `Roxygen: list(markdown
  = TRUE)`, `Config/testthat/edition: 3`.

## Tests

- `tests/testthat/test-se-input.R` — matrix vs SummarizedExperiment
  numerical equivalence for both `fit_stage_de` (DESeq2) and
  `fit_stage_de_protein` (limma), plus informative-error coverage
  on missing `colData` columns / `assay_name`.
- `tests/testthat/test-explain-compare-se.R` — `explain_score()`
  structure + missing-gene error, `compare_candidates()` row count +
  sort order + missing-gene padding, `as_summarized_experiment()`
  shape + reference-column invariant.
- `tests/testthat/test-marker-class-stability.R` — pins LGALS3BP
  (high_confidence), LTBP1 (supported_uncertain), ALB (penalized),
  GAPDH (excluded); verifies `compare_candidates()` ordering on the
  four case studies; round-trips `toy_counts` through
  `fit_stage_de()` + `classify_trajectory()`; checks
  `list_data_sources()` shape + filter and `atlas_provenance()`
  Zenodo-DOI invariants.
- Total testthat coverage at v0.99.0: **251 PASS**.

## Verification (R 4.5.2 conda, matched libs)

- `R CMD check --as-cran`: 0 ERROR, 1 WARN, 6 NOTE.
  The single WARN is `qpdf needed for PDF size checks` — env-only
  (the Bioconductor build farm has qpdf installed); no substantive
  WARN remains.
- `BiocCheck(new-package = TRUE)`: **0 ERROR, 0 WARN, 15 NOTE**.

# pdactrace 0.4.0

**User-data extension track + dual-classification cleanup** — additive,
backwards-compatible API that lets a new user run the v0.3.0 audit
framework on their own multi-omics cohort. The frozen 3-axis + 2-gate
rule, the audit labels, and all v0.3.0 evaluation results remain
unchanged.

## Removed (breaking schema change)

- v0.2.0 confidence-tier columns (`confidence_tier`,
  `early_onset_score`, `heterogeneity_factor`) are dropped from the
  bundled atlas. They were superseded by v0.3.0 `audit_class` and
  retained as legacy in v0.3.0 for ablation; the dual-classification
  was creating reviewer-grade ambiguity (e.g., LGALS3BP scoring
  `audit_class = high_confidence` for strong multi-layer convergence
  while the v0.2 tier flagged it `ARTIFACT` for low single-layer
  effect size). v0.3.0 `audit_class` is now the canonical
  per-gene classification.
- `list_candidates()` no longer accepts `min_tier`. Use the new
  `min_audit_class` parameter instead. Default changed from
  `"SILVER"` to `"ALL"` (no audit-class filter applied) — explicit
  filtering by audit_class is the intended workflow.
- The audit-framework numerical inputs `max_abs_beta_meta` and
  `max_I2_meta` are retained because they feed the
  `heterogeneity_gate` and rescue-eligibility logic.
- Atlas dimensions: 10,113 × 116 → 10,113 × **113** columns.

## New

- `plot_gene_hexagon(gene_symbol, comparison)` — compact 6-axis
  hexagonal radar chart of a gene's audit-feature profile
  (Multi-layer / Direction / Stage-onset / Serum bridge / Leakage
  safety / Cohort consistency). Multi-gene overlay and optional
  comparison polygon (e.g., `comparison = "high_confidence_mean"`)
  are supported; designed as a high-quality NCS-grade per-gene
  summary card.
- `fit_stage_de_protein(intensity, stage, cohort)` — limma-based
  parallel of `fit_stage_de` for user-supplied tissue-protein
  intensity matrices. Returns the same column schema so downstream
  trajectory classification works for either layer.
- `classify_protein_trajectory(fit, rho_cutoff)` — protein-side
  wrapper of `classify_trajectory` that renames the output column
  from `rna_pattern` to `prot_pattern` for assembly downstream.
- `assemble_user_evidence(rna_fit, prot_fit, scrna_summary,
  serum_summary, signal_peptide, cross_cohort_agreement, max_I2)` —
  combines optional per-layer user inputs into a per-gene evidence
  table conforming to the schema columns the bundled atlas exposes;
  layers the user does not supply are filled with `NA` and degrade
  gracefully through the existing feature constructor.

## Modified

- **12-template trajectory catalog** (was 8-template in v0.3.0). Adds
  Late × 2 (`Late_Burst_Up`, `Late_Loss_Down` — onset at Stage III/IV;
  Peak/Trough/Plateau are structurally degenerate at the last stage)
  and pan-stage Monotonic × 2 (`Monotonic_Up`, `Monotonic_Down` —
  linear progression with maximum at Late, distinct from
  `Mid_Plateau_*` which has its largest jump at Stage II). Atlas
  surface (`rna_pattern`) still exposes only Early × 4 by design; the
  other 8 templates compete in the matching step and are flagged via
  `excluded_mid_pattern`, `excluded_late_pattern`,
  `excluded_monotonic_pattern`. Effect on numbers: GAPDH and CDH13
  correctly reclassified as `Monotonic_Up` and dropped from the Early
  surface (housekeeping/leaky-plasma genes that previously contaminated
  Early × 4); RNA × Protein 4-Early-concordant set 1,424 → 1,356; full
  Tier 1 gold (any 12-template match) 1,449 → 1,375. Headline anchor
  enrichment numbers are invariant (deterministic top-100 secondary:
  7 hits / 39.3× / p = 2.18e-10 unchanged from v0.3.0; LOO median 41.6,
  bootstrap 95% CI [20.2, 56.3]).
- Atlas schema gains three columns `rna_wald_padj_E/M/L` — genuine
  per-contrast Wald padj (Early/Mid/Late vs Normal) from a
  `nbinomWaldTest()` refit of the same `~ dataset + stage_group` model
  used for the LRT. `rna_padj_E/M/L` remain in the schema for
  back-compat but are LRT padj copies (DESeq2 LRT mode returns the
  same omnibus padj for every per-contrast `results()` call). The
  audit framework continues to use `rna_lrt_padj` as its significance
  gate, so no v0.3.0 audit numbers change. `query_gene_detailed()`
  $per_stage now exposes both `lrt_padj` and `wald_padj` columns;
  `plot_stage_effect()` uses Wald padj for the per-stage `*` marks.
- `compute_audit_score(gene_symbol, evidence = NULL)` — new
  `evidence` parameter accepts a user-supplied evidence table (e.g.,
  the output of `assemble_user_evidence`). When `NULL` (default),
  the bundled v0.3.0 reference atlas is used; behaviour and numbers
  are bit-for-bit identical to v0.3.0.
- `propagate_uncertainty(gene_symbol, n_mc, seed, evidence = NULL)` —
  same evidence override; stored MC summaries are unavailable for
  user evidence and `n_mc` must be supplied in that case.
- `evaluate_anchor_enrichment(top_n, tier, score_col, evidence = NULL,
  anchors = NULL)` — both the evidence table and the anchor table
  can be overridden so a user can run anchor enrichment against
  their own evidence and their own anchor reference set.
- `.get_reference(reference = NULL)` — internal helper now accepts
  an optional reference override; default behaviour unchanged.
- `.audit_genes(genes, reference = NULL)` — internal helper takes an
  optional reference for membership checks against user evidence.

## New vignette

- `user_cohort_extension.Rmd` — five-step walkthrough on a small
  synthetic example showing user RNA + protein → assemble_user_evidence
  → compute_audit_score → evaluate_anchor_enrichment with a
  user-supplied anchor list. Runs in <60 s under
  `devtools::check()`.

## New tests

- `tests/testthat/test-user-extension.R` — round-trip on a fixed
  synthetic example: schema sanity for both fit functions, full
  audit-score column schema after `assemble_user_evidence`, and a
  backwards-compat assertion that
  `compute_audit_score(evidence = NULL)` matches the bundled v0.3.0
  result.

## Scope boundary

The PDAC reference atlas remains the default reference; v0.4 does
not provide a `build_user_atlas()` function or layer plugin
architecture. Disease-agnostic atlas building is planned for v0.5
(see Discussion §6 of the BiB manuscript).

## Suggests

- `limma` is added to `Suggests` so `fit_stage_de_protein()` can be
  loaded conditionally; users who do not need protein input do not
  need to install limma.

# pdactrace 0.3.0

**Audit framework track** — replaces opaque marker tiers with a
transparent, pre-specified evidence aggregation framework for PDAC
early-marker hypothesis prioritization. The deterministic rule
summarises each gene along **three evidence axes**
(`evidence_strength`, `biological_coherence`,
`translational_relevance`), applies **two reliability gates**
(`leakage_gate`, `heterogeneity_gate`), and assigns it to one of
interpretable audit labels (`high_confidence`, `supported_uncertain`,
`penalized`, `excluded`, `low`). This release does **not** train a
supervised biomarker predictor; external anchors are used only for
post-freeze enrichment evaluation.

## New

- `build_evidence_graph(gene)` — returns a descriptive evidence graph
  for a gene, including layer presence, direction agreement,
  early-pattern support, serum bridge, leakage gate, heterogeneity
  gate, and final audit score.
- `extract_graph_features(gene)` — exposes the deterministic feature
  components that compose the three evidence axes (mapping in
  supplement Table S1).
- `compute_audit_score(gene)` — computes the frozen v0.3.0 score
  `audit_score = positive_score × leakage_gate × heterogeneity_gate`
  with `positive_score = 0.40 × evidence_strength + 0.35 ×
  biological_coherence + 0.25 × translational_relevance`,
  normalized to the atlas-wide maximum, and assigns an `audit_class`
  label per gene.
- `propagate_uncertainty(gene)` — returns Monte Carlo score/rank
  intervals and a `confidence_class` (stable_high, high_uncertain,
  medium, low, excluded) reported alongside the deterministic
  `audit_class`.
- `evaluate_anchor_enrichment(top_n, tier)` — evaluates top-N
  enrichment of external PDAC biomarker anchors against the frozen
  audit ranking; anchors are never used for training or weight tuning.
- `case_study(genes)` — combines v0.2 tier, deterministic
  `audit_class`, MC uncertainty, anchor provenance, and the
  composite `pdactrace_call` label for manuscript panels.

## Reference DB schema (84 → 111 columns)

27 populated `audit_*` columns:

- **3 evidence axes (new in 3+2 reparameterization)**:
  `audit_evidence_strength`, `audit_biological_coherence`,
  `audit_translational_relevance`.
- **2 reliability gates (new)**: `audit_leakage_gate`,
  `audit_heterogeneity_gate`.
- **Audit-label assignment (new)**: `audit_class`
  (`high_confidence` / `supported_uncertain` / `penalized` /
  `excluded` / `low`).
- **Internal feature components** (kept for supplement
  reproducibility): `audit_score_layer`, `audit_score_direction`,
  `audit_score_early`, `audit_score_serum`, `audit_score_rescue`,
  `audit_positive_score`, `audit_leakage_mult`,
  `audit_heterogeneity_mult`, `audit_score_raw`, `audit_score`,
  `audit_is_housekeeping`, `audit_is_plasma_high_abundance`,
  `audit_rescue_eligible`.
- **Monte Carlo uncertainty**: `audit_score_median`,
  `audit_score_lo95`, `audit_score_hi95`, `audit_uncertainty_width`,
  `audit_rank_median`, `audit_rank_lo95`, `audit_rank_hi95`,
  `audit_confidence_class`.

## Audit validation

- External anchor enrichment is strong under the frozen 3+2 rule
  with deterministic `audit_score` ranking: secondary Tier 1+2
  anchors show **7 hits in the top 100** (`39.3×`, hypergeometric
  `p = 2.2e-10`; hits = MSLN, THBS2, CEACAM1, TIMP1, MUC16,
  SERPINA1, ANPEP) — strengthened from the prototype 7-feature
  flat-weight implementation (5 hits / 28.1× / `p = 6.6e-07`).
- Robustness under the v0.3.0 rule: leave-one-anchor-out median
  fold 41.6× (range 35.7–41.6×); bootstrap median 40.5× with 95% CI
  [20.2, 56.3] (1,000 iterations). MC-median ranking
  (`audit_score_median`) gives 5 hits / 28.1× and is reported as a
  sensitivity ranking (Methods §7).
- Negative rejection is explicit: housekeeping genes receive a hard
  `leakage_gate = 0` (audit_class = `excluded`), and high-abundance
  plasma proteins receive `leakage_gate = 0.5`
  (audit_class = `penalized`).
- Case-study behavior matches the intended framing:
  `LGALS3BP` and `THBS2` are `high_confidence`, `LTBP1` is
  `supported_uncertain` (heterogeneity-gated), `ALB` is
  `penalized`, and `GAPDH` is `excluded`.

## Manuscript figures (PDAC_졸업final5/)

- F1 — overview infographic (3 axes + 2 gates + audit labels).
- F2 — Sankey from v0.2 tier to audit labels + anchor target
  board (39.3×, p = 2.2e-10) + housekeeping rejection traffic
  light + audit-label showcase + robustness badges.
- F3 — multi-layer evidence aggregation (per-cohort RNA forest,
  RNA × Protein direction, scRNA cell origin, serum bridge,
  cross-layer support).
- F4 — top candidate landscape (4 mini-cards, top-30 ranked
  table, mechanism distribution, rescue showcase, layer
  coverage matrix).
- F5 — clinical demonstration (LTBP1+SERPINA1 LOOCV ROC,
  illustrative) + R API quick-start + anchor provenance + tool
  ecosystem badge.

# pdactrace 0.1.1

**Granularity track** — exposes hidden per-stage / per-cohort /
per-celltype / per-filter-step detail that v0.1.0 only summarized
at the headline level. No breaking changes: all v0.1.0 functions
preserved; new functions are *additive*.

## New

- `query_gene_detailed(gene)` — returns 5-slot `data.table` list:
  `per_stage`, `per_cohort`, `per_celltype` (all 11), `filter_diag`,
  `serum_per_cohort`. Stouffer p, agreement %, specificity tau,
  PDAC-vs-Pan t-pval exposed as table attributes.
- `summarize_gene_evidence(detail = TRUE)` — adds per-stage,
  per-cohort, pancreatitis raw means, and filter underlying-metric
  block to text summary.
- 4 new plot functions:
  * `plot_stage_effect()` — log2FC ± 95% CI forest with significance
    asterisks per stage.
  * `plot_per_cohort()` — TCGA / CPTAC / GSE224564 / GSE79668 trend
    + monotonic indicator + Stouffer summary annotation.
  * `plot_filter_diagnostics()` — 7-step filter bar with underlying
    metric annotation (e.g., `pool padj = NA`, `pan-vs-HC pval = 0.017`).
  * `plot_celltype_full()` — full 11-celltype expression bar with
    PDAC-relevant cell types highlighted in red.

## Reference DB schema (45 → 62 columns)

13 new populated columns:
- `rna_padj_E/M/L`, `rna_lfcSE_E/M/L` (per-stage)
- `rna_stouffer_p`, `rna_stouffer_padj`
- `rna_per_cohort_trend`, `rna_per_cohort_monotonic` (list-cols)
- `cell_specificity_tau`
- `ann_pan_vs_hc_logfc`, `ann_pan_vs_hc_pval`,
  `ann_pan_excluded_phase60`
- `ann_pdac_mean`, `ann_pan_mean`, `ann_hc_mean`,
  `ann_pdac_vs_pan_pval`

## Audit findings exposed by v0.1.1

- **LTBP1**: 4-cohort RNA divergence (TCGA Increasing + monotonic vs
  CPTAC Decreasing + monotonic vs 2 non-monotonic; Stouffer p=0.51,
  50% agreement). The headline `Early_Burst_Up` label hides this.
- **SERPINA1**: passes all 7 phase60 filters (only gene to do so
  among the 4 panel members). Per-stage ↑ +1.05/+0.76/+0.72.

## Tests

- 100 `testthat` passes (current 59 + 41 new).

# pdactrace 0.1.0

Initial public release of **PDAC-TRACE** — queryable PDAC stage-aware
multi-omics gene-level lookup atlas.

> **Central message**: *PDAC tissue biomarker ≠ always serum-up
> biomarker.* Tissue-derived candidates can preserve, invert, or
> decouple when projected into serum.

## Framework

Originally introduced as a 4-Early-pattern trajectory framework
(v0.1.0). Widened in v0.4.0 to a **12-template competitive catalog**
(Early × 4 + Mid × 4 + Late × 2 + Monotonic × 2) where each gene
must out-compete 11 alternatives before being surfaced as Early. The
atlas surface (`rna_pattern`) restricts visible calls to:

* `Early_Burst_Up` — Normal low → Early up → sustained
* `Early_Loss_Down` — Normal high → Early down → sustained
* `Early_Peak` — Early peak → Mid/Late decline
* `Early_Trough` — Early trough → Mid/Late recovery

Mid / Late / Monotonic best-matches are flagged via
`excluded_mid_pattern` / `excluded_late_pattern` /
`excluded_monotonic_pattern` for transparent provenance.

## Features

- **Reference atlas** — `pdactrace_reference` (10,113 genes × 45
  cols) bundling 8 source phase result CSVs across bulk RNA-seq,
  tissue proteomics, scRNA cell origin, serum proteomics,
  pancreatitis context, and `phase60` 7-step filter audit trail.
- **Lookup API** — `query_gene()`, `query_panel()`,
  `list_candidates()` (multi-param), `trace_filters()`,
  `summarize_gene_evidence()`, `list_atlas_metadata()`.
- **Trajectory scoring** — `fit_stage_de()`, `classify_trajectory()`,
  `score_trajectory()` for users with their own bulk RNA-seq.
- **Visualization** — `plot_gene_evidence()`, `plot_filter_trace()`,
  `plot_panel_heatmap()`, `plot_candidate_landscape()`.
- **Vignettes** — `01_lookup_basics`, `02_ltbp1_case_study`.

## LTBP1 narrative (inverse stromal exemplar)

* RNA: `Early_Burst_Up` (rho ≈ 1.0)
* Tissue protein: `Early_Burst_Up` (Tier 1 gold concordant)
* scRNA cell origin: myCAF / iCAF dominant
* Serum: detected, **Class B (inverse)** — PDAC down vs HC,
  pancreatitis-elevated
* Filter audit: `phase77_classB` route — passes phase77 strict +
  Class B + manual curation, while `phase60` strict 7-step funnel
  only captures SignalP

## Design

- **4-Early-only surface, by design** (v0.1.0; in v0.4.0 the catalog
  itself widens to 12 templates and Mid / Late / Monotonic best-
  matches are flagged via `excluded_*_pattern` rather than surfaced —
  see v0.4.0 entry above).
- **Class B is manually curated**, not predicted.
- **Two parallel narrowing funnels** (`phase60` strict vs `phase77`
  strict) can disagree — `trace_filters()` shows both routes
  honestly.

## Deferred

- v0.2 (optional): calibrated detectability classifier
  (`predict_detectability()`).
- Future paper: direction predictor (Class A vs B), inverse marker
  risk modeling, web/Shiny version, cross-cancer generalization.

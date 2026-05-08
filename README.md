# pdactrace

> Stage-aware PDAC multi-omics atlas and transparent evidence-audit
> framework for tissue-to-serum biomarker prioritization.

[![R-CMD-check](https://github.com/jibeomko/pdactrace/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jibeomko/pdactrace/actions/workflows/R-CMD-check.yaml)
[![Version](https://img.shields.io/github/v/release/jibeomko/pdactrace?include_prereleases&sort=semver&label=version&color=blue)](https://github.com/jibeomko/pdactrace/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-submission%20in%20preparation-lightgrey.svg)](https://www.bioconductor.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20076698.svg)](https://doi.org/10.5281/zenodo.20076698)

<p align="center">
  <img src="man/figures/pdactrace_overview.jpg" width="900"
       alt="pdactrace workflow: staged omics evidence (RNA, protein, scRNA, serum, user cohort) → 12-template trajectory matching → Early-onset atlas surface → multi-layer evidence integration → 3-axis + 2-gate audit scoring → user outputs (query_gene, explain_score, compare_candidates, trace_filters, project_user_cohort). Transparent prioritization, not a supervised diagnostic classifier.">
</p>

`pdactrace` is an R package for querying and prioritizing pancreatic
ductal adenocarcinoma (PDAC) biomarker candidates across bulk RNA-seq,
tissue proteomics, single-cell RNA-seq, serum proteomics, and
pancreatitis context.

It does **not** train a supervised biomarker classifier. Instead, it
uses a frozen, interpretable scoring rule and reports uncertainty
because PDAC early detection lacks robust gene-level ground truth.

## What You Get

- A bundled reference atlas: **10,113 genes x 113 columns**
- A **12-template competitive trajectory catalog** with an Early x 4
  atlas surface
- A transparent **3-axis + 2-gate audit score**
- Per-gene lookup, panel lookup, candidate listing, filter tracing,
  and visualization APIs
- User-cohort helpers for applying the same framework to new RNA or
  protein evidence

Core message: **a PDAC tissue biomarker is not always a serum-up
biomarker.** Tissue signals can preserve, invert, or decouple when
projected into serum.

## Install

```r
# install.packages("remotes")
remotes::install_github("jibeomko/pdactrace")
```

A Bioconductor submission is in preparation.

## Quick Start

```r
library(pdactrace)

# Look up one gene
query_gene("LGALS3BP")

# Compare a small panel
query_panel(c("LGALS3BP", "SERPINA1", "ALB", "GAPDH"))

# Inspect stage and cohort detail
query_gene_detailed("SERPINA1")$per_stage
query_gene_detailed("SERPINA1")$per_cohort

# List candidates
list_candidates(translation_class = "inverse")
list_candidates(onset = "Early", tissue_direction = "Up")

# Trace why selected genes pass or fail the tissue-to-serum funnel
trace_filters(c("SERPINA1", "SPARC", "CDH13", "GAPDH"))
```

## Audit Scoring

Each gene is scored with three evidence axes and two reliability gates:

```text
positive_score = 0.40 * evidence_strength
               + 0.35 * biological_coherence
               + 0.25 * translational_relevance

audit_score = positive_score * leakage_gate * heterogeneity_gate
```

The score is deterministic; Monte Carlo summaries report rank stability
under evidence perturbation. External anchors are used for post-freeze
evaluation only, not training.

| Audit label | Typical meaning |
|---|---|
| `high_confidence` | Strong multi-layer support with clean gates |
| `supported_uncertain` | Supported, but cohort-heterogeneous or moderate score |
| `penalized` | Plasma high-abundance penalty |
| `excluded` | Housekeeping or leakage artifact |
| `low` | Insufficient evidence |

```r
compute_audit_score(c("LGALS3BP", "SERPINA1", "ALB", "GAPDH"))
propagate_uncertainty(c("LGALS3BP", "SERPINA1", "GAPDH"))
anchor_enrichment(top_n = 100, tier = "secondary")  # alias of evaluate_anchor_enrichment
```

The frozen audit rule recovers **7 secondary-tier external anchors in
the top 100** (`39.3x`, hypergeometric `p = 2.18e-10`).

## Stage Harmonization

Every cohort is collapsed onto a single **4-level scale**
(`Normal / Early / Mid / Late`) before any per-gene trajectory or
audit step. The collapse is deterministic and applied identically
across bulk RNA-seq, tissue proteomics, and serum proteomics.

| Source label | Standardized | 4-level group |
|---|---|---|
| Histologically normal / matched-normal tissue | — | `Normal` |
| AJCC Stage I (IA, IB) | `Stage_I` | **Early** |
| AJCC Stage II (IIA, IIB) | `Stage_II` | **Mid** |
| AJCC Stage III | `Stage_III` | **Late** |
| AJCC Stage IV (M1) | `Stage_IV` | **Late** |
| TNM with M1 | → `Stage_IV` | **Late** |
| TNM with N1+ (no M1) | → `Stage_III` | **Late** |
| TNM T3 / T4 (M0, N0) | → `Stage_II` | **Mid** |
| TNM T1-T2 (M0, N0) | → `Stage_I` | **Early** |
| Cohort-specific category "resectable" | → `Stage_II` | **Mid** |
| Cohort-specific category "borderline" / "locAdvanced" | → `Stage_III` | **Late** |
| Cohort-specific category "metastatic" | → `Stage_IV` | **Late** |

Reference: AJCC Cancer Staging Manual (8th edition) for PDAC.
Source-of-truth code: [analysis/transcriptomics/scripts/figure1_wald_lrt_pipeline.R](analysis/transcriptomics/scripts/figure1_wald_lrt_pipeline.R)
in the [companion manuscript repo](https://github.com/jibeomko/PDAC_biomarker)
(lines 93–96 for the AJCC → 3-level mapping; the matching
`extract_all_clinical.R` carries the per-cohort TNM and free-text
normalisation).

**Why 3-level (Early / Mid / Late) instead of 4-level (I / II / III / IV)?**
Per-cohort AJCC IV cell counts are too small for stable per-stage
contrasts; collapsing III + IV into a single `Late` group preserves
statistical power without erasing the resectable (Early) vs
locally-advanced-or-metastatic (Late) distinction that drives the
clinical question. A 4-level **AJCC sensitivity analysis** is
preserved in parallel (`stage_ajcc4`, `stage_numeric` columns) so
reviewers can re-run the same per-stage tests at full granularity.

For applying the same harmonization to your own cohort, see
`vignette("user_cohort_extension")` — `fit_stage_de(stage = ...)`
expects the 4-level factor with `levels = c("Normal", "Early", "Mid",
"Late")`.

## Trajectory Framework

Each gene's Normal/Early/Mid/Late trajectory is matched against 12
pre-declared templates:

| Family | Templates | Atlas surface? |
|---|---|---|
| **Early x 4** | `Early_Burst_Up`, `Early_Loss_Down`, `Early_Peak`, `Early_Trough` | surfaced (`rna_pattern`) |
| **Mid x 4** | `Mid_Plateau_Up`, `Mid_Plateau_Down`, `Mid_Peak`, `Mid_Trough` | flagged via `excluded_mid_pattern` |
| **Late x 2** | `Late_Burst_Up`, `Late_Loss_Down` | flagged via `excluded_late_pattern` |
| **Monotonic x 2** | `Monotonic_Up`, `Monotonic_Down` | flagged via `excluded_monotonic_pattern` |

Only Early x 4 calls are surfaced in `rna_pattern` / `prot_pattern`.
Non-Early best matches remain visible as exclusion flags. This makes
Early calls stricter: a candidate must beat Mid, Late, and Monotonic
alternatives before being surfaced.

## Reference Atlas

The bundled `pdactrace_reference` object is a `data.table` with
**10,113 genes x 113 columns**. It includes RNA trajectory evidence,
tissue-protein support, scRNA cell-origin summaries, serum translation
features, pancreatitis context, audit scores, uncertainty summaries,
and provenance.

Important bundled data:

| Object | Role |
|---|---|
| `pdactrace_reference` | Main 10,113-gene atlas |
| `default_templates` | Canonical 12-template trajectory catalog |
| `pdactrace_external_anchors` | External anchor set for post-freeze evaluation |
| `meta_analysis` | Random-effects RNA meta-analysis summaries |

## Common Plots

```r
plot_gene_evidence("LGALS3BP")
plot_filter_trace(c("SERPINA1", "SPARC", "CDH13", "GAPDH"))
plot_candidate_landscape()
plot_stage_effect("SERPINA1")
plot_per_cohort("SERPINA1")
plot_meta_forest("SERPINA1", contrast = "Mid_vs_Early")
plot_gene_hexagon("LGALS3BP", comparison = "high_confidence_mean")

# Per-template trajectory atlas (12 PDFs per layer; v0.99.1)
plot_template_atlas("rna",     output_dir = tempdir())
plot_template_atlas("protein", output_dir = tempdir())

# Per-gene trajectory overlaid on its matched template panel
plot_gene_template("LTBP1",    layer = "rna")
plot_gene_template("LGALS3BP", layer = "protein")
```

## Use Your Own Cohort

The v0.4.0 API lets users apply the same trajectory and audit framework
to staged RNA or tissue-protein data:

```r
rna_fit <- fit_stage_de(counts, stage, cohort)
rna_pat <- classify_trajectory(rna_fit)

prot_fit <- fit_stage_de_protein(intensity, stage, cohort)
prot_pat <- classify_prot_trajectory(prot_fit)  # alias of classify_protein_trajectory

evidence <- assemble_user_evidence(rna_fit = rna_pat, prot_fit = prot_pat)
compute_audit_score(evidence = evidence)
```

## Function name aliases (v0.4.1)

The following short aliases are exported alongside their fully-spelled
originals; both forms refer to the same function.

| Long name | Short alias |
|---|---|
| `evaluate_anchor_enrichment()` | `anchor_enrichment()` |
| `early_pattern_names()` | `early_patterns()` |
| `mid_pattern_names_excluded()` | `mid_patterns()` |
| `align_patient_profile()` | `align_patient()` |
| `classify_protein_trajectory()` | `classify_prot_trajectory()` |

Use whichever you prefer — existing scripts using the long names continue
to work unchanged.

## Function reference

All 59 exported functions and constants, grouped by role. Long-name
aliases listed in the v0.4.1 alias table above are not repeated here.

### Lookup + summary

| Function | What it does |
|---|---|
| `query_gene(gene)` | Single-gene full evidence dump (RNA, protein, scRNA, serum, clinical, filter, annotation) |
| `query_gene_detailed(gene)` | Per-stage / per-cohort / per-cell-type / per-step breakdown |
| `query_panel(genes)` | Multi-gene join on the same atlas columns |
| `summarize_gene_evidence(gene, detail)` | Human-readable text summary |
| `list_candidates(...)` | Criterion-based candidate filter (onset, direction, translation_class, ...) |
| `list_atlas_metadata()` | Bundled `atlas_metadata` (version, build date, cohort manifest) |
| `case_study(name)` | Pre-configured case-study panel (e.g., LTBP1, LGALS3BP) |

### Audit framework

| Function | What it does |
|---|---|
| `compute_audit_score(genes, evidence)` | 3-axis × 2-gate audit score + class |
| `propagate_uncertainty(genes, n_mc, seed)` | Monte-Carlo rank stability summary |
| `evaluate_anchor_enrichment(top_n, tier, score_col)` | External-anchor enrichment (frozen post-freeze evaluation) |
| `anchor_enrichment(...)` | Short alias of the above |
| `explain_score(gene)` | Plain-English decomposition of the audit score (3 axes + 2 gates + reason) |
| `compare_candidates(genes)` | Side-by-side comparison table (audit class, pattern, cell origin, redundancy hint) |
| `format_provenance(provenance, style)` | Human-readable evidence labels for `phaseXX` tags |

### Trajectory framework

| Function | What it does |
|---|---|
| `fit_stage_de(object, ...)` | DESeq2 LRT wrapper for stage-aware DE (matrix or `SummarizedExperiment`) |
| `fit_stage_de_protein(object, ...)` | limma parallel for tissue-protein intensity |
| `classify_trajectory(fit)` | Best-match against the 12-template catalog |
| `classify_protein_trajectory(fit)` | Protein-side wrapper (`classify_prot_trajectory` alias) |
| `score_trajectory(pat, gene)` | 12-template Pearson rho vector per gene |
| `assemble_user_evidence(...)` | Combine optional per-layer user inputs into a per-gene evidence table |
| `align_patient_profile(rna_logfc, ...)` | Sample-level alignment of one patient's profile against atlas stage axes |
| `align_patient(...)` | Short alias |
| `project_user_cohort(rna, coldata, ...)` | End-to-end wrapper around `fit_stage_de → classify_trajectory → assemble_user_evidence → compute_audit_score` |
| `early_pattern_names()` / `early_patterns()` | Names of the 4 surfaced Early-onset templates |
| `mid_pattern_names_excluded()` / `mid_patterns()` | Names of the 4 Mid templates excluded from atlas surface |

### Filter audit

| Function | What it does |
|---|---|
| `trace_filters(genes)` | 7-step filter audit + class route per gene |

### Visualization

| Function | What it does |
|---|---|
| `plot_gene_evidence(gene)` | Multi-panel composite (trajectory + cell origin + serum + summary) |
| `plot_gene_hexagon(gene, comparison)` | 6-axis evidence radar |
| `plot_stage_effect(gene)` | Per-stage forest with log2FC ± 1.96·SE |
| `plot_per_cohort(gene)` | Cohort-by-cohort trend bar plot |
| `plot_meta_forest(gene, contrast)` | Random-effects meta-analysis forest |
| `plot_filter_trace(genes)` | Pass/fail step bar across the 7-step filter |
| `plot_filter_diagnostics(...)` | Per-step diagnostic counters across the atlas |
| `plot_panel_heatmap(genes)` | Gene × evidence-axis comparison heatmap |
| `plot_candidate_landscape(...)` | Tissue × serum scatter, Class A/B coloured |
| `plot_celltype_full(...)` | Full cell-type-of-origin overview from the scRNA atlas |
| `plot_template_atlas(layer, output_dir)` | 12 PDFs per layer (RNA + protein) showing each template's cohort |
| `plot_gene_template(gene, layer)` | One PDF: gene's matched template + its trajectory overlaid |

### Reporting

| Function | What it does |
|---|---|
| `report_gene(genes, output_dir)` | Self-contained HTML evidence report (single gene or panel) |

### Bioconductor-native helpers

| Function | What it does |
|---|---|
| `as_summarized_experiment(reference)` | Convert the bundled atlas to a `SummarizedExperiment` (2 assays + `rowData` + `colData` + `metadata`) |
| `atlas_provenance()` | One-call provenance dossier (version, repo URL, both Zenodo DOIs, cohort count, layers) |
| `list_data_sources(layer)` | 26 contributing public datasets with accession + URL |
| `schema_spec()` | Canonical column list for the atlas + T2.5/T3 status |
| `build_evidence_graph(gene)` | Audit-score evidence as a small node/edge graph |
| `extract_graph_features(genes)` | Per-gene 6-axis feature vector (used by `plot_gene_hexagon` and the audit rule) |

### Theme + plotting infrastructure

| Object / function | What it does |
|---|---|
| `pdactrace_axes_theme()` | NCS-grade `ggplot2` theme for plots with axes |
| `pdactrace_panel_theme()` | NCS-grade `ggplot2` theme for schematic panels |
| `pdactrace_save(p, dir, name, w, h)` | Cairo-PDF writer (BiocCheck-clean) |
| `pdactrace_pal_class` | Palette for translation classes (A / B / C / Other) |
| `pdactrace_pal_group` | Palette for HC / Pancreatitis / PDAC sample groups |
| `pdactrace_pal_dir` | Palette for UP / DOWN / NS direction |
| `pdactrace_pal_pattern` | Palette for the 4 Early-surfaced patterns |
| `NCS_W_SINGLE` / `NCS_W_15COL` / `NCS_W_DOUBLE` / `NCS_W_TRIPLE` | Width constants in inches (3.46 / 4.72 / 7.08 / 10.50) |

For every function, see `?function_name` for the full Rd page with
arguments, return shape, and an example.

## Vignettes

```r
vignette("lookup_basics", package = "pdactrace")
vignette("audit_case_studies", package = "pdactrace")
vignette("audit_framework", package = "pdactrace")
vignette("user_cohort_extension", package = "pdactrace")
```

## Reproducibility

The R package is **self-contained for atlas re-derivation.** It
ships with the reference atlas (`data/*.rda`), the canonical
12-template trajectory catalog, the external anchor set, all
testthat suites, the build scripts under `data-raw/`, and — new
in v0.99.3 — the **six small downstream phase tables** under
`inst/extdata/phase{2c,29,42,60,77,80}_*.csv.xz`. The build chain
runs end-to-end from these bundled inputs alone:

```r
# Rebuild the bundled atlas + auxiliary objects from data-raw/
source("data-raw/build_reference.R")
source("data-raw/build_templates.R")
source("data-raw/build_meta_analysis.R")
source("data-raw/build_atlas_metadata.R")
```

The two large upstream CSVs (`phase33_deseq2_coef_12template.csv`
RNA fit, `phase34_protein_pooled_12template.csv` protein fit) are
not bundled in `inst/extdata` because together they push the
tarball past Bioconductor's 5 MB ceiling. The build scripts
fall through a clean lookup chain — `inst/extdata/` →
`$PDAC_BASE_DIR/...` → an informative error pointing to the
manuscript Zenodo archive
([10.5281/zenodo.20067849](https://doi.org/10.5281/zenodo.20067849))
where `phase33_deseq2_coef_12template.csv` and
`phase34_protein_pooled_12template.csv` are available verbatim.

User cohorts can be projected end-to-end through the same audit
framework via [`project_user_cohort()`](#use-your-own-cohort)
without touching `data-raw/` at all — the atlas itself is the
prebuilt deliverable, and the rebuild path is for verifying it
or extending it to new cohorts.

Key bundled data:

| File | Role |
|---|---|
| `data/pdactrace_reference.rda` | Bundled 10,113-gene atlas |
| `data/pdactrace_protein_betas.rda` | Per-stage tissue-protein effect sizes |
| `data/default_templates.rda` | Canonical 12-template catalog |
| `data/pdactrace_external_anchors.rda` | External anchor set |
| `data/meta_analysis.rda` | Random-effects RNA meta-analysis summaries |
| `data/atlas_metadata.rda` | Atlas provenance + cohort manifest |
| `data/pdactrace_data_sources.rda` | 26 contributing public datasets |
| `inst/extdata/phase*.csv.xz` | Six bundled downstream phase tables |

## Citation

If you use `pdactrace`, please cite the software via its Zenodo DOI:

> Ko, J. (2026). *pdactrace: Queryable Stage-Aware PDAC Tissue-to-Serum
> Biomarker Reference Atlas* (v0.99.0) [Software]. Zenodo.
> [10.5281/zenodo.20076698](https://doi.org/10.5281/zenodo.20076698)

```bibtex
@software{pdactrace2026,
  author       = {Ko, Jibeom},
  title        = {{pdactrace: Queryable Stage-Aware PDAC
                   Tissue-to-Serum Biomarker Reference Atlas}},
  year         = 2026,
  version      = {v0.99.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20076698},
  url          = {https://doi.org/10.5281/zenodo.20076698}
}
```

The accompanying *Briefings in Bioinformatics* manuscript reference
will be added once the preprint or journal DOI is available.

## License

MIT. See [LICENSE](LICENSE).

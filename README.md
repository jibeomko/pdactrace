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

## Table of contents

- [What you get](#what-you-get)
- [Install](#install)
- [Step-by-step walkthrough](#step-by-step-walkthrough) — six end-to-end scenarios
- [Audit scoring rule](#audit-scoring-rule)
- [Stage harmonization](#stage-harmonization)
- [Trajectory framework](#trajectory-framework)
- [Reference atlas](#reference-atlas)
- [Function reference](#function-reference) — detailed blocks for major functions + table for the rest
- [Function name aliases](#function-name-aliases-v041)
- [Vignettes](#vignettes)
- [Reproducibility](#reproducibility)
- [Citation](#citation)

## What you get

- A bundled reference atlas: **10,113 genes × 113 columns**
- A **12-template competitive trajectory catalog** with an Early × 4
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

A Bioconductor submission is in preparation; once accepted you will
also be able to do:

```r
# BiocManager::install("BiocManager")
BiocManager::install("pdactrace")  # not yet available
```

To verify your install:

```r
library(pdactrace)
atlas_provenance()
#> $package_version  "0.99.4"
#> $atlas_version    "v0.4.0"
#> $package_doi      "10.5281/zenodo.20076698"
#> $manuscript_doi   "10.5281/zenodo.20067849"
#> ...
```

# Step-by-step walkthrough

Six end-to-end scenarios, in increasing order of complexity. Every
example is copy-paste runnable on the bundled atlas — no external
data download required for scenarios 1–4.

## Scenario 1 — Look up one gene

You have a candidate gene (e.g., `LTBP1`) and want every layer of
evidence the atlas has on it.

```r
library(pdactrace)

ev <- query_gene("LTBP1")
ev
```

What you'll see:

```
LTBP1: Early_Burst_Up (rho=1.00, padj=1.7e-02) | Tier1_gold | serum-detected (n=1 cohort, Class B).
Evidence:  RNA trajectory + Tissue protein trajectory + Multi-cohort RNA consistency
           + scRNA cell origin + 7-step serum filter audit
           + Strict RNA-protein-serum bridge + Predeclared panel member
Technical: phase33, phase34, stouffer_consistency, phase2c, phase60, phase77, phase80

Layers loaded: rna, protein, scrna, serum, clinical, filter_status, annotation
Use $rna / $protein / $scrna / $serum / $clinical /
$filter_status / $annotation for full evidence.
```

The first line is the **headline** (matched template, audit class,
serum status, translation class). The next two lines surface the
**evidence sources** in plain English plus the technical phase IDs
for traceability.

Drill into specific layers:

```r
ev$rna           # per-cohort RNA trajectory + Stouffer consistency
ev$protein       # tissue-protein 12-template match
ev$scrna         # cell-of-origin from the scRNA atlas
ev$serum         # serum direction + 3-class translation label
ev$filter_status # 7-step tissue-to-serum filter route
```

For the per-stage / per-cohort breakdown:

```r
detail <- query_gene_detailed("LTBP1")
detail$per_stage
detail$per_cohort
detail$per_celltype
detail$filter_diag       # 7 rows, one per filter step
detail$serum_per_cohort
```

## Scenario 2 — Understand *why* a gene got its score

The audit score is a frozen, deterministic weighted sum:

```text
positive_score = 0.40 * evidence_strength
               + 0.35 * biological_coherence
               + 0.25 * translational_relevance
audit_score    = positive_score * leakage_gate * heterogeneity_gate
```

[`explain_score()`](#explain_scoregene_symbol-verbose--true) breaks
that formula apart for one gene and tells you in plain English which
gate or axis is driving the class label:

```r
explain_score("LTBP1")
```

Output:

```
LTBP1 lands in `supported_uncertain` with audit_score = 0.577.
  positive_score = 0.40 * 0.85 + 0.35 * 0.82 + 0.25 * 0.66 = 0.792
  audit_score = 0.792 * leakage_gate(1.00) * heterogeneity_gate(0.70) = 0.554
  Positive score is solid (0.79) but the heterogeneity_gate = 0.70
  penalises cohort divergence (max meta I² = 75% in [70%, 90%)).
```

Compare cleanly across the four canonical case studies:

```r
explain_score("LGALS3BP")  # high_confidence — both gates clean
explain_score("LTBP1")     # supported_uncertain — heterogeneity gate penalty
explain_score("ALB")       # penalized — plasma-high-abundance leakage gate
explain_score("GAPDH")     # excluded — housekeeping leakage gate hard-zeroes
```

If you want the breakdown as a structured object (for downstream
filtering or custom plotting), capture the invisible return value:

```r
res <- explain_score("LTBP1", verbose = FALSE)
res$axes      # data.table: axis, weight, value, contribution
res$gates     # data.table: gate, value, triggered_by
res$audit_class
res$audit_score
```

## Scenario 3 — Compare a small panel side-by-side

You have several candidate genes and want to pick the strongest, or
spot redundancy:

```r
cmp <- compare_candidates(c("LGALS3BP", "LTBP1", "SERPINA1",
                              "ALB", "GAPDH"))
cmp
```

The result is a `data.table` sorted by `audit_score` descending:

```
   gene_symbol  audit_class           audit_score   rna_pattern
1: LGALS3BP     high_confidence       1.000         Early_Burst_Up
2: LTBP1        supported_uncertain   0.577         Early_Burst_Up
3: SERPINA1     supported_uncertain   0.430         Early_Loss_Down
4: ALB          penalized             0.240         Early_Trough
5: GAPDH        excluded              0.000         <NA>

   prot_pattern    translation_class   cell_origin_top   serum_detected
1: Early_Burst_Up  A                   Ductal            TRUE
2: Early_Burst_Up  B                   myCAF             TRUE
3: Early_Loss_Down A                   Acinar            TRUE
4: Early_Burst_Up  <NA>                <NA>              FALSE
5: Early_Burst_Up  B                   <NA>              TRUE

   max_I2_meta   redundancy_with
1: 25.5
2: 75.2
3: 32.1
4: 56.8
5: 51.8          (none)
```

`redundancy_with` flags genes that share **both** `rna_pattern` and
`cell_origin_top` — useful when you want a non-redundant 5-gene
panel (LGALS3BP and a hypothetical second Ductal Early_Burst_Up
gene would land here).

Visualize the panel as a single composite radar:

```r
plot_gene_hexagon(c("LGALS3BP", "LTBP1", "SERPINA1"))
```

## Scenario 4 — Generate a one-gene HTML evidence report

For a polished, share-with-collaborator artifact:

```r
fp <- report_gene("LTBP1", output_dir = tempdir())
fp
#> "/tmp/.../LTBP1_pdactrace_report.html"

browseURL(fp)   # open in browser
```

The HTML file (~400 KB, self-contained, no Bootstrap) includes:

- Audit-class colored tag + audit_score
- Audit-component table (3 axes × 2 gates)
- Six-axis evidence radar vs the `high_confidence_mean` reference
- Stage-trajectory plot (Normal / Early / Mid / Late)
- Per-cohort breakdown
- Per-stage detail table
- Tissue-to-serum filter trace
- Atlas-version provenance footer

Multi-gene panel report (one HTML for all genes):

```r
fp <- report_gene(c("LGALS3BP", "LTBP1", "GAPDH"),
                   output_dir = tempdir())
```

The panel template starts with a `compare_candidates()` table at
the top, then a multi-gene evidence radar, then one
`explain_score()`-derived rationale card per gene.

## Scenario 5 — Single-patient trajectory alignment

You have a single patient's tumor-vs-matched-normal log2 fold-change
profile and want to see which atlas stage axis (Early / Mid / Late)
it best matches.

```r
# Synthetic Late-stage profile for demonstration
ref <- pdactrace:::.get_reference()
high_audit <- ref[!is.na(rna_beta_L) & audit_score >= 0.3][seq_len(300)]

set.seed(7)
patient_rna <- setNames(
  high_audit$rna_beta_L + rnorm(nrow(high_audit), sd = 0.3),
  high_audit$gene_symbol)

aln <- align_patient_profile(patient_rna, top_n_genes = 200)
print(aln)
```

Output:

```
Best-match by rho: Late (rho=0.97, p=9.2e-204, 200 genes);
voting agrees: Late (vote_share=0.42).

   stage    cor_to_stage_axis   cor_pval     vote_share   weighted_dist
1: Normal   NA                  NA           0.108        1.057
2: Early    0.949               1.6e-144     0.213        0.319
3: Mid      0.964               2.5e-173     0.255        0.282
4: Late     0.974               9.2e-204     0.424        0.242
```

**This is an alignment readout, not a stage prediction.** No
supervised model is fit; the atlas is frozen and the function is
deterministic given (input, atlas). Misalignment with all three
stage axes is itself a finding (it does **not** imply Normal).

## Scenario 6 — Project your own staged cohort

The most common "use pdactrace on my own data" path. Inputs:

- A count matrix (genes × samples) — integer matrix or `data.frame`,
  or a `SummarizedExperiment` with the counts in its `assay()`.
- A `coldata` `data.frame` with at least `stage` (must contain levels
  from `c("Normal", "Early", "Mid", "Late")`) and optionally
  `cohort`.
- (optional) A matched tissue-protein log2 intensity matrix on the
  same samples.

The package ships a tiny synthetic toy cohort so you can test the
pipeline end-to-end:

```r
data(toy_counts);  data(toy_coldata);  data(toy_protein)
str(toy_counts)
#  int [1:50, 1:24] ...                — 50 genes × 24 samples
str(toy_coldata)
#  'data.frame': 24 obs. of 3 variables — sample, stage, cohort
```

Run the full chain in one call:

```r
res <- project_user_cohort(
  rna        = toy_counts,
  coldata    = toy_coldata,
  stage_col  = "stage",
  cohort_col = "cohort",
  protein    = toy_protein)

res
#> <pdactrace_user_projection>
#>   Input:   50 genes; protein layer: yes
#>   Audit:   50 genes scored
#>   Classes: high=12 supp_unc=8 penalized=2 excluded=0 low=28
```

Drill into intermediates:

```r
res$rna_fit       # data.table from fit_stage_de()  (DESeq2 LRT)
res$rna_pattern   # data.table from classify_trajectory()
res$prot_fit      # data.table from fit_stage_de_protein() (limma)
res$prot_pattern  # protein 12-template match
res$evidence      # assemble_user_evidence() output
res$audit         # compute_audit_score(evidence = ...) output
res$summary       # one-row class-count summary
```

Same chain with a `SummarizedExperiment` (Bioconductor-native):

```r
library(SummarizedExperiment)
se <- SummarizedExperiment(
  assays  = list(counts = toy_counts),
  colData = DataFrame(toy_coldata))

res <- project_user_cohort(rna = se,
                             stage_col = "stage",
                             cohort_col = "cohort")
```

# Audit scoring rule

The frozen v0.3.0 audit rule combines three weighted evidence axes
and two reliability gates:

```text
positive_score = 0.40 * evidence_strength
               + 0.35 * biological_coherence
               + 0.25 * translational_relevance

audit_score = positive_score * leakage_gate * heterogeneity_gate
```

The score is deterministic given (gene, atlas). Monte Carlo
summaries via [`propagate_uncertainty()`](#propagate_uncertaintygenes-n_mc-seed)
report rank stability under evidence perturbation. External anchors
are used for **post-freeze evaluation only**, never for training.

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
anchor_enrichment(top_n = 100, tier = "secondary")
```

The frozen audit rule recovers **7 secondary-tier external anchors
in the top 100** (`39.3×`, hypergeometric `p = 2.18e-10`).

# Stage harmonization

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

Reference: AJCC Cancer Staging Manual (8th edition).

**Why 3-level (Early / Mid / Late) instead of 4-level (I / II / III / IV)?**
Per-cohort AJCC IV cell counts are too small for stable per-stage
contrasts; collapsing III + IV into a single `Late` group preserves
statistical power without erasing the resectable (Early) vs
locally-advanced-or-metastatic (Late) distinction that drives the
clinical question. A 4-level **AJCC sensitivity analysis** is
preserved in parallel (`stage_ajcc4`, `stage_numeric` columns).

For applying the same harmonization to your own cohort: pass the
4-level factor with `levels = c("Normal", "Early", "Mid", "Late")`
to `fit_stage_de()` (or set the corresponding column in the
`coldata` you pass to `project_user_cohort()`).

# Trajectory framework

Each gene's Normal/Early/Mid/Late trajectory is matched against
12 pre-declared templates:

| Family | Templates | Atlas surface? |
|---|---|---|
| **Early × 4** | `Early_Burst_Up`, `Early_Loss_Down`, `Early_Peak`, `Early_Trough` | surfaced (`rna_pattern`) |
| **Mid × 4** | `Mid_Plateau_Up`, `Mid_Plateau_Down`, `Mid_Peak`, `Mid_Trough` | flagged via `excluded_mid_pattern` |
| **Late × 2** | `Late_Burst_Up`, `Late_Loss_Down` | flagged via `excluded_late_pattern` |
| **Monotonic × 2** | `Monotonic_Up`, `Monotonic_Down` | flagged via `excluded_monotonic_pattern` |

Only Early × 4 calls are surfaced in `rna_pattern` / `prot_pattern`.
Non-Early best matches remain visible as exclusion flags. This makes
Early calls stricter: a candidate must beat Mid, Late, and Monotonic
alternatives before being surfaced.

To visualise the cohort that matches each of the 12 templates:

```r
plot_template_atlas("rna",     output_dir = tempdir())
plot_template_atlas("protein", output_dir = tempdir())
```

Each call writes 12 PDFs (one per template) showing the cohort's
z-scored Normal/Early/Mid/Late trajectories as thin lines plus mean
and ±1 SD ribbon. To overlay one specific gene on its matched
template panel:

```r
plot_gene_template("LTBP1", layer = "rna")
```

# Reference atlas

The bundled `pdactrace_reference` object is a `data.table` with
**10,113 genes × 113 columns**: RNA trajectory evidence,
tissue-protein support, scRNA cell-origin summaries, serum
translation features, pancreatitis context, audit scores,
uncertainty summaries, and provenance.

| Object | Role |
|---|---|
| `pdactrace_reference` | Main 10,113-gene atlas |
| `pdactrace_protein_betas` | Per-stage tissue-protein effect sizes (5,917 proteins) |
| `default_templates` | Canonical 12-template catalog |
| `pdactrace_external_anchors` | External anchor set for post-freeze evaluation |
| `meta_analysis` | Random-effects RNA meta-analysis summaries |
| `atlas_metadata` | Atlas provenance + cohort manifest |
| `pdactrace_data_sources` | 26 contributing public datasets |
| `inst/extdata/phase*.csv.xz` | Six bundled downstream phase tables |

# Function reference

Detailed parameter / return / example blocks for the 20 most-used
functions. The full 59-export catalogue is in the
[summary table](#summary-of-all-59-exported-objects) at the end.

## Lookup

### `query_gene(gene_symbol)`

Single-gene full evidence dump across every layer of the atlas.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character(1) | required | HGNC gene symbol |

**Returns:** an `S3` list of class `pdactrace_gene_evidence` with
slots `$rna`, `$protein`, `$scrna`, `$serum`, `$clinical`,
`$filter_status`, `$annotation`, `$summary`, `$provenance`.
`print()` renders the human-readable block shown in
[Scenario 1](#scenario-1--look-up-one-gene).

```r
ev <- query_gene("LTBP1")
ev$rna$rna_pattern
ev$serum$translation_class
```

### `query_gene_detailed(gene_symbol)`

Per-stage, per-cohort, per-cell-type, and per-filter-step breakdown
of one gene. Use this when `query_gene()` has shown you the
high-level call and you want to know exactly which cohort / stage /
filter step is driving it.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character(1) | required | HGNC gene symbol |

**Returns:** a list with `$per_stage`, `$per_cohort`,
`$per_celltype`, `$filter_diag` (7 rows, one per filter step), and
`$serum_per_cohort` `data.table`s.

```r
detail <- query_gene_detailed("LTBP1")
detail$per_stage          # 4 rows: Normal/Early/Mid/Late
detail$filter_diag        # which of the 7 filter steps passed/failed
```

### `query_panel(gene_symbols)`

Multi-gene join on the atlas's main columns. Faster than calling
`query_gene()` in a loop.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbols` | character | required | One or more HGNC symbols |

**Returns:** `data.table` with 30 columns covering RNA pattern,
protein pattern, translation class, cell origin, serum direction,
filter route, and pooled annotation. One row per input gene.

```r
qp <- query_panel(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
qp[, .(gene_symbol, rna_pattern, prot_pattern, translation_class)]
```

### `list_candidates(onset, tissue_direction, translation_class, min_audit_class, ...)`

Filter the atlas by criterion. Most commonly used to "give me all
Early-onset, Up genes that are at least supported_uncertain".

| Argument | Type | Default | Description |
|---|---|---|---|
| `onset` | character | `NULL` | `"Early"` / `"Mid"` / `NULL` (any) |
| `tissue_direction` | character | `NULL` | `"Up"` / `"Down"` / `NULL` |
| `translation_class` | character | `NULL` | `"A"` / `"B"` / `"C"` / `"inverse"` / `NULL` |
| `min_audit_class` | character | `"ALL"` | `"high_confidence"` / `"supported_uncertain"` / `"penalized"` / `"ALL"` |
| `serum_detected` | logical | `NULL` | filter on serum detectability |
| `min_score` | numeric | `0` | minimum `audit_score` |

**Returns:** `data.table` of matching genes, sorted by audit score.

```r
list_candidates(onset = "Early", tissue_direction = "Up",
                  min_audit_class = "supported_uncertain")
list_candidates(translation_class = "inverse")
```

## Audit framework

### `compute_audit_score(gene_symbols, evidence = NULL, ...)`

Deterministic 3-axis × 2-gate audit score and class label.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbols` | character | `NULL` | atlas genes to score; if `NULL`, scores all |
| `evidence` | data.table | `NULL` | optional user-supplied evidence (from `assemble_user_evidence()`); when supplied, scores user genes instead of atlas |

**Returns:** `data.table` with `gene_symbol`, `evidence_strength`,
`biological_coherence`, `translational_relevance`, `leakage_gate`,
`heterogeneity_gate`, `positive_score`, `audit_score`, `audit_class`.

```r
compute_audit_score(c("LGALS3BP", "LTBP1"))
# Or score user evidence (see project_user_cohort)
compute_audit_score(evidence = my_evidence)
```

### `explain_score(gene_symbol, verbose = TRUE)`

Plain-English decomposition of one gene's audit score: which axis
contributed how much, which gate (if any) penalised it, and the
class-label rationale.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character(1) | required | HGNC symbol |
| `verbose` | logical | `TRUE` | if `TRUE`, prints to console |

**Returns (invisibly):** list with `$audit_class`, `$audit_score`,
`$positive_score`, `$axes` (3-row data.table: axis, weight, value,
contribution), `$gates` (2-row data.table: gate, value,
triggered_by), `$explanation` (the printed paragraph).

```r
explain_score("LTBP1")     # prints to console
res <- explain_score("LTBP1", verbose = FALSE)  # capture only
res$gates
```

### `compare_candidates(gene_symbols)`

Side-by-side comparison of several genes with a redundancy hint.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbols` | character | required | 2 or more HGNC symbols |

**Returns:** `data.table` sorted by `audit_score` descending with
columns `gene_symbol`, `audit_class`, `audit_score`, `rna_pattern`,
`prot_pattern`, `translation_class`, `cell_origin_top`,
`serum_detected`, `serum_log2fc_PDAC_vs_HC`, `max_I2_meta`,
`redundancy_with`. Genes outside the atlas appear with `audit_class
= NA`.

```r
compare_candidates(c("LGALS3BP", "LTBP1", "SERPINA1", "ALB", "GAPDH"))
```

### `propagate_uncertainty(gene_symbols, n_mc = 200, seed = 1)`

Monte Carlo rank stability: how often does this gene stay in the
top-N under bootstrapped evidence perturbation?

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbols` | character | required | one or more HGNC symbols |
| `n_mc` | integer | `200` | number of Monte Carlo iterations |
| `seed` | integer | `1` | RNG seed (scoped via `withr::local_seed()`) |

**Returns:** `data.table` with median rank, 95% CI of rank, and
median audit_score per gene.

### `evaluate_anchor_enrichment(top_n, tier, score_col = "audit_score")`

Hypergeometric enrichment of the bundled external anchor set in the
top-N atlas ranking. **Post-freeze evaluation only** — the anchors
are not used during scoring.

| Argument | Type | Default | Description |
|---|---|---|---|
| `top_n` | integer | `c(50, 100, 200, 500, 1000)` | top-N cutoffs |
| `tier` | character | `c("primary", "secondary", "exploratory", "all")` | anchor tier |
| `score_col` | character | `"audit_score"` | ranking column |

Short alias: `anchor_enrichment()`.

```r
anchor_enrichment(top_n = 100, tier = "secondary")
#> 7 hits in top 100, 39.3× enrichment, hypergeometric p = 2.18e-10
```

### `format_provenance(provenance, style)`

Human-readable evidence labels for the bundled phase tags. Used
internally by `query_gene()` print methods; exported so users who
read the raw `provenance` column can relabel it the same way.

| Argument | Type | Default | Description |
|---|---|---|---|
| `provenance` | character | required | comma-separated tag string or vector |
| `style` | character | `"compact"` | `"compact"` / `"verbose"` / `"raw"` |

**Returns:** character(1).

| `style` | Output example |
|---|---|
| `"compact"` | `"RNA trajectory + Tissue protein + ..."` |
| `"verbose"` | `"- RNA trajectory: matched in bulk RNA-seq stage model\n- ..."` |
| `"raw"` | `"phase33, phase34, ..."` |

```r
format_provenance("phase33,phase34,phase60", "verbose")
```

## Trajectory framework

### `fit_stage_de(object, ...)`

S4 generic — DESeq2 LRT wrapper. Two interfaces:

**Matrix interface** (original):

| Argument | Type | Default | Description |
|---|---|---|---|
| `object` | matrix / data.frame | required | counts (genes × samples) |
| `stage` | factor / character | required | length = `ncol(object)`, levels include `Normal/Early/Mid/Late` |
| `cohort` | factor / character | `NULL` | optional cohort labels (used as fixed effect) |
| `min_count` | numeric | `10` | row-sum filter |
| `padj_cutoff` | numeric | `0.05` | LRT padj cutoff for `lrt_significant` |

**SummarizedExperiment interface**:

| Argument | Type | Default | Description |
|---|---|---|---|
| `object` | SummarizedExperiment | required | `assay()` is the count matrix |
| `stage_col` | character(1) | required | column in `colData(object)` with stage labels |
| `cohort_col` | character(1) | `NULL` | optional cohort column |
| `assay_name` | character(1) | `"counts"` | which assay slot to pull |

**Returns:** `data.table` with one row per gene and columns
`gene_symbol`, `beta_N`, `beta_E`, `beta_M`, `beta_L`, `lfcSE_*`,
`lrt_padj`, `lrt_significant`.

```r
fit <- fit_stage_de(my_counts, my_stage, my_cohort)
# OR
fit <- fit_stage_de(my_se, stage_col = "stage", cohort_col = "cohort")
```

### `fit_stage_de_protein(object, ...)`

limma parallel of `fit_stage_de()` for log2 protein intensity.
Same S4 dispatch on matrix / `SummarizedExperiment`. The SE method
defaults to `assay_name = "intensity"` instead of `"counts"`.

### `classify_trajectory(fit, rho_cutoff = 0.85, sig_only = TRUE)`

Match each LRT-significant gene's z-scored 4-point profile against
the 12-template catalog by Pearson rho.

| Argument | Type | Default | Description |
|---|---|---|---|
| `fit` | data.table | required | output of `fit_stage_de()` |
| `rho_cutoff` | numeric | `0.85` | minimum rho for a template assignment |
| `sig_only` | logical | `TRUE` | restrict to `lrt_significant == TRUE` |

**Returns:** input `fit` augmented with `rna_pattern` (Early × 4
only — non-Early best-matches set to `NA`),
`rna_pattern_rho`, `rna_pattern_rho_runner_up`.

Protein-side wrapper: `classify_protein_trajectory()` (alias
`classify_prot_trajectory()`) — same arguments, but renames output
columns to `prot_pattern*`.

### `score_trajectory(fit, gene)`

Per-gene 12-template Pearson rho vector — useful for inspecting
ambiguous calls.

| Argument | Type | Default | Description |
|---|---|---|---|
| `fit` | data.table | required | `classify_trajectory()` output (or fit) |
| `gene` | character | required | one or more HGNC symbols |

**Returns:** `data.table` with `gene_symbol` + 12 rho columns:
`rho_Early_Burst_Up`, `rho_Early_Loss_Down`, ..., `rho_Monotonic_Down`.

### `align_patient_profile(rna_logfc, prot_logfc, ...)`

Sample-level alignment of one patient's tumor-vs-matched-normal
log2FC profile against the atlas's frozen stage axes (Early / Mid
/ Late). **An alignment readout, not a stage prediction.**

| Argument | Type | Default | Description |
|---|---|---|---|
| `rna_logfc` | named numeric | required | gene_symbol → log2(tumor/normal) |
| `prot_logfc` | named numeric | `NULL` | optional protein layer |
| `weight_by` | character | `"audit_score"` | per-gene weight in rho + vote share |
| `min_audit_score` | numeric | `0.3` | minimum atlas audit_score for inclusion |
| `top_n_genes` | integer | `500` | gene-dictionary cap |
| `min_genes` | integer | `50` | warn-only floor; hard-error below 10 |

**Returns:** `pdactrace_patient_alignment` list with `$rna` (4-row
data.table: stage / cor_to_stage_axis / cor_pval / cor_lo95 /
cor_hi95 / vote_share / weighted_dist / n_genes_used),
`$prot` (categorical concordance, optional), `$summary` (one-line
text), `$attrs` (parameters + provenance).

Short alias: `align_patient()`.

### `project_user_cohort(rna, coldata, stage_col, cohort_col, protein, ...)`

End-to-end wrapper: `fit_stage_de()` → `classify_trajectory()` →
`assemble_user_evidence()` → `compute_audit_score()`. Optional
protein layer.

| Argument | Type | Default | Description |
|---|---|---|---|
| `rna` | matrix / data.frame / SE | required | count matrix or `SummarizedExperiment` |
| `coldata` | data.frame / DataFrame | `NULL` | required when `rna` is a matrix |
| `stage_col` | character(1) | required | name of stage column |
| `cohort_col` | character(1) | `NULL` | optional cohort column |
| `protein` | matrix / SE | `NULL` | optional log2 intensity matrix or SE |
| `protein_assay_name` | character(1) | `"intensity"` | for SE protein input |
| `signal_peptide` | character | `NULL` | gene symbols carrying SP-positive |
| `sig_only` | logical | `FALSE` | for the protein-side classification |

**Returns:** `pdactrace_user_projection` list with `$rna_fit`,
`$rna_pattern`, `$prot_fit`, `$prot_pattern`, `$evidence`,
`$audit`, `$summary`.

## Visualization

### `plot_gene_evidence(gene_symbol, layers = c("rna","protein","scrna","serum","summary"))`

Multi-panel composite figure for one gene. Returns a single ggplot
(actually a `patchwork::wrap_plots()` composite) you can save via
`pdactrace_save()` or `ggplot2::ggsave()`.

### `plot_gene_hexagon(gene_symbol, comparison = NULL)`

6-axis evidence radar for one or more genes. The six axes are
`Multi-layer`, `Direction`, `Stage-onset`, `Serum bridge`,
`Leakage safety`, `Cohort consistency`.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character | required | one or more HGNC symbols |
| `comparison` | character(1) | `NULL` | reference polygon: `"high_confidence_mean"`, `"supported_uncertain_mean"`, `"penalized_mean"`, `"excluded_mean"`, or a single gene symbol |

```r
plot_gene_hexagon("LTBP1", comparison = "high_confidence_mean")
plot_gene_hexagon(c("LGALS3BP", "LTBP1"))   # multi-gene overlay
```

### `plot_template_atlas(layer, output_dir, templates, ...)`

12 PDFs per layer — one per template. Each panel shows the cohort
of genes that match the template best.

| Argument | Type | Default | Description |
|---|---|---|---|
| `layer` | character(1) | `"rna"` | `"rna"` or `"protein"` |
| `templates` | character | `NULL` (= all 12) | subset of templates |
| `output_dir` | character(1) | `NULL` | if non-NULL, writes one PDF per template |
| `width`, `height` | numeric | `1.55`, `1.40` | inches (fig2C-compact) |

**Returns:** named `list` of 12 `ggplot` objects.

### `plot_gene_template(gene_symbol, layer, output_file, ...)`

One PDF showing the gene's matched template panel with the gene's
own trajectory overlaid in highlight colour.

### Other plot helpers

| Function | What it shows |
|---|---|
| `plot_stage_effect(gene)` | Per-stage forest with log2FC ± 1.96·SE |
| `plot_per_cohort(gene)` | Cohort-by-cohort trend bar plot |
| `plot_meta_forest(gene, contrast)` | Random-effects meta-analysis forest |
| `plot_filter_trace(genes)` | Pass/fail step bar across the 7-step filter |
| `plot_filter_diagnostics(...)` | Per-step diagnostic counters across the atlas |
| `plot_panel_heatmap(genes)` | Gene × evidence-axis comparison heatmap |
| `plot_candidate_landscape(...)` | Tissue × serum scatter, Class A/B coloured |
| `plot_celltype_full(...)` | Full cell-type-of-origin overview from scRNA |

All plot functions return a `ggplot` object. Save via:

```r
pdactrace_save(p, dir = "fig", name = "ltbp1_hex",
                w = NCS_W_SINGLE, h = 2.5)   # cairo_pdf, NCS-grade
# OR
ggplot2::ggsave("ltbp1_hex.png", p, width = 3.5, height = 2.5,
                  dpi = 300)
```

## Reporting

### `report_gene(gene_symbol, output_dir = tempdir(), ...)`

Self-contained HTML evidence report (single gene or multi-gene
panel). Renders [`inst/rmd/gene_report.Rmd`](inst/rmd/gene_report.Rmd)
or [`inst/rmd/panel_report.Rmd`](inst/rmd/panel_report.Rmd) via
`rmarkdown::render()`.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character | required | length 1 = single-gene; length 2+ = panel |
| `output_dir` | character(1) | `tempdir()` | output directory |
| `output_file` | character(1) | auto | `<gene>_pdactrace_report.html` |
| `quiet` | logical | `TRUE` | suppress rmarkdown progress |

**Returns (invisibly):** absolute path of the rendered HTML.

```r
report_gene("LTBP1")                         # single-gene HTML
report_gene(c("LGALS3BP", "LTBP1", "GAPDH")) # 3-gene panel HTML
```

## Bioconductor-native helpers

### `as_summarized_experiment(reference = NULL)`

Convert the bundled atlas to a `SummarizedExperiment`:

- Two assays: `rna_beta` (10,113 × 4 stages), `rna_lfcSE`.
- 4-row `colData`: `stage` ∈ {Normal, Early, Mid, Late},
  `reference_level` (TRUE for Normal).
- ~109-column `rowData`: every non-stage-axis atlas column.
- `metadata`: atlas version + cohort count + Zenodo DOIs.

```r
se <- as_summarized_experiment()
SummarizedExperiment::assay(se, "rna_beta")[1:5, ]
SummarizedExperiment::rowData(se)$audit_class[1:5]
```

### `atlas_provenance()`

One-call provenance dossier — package version, atlas snapshot, both
Zenodo DOIs, both repo URLs, cohort count, per-layer dataset count.

### `list_data_sources(layer = NULL)`

26 contributing public datasets with accession + URL + role.

| Argument | Type | Default | Description |
|---|---|---|---|
| `layer` | character | `NULL` | filter: `"RNA"`, `"Protein"`, `"scRNA"`, `"Serum"`, `"Pancreatitis"`, `"Validation"` |

```r
list_data_sources()              # all 26 rows
list_data_sources(layer = "RNA") # 6 rows
```

### `download_phase_csvs(target = "both", ref = "main", cache = NULL)`

Fetch the two large upstream phase CSVs (`phase33` RNA fit, `phase34`
protein fit) from `raw.githubusercontent.com/jibeomko/PDAC_biomarker`
and cache them via `BiocFileCache::BiocFileCache()`. Required only
for users who want to fully re-run `data-raw/build_reference.R`
without cloning the manuscript-monorepo (see
[Reproducibility](#reproducibility)).

## Summary of all 59 exported objects

The above blocks cover the most-used 20. Here are the rest grouped
by role:

| Role | Object |
|---|---|
| Lookup + summary | `summarize_gene_evidence`, `list_atlas_metadata`, `case_study` |
| Audit | `build_evidence_graph`, `extract_graph_features` |
| Trajectory | `assemble_user_evidence`, `early_pattern_names` (alias `early_patterns`), `mid_pattern_names_excluded` (alias `mid_patterns`) |
| Filter | `trace_filters` |
| Schema | `schema_spec` |
| Theme | `pdactrace_axes_theme`, `pdactrace_panel_theme`, `pdactrace_save` |
| Palettes | `pdactrace_pal_class`, `pdactrace_pal_group`, `pdactrace_pal_dir`, `pdactrace_pal_pattern` |
| Width constants | `NCS_W_SINGLE` (3.46 in), `NCS_W_15COL` (4.72 in), `NCS_W_DOUBLE` (7.08 in), `NCS_W_TRIPLE` (10.50 in) |

For every function, `?function_name` opens the full Rd page with
arguments, return shape, and a runnable example.

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

Use whichever you prefer — existing scripts using the long names
continue to work unchanged.

## Vignettes

```r
vignette("lookup_basics",          package = "pdactrace")
vignette("audit_case_studies",     package = "pdactrace")
vignette("audit_framework",        package = "pdactrace")
vignette("user_cohort_extension",  package = "pdactrace")
vignette("reproducibility",        package = "pdactrace")
```

## Reproducibility

The R package is **self-contained for atlas re-derivation.** It
ships with the reference atlas (`data/*.rda`), the canonical
12-template trajectory catalog, the external anchor set, all
testthat suites, the build scripts under `data-raw/`, and the
**six small downstream phase tables** under
`inst/extdata/phase{2c,29,42,60,77,80}_*.csv.xz`. The build chain
runs end-to-end from these bundled inputs alone, plus the two
upstream CSVs fetched via [`download_phase_csvs()`](#download_phase_csvstarget--both-ref--main-cache--null).

Four reproducibility layers (see
`vignette("reproducibility")` for the full walkthrough):

1. **Layer 1 — offline, bundled** — use `data/*.rda` directly.
   Most users start and stop here.
2. **Layer 2 — user cohort** — project a count / intensity matrix
   through `project_user_cohort()`.
3. **Layer 3 — re-derive the atlas** — `download_phase_csvs()` +
   `data-raw/build_reference.R`. Self-contained from public inputs.
4. **Layer 4 — full FASTQ → counts → fits pipeline** — out of
   scope for the package; lives in the manuscript-monorepo
   ([github.com/jibeomko/PDAC_biomarker](https://github.com/jibeomko/PDAC_biomarker),
   Zenodo
   [10.5281/zenodo.20067849](https://doi.org/10.5281/zenodo.20067849)).

```r
# Layer 1 — direct atlas use
library(pdactrace)
query_gene("LTBP1")

# Layer 3 — re-derive the atlas
download_phase_csvs("both")
# ...then run data-raw/build_reference.R as the vignette describes
```

## Citation

If you use `pdactrace`, please cite the software via its Zenodo DOI:

> Ko, J. (2026). *pdactrace: Queryable Stage-Aware PDAC Tissue-to-Serum
> Biomarker Reference Atlas* (v0.99.4) [Software]. Zenodo.
> [10.5281/zenodo.20076698](https://doi.org/10.5281/zenodo.20076698)

```bibtex
@software{pdactrace2026,
  author       = {Ko, Jibeom},
  title        = {{pdactrace: Queryable Stage-Aware PDAC
                   Tissue-to-Serum Biomarker Reference Atlas}},
  year         = 2026,
  version      = {v0.99.4},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20076698},
  url          = {https://doi.org/10.5281/zenodo.20076698}
}
```

The accompanying *Briefings in Bioinformatics* manuscript reference
will be added once the preprint or journal DOI is available.

## License

MIT. See [LICENSE](LICENSE).

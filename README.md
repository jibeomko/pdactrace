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

`pdactrace` is a **transparent, deterministic multi-omics
framework** for prioritising tissue-to-serum biomarker candidates
by integrating per-stage trajectory matching across bulk RNA-seq,
tissue proteomics, single-cell origin, and serum proteomics, with
explicit anchor-enrichment evaluation against a curated external
reference set.

The framework is **methodology-first**: a 12-template competitive
trajectory catalog, a 3-axis + 2-gate audit score, an Evidence
Math layer, an optional interpretable elastic-net prioritisation
layer, and a tissue-to-serum *translation discipline*
(Class A / B / C) that distinguishes direction-preserved,
direction-inverted, and decoupled candidates. The bundled PDAC
reference atlas is the **demonstration cohort** -- the same
framework applies unchanged to other cancers given a per-stage
expression model.

It does **not** train a supervised biomarker classifier and ships
no pretrained predictor. Instead, it uses a frozen, interpretable
scoring rule and reports uncertainty, because validated
non-circular early-detection ground truth is unavailable in PDAC
and several adjacent cancers. As a post-freeze evaluation-only
sanity check, the frozen audit rule preferentially ranks
**7 of the curated secondary-tier external anchor biomarkers in
the top 100** candidates (39.3x hypergeometric enrichment,
p = 2.18e-10; LOO median 41.6x; bootstrap 95% CI [20.2, 56.3]) --
this is a sanity check, not a fully-blinded external validation.

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

**Methodology (transferable to any cancer):**

- A **12-template competitive trajectory catalog** that classifies
  per-gene Normal->Early->Mid->Late profiles into one of 12
  pre-declared shapes (Early x 4 + Mid x 4 + Late x 2 + Monotonic x 2)
  by Pearson correlation against z-scored templates.
- A transparent **3-axis + 2-gate audit score** with closed-form
  decomposition (`explain_score()`), Monte Carlo uncertainty
  (`propagate_uncertainty()`), and a post-freeze evaluation-only
  sanity check showing **39.3x anchor enrichment**
  (hypergeometric p = 2.18e-10; LOO median 41.6x; bootstrap 95%
  CI [20.2, 56.3]) -- the curated anchor set was frozen before
  audit parameters were finalised, but anchor selection and
  framework design share domain expert input, so this is a
  sanity check rather than fully-blinded external validation.
- An **Evidence Math layer** (`evidence_math()`,
  `compare_genes()`) that exposes the per-axis math values
  (delta_rho, ‖beta‖2, RNA-protein cosine, Stouffer Z, tau
  specificity, 7-step filter pass count) without folding them
  into a black-box composite.
- An **optional interpretable ML layer** (`make_evidence_features()`,
  `score_anchor_similarity()`, `fit_user_evidence_model()` over
  elastic net) that runs only on user-supplied labels -- no
  pretrained classifier is shipped.
- A **tissue-to-serum translation discipline** that classifies
  candidates as Class A (same-direction tissue<->serum), Class B
  (direction-inverted) or Class C (decoupled), surfaced as
  `translation_class` in the bundled atlas.

**PDAC demonstration cohort** (bundled with the package):

- A pre-built reference atlas: **10,113 genes x 113 columns**
  across 11 RNA-seq cohorts, pooled tissue proteomics, scRNA cell
  origin, and 3 serum proteomics cohorts.
- Per-gene lookup, panel lookup, candidate listing, filter tracing,
  and a single-call visual canvas (`viz_gene()`).
- A user-cohort wrapper (`project_user_cohort()`) applies the same
  framework to a user-supplied count or intensity matrix.

Core finding (PDAC demonstration): **a tissue biomarker is not
always a serum-up biomarker.** Tissue signals can preserve,
invert, or decouple when projected into serum. Four canonical
case studies illustrate the four audit classes -- LGALS3BP
(`high_confidence`), LTBP1 (`supported_uncertain`, Class B
inverse), ALB (`penalized`, plasma-high-abundance gate), and
GAPDH (`excluded`, housekeeping leakage gate) -- with no single
gene treated as a flagship; each is one demonstration of how the
framework discriminates clean signal from leaky signal.

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

**Visual entry point** (recommended for clinicians and biologists --
one call, the whole evidence picture, no per-axis function names
required):

```r
library(pdactrace)

viz_gene("LTBP1")
```

This produces a 2x2 patchwork canvas with: per-stage trajectory
forest (top-left), per-cohort sign-vote bar (top-right), scRNA
cell-of-origin distribution (bottom-left), and 7-step
tissue-to-serum filter trace (bottom-right), under a one-line
title strip naming the matched template, audit class, and
translation class. Pass `ncol = 1` for a vertical strip suitable
for narrow embedding.

**Text view** (for scripting, knitr reports, or when you need
specific column values):

```r
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

**Visualize.** All four evidence layers in a single composite figure:

```r
plot_gene_evidence("LTBP1")
# Default panels: trajectory + cell_origin + serum + summary.
# Subset with layers = c("trajectory", "serum") for a smaller figure.
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

### Going deeper: per-axis Evidence Math (v0.99.5)

`explain_score()` rolls everything up into one weighted-sum number
plus two gates. If you want the actual *math* behind each evidence
axis, use `evidence_math()` and `explain_gene()`. These show the
underlying values per axis (rho_best, delta_rho, ‖β‖₂, Stouffer Z,
RNA-protein cosine, tau, A/B/C class, 7/7 filters) without folding
them into a composite score:

```r
explain_gene("LTBP1", view = "math")
```

Truncated output:

```
LTBP1 — Evidence Math Summary  [supported_uncertain]

Trajectory fit
  pattern        : Early_Burst_Up
  rho_best       : 0.999
  rho_runner_up  : 0.810
  delta_rho      : 0.188    (specificity margin; >0.10 = clean)

Effect magnitude
  ||beta_RNA||_2 : 0.848    (target stage: E)
  ||beta_prot||_2: 1.828

Cohort consistency
  Stouffer Z     : 0.658    (padj = 1.0000)
  agreement      : 0.50  fraction of cohorts agreeing
  max meta I2    : 75%

RNA-protein coupling
  cosine(beta_RNA, beta_prot) : 0.997
  prot_pattern   : Early_Burst_Up
  prot_tier      : Tier1_gold

Serum bridge
  translation_class      : B
  ...
Cell specificity
  cell_origin_top: myCAF
  tau index      : 0.883

Filter survival
  passed 1 / 7 steps
```

Use `view = "evidence"` for the plain-English provenance only,
`view = "math"` for math only, `view = "both"` for both. Pass
`verbose = FALSE` to suppress printing and just capture the
underlying list.

The pure-data accessor is `evidence_math()` — same content, no
formatting:

```r
m <- evidence_math("LTBP1")
m$trajectory_fit$delta_rho        # 0.188
m$rna_protein_coupling$cosine     # 0.997
m$filter_survival$passed          # 1
```

To compare a panel of genes side-by-side on the math layer, use
`compare_genes()`:

```r
compare_genes(c("LGALS3BP", "LTBP1", "TIMP1"),
              axes = c("trajectory_fit", "rna_protein_coupling"),
              wide = TRUE)
#>     gene  trajectory_fit.delta_rho  rna_protein_coupling.cosine  ...
#> 1: LGALS3BP                  0.093                       0.999
#> 2:    LTBP1                  0.188                       0.997
#> 3:    TIMP1                  0.131                       0.961
```

Default `wide = FALSE` returns a long table with one row per
(gene × axis × metric) — convenient for `data.table::dcast()` into
custom layouts. Note that `compare_candidates()` is the audit-score
ranking layer; `compare_genes()` is the Evidence Math layer. Both
read the same atlas, neither replaces the other.

**Visualize.** The 7-step tissue-to-serum filter trace (matches
`filter_survival` from `evidence_math()`):

```r
plot_filter_trace("LTBP1")
plot_filter_trace(c("LGALS3BP", "LTBP1", "ALB"))   # multi-gene comparison
```

Multi-cohort RNA forest plot (the data behind the
`heterogeneity_gate`):

```r
plot_per_cohort("LTBP1")
```

### Optional: interpretable ML prioritization (v0.99.6)

The deterministic `audit_score` and the per-axis `evidence_math()`
layer cover most use cases. For callers who want a continuous
ML-flavoured prioritisation signal — feature-space proximity to
known biomarkers, or a supervised model trained on their *own*
positive set — v0.99.6 ships an opt-in interpretable ML layer. It
follows three strict design choices:

1. **No deep learning** in the core package — DL belongs in the
   manuscript-monorepo.
2. **No pretrained classifier** is shipped — every supervised fit
   is owned by the user and trained on user-supplied labels.
3. **Per-axis interpretability is preserved** — the user-facing
   output names which evidence features pushed a gene up or down,
   not just a black-box probability.

Step 1 — build a flat feature matrix:

```r
feats <- make_evidence_features(scale = "z", impute = "mean")
dim(feats)
#> [1] 10113    22
```

`make_evidence_features()` returns one row per gene and ~21 numeric
feature columns (trajectory_delta_rho, RNA / protein ‖β‖₂,
RNA-protein cosine, Stouffer Z, cohort agreement, cell-specificity
tau, filter pass fraction, etc.). `impute = "mean"` fills the NA
gaps from the protein and serum layers — necessary because most
genes are missing at least one layer.

Step 2 — descriptive similarity to bundled anchors (no training):

```r
sim <- score_anchor_similarity(tier = "primary")
head(sim, 5)
#>    gene_symbol anchor_similarity anchor_n anchor_tier
#> 1:      IGFBP3             0.957        7     primary
#> 2:      IGFBP2             0.957        7     primary
#> 3:         FN1             0.955        7     primary
#> 4:        ECM1             0.953        7     primary
#> 5:       ANPEP             0.951        7     primary
```

The bundled `pdactrace_external_anchors` set defines a centroid in
the z-scored feature space; `anchor_similarity` is each gene's
cosine to that centroid. **No supervised classifier is trained** —
this is descriptive feature-space distance, defensible to
reviewers as the same evaluation-only discipline already documented
for `evaluate_anchor_enrichment()`.

The new `score_col` value plugs straight into the existing
evaluation harness:

```r
ev <- merge(sim, feats[, .(gene_symbol)], by = "gene_symbol")
evaluate_anchor_enrichment(score_col = "anchor_similarity",
                            evidence = ev, tier = "primary")
```

Step 3 — supervised user-fitted model (only with user labels):

```r
# User supplies their own positive set (their lab's validated panel,
# their literature shortlist, etc.). The package ships nothing.
my_positives <- c("LTBP1", "TIMP1", "LGALS3BP", "LRG1", "CEACAM5")
y <- as.integer(feats$gene_symbol %in% my_positives)

fit <- fit_user_evidence_model(feats, y, alpha = 0.5)
fit
#> <pdactrace_user_model>
#>   method:      elastic_net
#>   n_train:     10113  (positives = 5)
#>   n_features:  21
#>   CV AUC:      0.96 +/- 0.01 at lambda.min = 0.0013
#>   alpha:       0.5  seed: 1

explain_user_evidence_model(fit, top_n = 5)
#> Main positive contributors (push score UP):
#>   +0.276  filter_pass_fraction
#>   +0.255  prot_beta_max_abs
#>   +0.134  serum_n_cohorts_detected
#>   +0.131  serum_abs_log2fc
#>   +0.087  prot_beta_norm
#>
#> Main negative contributors (push score DOWN):
#>   ...

predict_user_evidence_model(fit, feats)   # rank the rest of the atlas
```

Step 4 — model card for reviewers:

```r
model_card("anchor_similarity")           # describes the descriptive layer
model_card("user_model", model = fit)     # describes the user's fit
```

The model card is the single point of reviewer-facing
documentation: feature set version, anchor count, what labels were
or were not used, leakage controls, intended use, limitations.

---

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

> **System requirement.** `report_gene()` renders an HTML report
> via `rmarkdown::render()`, which needs the system tool **pandoc**
> (>= 1.12.3). RStudio bundles pandoc, so simply running R inside
> RStudio is enough. On a bare R install: `sudo apt install pandoc`
> (Debian / Ubuntu), `brew install pandoc` (macOS), or see
> <https://pandoc.org/installing.html>. All other pdactrace
> functions (`query_gene`, `evidence_math`, `plot_*`, etc.) work
> without pandoc.

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

**Visualize the report's component panels standalone.** Useful when
you want to embed individual panels in a slide or paper without the
HTML wrapper:

```r
plot_celltype_full("LTBP1") # full cell-type origin distribution
plot_panel_heatmap(c("LGALS3BP", "LTBP1", "GAPDH"))  # gene x evidence heatmap
plot_candidate_landscape()  # tissue x serum scatter, Class A/B colored
```

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

**Visualize.** `align_patient_profile()` returns a structured list
(`$rna`, `$prot`, `$summary`); render the four-stage rho profile
with base barplot:

```r
barplot(setNames(aln$rna$cor_to_stage_axis, aln$rna$stage),
         ylim = c(0, 1), ylab = "rho vs stage axis",
         main = "Patient stage alignment")
```

## Scenario 6 — Project your own data through the framework

The user-data side of pdactrace is **asymmetric**, mirroring the
data availability in the field: the **tissue layer** (RNA / tissue
protein) requires stage-aware input (Normal / Early / Mid / Late),
while the **serum layer** only needs binary group contrasts
(PDAC vs HC, optionally Pancreatitis vs HC). The two paths share
the same downstream evidence-assembly + audit-score machinery.

> ⚠️ **Important data requirement (tissue side).** The 12-template
> trajectory framework matches per-gene **N → E → M → L curve
> shapes** to one of 12 pre-declared templates. **Binary
> Normal vs Tumor input does NOT apply at the tissue layer** —
> the framework needs at least 2 of the 4 stage levels with
> samples to discriminate trajectory shape. If you only have
> N vs T for your tissue cohort, you have three options:
>
> 1. Add stage info from clinical metadata if available
>    (TNM / AJCC mapping table is in the `user_cohort_extension`
>    vignette).
> 2. Use a different tool for tissue analysis (e.g., direct
>    DESeq2 / limma) and bring the **serum** side through
>    `project_user_serum_cohort()` (Scenario 6b below) which
>    needs only N vs T.
> 3. Re-cast all tumor samples to one stage (degenerate; loses
>    framework value).

### Scenario 6a — Tissue projection (stage-aware)

**Inputs:**

- A count matrix (genes × samples) — integer matrix or `data.frame`,
  or a `SummarizedExperiment` with the counts in its `assay()`.
- A `coldata` `data.frame` with at least `stage` (must contain
  levels from `c("Normal", "Early", "Mid", "Late")`) and
  optionally `cohort`. Minimum 8 samples across stages.
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

**Visualize.** Atlas-wide reference panels for sanity-checking your
cohort's trajectories, plus a single-gene overlay onto the matched
template:

```r
plot_template_atlas("rna",     output_dir = tempdir())
plot_template_atlas("protein", output_dir = tempdir())
# 12 PDFs per layer (Early × 4 + Mid × 4 + Late × 2 + Monotonic × 2).

plot_gene_template("LTBP1", layer = "rna")
plot_gene_template("LTBP1", layer = "protein")
```

### Scenario 6b — Serum projection (binary group contrast)

**Inputs:** much simpler than 6a, no stage labels required.

- A **log2 protein abundance matrix** (genes × samples) — `matrix`
  or `data.frame`. Pre-log2 transformation is required; raw MS
  intensities will produce nonsensical log2FC values.
- A `coldata` `data.frame` with at least one column carrying
  group labels matching `pdac_label`, `hc_label`, and optionally
  `pan_label`. `nrow(coldata)` must equal `ncol(intensity)`,
  same order. Minimum 3 samples per group recommended.

```r
# Example: 50 genes × 24 samples (8 PDAC, 8 HC, 8 Pancreatitis)
set.seed(1)
g <- 50; n <- 24
intensity <- matrix(rnorm(g * n, mean = 5),
                     nrow = g, ncol = n,
                     dimnames = list(paste0("GENE", seq_len(g)),
                                      paste0("S", seq_len(n))))
intensity[1:10, 1:8] <- intensity[1:10, 1:8] + 2   # PDAC-up signal

coldata <- data.frame(
  sample = colnames(intensity),
  group  = rep(c("PDAC", "HC", "Pancreatitis"), each = 8))

serum_summary <- project_user_serum_cohort(
  intensity, coldata,
  group_col  = "group",
  pdac_label = "PDAC",
  hc_label   = "HC",
  pan_label  = "Pancreatitis")    # optional

head(serum_summary)
#>    gene_symbol serum_log2fc_PDAC_vs_HC serum_padj_PDAC_vs_HC
#> 1:       GENE1                    2.05                3.4e-04
#> 2:       GENE2                    1.92                4.1e-04
#> ...
#>    serum_log2fc_Pan_vs_HC translation_class serum_detected
#> 1:                  -0.31                 B           TRUE
```

The function:

- Computes per-gene `serum_log2fc_PDAC_vs_HC` and (if `pan_label`
  is set) `serum_log2fc_Pan_vs_HC`.
- Runs a per-gene two-group test
  (`test = "limma_eBayes"` default; falls back to `"wilcox"` or
  `"t_test"`) and BH-adjusts.
- Joins the bundled `pdactrace_reference`'s tissue direction to
  assign `translation_class` (A = same direction tissue ↔ serum;
  B = opposite, the rare inverse-translation case;
  C = decoupled).
- Flags `serum_detected = TRUE` for genes with
  `padj_PDAC_vs_HC <= padj_cutoff` (default 0.05).

The output table plugs directly into
`assemble_user_evidence(serum_summary = ...)`.

### Scenario 6c — Combined tissue + serum into final audit score

If the user has both tissue and serum cohorts (rare, but ideal):

```r
# Step 1: tissue side (stage-aware)
tis <- project_user_cohort(rna = my_counts, coldata = my_tissue_cd,
                            stage_col = "stage", cohort_col = "cohort",
                            protein   = my_tissue_protein)

# Step 2: serum side (binary)
ser <- project_user_serum_cohort(intensity = my_serum,
                                   coldata = my_serum_cd,
                                   pan_label = "Pancreatitis")

# Step 3: combine into one evidence frame + score
ev <- assemble_user_evidence(rna_fit       = tis$rna_pattern,
                              prot_fit      = tis$prot_pattern,
                              serum_summary = ser)
audit <- compute_audit_score(evidence = ev)
head(audit[order(-audit_score)], 20)
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

As an evaluation-only sanity check, the frozen audit rule
preferentially ranks **7 of the curated secondary-tier external
anchors into the top 100** (`39.3x`, hypergeometric
`p = 2.18e-10`). This is post-freeze (the anchor set was
finalised before audit parameters), but the curators of the
anchor set and the framework overlap, so the result is best
interpreted as a sanity check rather than fully-blinded
external validation. A held-out PDAC cohort run is planned for
v1.0.

## Why these specific cutoffs?

This subsection makes the framework's numeric choices explicit so
they can be audited or perturbed. Three groups of cutoffs:

**Trajectory cutoff (`rho_cutoff = 0.85`).**
A z-scored 4-stage profile (Normal / Early / Mid / Late) is
matched against 12 pre-declared templates by Pearson correlation.
The `rho >= 0.85` threshold is **conservative** -- four-point
correlations crowd toward 1, so 0.85 still requires the gene's
trajectory shape to closely follow one specific template. The
`rho` value used by the surface call is also recorded
(`rna_pattern_rho`); users wanting a different stringency can
re-run [classify_trajectory()] with their own cutoff. The
companion `methodology_validation` vignette includes a sweep
across `rho_cutoff` in `{0.80, 0.85, 0.90}` showing that
top-100 anchor enrichment is robust across this range.

**Three-axis weights (`0.40 / 0.35 / 0.25`).**
The weights map onto a simple priority order:
**evidence_strength** (which layers are present at all)
> **biological_coherence** (do they agree on direction across
RNA cohorts and concord with protein) > **translational_relevance**
(does the gene reach serum, with healthy/pancreatitis discrimination).
This ordering reflects the practical bottleneck in PDAC biomarker
discovery: **layer presence is the limiting factor** -- many genes
have RNA but no serum data, so weighting that axis highest
reflects how much it constrains the prioritisation. The weights
sum to 1.0 by design so `positive_score in [0, 1]`. Sensitivity
to the exact weight values is documented in the
`methodology_validation` vignette (Section B).

**Two-gate multipliers.**
The leakage gate enforces hard discounts for known artifact
classes:

| Trigger | Multiplier | Rationale |
|---|---|---|
| housekeeping flag | `0.00` | "Housekeeping" by definition is invariant; trajectory signal is presumed artifactual. Hard zero. |
| plasma_high_abundance flag | `0.50` | Top-decile plasma proteins (e.g. ALB, IGHG1) elevated in any inflammatory state, not PDAC-specific. Half-discount, not zero, because some are still meaningful (ALB drops in advanced disease). |
| neither | `1.00` | Pass through. |

The heterogeneity gate uses the **Higgins-Thompson I-squared
convention** (low / moderate / high boundaries at 25 / 50 / 75)
adapted to a pdactrace-specific binning that reflects the
v0.3.0 audit experiment:

| max meta I^2 across N-vs-E / M-vs-E / L-vs-E | Multiplier |
|---|---|
| `< 70%` (low / moderate) | `1.00` |
| `70-90%` (high) | `0.70` |
| `>= 90%` (very high) | `0.30` |

Both gates emit deterministic per-gene values; their stability
under evidence perturbation is reported by
[`propagate_uncertainty()`](#propagate_uncertaintygenes-n_mc-seed).

**Audit class boundaries (`0.5` / `0.3`).**
Genes with `audit_score >= 0.5` are `high_confidence`;
`>= 0.3` is `supported_uncertain`; below `0.3` and not gate-zeroed
is `low`. Gate-zeroed genes are `excluded` (leakage = 0) or
`penalized` (leakage = 0.5).

**Trajectory effect-size threshold (`0.585 = log2(1.5)`).**
Used internally in the rescue-eligibility check
(`max_abs_beta_meta < 0.585`). 1.5-fold change is the default
"meaningfully detectable" effect size for log2 fold-change
discussions in the bulk RNA-seq + proteomics literature.

**Cohort independence assumption (Stouffer meta-analysis).**
The cross-cohort consistency Z-statistic (`rna_stouffer_z`)
combines per-cohort one-sided p-values via Stouffer's method,
which **assumes independence between cohorts**. In the bundled
atlas, this assumption is partially violated: TCGA-PAAD and
CPTAC PDAC have overlapping donor recruitment (CPTAC re-quantified
a subset of TCGA samples for proteomics). The framework absorbs
this via the **heterogeneity gate** -- genes with high I^2 across
the four contributing cohorts (TCGA / CPTAC / GSE224564 /
GSE79668) get a multiplier of `0.7` or `0.3`, capping the
inflated significance that the Stouffer assumption would
otherwise produce. Users running [`project_user_cohort()`] on
their own data with non-overlapping donors are not affected by
this caveat.

## Projection stress test and evidence-layer dependence

We performed a three-way evaluation comparing the bundled full
multi-layer atlas, a bundled RNA-only baseline, and an external
RNA-only projection. The full atlas showed the strongest anchor
enrichment, the bundled RNA-only baseline retained partial
enrichment, and the held-out RNA-only projection did not recover
anchor enrichment. In contrast, housekeeping negative-control
discipline was reproduced in the held-out projection after
transfer of atlas-defined leakage annotations.

| Source | Top-50 fold | Top-100 fold | Top-200 fold | Top-500 fold |
|---|---:|---:|---:|---:|
| Bundled FULL, multi-layer | 44.9x | 39.3x | 22.5x | 10.1x |
| Bundled RNA-only | 33.7x | 16.9x | 8.4x | 4.5x |
| Held-out RNA-only | 0 | 0 | 0 | 0 |

Negative-control housekeeping bottom-500 enrichment:

| Source | hits / 31 | fold | p-value |
|---|---:|---:|---:|
| Bundled FULL | 31 | 20.2x | 1.3e-41 |
| Bundled RNA-only (audit_* recomputed) | 1 | 0.65 | 0.79 |
| Held-out RNA-only (atlas leakage annotations applied) | 28 | 17.8x | 1.1e-33 |

These results indicate that the frozen audit score gains
specificity from multi-layer evidence integration rather than
stage-aware RNA trajectory alone. The held-out RNA-only
projection is therefore reported as a **projection stress test**,
not as a definitive external multi-omics validation. Full
methodology and code is in `vignette("methodology_validation")`
Sections D and E.

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

### `evidence_math(gene_symbol, reference = NULL)`

Returns the per-axis mathematical evidence values that fed the audit
decisions, organised by axis. The pure-data accessor underneath
`explain_gene()` and `compare_genes()`.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character(1) | required | HGNC symbol |
| `reference` | data.table | `NULL` | optional atlas override (testing) |

**Returns:** a list with one element per axis:

- `trajectory_fit` — `rna_pattern`, `rho_best`, `rho_runner_up`,
  `delta_rho`, `note`.
- `effect_magnitude` — `rna_beta_norm` (‖β‖₂ over E,M,L),
  `rna_beta_max_abs`, `rna_max_at_stage`, `rna_target_stage`,
  protein-side counterparts.
- `cohort_consistency` — `stouffer_z`, `stouffer_p`,
  `stouffer_padj`, `cohort_agreement`, `max_meta_I2`.
- `rna_protein_coupling` — `cosine` (E,M,L), `prot_pattern`,
  `prot_pattern_rho`, `prot_tier`, `rnaprot_concordant`,
  `prot_in_atlas`, `note`.
- `serum_bridge` — `translation_class` (A/B/C),
  `serum_log2fc_PDAC_vs_HC`, `serum_log2fc_Pan_vs_HC`,
  `serum_n_cohorts_detected`, `phase77_strict`.
- `cell_specificity` — `cell_origin_top`, `tau`, `cell_origin_padj`.
- `filter_survival` — `passed`, `total`, `per_step` (named logical).
- `clinical_role` — `resectable_marker`, `panel_member`.

```r
m <- evidence_math("LTBP1")
m$trajectory_fit$delta_rho
m$rna_protein_coupling$cosine
m$filter_survival$per_step
```

### `explain_gene(gene_symbol, view, verbose, reference)`

Console-friendly text formatter wrapping `format_provenance()` (the
plain-English provenance summary) and `evidence_math()` (the
per-axis math). Mirrors `explain_score()` in pattern: prints when
`verbose = TRUE`, returns the structured data invisibly either way.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbol` | character(1) | required | HGNC symbol |
| `view` | character(1) | `"evidence"` | `"evidence"`, `"math"`, or `"both"` |
| `verbose` | logical | `TRUE` | if `TRUE`, prints sections to console |
| `reference` | data.table | `NULL` | optional atlas override (testing) |

**Returns (invisibly):** list with `$gene`, `$view`, `$audit_class`,
`$provenance` (compact one-liner), `$math` (the `evidence_math()`
list, or `NULL` for `view = "evidence"`).

```r
explain_gene("LTBP1", view = "math")    # math only
explain_gene("LTBP1", view = "both")    # provenance + math
res <- explain_gene("LTBP1", view = "math", verbose = FALSE)
res$math$trajectory_fit$delta_rho
```

### `compare_genes(gene_symbols, axes, wide, reference)`

Multi-gene tidy-table pivot of `evidence_math()`. Long form by
default (`gene, axis, metric, value`); `wide = TRUE` returns one row
per gene with `axis.metric` columns suitable for manuscript tables.

| Argument | Type | Default | Description |
|---|---|---|---|
| `gene_symbols` | character | required | HGNC symbols |
| `axes` | character | `NULL` (all) | axis filter — any of the 8 axis names |
| `wide` | logical | `FALSE` | one row per gene with `axis.metric` columns |
| `reference` | data.table | `NULL` | optional atlas override (testing) |

**Returns:** a `data.table`. Long form columns: `gene, axis,
metric, value` (numeric values formatted to 4 sig figs as character
to keep one column type). Wide form: `gene` plus one column per
`axis.metric` pair.

```r
compare_genes(c("LGALS3BP", "LTBP1", "TIMP1"))
compare_genes(c("LGALS3BP", "LTBP1"),
              axes = c("trajectory_fit", "rna_protein_coupling"),
              wide = TRUE)
```

Note: `compare_candidates()` is the audit-score ranking layer (with
redundancy grouping); `compare_genes()` is the Evidence Math layer.
Both read the same atlas, neither replaces the other.

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

**Advanced / manual rebuild only — not part of the default user
path.** Helper for users who want to rebuild
`data/pdactrace_reference.rda` from the documented public processed
inputs (Layer 3 in the [Reproducibility](#reproducibility) section).
Caches the two large upstream phase tables via
`BiocFileCache::BiocFileCache()` so subsequent rebuilds reuse the
cache.

Requires network access at first call and is intentionally **not
evaluated** by the package's vignettes, examples, or unit tests.
Only use when manually re-running `data-raw/build_reference.R`.

```r
## Not evaluated during package checks.
## Requires network access; intended only for manual rebuilds.
# download_phase_csvs("both")
```

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

The package provides a **pre-built PDAC trajectory reference atlas
for offline use**, together with the 12-template catalog, the
external anchor set, the unit-test suite, and documented build
scripts (`data-raw/`) used to assemble the distributed reference
object.

For ordinary use, **no network access is required**:

```r
library(pdactrace)
query_gene("LTBP1")
```

The package supports several reproducibility levels (see
`vignette("reproducibility")` for the full walkthrough):

- **Layer 1 — bundled reference atlas.** Use the distributed
  `data/*.rda` objects directly. This is the default user path and
  the only path required for normal evidence lookup, scoring,
  reporting, and visualization.
- **Layer 2 — user cohort projection.** Project a user-supplied
  count or intensity matrix onto the canonical trajectory templates
  with `project_user_cohort()`. Runs entirely on local inputs.

The two layers above are the package's primary scope. The two
layers below are documented as **advanced reproducibility / data
provenance** for users who want to rebuild the distributed object
from public inputs:

- **Layer 3 — processed-input atlas rebuild.** Rebuild
  `data/pdactrace_reference.rda` from the documented processed
  phase tables and `data-raw/build_reference.R`. As of v0.99.6,
  the seven processed inputs the build script reads
  (`multi_cohort_consistency.csv` plus six phase tables) are all
  bundled in `inst/extdata/*.csv.xz`, so the rebuild runs without
  the manuscript-monorepo. Input file provenance and build
  determinism are described in the reproducibility vignette.
- **Layer 4 — raw-data reanalysis.** FASTQ / raw proteomics
  processing, count generation, model fitting, and large
  intermediate files are **outside the scope of this software
  package** and are documented in the associated manuscript
  workflow archive.

The distributed atlas is **self-contained for offline use**.
Re-derivation from processed public inputs (Layer 3) is documented
separately in the vignette.

### Where the large quantification matrices live

The raw quantification matrices (per-cohort RNA count tables,
FragPipe protein intensity tables, the 372k-cell scVI atlas
embedding) are intentionally **not** bundled in the package
tarball. The appropriate Bioconductor mechanism for hosting these
artefacts is a companion **`pdactraceData` ExperimentHub package**
(`biocViews: ExperimentData, ExperimentHub`); that package is
planned for a future release and will allow lazy on-demand access
via `ExperimentHub::ExperimentHub()`.

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

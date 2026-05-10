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

`pdactrace` is a **transparent, deterministic multi-omics framework**
for prioritising tissue-to-serum biomarker candidates. It integrates
per-stage trajectory matching across bulk RNA-seq, tissue
proteomics, single-cell origin, and serum proteomics under a
**12-template competitive trajectory catalog**, a **3-axis + 2-gate
audit score** (closed-form, deterministic), an **Evidence Math
layer** that exposes the per-axis math values rather than folding
them into a black box, and a **tissue-to-serum translation
discipline** (Class A/B/C) that distinguishes direction-preserved,
direction-inverted, and decoupled candidates.

It does **not** train a supervised biomarker classifier and ships
no pretrained predictor — validated non-circular early-detection
ground truth is unavailable in PDAC. As a post-freeze
evaluation-only sanity check, the frozen audit rule preferentially
ranks **7 of the curated secondary-tier external anchor biomarkers
in the top 100** candidates (39.3x hypergeometric enrichment,
p = 2.18e-10; LOO median 41.6x; bootstrap 95% CI [20.2, 56.3]).
This is a sanity check, not a fully-blinded external validation.

The bundled PDAC reference atlas is the focus of this release.
The API is designed so the same evidence-aggregation and audit
logic could be reused for other cancers given an analogous
per-stage expression atlas, but the demonstrated value, the
curated external anchor set, and the validation work are all
PDAC-specific.

## Install

```r
# install.packages("remotes")
remotes::install_github("jibeomko/pdactrace")
library(pdactrace)
atlas_provenance()
```

A Bioconductor submission is in preparation; once accepted,
`BiocManager::install("pdactrace")` will work.

## Quick start

```r
library(pdactrace)

# 1. One-call visual canvas — every evidence layer for a gene
viz_gene("LTBP1")

# 2. Plain-English audit-score decomposition
explain_score("LTBP1")

# 3. Per-axis math (delta_rho, ‖β‖₂, RNA-protein cosine, Stouffer Z, ...)
explain_gene("LTBP1", view = "math")

# 4. Side-by-side panel comparison
compare_candidates(c("LGALS3BP", "LTBP1", "SERPINA1", "ALB", "GAPDH"))

# 5. Project your own stage-aware cohort through the framework
res <- project_user_cohort(rna = my_counts, coldata = my_cd,
                            stage_col = "stage", cohort_col = "cohort")
```

# What you get

**Methodology:**

- **12-template competitive trajectory catalog** — Early × 4 +
  Mid × 4 + Late × 2 + Monotonic × 2; only Early × 4 are surfaced
  in `rna_pattern`, the other 8 act as negative-evidence pool.
- **3-axis + 2-gate audit score** — `evidence_strength` +
  `biological_coherence` + `translational_relevance`, multiplied
  by `leakage_gate × heterogeneity_gate`. Closed-form
  decomposition via `explain_score()`; Monte Carlo rank stability
  via `propagate_uncertainty()`.
- **Evidence Math layer** — `evidence_math()` and `compare_genes()`
  expose per-axis values (delta_rho, ‖β‖₂, RNA-protein cosine,
  Stouffer Z, tau specificity, 7-step filter pass count) without
  folding them into a composite.
- **Optional interpretable ML layer** —
  `make_evidence_features()`, `score_anchor_similarity()`,
  `fit_user_evidence_model()` (elastic net). Runs only on
  user-supplied labels; no pretrained classifier shipped.
- **Tissue-to-serum translation discipline** — Class A
  (same-direction tissue↔serum), Class B (inverted), Class C
  (decoupled), surfaced as `translation_class`.

**PDAC demonstration cohort** (bundled with the package):

- 10,113 genes × 113 columns across 11 RNA-seq cohorts, pooled
  tissue proteomics, scRNA cell origin, and 3 serum proteomics
  cohorts.
- Per-gene lookup, panel comparison, filter tracing, single-call
  visual canvas (`viz_gene()`), and a user-cohort wrapper
  (`project_user_cohort()`).
- Four canonical case studies illustrate the four audit classes:
  LGALS3BP (`high_confidence`), LTBP1 (`supported_uncertain`,
  Class B inverse), ALB (`penalized`, plasma-high gate), GAPDH
  (`excluded`, housekeeping gate). No single gene is treated as
  a flagship.

**Core finding:** *a tissue biomarker is not always a serum-up
biomarker.* Tissue signals can preserve, invert, or decouple when
projected into serum.

# Step-by-step walkthrough

Six end-to-end scenarios. Detailed walkthroughs live in the
listed vignettes; the calls below are copy-paste runnable on the
bundled atlas.

## Scenario 1 — Look up one gene

```r
viz_gene("LTBP1")            # 6-panel visual canvas
query_gene("LTBP1")          # text view (full evidence dump)
query_gene_detailed("LTBP1") # per-stage / per-cohort / per-celltype
```

`viz_gene()` produces a 2×3 patchwork canvas with bulk RNA,
tissue protein, scRNA cell-origin, serum direction, 7-step filter
trace, and 6-axis audit hexagon. Pass `layout = "split"` to get
each panel as a full-size figure. See
`vignette("lookup_basics")`.

## Scenario 2 — Understand *why* a gene got its score

The audit score is a frozen, deterministic weighted sum:

```text
positive_score = 0.40 * evidence_strength
               + 0.35 * biological_coherence
               + 0.25 * translational_relevance
audit_score    = positive_score * leakage_gate * heterogeneity_gate
```

```r
explain_score("LTBP1")       # plain-English decomposition
evidence_math("LTBP1")       # per-axis math values
explain_gene("LTBP1", view = "math")     # math-only print
compare_genes(c("LGALS3BP", "LTBP1"), wide = TRUE)
plot_filter_trace("LTBP1")   # 7-step tissue-to-serum filter
```

`explain_score()` reports the audit class with the contributing
axes and gates. `evidence_math()` is the pure-data accessor for
the per-axis math (no formatting). See
`vignette("audit_case_studies")`.

### Optional: interpretable ML prioritization

For callers who want a continuous ML-flavoured prioritisation
signal, v0.99.6 ships an opt-in interpretable layer. Three design
constraints: no deep learning, no pretrained classifier, per-axis
interpretability preserved.

```r
feats <- make_evidence_features(scale = "z", impute = "mean")
sim   <- score_anchor_similarity(tier = "primary")
fit   <- fit_user_evidence_model(feats, my_user_labels, alpha = 0.5)
explain_user_evidence_model(fit, top_n = 5)
model_card("user_model", model = fit)
```

`score_anchor_similarity()` is descriptive cosine to the bundled
anchor centroid — no training. `fit_user_evidence_model()`
trains elastic net on the user's *own* positive set; the package
ships no pretrained predictor. See `vignette("audit_framework")`.

## Scenario 3 — Compare a panel side-by-side

```r
compare_candidates(c("LGALS3BP", "LTBP1", "SERPINA1", "ALB", "GAPDH"))
plot_gene_hexagon(c("LGALS3BP", "LTBP1", "SERPINA1"))
plot_panel_heatmap(c("LGALS3BP", "LTBP1", "GAPDH"))
```

`compare_candidates()` returns a `data.table` sorted by
`audit_score` descending, with a `redundancy_with` column that
flags genes sharing both `rna_pattern` and `cell_origin_top` —
useful for picking a non-redundant panel.

## Scenario 4 — One-gene HTML evidence report

Requires the system tool `pandoc` (RStudio bundles it; on a bare
R install: `sudo apt install pandoc` / `brew install pandoc`).

```r
fp <- report_gene("LTBP1", output_dir = tempdir())
browseURL(fp)

report_gene(c("LGALS3BP", "LTBP1", "GAPDH"))    # multi-gene panel
```

The HTML (~400 KB, self-contained) includes audit components,
6-axis radar, stage trajectory, per-cohort breakdown, filter
trace, and an atlas-version provenance footer.

## Scenario 5 — Single-patient trajectory alignment

```r
aln <- align_patient_profile(patient_rna, top_n_genes = 200)
print(aln)
```

Aligns one patient's tumor-vs-matched-normal log2FC profile
against the atlas's frozen Early/Mid/Late stage axes. **An
alignment readout, not a stage prediction**: misalignment with
all three axes is itself a finding (it does *not* imply Normal).
Deterministic given (input, atlas).

## Scenario 6 — Project your own data through the framework

The user-data side is **asymmetric**, mirroring data availability
in the field. The tissue layer needs stage-aware input
(N/E/M/L); the serum layer only needs binary group contrasts
(PDAC vs HC, optionally Pancreatitis vs HC).

> **Tissue requirement.** The 12-template framework matches per-gene
> N→E→M→L curve shapes. Binary Normal vs Tumor input does **not**
> apply at the tissue layer — the framework needs at least 2 of
> the 4 stage levels to discriminate trajectory shape. If you only
> have N vs T, either add stage from clinical metadata, or use
> only the serum-side wrapper (Scenario 6b).

```r
# 6a. Tissue projection (stage-aware)
data(toy_counts);  data(toy_coldata);  data(toy_protein)
res <- project_user_cohort(rna = toy_counts, coldata = toy_coldata,
                            stage_col = "stage", cohort_col = "cohort",
                            protein   = toy_protein)
res$audit       # data.table from compute_audit_score()

# 6b. Serum projection (binary contrast only)
ser <- project_user_serum_cohort(intensity = my_serum,
                                  coldata = my_serum_cd,
                                  pdac_label = "PDAC",
                                  hc_label   = "HC",
                                  pan_label  = "Pancreatitis")

# 6c. Combine into one evidence frame + score
ev    <- assemble_user_evidence(rna_fit       = tis$rna_pattern,
                                 prot_fit      = tis$prot_pattern,
                                 serum_summary = ser)
audit <- compute_audit_score(evidence = ev)
```

`project_user_cohort()` accepts both matrix + coldata and a
`SummarizedExperiment`. See `vignette("user_cohort_extension")`
for the TNM/AJCC mapping table and worked examples.

# Audit scoring rule

The frozen v0.3.0 audit rule combines three weighted evidence
axes and two reliability gates:

```text
positive_score = 0.40 * evidence_strength
               + 0.35 * biological_coherence
               + 0.25 * translational_relevance
audit_score    = positive_score * leakage_gate * heterogeneity_gate
```

Deterministic given (gene, atlas). Monte Carlo summaries via
`propagate_uncertainty()` report rank stability. External anchors
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
#> 7 hits in top 100, 39.3× enrichment, hypergeometric p = 2.18e-10
```

## Why these specific cutoffs?

- **Trajectory cutoff (`rho_cutoff = 0.85`).** Four-point
  z-scored profiles crowd toward 1, so 0.85 still requires the
  gene's shape to closely follow one specific template.
  `rna_pattern_rho` is preserved for users who want a different
  stringency. The `methodology_validation` vignette sweeps
  `rho_cutoff ∈ {0.80, 0.85, 0.90}` and confirms top-100 anchor
  enrichment is robust.
- **Three-axis weights (`0.40 / 0.35 / 0.25`).** Reflects the
  practical bottleneck: layer presence is the limiting factor in
  PDAC biomarker discovery (many genes have RNA but no serum
  data). Weights sum to 1.0 so `positive_score ∈ [0, 1]`.
  Sensitivity sweep is in the validation vignette (Section B).
- **Leakage gate.** housekeeping flag → `0.00` (housekeeping is
  invariant by definition; trajectory signal is artefactual);
  plasma_high_abundance flag → `0.50` (top-decile plasma proteins
  inflate in any inflammation, not PDAC-specific); neither →
  `1.00`.
- **Heterogeneity gate.** Higgins-Thompson I² adapted to the
  v0.3.0 audit experiment: `< 70%` → `1.00`; `70–90%` → `0.70`;
  `≥ 90%` → `0.30`. Caps the Stouffer significance inflation
  that arises from partially overlapping cohorts (e.g.
  TCGA-PAAD and CPTAC-PDAC share donors).
- **Audit class boundaries.** `audit_score ≥ 0.5` →
  `high_confidence`; `≥ 0.3` → `supported_uncertain`; below 0.3
  and not gate-zeroed → `low`. Gate-zeroed genes are `excluded`
  (leakage = 0) or `penalized` (leakage = 0.5).
- **Effect-size threshold (`0.585 = log2(1.5)`).** 1.5-fold
  change is the default "meaningfully detectable" effect size in
  bulk RNA-seq + proteomics literature; used by the
  rescue-eligibility check (`max_abs_beta_meta < 0.585`).
- **Cohort independence assumption.** The Stouffer Z combiner
  assumes independence between cohorts. The bundled atlas
  partially violates this (TCGA / CPTAC overlap); the
  heterogeneity gate is the absorber. Users running
  `project_user_cohort()` on non-overlapping donors are
  unaffected.

## Projection stress test and evidence-layer dependence

Three-way evaluation: bundled full multi-layer atlas, bundled
RNA-only baseline, and an external RNA-only projection (held-out
PDAC cohort, GSE253260 + GTEx normals).

| Source | Top-50 fold | Top-100 fold | Top-200 fold | Top-500 fold |
|---|---:|---:|---:|---:|
| Bundled FULL, multi-layer | 44.9x | 39.3x | 22.5x | 10.1x |
| Bundled RNA-only | 33.7x | 16.9x | 8.4x | 4.5x |
| Held-out RNA-only | 0 | 0 | 0 | 0 |

Negative-control (housekeeping bottom-500):

| Source | hits / 31 | fold | p-value |
|---|---:|---:|---:|
| Bundled FULL | 31 | 20.2x | 1.3e-41 |
| Bundled RNA-only (audit_* recomputed) | 1 | 0.65 | 0.79 |
| Held-out RNA-only (atlas leakage annotations applied) | 28 | 17.8x | 1.1e-33 |

The frozen audit score gains specificity from **multi-layer
evidence integration** rather than stage-aware RNA trajectory
alone. The held-out RNA-only run is reported as a **projection
stress test**, not as a definitive external multi-omics
validation. Full methodology in
`vignette("methodology_validation")` Sections D and E.

# Stage harmonization

Every cohort is collapsed onto a single 4-level scale
(`Normal / Early / Mid / Late`) before any per-gene trajectory or
audit step. Mapping is deterministic and applied identically
across RNA, tissue protein, and serum.

| AJCC | 4-level group |
|---|---|
| Stage I (IA, IB) | Early |
| Stage II (IIA, IIB) | Mid |
| Stage III, Stage IV | Late |

TNM with M1 → Stage IV → Late; T3/T4 (M0, N0) → Stage II → Mid;
T1–T2 (M0, N0) → Stage I → Early. Cohort-specific labels
("resectable" → Mid; "borderline"/"locAdvanced" → Late;
"metastatic" → Late). Reference: AJCC 8th ed.

III + IV are collapsed into `Late` because per-cohort AJCC IV
counts are too small for stable per-stage contrasts; a 4-level
sensitivity copy is preserved as `stage_ajcc4` /
`stage_numeric`.

# Trajectory framework

| Family | Templates | Atlas surface? |
|---|---|---|
| Early × 4 | `Early_Burst_Up`, `Early_Loss_Down`, `Early_Peak`, `Early_Trough` | surfaced (`rna_pattern`) |
| Mid × 4 | `Mid_Plateau_Up`, `Mid_Plateau_Down`, `Mid_Peak`, `Mid_Trough` | flagged via `excluded_mid_pattern` |
| Late × 2 | `Late_Burst_Up`, `Late_Loss_Down` | flagged via `excluded_late_pattern` |
| Monotonic × 2 | `Monotonic_Up`, `Monotonic_Down` | flagged via `excluded_monotonic_pattern` |

A gene's z-scored 4-point trajectory is matched against all 12
templates by Pearson rho. A candidate must beat Mid, Late, and
Monotonic alternatives before being surfaced as Early-onset.

```r
plot_template_atlas("rna",     output_dir = tempdir())
plot_gene_template("LTBP1",    layer = "rna")
```

# Reference atlas

`pdactrace_reference` is a `data.table`: 10,113 genes × 113
columns (RNA trajectory, tissue protein, scRNA cell-origin, serum
translation, pancreatitis context, audit scores, uncertainty
summaries, provenance).

| Object | Role |
|---|---|
| `pdactrace_reference` | Main 10,113-gene atlas |
| `pdactrace_protein_betas` | Per-stage tissue-protein effect sizes (5,917 proteins) |
| `default_templates` | Canonical 12-template catalog |
| `pdactrace_external_anchors` | External anchor set for post-freeze evaluation |
| `meta_analysis` | Random-effects RNA meta-analysis summaries |
| `atlas_metadata` | Atlas provenance + cohort manifest |
| `pdactrace_data_sources` | 26 contributing public datasets |
| `inst/extdata/phase*.csv.xz` | Bundled downstream phase tables |

# Function reference

For every function, `?function_name` opens the full Rd page with
arguments, return shape, and a runnable example. Grouped summary
of the 59 exports:

| Role | Functions |
|---|---|
| **Lookup** | `query_gene`, `query_gene_detailed`, `query_panel`, `list_candidates`, `summarize_gene_evidence`, `case_study`, `format_provenance`, `list_atlas_metadata`, `atlas_provenance`, `list_data_sources` |
| **Audit score** | `compute_audit_score`, `explain_score`, `propagate_uncertainty`, `evaluate_anchor_enrichment` (alias `anchor_enrichment`), `evidence_math`, `explain_gene`, `compare_genes`, `compare_candidates`, `build_evidence_graph`, `extract_graph_features` |
| **Optional ML layer** | `make_evidence_features`, `score_anchor_similarity`, `fit_user_evidence_model`, `predict_user_evidence_model`, `explain_user_evidence_model`, `model_card` |
| **Trajectory + user cohort** | `fit_stage_de`, `fit_stage_de_protein`, `classify_trajectory`, `classify_protein_trajectory` (alias `classify_prot_trajectory`), `score_trajectory`, `align_patient_profile` (alias `align_patient`), `assemble_user_evidence`, `project_user_cohort`, `project_user_serum_cohort`, `early_pattern_names` (alias `early_patterns`), `mid_pattern_names_excluded` (alias `mid_patterns`) |
| **Filter audit** | `trace_filters`, `plot_filter_trace`, `plot_filter_diagnostics` |
| **Visualization** | `viz_gene`, `plot_gene_evidence`, `plot_gene_hexagon`, `plot_stage_effect`, `plot_per_cohort`, `plot_meta_forest`, `plot_serum_direction`, `plot_template_atlas`, `plot_gene_template`, `plot_panel_heatmap`, `plot_candidate_landscape`, `plot_celltype_full` |
| **Reporting** | `report_gene` |
| **Bioconductor-native** | `as_summarized_experiment`, `download_phase_csvs` |
| **Schema + theme** | `schema_spec`, `pdactrace_axes_theme`, `pdactrace_panel_theme`, `pdactrace_save`, `pdactrace_pal_class`, `pdactrace_pal_group`, `pdactrace_pal_dir`, `pdactrace_pal_pattern`, `NCS_W_SINGLE`, `NCS_W_15COL`, `NCS_W_DOUBLE`, `NCS_W_TRIPLE` |

All plot functions return a `ggplot` object; save via
`pdactrace_save(p, dir = "fig", name = "...", w = NCS_W_SINGLE)`
or `ggplot2::ggsave()`. `viz_gene(layout = "split")` returns a
named list of six ggplots so each panel can be saved at full
size.

# Limitations

- **PDAC-focused tool.** The bundled atlas, the curated
  external anchor set, and the validation work are PDAC-specific.
  The API is designed so the same evidence-aggregation and audit
  logic could be reused for other cancers given an analogous
  per-stage expression atlas, but no other-cancer atlas is
  shipped or evaluated.
- **Cohort independence assumption.** TCGA-PAAD and CPTAC-PDAC
  share donor recruitment; the heterogeneity gate caps the
  resulting Stouffer-Z inflation but does not eliminate it.
- **Held-out RNA-only baseline does not recover anchor
  enrichment** — multi-layer evidence is required (see
  Projection stress test above).
- **No supervised classifier shipped.** All ML is opt-in and
  user-trained on user-supplied labels; the package distributes
  no pretrained predictor.

# Vignettes

```r
vignette("lookup_basics",            package = "pdactrace")
vignette("audit_case_studies",       package = "pdactrace")
vignette("audit_framework",          package = "pdactrace")
vignette("user_cohort_extension",    package = "pdactrace")
vignette("methodology_validation",   package = "pdactrace")
vignette("reproducibility",          package = "pdactrace")
```

# Reproducibility

For ordinary use, **no network access is required**:

```r
library(pdactrace)
query_gene("LTBP1")
```

Four reproducibility layers (full walkthrough in
`vignette("reproducibility")`):

- **Layer 1 — bundled reference atlas.** Use the distributed
  `data/*.rda` directly. The default user path; covers all
  evidence lookup, scoring, reporting, visualization.
- **Layer 2 — user cohort projection.**
  `project_user_cohort()` runs entirely on local inputs.
- **Layer 3 — processed-input atlas rebuild.** Rebuild
  `data/pdactrace_reference.rda` from
  `data-raw/build_reference.R`. As of v0.99.6 the seven
  processed inputs (`multi_cohort_consistency.csv` plus six
  phase tables) are bundled in `inst/extdata/*.csv.xz`, so the
  rebuild runs without the manuscript-monorepo.
- **Layer 4 — raw-data reanalysis.** FASTQ / raw proteomics
  processing is **outside the scope of this software package**
  and is documented in the associated manuscript workflow
  archive.

The raw quantification matrices (per-cohort RNA count tables,
FragPipe protein intensity, scVI cell embedding) are
intentionally **not** bundled in the package tarball. A
companion **`pdactraceData` ExperimentHub package** (`biocViews:
ExperimentData, ExperimentHub`) is planned for a future release
and will allow lazy on-demand access via
`ExperimentHub::ExperimentHub()`.

# Citation

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

The accompanying *Briefings in Bioinformatics* manuscript
reference will be added once the preprint or journal DOI is
available.

# License

MIT. See [LICENSE](LICENSE).

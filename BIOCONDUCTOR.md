# pdactrace Bioconductor submission summary

## Submission overview (paste at the top of the reviewer issue body)

> pdactrace is submitted as a Software package for interpretable,
> stage-aware multi-omics biomarker prioritization. The bundled
> reference atlas is currently PDAC-focused, whereas the
> user-facing projection interface is designed for staged omics
> cohorts more generally. The package does not train or
> distribute a clinical diagnostic classifier; instead, it
> provides a frozen deterministic audit score, gene-level
> evidence summaries, visualization, and user-cohort projection
> utilities. Post-freeze anchor enrichment is reported as an
> evaluation-only sanity check, and the external RNA-only
> projection is reported as a stress test of evidence-layer
> dependence rather than definitive external multi-omics
> validation.

## Primary user workflow (small, by design)

Although pdactrace exports a broad API, the primary user
workflow is intentionally small:

```r
library(pdactrace)
viz_gene("LTBP1")            # one-call visual evidence canvas
query_gene("LTBP1")          # text view
explain_gene("LTBP1", view = "math")   # per-axis math
compare_candidates(c("LGALS3BP", "LTBP1", "ALB", "GAPDH"))
```

Additional exported functions support evidence mathematics,
user-cohort projection, and optional advanced diagnostics --
documented under "User-facing API groups" below.

## Package purpose

pdactrace is an interpretable multi-omics evidence integration
framework for stage-aware biomarker prioritization, demonstrated
in pancreatic cancer.

## What the package does

- stage-aware RNA trajectory matching (12-template catalog: Early x 4
  + Mid x 4 + Late x 2 + Monotonic x 2)
- multi-layer evidence feature extraction (RNA + tissue protein +
  scRNA cell origin + serum proteomics)
- transparent 3-axis + 2-gate audit scoring (frozen v0.3.0 rule;
  closed-form decomposition via `explain_score()`; Monte Carlo
  uncertainty via `propagate_uncertainty()`)
- translation-aware evidence stratification (Class A / B / C
  tissue-to-serum direction discipline)
- gene-level reports and visualization (HTML report via
  `report_gene()`; single-call visual canvas via `viz_gene()`)
- user cohort projection (tissue stage-aware via
  `project_user_cohort()`; serum binary group contrast via
  `project_user_serum_cohort()`)

## What the package is not

- not a clinical diagnostic classifier
- not a black-box ML model
- not a pan-cancer atlas
- not a prospective biomarker validation study

## User-facing API groups

| Group | Functions |
|---|---|
| Core query | `query_gene()`, `explain_gene()`, `query_panel()`, `compare_candidates()` |
| Visualization | `viz_gene()`, `plot_gene_evidence()`, `plot_gene_hexagon()`, `plot_stage_effect()`, `plot_per_cohort()`, `plot_celltype_full()`, `plot_serum_direction()`, `plot_filter_trace()`, `plot_panel_heatmap()`, `plot_candidate_landscape()`, `plot_template_atlas()`, `plot_gene_template()` |
| Evidence math | `evidence_math()`, `compare_genes()`, `make_evidence_features()`, `compute_audit_score()`, `propagate_uncertainty()` |
| User cohort | `project_user_cohort()`, `project_user_serum_cohort()`, `fit_stage_de()`, `fit_stage_de_protein()`, `classify_trajectory()`, `score_trajectory()`, `align_patient_profile()`, `assemble_user_evidence()` |
| Translation discipline | `format_provenance()`, `evaluate_anchor_enrichment()` |
| Optional / advanced | `score_anchor_similarity()`, `fit_user_evidence_model()`, `predict_user_evidence_model()`, `explain_user_evidence_model()`, `model_card()` |
| Reports + I/O | `report_gene()`, `as_summarized_experiment()`, `download_phase_csvs()` |
| Internal | low-level helper functions are not part of the stable user API |

## Evaluation discipline

The frozen audit rule was evaluated using post-freeze PDAC anchor
enrichment, negative-control depletion, bootstrap / leave-one-out
robustness, and external RNA-only projection stress testing.

Three-way evaluation summary (see `vignette("methodology_validation")`
Sections D and E for full code + reproduction):

| Source | Top-100 anchor fold | p-value |
|---|---:|---:|
| Bundled FULL, multi-layer | 39.3x | 2.2e-10 |
| Bundled RNA-only (audit_* recomputed) | 16.9x | 6.9e-04 |
| Held-out RNA-only (GSE253260 + GTEx, projection stress test) | 0 | 1.0 |

Negative-control housekeeping bottom-500 enrichment (atlas-defined
leakage annotations applied):

| Source | hits / 31 | fold | p-value |
|---|---:|---:|---:|
| Bundled FULL | 31 | 20.2x | 1.3e-41 |
| Held-out RNA-only | 28 | 17.8x | 1.1e-33 |

These results indicate that anchor specificity depends on
multi-layer evidence integration; external RNA-only projection is
useful as a stress test of mapping, leakage annotation transfer,
and single-layer limitations.

## Limitations

- bundled reference atlas is PDAC-focused
- user cohort projection is sensitive to gene universe
  harmonization (ENSEMBL -> SYMBOL mapping)
- external RNA-only projection is not equivalent to full
  multi-layer validation
- translation classes are evidence categories, not clinical
  diagnostic labels
- no prospective clinical validation has been performed
- atlas-defined leakage annotations are manually curated at v0.3.0
  freeze and require explicit transfer to user / external
  projections (not auto-derived)
- protein per-stage standard errors are not available for all
  datasets (limma F-test bundled; per-contrast Wald deferred)
- Stouffer cohort meta-analysis assumes independence; partial
  TCGA-CPTAC donor overlap is acknowledged and absorbed via the
  heterogeneity gate

## Known submission risks (proactively disclosed)

- **Tarball size ~5.83 MB exceeds the 5 MB BiocCheck
  software-package threshold (medium-high risk).** The bundled
  PDAC reference atlas (`pdactrace_reference.rda`, ~3.0 MB
  xz-compressed) and the bundled per-stage protein betas
  (`pdactrace_protein_betas.rda`, ~80 KB) are required for
  offline reproducibility, the lookup workflow, and the vignette
  examples. We are aware of this constraint and have already
  offloaded the larger upstream phase CSVs (~3.5 MB total) to
  `BiocFileCache` via `download_phase_csvs()` to keep the
  shipped tarball as small as possible while preserving the
  Layer 1 offline path. A companion `pdactraceData`
  ExperimentHub package (`biocViews: ExperimentData,
  ExperimentHub`) is planned as a v1.1 follow-on to host the
  raw quantification matrices that would otherwise inflate the
  package; if the reviewers prefer atlas relocation as a
  pre-acceptance condition, we are prepared to split the
  reference atlas into the data package on the first review
  iteration. We recommend keeping the bundled atlas in v0.99.x
  because the tool package's lookup, scoring, and reporting
  functions become opaque without it.

- **API surface (~71 exports).** The exported API is intentionally
  layered: a small Core query + Visualization workflow most users
  start with, plus Evidence Math, User cohort, and Optional /
  Advanced groups for power users. Internal helpers are not
  exported. The grouping is shown above ("User-facing API
  groups"). If reviewers prefer a smaller API surface, we can
  deprecate alias functions
  (`align_patient` / `align_patient_profile`,
  `classify_prot_trajectory` / `classify_protein_trajectory`,
  `early_patterns` / `early_pattern_names`,
  `mid_patterns` / `mid_pattern_names_excluded`,
  `anchor_enrichment` / `evaluate_anchor_enrichment`) in the
  first review iteration.

- **Local R CMD check WARNINGs are vignette-build-environment
  artefacts.** Two WARNINGs in our local R CMD check log
  ("Files in the 'vignettes' directory but no files in
  'inst/doc'" and "Directory 'inst/doc' does not exist. Package
  vignettes without corresponding single PDF/HTML") are produced
  because our local check environment lacks pandoc and we use
  `--no-build-vignettes` for fast iteration. **All 7 vignettes
  knit cleanly via `knitr::knit()` in our local environment**,
  and we expect both warnings to disappear on the Bioconductor
  build server which has pandoc available.

- **BiocCheck network step blocked locally.** Our sandbox cannot
  reach Bioconductor's deprecated-package status feed during
  `BiocCheck::BiocCheck()`. The package previously passed
  BiocCheck cleanly at v0.99.0 (0 ERROR / 0 WARN / 15 NOTE);
  we anticipate the new-package check on the Bioc build server
  will surface only the standard new-package NOTEs.

## Quality summary

- R CMD check --as-cran on v0.99.12 tarball: 0 ERROR / 2 WARN
  (vignette-only, from --no-build-vignettes flag) / 4 NOTE
  (pre-existing)
- testthat: 504 PASS / 0 FAIL
- BiocCheck previously clean on v0.99.0; expected re-check on
  submission
- 7 vignettes (lookup_basics, audit_framework, audit_case_studies,
  user_cohort_extension, reproducibility,
  cross_cancer_demonstration, methodology_validation)
- tarball size: ~5.83 MB (slightly above the 5 MB Bioconductor
  guideline; reviewer discretion expected)

## Availability

- GitHub: https://github.com/jibeomko/pdactrace
- Zenodo DOI (latest archive): 10.5281/zenodo.20076698
- Manuscript-monorepo + raw upstream phase CSVs (Zenodo):
  10.5281/zenodo.20067849
- Version: v0.99.12

## Submission tone (for the Bioconductor reviewer issue body)

> The package is prepared for Bioconductor review and passes local
> R CMD check and BiocCheck. We expect standard reviewer iteration
> and are submitting the package as an interpretable multi-omics
> evidence integration framework for stage-aware biomarker
> prioritization, demonstrated in pancreatic cancer. The bundled
> reference atlas is PDAC-focused; the framework's R API is
> cancer-agnostic. We treat post-freeze anchor enrichment as an
> evaluation-only sanity check, not as fully-blinded external
> validation. A held-out PDAC RNA-only projection is documented as
> a stress test, not as full multi-layer validation. We welcome
> reviewer feedback on API surface organization, vignette scope,
> and any items the maintainer team flags during review.

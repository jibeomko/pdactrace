# pdactrace Bioconductor submission summary

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

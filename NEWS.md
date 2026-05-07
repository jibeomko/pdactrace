# pdactrace 0.99.1

**Single-gene evidence report + plot polish.** Adds a one-page
self-contained HTML report for a single gene, and migrates
`plot_gene_hexagon()` onto the NCS-grade theme so all 10 plot
functions now share the same publication-grade typography.

## New

- `report_gene(gene_symbol, output_dir = tempdir())` — renders
  `inst/rmd/gene_report.Rmd` for one gene. Output is a
  self-contained HTML file with: audit-class tag, audit-component
  table, six-axis evidence hexagon (versus the
  `high_confidence_mean` reference), stage-trajectory plot,
  per-cohort breakdown, per-stage detail table, and the
  tissue-to-serum filter trace. `rmarkdown` and `knitr` are picked
  up at runtime from `Suggests` via `requireNamespace()` — no
  hard dependency added.

## Changed

- `plot_gene_hexagon()` now uses `pdactrace_panel_theme()` instead
  of bare `theme_void()`, and the four `geom_text` size settings
  (axis labels / score lines / comparison label / ring scale) are
  reduced to the NCS-grade range (1.7–2.2 mm ≈ 5–6.5 pt) so the
  hexagon embeds cleanly at single-column width (88 mm) without
  visible label clipping.

## Documentation

- `inst/rmd/gene_report.Rmd` — bundled report template, covered by
  R CMD check via `requireNamespace()` guards in `report_gene()`
  rather than as a vignette (it is parameterised on `gene` and
  intended for end-user invocation, not pre-render).

# pdactrace 0.99.0

**First Bioconductor submission version.** Per Bioconductor convention,
new packages submit at 0.99.0; subsequent versions in the submission
review queue go 0.99.1, 0.99.2, ... and on acceptance into the release
branch the version becomes 1.0.0. Functionality is unchanged from 0.4.1
beyond the Bioconductor pre-submission cleanup listed below.

## Bioconductor pre-submission cleanup

- `DESCRIPTION` — `Version: 0.99.0`, dropped redundant `Author:` /
  `Maintainer:` fields (kept `Authors@R` only, per BiocCheck), bumped
  `Depends: R (>= 4.5.0)`, added `withr` to `Imports`.
- `R/audit_score.R` — replaced the manual `set.seed()` +
  `on.exit()` save/restore block in `.audit_mc_table()` with a single
  `withr::local_seed(seed)` call (BiocCheck-compliant; functionally
  identical scope).
- `R/theme.R` — replaced `cat()` in `pdactrace_save()` with
  `message()` (BiocCheck disallows `cat` outside `show` methods).
- `R/pdactrace-package.R` — extended `globalVariables()` to cover
  `plot_gene_hexagon()`'s data.table NSE expressions (`y`, `x_end`,
  `y_end`) and silence the remaining R CMD check NOTE.

## Single-patient trajectory alignment

Adds one new exported function that takes a single patient's
tumor-vs-matched-normal log2 fold-change
profile and reports how it aligns with the frozen atlas's stage-trajectory
axes. The result is an *alignment readout*, not a stage prediction; no
supervised model is fit and the function is deterministic given (input,
atlas).

## New

- `align_patient_profile(rna_logfc, prot_logfc = NULL, ...)` — sample-level
  analog of `score_trajectory()`. Computes audit-weighted Pearson rho
  between the patient's per-gene log2FC vector and each of the three atlas
  stage axes (`rna_beta_E`, `rna_beta_M`, `rna_beta_L`), with Fisher-z 95%
  CI, plus a secondary audit-weighted vote share that surfaces
  disagreement-with-rho as a diagnostic. Returns a `pdactrace_patient_alignment`
  list with `$rna` (4-row data.table), `$prot` (categorical concordance,
  optional), `$summary` (one-line), and `$attrs` (provenance). Misalignment
  with all three stage axes is a legitimate output and does not imply
  Normal. Discipline-aligned: deliberately not a supervised classifier (a
  RF / GLM trained on TCGA stage labels would break the package's "no
  supervised fitting under limited ground truth" core principle).

## Documentation

- New "Step 6 — Single-patient trajectory alignment" section in
  `vignette("user_cohort_extension")` with a synthetic Late-stage
  example. Scope boundary block updated to enumerate the limitations
  (rho unstable below ~50 genes, protein arm categorical-only, etc.).

## Function name aliases

Five exported aliases were added alongside their fully-spelled originals
to reduce typing friction. Both forms refer to the same function and are
documented in the same Rd page; existing scripts using the long names
continue to work unchanged.

| Long name (kept) | Short alias (new) |
|---|---|
| `evaluate_anchor_enrichment()` | `anchor_enrichment()` |
| `early_pattern_names()` | `early_patterns()` |
| `mid_pattern_names_excluded()` | `mid_patterns()` |
| `align_patient_profile()` | `align_patient()` |
| `classify_protein_trajectory()` | `classify_prot_trajectory()` |

Total exports: 43 → 49.

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

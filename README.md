# Influenza A virus drives a strain-dependent activation programme and α-synuclein accumulation in human microglia

![License](https://img.shields.io/badge/license-GPL--3.0-blue)
![R](https://img.shields.io/badge/R-4.4.2-276DC3)
![Python](https://img.shields.io/badge/Python-3.9-3776AB)

Analysis code for the study of how human microglia (HMC3) respond to two H1N1 influenza A
strains (**A/PR/8/34** and **A/NWS/33**), with and without bacterial (LPS) priming and a
poly(I:C) comparator, profiled across **flow cytometry, ELISA, RT-qPCR, Western blot,
proteasome-activity assays and confocal microscopy**, and integrated at the single-culture
level. This repository contains the cleaned, reproducible scripts behind the statistics and
figures of the manuscript.

---

## Authors

Thi Len Ho<sup>1,#</sup>, Huu Dat Nguyen<sup>2,3,#</sup>, Eun-Ju Ko<sup>1,\*</sup>, Young Eun Kim<sup>3,\*</sup>

<sup>#</sup> These authors contributed equally. <sup>\*</sup> Corresponding authors.

## Affiliations

1. Interdisciplinary Graduate Program in Advanced Convergence Technology & Science, Jeju National University, Jeju, Republic of Korea
2. Laboratory of Parkinson’s Disease and Neurodegenerative disease, Hallym Institute for Translational Medicine, Anyang, Gyeonggi, Republic of Korea
3. Department of Neurology, Hallym University Sacred Heart Hospital, Hallym University, Anyang, Gyeonggi 14068, Republic of Korea

## Correspondence

- Prof. Eun-Ju Ko DVM., Ph.D.— [Email:]
- Prof. Young Eun Kim M.D., Ph.D. — [Email:]

---

## Repository structure

Each top-level folder is a self-contained analysis module; scripts are numbered in run order.

| Folder | Platform / analysis |
|---|---|
| `01_confocal/` | Confocal imaging of α-synuclein localisation relative to viral antigen |
| `02_proteasome/` | Chymotrypsin-like proteasome activity (Ac-ANW-AMC, Suc-LLVY-AMC) |
| `03_flow_elisa/` | Surface markers (flow) + secreted cytokines (ELISA): univariate → composite → multivariate, and single-cell infected-vs-bystander |
| `04_rtqpcr/` | RT-qPCR transcriptional profiling (8 genes) |
| `05_western_blot/` | Bayesian densitometry of α-synuclein species and signalling |
| `06_crossplatform_integration/` | Single-culture integration of α-synuclein with immune/transcriptional markers |

### Module detail

**`01_confocal/`**
`01_lif_to_tiff_conversion.py` → `02_cellpose_segmentation_quantification.py` →
`03_colocalization_analysis.py` (Costes/Manders) → `04_single_cell_mixed_models.R`
(linear mixed-effects models, single cells nested in fields of view; cluster bootstrap) →
`05_viral_load_correlations.R`.

**`02_proteasome/`**
`01_proteasome_activity_factorial.R` (per-well rate fitting, specific activity, 2 × 3 factorial
with HC3 robust SE and planned BH-adjusted contrasts); `02_rate_estimation_qc.R`.

**`03_flow_elisa/`**
`prep/` assembles tidy tables from FlowJo/ELISA exports. Whole-population pipeline:
`00_data_ingestion.R` → `01_univariate_factorial.R` (HC3, Holm) → `02_composite_indices.R`
(control-referenced effect sizes, surface/secreted composites, Decoupling & Polarization indices,
leave-one-marker-out, BCa bootstrap) → `03_multivariate_pca_permanova.R` (correlation-PCA +
PERMANOVA + betadisper). The `infected_vs_bystander/` sub-pipeline models raw event counts with
beta-binomial GLMMs (NP⁺ vs NP⁻): `00_count_ingestion.R`, `01_infection_kinetics.R`,
`02_infected_vs_bystander_glmm.R`, `03_paracrine_activation.R`, `04_dose_response.R`,
`05_bayesian_4h_sensing.R` (+ `_model_helpers.R`, `run_all.R`).

**`04_rtqpcr/`**
`01_preprocessing_pipeline.R` (.eds/.xls → per-replicate fold-change, SDM) →
`02_efficiency_calibration.R` (standard curves, geNORM with GAPDH+RPS11) →
`03_factorial_statistics.R` (RCBD `lm(log2 FC ~ Timepoint*LPS*Virus + Bio_Rep)`, Type III, Holm).

**`05_western_blot/`**
`00_config_and_ingestion.R` (config, ingestion, priors, MCMC/QC helpers) →
`01_membrane_diagnostic.R` → `02_bayesian_densitometry.R` (brms Phase-1 factorial + Phase-2 kinetic,
planned/polynomial contrasts, HC3 cross-check) → `03_sensitivity_analysis.R` (prior + ROPE) →
`04_supplementary_tables.R`.

**`06_crossplatform_integration/`**
`_labels.R` / `_bayes_common.R` (shared labels, MCMC config, convergence gate) →
`00_build_matched_dataset.R` → `01_within_condition_correlation.R` (Bayesian partial correlation,
one-step multivariate `rescor`) → `02_temporal_trajectories.R` → `03_same_lysate_correlation.R` →
`04_integrated_pca_permanova.R` → `05_sensitivity_analysis.R` → `06_supplementary_tables.R`.


## Statistical approach (summary)

- **Frequentist factorial** (flow, ELISA, RT-qPCR, proteasome): linear models with Type III ANOVA and
  **HC3 heteroscedasticity-consistent** robust standard errors; pre-specified `emmeans` contrasts;
  **Holm** within each marker/gene family (Benjamini–Hochberg for proteasome). Effect sizes: Cohen's *d*,
  partial η²ₚ, partial *r²*.
- **Clustered / count data**: linear mixed-effects models (confocal, field-of-view random effect) and
  **beta-binomial GLMMs** on raw event counts (infected-vs-bystander), with cluster bootstrap and a
  Bayesian cross-check.
- **Bayesian** (Western blot, cross-platform): `brms`/Stan with robust Student-*t* likelihoods and
  weakly-informative priors; reported as posterior **median, 95% HDI, probability of direction (pₐ)**
  and **ROPE**. All models met a strict convergence gate (**0 divergent transitions, R̂ ≤ 1.01, and
  effective sample size > 10,000** for the reported parameters), with **prior and ROPE-width sensitivity**
  analyses provided; reporting follows the **BARG** guidelines.
- Transforms: log₂ (densitometry, qPCR fold-change, cytokines), natural-log (confocal ratios), logit
  (flow proportions). Cross-platform and temporal links are interpreted as **associations, not causal effects**.

## Software / dependencies

- **R 4.4.2** — `lme4`, `lmerTest`, `emmeans`, `car`, `sandwich`, `effectsize`, `glmmTMB`, `brms`,
  `bayestestR`, `posterior`, `correlation`, `vegan`, `tidyverse`, `openxlsx`, `ggplot2`, `patchwork`,
  `ggdist`, `ggtext`.
- **Python 3.9** — `cellpose`, `python-bioformats`, `numpy`, `pandas`, `scikit-image`, `matplotlib`;
  (dry-lab) `scvi-tools`, `scanpy`, `anndata`.

## Data availability

Raw and processed data are not included here. Paths in the scripts are relative placeholders (`.`) and
reflect the original module layout; set the input paths to your own copy of the data. Source data and
per-test statistical outputs accompany the manuscript.

## License

Released under the **GNU General Public License v3.0 (GPL-3.0)**; see `LICENSE`.

## Citation

If you use this code, please cite the associated manuscript (details to be added on publication).

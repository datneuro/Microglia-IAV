# Influenza A virus infection drives strain-dependent α-synuclein remodeling with impaired proteasomal function in microglia

![License](https://img.shields.io/badge/license-GPL--3.0-blue)
![R](https://img.shields.io/badge/R-4.4.2-276DC3)
![Python](https://img.shields.io/badge/Python-3.9+-3776AB)

Analysis code for the manuscript. Human HMC3 microglia and mouse brain microglia were studied after
infection with two H1N1 influenza A strains, **A/PR/8/34** and **A/NWS/33**.

**Authors:** Thi Len Ho<sup>1,#</sup>, Huu Dat Nguyen<sup>2–4,#</sup>, Eun-Ju Ko<sup>1,5,\*</sup>, Young Eun Kim<sup>2–4,\*</sup>  


<sup>1</sup> Interdisciplinary Graduate Program in Advanced Convergence Technology & Science, Jeju National University, Republic of Korea ·  
<sup>2</sup> Laboratory of Parkinson’s Disease and Neurodegenerative Disease, Hallym Institute for Translational Medicine, Anyang, Republic of Korea ·  
<sup>3</sup> Department of Neurology, Hallym University Sacred Heart Hospital, Anyang, Republic of Korea ·  
<sup>4</sup> Hallym Neurological Institute, Hallym University, Anyang, Republic of Korea ·  
<sup>5</sup> Department of Veterinary Medicine, College of Veterinary Medicine, Jeju National University, Republic of Korea

  
<sup>#</sup> equal contribution; <sup>\*</sup> corresponding authors: eunju@jejunu.ac.kr, yekneurology@hallym.or.kr
## Repository structure

| Folder | Analysis | Manuscript |
|---|---|---|
| `01_Flow_ELISA/` | Flow cytometry + ELISA factorial models; infected-vs-bystander GLMMs | Fig. 1; Suppl. Figs. S2–S4; Suppl. Tables S3–S6 |
| `02_RT-qPCR/` | Standard curves, geNORM, SDM Cq, Pfaffl fold change, factorial models | Fig. 2a; Suppl. Table S7 |
| `03_Proteasome/` | Chymotrypsin-like activity (rate estimation, factorial model) | Fig. 2b–c |
| `04_Western_Blot/` | Bayesian densitometry (brms) | Figs. 3–4; Suppl. Tables S8–S9 |
| `05_Confocal/` | LIF conversion, Cellpose-SAM segmentation, colocalisation, mixed models | Fig. 5; Suppl. Fig. S6 |
| `06_CrossPlatform_Integration/` | Bayesian partial correlations, PCA/PERMANOVA | Fig. 6; Suppl. Fig. S7; Suppl. Table S10 |
| `07_In_Vivo/` | Body weight, lung/brain ELISA, lung/brain flow cytometry | Figs. 7–8; Suppl. Figs. S8–S10 |
| `Figures/` | Figure assembly and frequentist supplementary tables | Figs. 1–6 |

Scripts are numbered in run order; files starting with `_` are helpers sourced by the module.
Run every script **from the repository root**.

| Module | Run order |
|---|---|
| `01_Flow_ELISA/` | `Prep/01–03` → `00–03`; `Infected_vs_Bystander/Run-All.R` |
| `02_RT-qPCR/` | `00` → `01` → `02` → `03` |
| `03_Proteasome/` | `02` (rate windows) → `01` |
| `04_Western_Blot/` | `00` → `01` → `02` → `03` → `04` |
| `05_Confocal/` | `01` → `02` → `03` (Python) → `04` → `05` → `06` (R) |
| `06_CrossPlatform_Integration/` | `00` → `06` |
| `07_In_Vivo/` | `Body_Weight/01–03`, `ELISA/01–02`; `Flow/`: brain `01`–`73`, lung `79`–`89`, integration `90`–`94` (settings in `Flow/Config.R`) |
| `Figures/` | after the modules; Fig. 5a–b: `Figure-05_N2V-Prepare-Data.py` → `Figure-05_N2V-Train.py` → `Figure-05_Representative-Images.py` |


## Data layout

Data are not included. Place (or symlink) the data folders in the repository root:

| Folder | Used by |
|---|---|
| `Flow-ELISA/` | `01_Flow_ELISA/`, `06_`, `Figures/` |
| `2606_PCR_VirusHMC3/` | `02_RT-qPCR/`, `06_`, `Figures/` |
| `Proteasome/` | `03_Proteasome/`, `Figures/` |
| `Western blot/` | `04_Western_Blot/`, `06_`, `Figures/` |
| `Confocal/` | `05_Confocal/`, `Figures/` |
| `CrossPlatform_aSyn/` | `06_CrossPlatform_Integration/`, `Figures/` |
| `ELISA/`, `Body weight/` | `07_In_Vivo/` |
| FCS files (`INVIVO_DATA`) | `07_In_Vivo/Flow/` |

Source data are available from the corresponding authors upon reasonable request.

## Software

- **R 4.4.2**: tidyverse, readxl, openxlsx, lme4, lmerTest, emmeans, car, sandwich, lmtest, effectsize,
  glmmTMB, brms, bayestestR, posterior, correlation, vegan, patchwork, cowplot, ggdist, ggtext,
  ggrepel, ggsignif, ctrlGene; in vivo: flowCore, nlme, survival and the packages in `07_In_Vivo/Flow/Common.R`.
- **Python 3.9+**: aicsimageio[bioformats], tifffile, cellpose, csbdeep, numpy, pandas, scipy,
  scikit-image, matplotlib, openpyxl; careamics 0.0.15 (Fig. 5a–b only).

## License and citation

GPL-3.0 (see `LICENSE`). Please cite the article and this code (see `CITATION.cff`).

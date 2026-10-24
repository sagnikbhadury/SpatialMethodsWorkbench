# Reproducibility materials for the JASA ACS submission

Manuscript: "Spatially Varying Conditional-Dependence Networks for Spatial Transcriptomics and Multiplexed Imaging" (GP-GHS).

JASA requires the Author Contributions Checklist (ACC) form, code, data and the workflow with the
*revised* manuscript (reproducibility review starts when a revision is invited); they are encouraged
at initial submission. This folder documents what exists and how to package it. The ACC form itself is
filled in on the journal's web form (https://jasa-acs.github.io/repro-guide/pages/acc.html); a draft of
every answer is in `ACC_form_draft.md`.

## What to deposit

Suggested layout of `reproducibility_materials.zip` (the repro-guide template at
https://github.com/jasa-acs/repro-template is optional):

```
reproducibility_materials/
  README.md                     this file, with the run order below
  code/
    gp_group_horseshoe_spatial_graph.R   sampler, Hilbert-space basis, summaries (sourced by everything)
    tests/test_sampler_quadrature.R      sampler check against numerical integration (Section S1)
    simulation_study_v4.R                main simulation (Section 3, Tables 1 and S1, Figures S1-S2)
    simulation_figures_v4.R              simulation figures and tables
    simulation_field_detail_v4.R         recovery of the coefficient functions (Table S2)
    simulation_gene_regime_v4.R          setting with the dimensions of the gene application (Table S3)
    simulation_counts_v4.R               count-valued experiment (Section S2.6)
    st_prepare.R, st_fit.R, st_network_summary.R, st_zone_figure.R, st_checks.R
                                         index section: gene selection, fits, Table 2, Figures 2-3, S3, permutation control
    st_prepare_panel.R, run_st_all_slides.ps1, st_replication_summary.R, st_replication_figure.R,
    st_replication_intervals.R           the other 13 sections, Table 3 (bootstrap intervals), Figure 4, Table S4
    crc_analysis.R, crc_network_summary.R, crc_sensitivity.R, crc_group_gee.R, crc_grid_sensitivity.R
                                         CODEX fits, Table 4, Figures S4-S5, Sections S4.3-S4.5
    make_workflow_figure.R               Figure 1
  data/
    CRC_data_cleaned_corrected.csv       CODEX cell table (Schurch et al. 2020; 140 images, phenotypes supplied with the data)
    visium/                              not included: download from the CZ CELLxGENE collection of Valdeolivas et al. (2024);
                                         st_prepare.R documents the h5ad files used and the pathologist annotation field
  output/                                the saved summaries that the manuscript numbers come from (csv and log files)
```

## Software

R 4.5.2 for the simulations (packages MASS, glmnet, glasso, parallel, geepack) and R 4.4.1 for the
data analyses (lme4, lmerTest, hdf5r, posterior, igraph, ggplot2). All runs used a 20-core workstation;
seeds are fixed inside each script (SEED_BASE constants).

## Run order and timings

```
Rscript tests/test_sampler_quadrature.R                       # R 4.5.2, seconds
Rscript simulation_study_v4.R ; Rscript simulation_figures_v4.R   # about 2.5 h on 18 cores
Rscript simulation_field_detail_v4.R                          # about 25 min on 8 cores
Rscript simulation_gene_regime_v4.R                           # 96 min on 18 cores
Rscript simulation_counts_v4.R
# R 4.4.1 from here
Rscript st_prepare.R ; Rscript st_fit.R                       # index section, 14 min on 18 cores
ST_PERMUTE=1 ST_TAG=_permuted Rscript st_fit.R                # permutation control
Rscript st_network_summary.R ; Rscript st_zone_figure.R ; Rscript st_checks.R
powershell -File run_st_all_slides.ps1                        # other 13 sections, about 1 h 20 min
Rscript st_replication_summary.R ; Rscript st_replication_figure.R ; Rscript st_replication_intervals.R
CRC_FEATURE=count CRC_N_GRID=15 CRC_M_BASIS=4 CRC_OUT_DIR=crc_results_counts Rscript crc_analysis.R   # 12 min, 15 cores
CRC_OUT_DIR=crc_results_counts Rscript crc_network_summary.R
CRC_FEATURE=count CRC_N_GRID=15 CRC_M_BASIS=4 CRC_DETREND=1 CRC_FIT_ONLY=1 CRC_OUT_DIR=crc_results_counts_detrend Rscript crc_analysis.R
CRC_FEATURE=count CRC_N_GRID=15 CRC_M_BASIS=4 CRC_MIN_CELLS=10 CRC_FIT_ONLY=1 CRC_OUT_DIR=crc_results_counts_min10 Rscript crc_analysis.R
CRC_FEATURE=count CRC_N_GRID=20 CRC_M_BASIS=4 CRC_FIT_ONLY=1 CRC_OUT_DIR=crc_results_counts_grid20 Rscript crc_analysis.R
CRC_FEATURE=count CRC_N_GRID=10 CRC_M_BASIS=4 CRC_FIT_ONLY=1 CRC_OUT_DIR=crc_results_counts_grid10 Rscript crc_analysis.R
Rscript crc_sensitivity.R ; Rscript crc_grid_sensitivity.R
Rscript crc_group_gee.R                                       # R 4.5.2 (geepack)
Rscript make_workflow_figure.R                                # needs the Visium and CODEX fits
```

## Anonymization

For review, deposit the materials without author names (the scripts carry no names) or through an
anonymous repository (for example Anonymous GitHub); the public repository named in the manuscript
replaces it after acceptance.

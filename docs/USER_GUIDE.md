# Usage, restrictions, and local installation

The Spatial Methods Workbench is a research application for guided analysis of spatial, imaging, network, predictive, causal, and structured biomedical data. It validates an uploaded table, shows which pipelines satisfy their minimum input contract, runs a selected workflow, and exports results with settings, software objects, provenance, and citations.

> The application computes statistical outputs but does not explain results, assign biological mechanisms, or make scientific or clinical conclusions. Contact Sagnik Bhadury at [bhadury@umich.edu](mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration) for study-specific interpretation or collaboration.

## Restrictions and responsibilities

1. **Research use only.** The suite is not a medical device and must not be used for diagnosis, treatment selection, patient management, or another clinical decision.
2. **No protected health information on the public service.** Do not upload names, medical-record numbers, dates of birth, addresses, detailed dates, free text, rare identifiers, raw DICOM headers, or other protected/re-identifiable information. Use a locally controlled installation for sensitive data and follow institutional policy and IRB/data-use requirements.
3. **Uploads are session-scoped, not a durable data store.** The application code does not intentionally save uploaded data, but hosting infrastructure may maintain operational logs. A sleeping, restarted, or terminated instance loses session state.
4. **Users remain responsible for preprocessing.** Registration, segmentation, normalization, batch correction, quality control, coordinate units, zone definitions, feature identity, and subject independence must be established before analysis.
5. **Compatibility is necessary, not sufficient.** A green readiness indicator means required columns and software are present. It does not establish that the sampling design, estimand, sample size, model assumptions, or biological question are appropriate.
6. **Association is not causation.** Network edges are conditional-association summaries, not automatic evidence of physical contact, signaling, intervention effects, or causal mechanisms.
7. **Prediction requires honest validation.** Random held-out splits may be optimistic when neighboring observations, sections, fields of view, or repeated measures leak across folds. Use patient-level, slide-level, spatially blocked, temporal, or external validation as required by the design.
8. **Imaging analyses require alignment.** Wide-image columns must represent the same registered pixels/voxels/features in the same order for every subject. Do not mix masks, resolutions, coordinate systems, or preprocessing versions.
9. **Causal analyses require design assumptions.** Mediation output is model-based. Temporal ordering, consistency, positivity, correct specification, and control of exposure–mediator, exposure–outcome, and mediator–outcome confounding must be defended outside the app.
10. **Advanced engines are compute intensive.** ISPat, GP–GHS, and ISPat 3D may exceed public-host memory, CPU, connection, or active-hour limits for real studies. Large analyses should run locally, on an institutional cluster, or on a controlled container/Posit Connect deployment.
11. **Review all outputs.** Check convergence, sensitivity to tuning parameters, missing-data handling, uncertainty, multiplicity, robustness, and domain plausibility before reporting results.
12. **Cite the software and method.** Acknowledgement is required before execution. Every download includes the workbench citation and applicable method/publication references.

## Quick-start workflow

1. Open the application and begin with the synthetic demonstration.
2. Upload a comma-separated file (`.csv`) with one row per observational unit for the intended pipeline.
3. Map coordinate, zone, ID, outcome, exposure, mediator, and feature columns. Unused roles may remain unselected.
4. Review **Readiness**. A blocked method states the missing columns or server engine.
5. Select only the features intended for the analysis. Avoid identifiers, leakage variables, and post-outcome measurements.
6. Configure the method-specific controls.
7. Acknowledge the citation requirement and run the analysis.
8. Inspect numerical tables and figures. The app does not provide a scientific explanation.
9. Download and retain the reproducibility bundle with the analysis record.
10. Seek statistical/domain review before publication or decision-making.

## Input contracts

### Cell-, spot-, or location-level spatial table

Use one row per cell, spot, grid location, or field-level observation. Typical columns are:

- numeric `x` and `y` coordinates;
- optional numeric `z`/section/depth coordinate;
- a region, density zone, compartment, or cluster label;
- a sample, patient, field, contour, or volume ID;
- numeric marker, abundance, expression, or engineered feature columns;
- an outcome/exposure/mediator only where the corresponding analysis requires it.

Coordinates must use a consistent unit and reference frame. Repeated subjects or slides must be represented explicitly and handled in validation outside workflows that assume independent rows.

### Tumor-contour table

Use one row per ordered boundary landmark with a contour ID and numeric x/y coordinates. Every contour needs at least six ordered points. The workflow closes, centers, size-normalizes, and resamples contours; it does not repair self-intersections or infer boundary order.

### Scalar-on-image wide table

Use one row per subject and one selected numeric column per aligned pixel, voxel, or image-derived feature. Map a numeric scalar outcome. The public interactive workflow limits selection to 5,000 image features and is a fast ridge screen, not the full SV-NN or ST-CAR posterior algorithm.

### Image-to-image wide table

Use one row per subject. Prefix aligned predictor-image columns with `input__` and outcome-image columns with `output__`; select both groups as features. At least three columns in each group and 30 complete subjects are required. This is a low-rank screening workflow, not the full SBLF posterior engine.

### ISPat 3D table

Use one row per volumetric location/cell with x/y/z coordinates, a zone, a volume or spot ID, and multiple numeric markers. Serial sections must already be registered into a common 3D coordinate system. Validate section spacing, missing sections, tissue masks, orientation, and zone consistency before fitting.

## Pipeline requirements

| Pipeline | Minimum mapped input | Main restriction |
|---|---|---|
| Spatial exploration and Moran's I | x, y, one numeric feature | Descriptive under a selected k-nearest-neighbor graph |
| Region-aware conditional network | at least three numeric features; zone optional | Ridge exploratory network; edges are not causal |
| ISPat 2D | x, y, zone, at least three markers, installed `ISPAT` | Bayesian/variational engine; intensive |
| GP–GHS | x, y, at least three markers, installed `GPGHS` | MCMC/shrinkage settings and runtime require review |
| ISPat 3D | x, y, z, zone, at least three markers, installed `ISPAT3D` | Registered volumetric data; intensive |
| Spatial machine-learning prediction | outcome and at least three features | Random held-out split is exploratory |
| Spatial phenotype clustering | x, y, at least three features | Cluster count and coordinate influence are user choices |
| Shallow-neural prediction | outcome and at least three features | Binary classification or numeric regression; max 250 predictors |
| Wide-image scalar regression | numeric outcome and at least five image features | Max 5,000 aligned features; ridge screen |
| Latent image-to-image regression | at least three `input__*` and three `output__*` features | Aligned wide images; low-rank screen |
| Mediation | numeric exposure, mediator, and outcome | Model output is not automatically causal |
| Tumor shape PCA | contour ID plus ordered x/y landmarks | Five or more valid contours |

## Preparing real-life data

Before upload, create a frozen analysis-ready dataset and data dictionary. At minimum:

- define the observational unit and independence/repeated-measures structure;
- remove direct and indirect identifiers;
- preserve a read-only raw-data copy outside the app;
- document assay, segmentation, registration, normalization, transformation, and filtering versions;
- confirm coordinate axes, units, origins, orientation, and section spacing;
- harmonize feature names and biological identities across samples;
- inspect missingness, zero variance, extreme values, sparsity, and class imbalance;
- define training/test grouping before predictive analysis;
- prespecify zones, outcomes, exposures, mediators, tuning choices, and primary contrasts where possible;
- keep a row/sample exclusion log and a checksum/version for the uploaded file.

## Install the complete suite locally

Local installation is recommended for real data, long-running Bayesian analyses, and controlled information environments.

### Prerequisites

- Git
- R 4.5.x (or a compatible current R release)
- RStudio is optional but convenient
- Windows users should install the Rtools version matching their R version when compiling packages
- Adequate RAM and CPU for the selected model; Bayesian and 3D workloads may require a workstation or cluster

### Clone the application

```text
git clone https://github.com/sagnikbhadury/SpatialMethodsWorkbench.git
cd SpatialMethodsWorkbench
```

### Install from a Zenodo release instead

Download the newest source ZIP from the [Zenodo concept record](https://doi.org/10.5281/zenodo.21763606), extract it, and open a terminal in the extracted directory. The complete distribution includes the exact public ISPAT, GP-GHS, and ISPAT-3D package sources, so Zenodo users can run:

```text
Rscript install-complete.R
Rscript verify-installation.R
Rscript tests/testthat.R
```

Use `Rscript install-complete.R --include-related` to compile the bundled STCAR package and, on Linux, SBLF. SV-NN and the selected public AI/ML/deep-learning sources are archived under `vendor/reference-workflows` for reproducibility and attribution; they are not silently substituted for the explicitly labelled built-in screening and adaptation paths. See `vendor/manifest.json` for exact commits, licenses, retained scope, and role.

The installer may obtain missing CRAN dependencies, but it does not clone the three direct advanced-engine repositories. A source archive cannot provide one portable compiled R library for Windows, macOS, and Linux.

### Install the interface and built-in workflow packages

Run in R:

```r
install.packages(c(
  "shiny", "bslib", "DT", "ggplot2", "jsonlite",
  "markdown", "zip", "testthat", "nnet", "remotes"
))
```

### Install all public advanced network engines

Run in R:

```r
remotes::install_github("sagnikbhadury/GP-GHS", upgrade = "never")
remotes::install_github("sagnikbhadury/ISPAT-3D", upgrade = "never")
remotes::install_github("sagnikbhadury/ISPAT", upgrade = "never")
```

These install `GPGHS`, `ISPAT3D`, and `ISPAT` plus their declared dependencies. ISPAT uses Stan; a first source installation or model compilation can take substantial time. Keep R, compiler, Stan, and package versions fixed for a formal analysis.

### Verify the installation

From the repository directory:

```text
Rscript tests/testthat.R
```

In R:

```r
stopifnot(
  requireNamespace("ISPAT", quietly = TRUE),
  requireNamespace("ISPAT3D", quietly = TRUE),
  requireNamespace("GPGHS", quietly = TRUE)
)
```

### Run locally

```r
shiny::runApp(".", host = "127.0.0.1", port = 3838)
```

Open `http://127.0.0.1:3838`. Binding to `127.0.0.1` keeps the service on the local machine. Do not bind to a public interface without authentication, TLS, firewall rules, upload limits, monitoring, and a reviewed data-retention policy.

### Run with Docker

For the built-in workflows:

```text
docker compose up --build
```

To build the image with public advanced engines:

```text
docker build --build-arg INSTALL_ADVANCED_ENGINES=true -t spatial-methods-workbench:advanced .
docker run --rm -p 127.0.0.1:3838:3838 spatial-methods-workbench:advanced
```

Pin image and package revisions before a production analysis. Do not mount sensitive directories read-write unless explicitly required.

## Reproducibility bundle

The downloaded ZIP contains the results table, figure where applicable, `manifest.json`, an R result object, `CITATION.txt`, and `COLLABORATION.txt`. It deliberately contains no automated scientific interpretation. Preserve the original input separately under the study's approved governance plan; the bundle does not duplicate the uploaded dataset.

## Public-host limitations

The shinyapps.io service is appropriate for demonstrations and moderate public research workloads. It is not the preferred location for protected data or large Bayesian/3D jobs. Session timeouts, active-hour quotas, instance memory, concurrent-worker limits, and package compilation affect availability. For operational research use, prefer a controlled local installation, institutional HPC batch workflow, Posit Connect, or a dedicated container deployment with authentication and monitoring.

## Citation and collaboration

Use the workbench citation included in the bundle and the method-specific citations added for the selected pipeline. Cite the version-specific DOI shown on the downloaded record. The all-versions DOI is [10.5281/zenodo.21763606](https://doi.org/10.5281/zenodo.21763606).

For method selection, scientific interpretation, a production deployment, or collaboration, email [bhadury@umich.edu](mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration).

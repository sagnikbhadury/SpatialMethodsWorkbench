# Spatial Methods Workbench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21763606.svg)](https://doi.org/10.5281/zenodo.21763606)

A functional Shiny application for guided analysis of spatial and structured biomedical data. Users upload a CSV, map its columns, see which analyses are compatible, configure a method, run it, and download a reproducibility bundle containing results, a figure, settings, citations, and the R result object.

The public application is intentionally separated from private research code. Its method registry contains only workflows intentionally implemented here or adapters to already-public packages.

## What works now

- Spatial feature maps and permutation-based Moran's I.
- Region-aware ridge-regularized partial-correlation networks.
- Spatial prediction with a held-out test set.
- Spatially weighted phenotype clustering adapted from public machine-learning curricula.
- Shallow-neural spatial prediction with held-out validation.
- Wide-image scalar regression screening for aligned pixel or voxel columns.
- Latent image-to-image regression using `input__*` and `output__*` feature groups.
- Bootstrap mediation analysis with explicit causal caveats.
- Contour alignment, resampling, and shape PCA.
- Optional adapters for the public `ISPAT`, `GPGHS`, and `ISPAT3D` R packages. The UI exposes these only when the corresponding package is installed.
- Synthetic demonstration data, automatic column suggestions, structural validation, compatibility gating, and downloadable result bundles.
- No automated interpretation of results. Researchers are directed to contact Sagnik Bhadury for study-specific guidance or collaboration.

The advanced Bayesian engines are computationally intensive. Production deployments should enforce upload, memory, worker, and runtime limits and should use a job queue for large studies.

## Expected tabular input

Use one row per observation (cell, spot, location, landmark, or subject-feature record). A single CSV may contain:

- numeric `x`, `y`, and optional `z`/section coordinates;
- a region or tumor-density zone;
- a sample, patient, contour, spot, or volume identifier;
- an outcome, exposure, and mediator when relevant;
- three or more numeric marker/feature columns for network analyses.

The application never assumes column names: users confirm every mapping. The bundled synthetic dataset opens by default and demonstrates the required structure.

See the complete [Usage, Restrictions, and Local Installation Guide](docs/USER_GUIDE.md) before applying the suite to real data. The same guide is available inside the application.

For wide-image workflows, use one row per subject. Scalar-on-image screening treats the selected numeric columns as aligned image features. Image-to-image screening requires selected predictor columns named `input__*` and outcome-image columns named `output__*`. Registration, mask, resolution, and feature order must already be harmonized.

The spatial AI/ML workflows are deployable adaptations informed by the public Microsoft AI/ML curricula and deep-learning course repositories; they do not claim that a curriculum repository is a statistical package. The image-regression screens cite the related SV-NN, ST-CAR, and SBLF publications and state explicitly when the interactive implementation is a faster screening model rather than the full Bayesian posterior engine.

No uploaded values or analysis results are sent to OmniRoute or an LLM. Public LLM engineering repositories inform deployment design only; automated result explanation remains outside the application by policy.

## Run locally

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "jsonlite", "markdown", "zip"))
shiny::runApp(".")
```

Or from a terminal with this repository as the working directory:

```text
Rscript run-local.R
```

Open `http://127.0.0.1:3838`.

## Install public advanced engines

```r
install.packages("remotes")
remotes::install_github("sagnikbhadury/ISPAT")
remotes::install_github("sagnikbhadury/GP-GHS")
remotes::install_github("sagnikbhadury/ISPAT-3D")
```

The app checks availability with `requireNamespace()` and never substitutes a different algorithm while displaying an advanced method's name.

## Test

```text
Rscript tests/testthat.R
```

The suite covers every built-in analysis, data compatibility, result bundles, and a Shiny server execution path.

## Deploy

The repository is suitable for:

- **shinyapps.io** for an initial public deployment of the fast modules;
- **Posit Connect** for authenticated access, process controls, and scheduled content;
- **Docker/Shiny Server** for a dedicated host with compiled advanced engines.

GitHub Pages cannot execute Shiny. Link the deployed application from `sagnikbhadury.github.io`, while keeping this repository as the source of record.

## Privacy and interpretation

Do not upload protected health information to a public instance. Uploaded data live in the active Shiny process and are not intentionally persisted by this application, but the host's infrastructure and logs must also be configured appropriately. See [SECURITY.md](SECURITY.md).

All network edges are conditional associations unless a method explicitly establishes something stronger. They are not automatic evidence of physical contact, molecular signaling, intervention effects, or causation. This software is for research, not clinical decision-making. The application deliberately does not interpret results; contact [Sagnik Bhadury](mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration) for study-specific interpretation or collaboration.

## Citation

Running an analysis requires acknowledgement of the citation condition. Every result bundle includes a `CITATION.txt` file and records the applicable citations in `manifest.json`. GitHub also renders the repository-level `CITATION.cff` through its **Cite this repository** interface.

When using ISPat, also cite:

> Bhadury S, et al. Informed spatially aware patterns for multiplexed immunofluorescence data. *Scientific Reports*. 2026;16:5015. https://doi.org/10.1038/s41598-026-35341-8

Each advanced engine retains its own citation and license requirements.

Citation acknowledgement is a strong research-norm and provenance mechanism, not a technical guarantee about a later publication's bibliography. A tagged release can be archived with Zenodo to add a persistent DOI.

- Cite the current `v0.2.0` release using [DOI 10.5281/zenodo.21763796](https://doi.org/10.5281/zenodo.21763796).
- Use the [concept DOI 10.5281/zenodo.21763606](https://doi.org/10.5281/zenodo.21763606) when referring to the workbench across all versions.

## Design evidence

The product and literature rationale, scope boundaries, and deployment recommendation are documented in [docs/LITERATURE_REVIEW.md](docs/LITERATURE_REVIEW.md).

# Spatial Methods Workbench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21763606.svg)](https://doi.org/10.5281/zenodo.21763606)
[![Version](https://img.shields.io/badge/release-v0.3.0-176b68)](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/releases/tag/v0.3.0)
[![License: MIT](https://img.shields.io/badge/License-MIT-d46245.svg)](LICENSE.md)
[![R tests](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/actions/workflows/test.yml/badge.svg)](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/actions/workflows/test.yml)

A functional Shiny application for guided analysis of spatial and structured biomedical data. Users upload a CSV, map its columns, see which analyses are compatible, configure a method, run it, and download a reproducibility bundle containing results, a figure, settings, citations, and the R result object.

**[Launch the live application](https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/)** · **[Read the complete usage and installation guide](docs/USER_GUIDE.md)** · **[Open the function-by-function reference](docs/PIPELINE_REFERENCE.md)** · **[Cite the software](https://doi.org/10.5281/zenodo.21763606)**

## Video tutorial

[![Watch the Spatial Methods Workbench tutorial on YouTube](https://img.youtube.com/vi/Ft9groM2ilI/maxresdefault.jpg)](https://www.youtube.com/watch?v=Ft9groM2ilI)

**[Watch the complete Spatial Methods Workbench usage tutorial on YouTube](https://www.youtube.com/watch?v=Ft9groM2ilI).** It demonstrates the live application, synthetic-data workflow, column mapping, compatibility checks, analysis selection, controls, citation requirements, reproducibility downloads, and local installation. The tutorial explains operation and statistical purpose without interpreting a user's scientific results.

## Reproducible vignettes

**[Open the simulated-data vignette index](docs/VIGNETTES.md).** Five executable vignettes reproduce the Shiny workflow and every analysis family using small deterministic datasets. Upload-ready CSV files are included under [`inst/extdata`](inst/extdata), and every dataset can be regenerated from a fixed seed with the exported simulation functions.

The public application is intentionally separated from private research code. Its method registry contains only workflows intentionally implemented here or adapters to already-public packages.

## Analyses available in the application

The Workbench currently hosts 12 executable analysis paths. It first checks the uploaded data against each method's input contract and exposes only compatible choices.

| Analysis | Designed for | Required data | Principal output | Implementation scope |
|---|---|---|---|---|
| Spatial exploration and autocorrelation | Detecting whether a measured feature is spatially patterned | Numeric x/y coordinates and one numeric feature | Spatial feature map, Moran's I, permutation p-value | Built-in executable diagnostic |
| Region-aware conditional network | Exploring adjusted relationships overall and within tissue regions | At least three numeric features; region is optional | Ridge-regularized partial-correlation edge tables and heat maps | Built-in exploratory network workflow |
| ISPat 2D Bayesian network | Shared and region-specific network estimation after spatial adjustment | x/y coordinates, tissue zone, and at least three numeric markers | Shared and zone-specific conditional-association networks and fitted object | Direct adapter to the public `ISPAT` package |
| Spatially varying GP–horseshoe network | Allowing network edges to vary continuously over tissue space | x/y coordinates and at least three numeric features | Posterior-selected adjacency, edge table, and full fit | Direct adapter to the public `GPGHS` package |
| Volumetric ISPat 3D network | Network estimation from registered 3D images or serial sections | x/y/z coordinates, zone, volume/spot ID, and numeric markers | Shared and zone-specific 3D networks and fitted objects | Direct adapter to the public `ISPAT3D` package |
| Spatial machine-learning prediction | Predicting a numeric outcome from measured covariates and optional coordinates | Numeric outcome and at least three numeric predictors | Held-out predictions, error metrics, and variable importance | Built-in held-out prediction workflow |
| Spatially weighted phenotype clustering | Discovering multivariate phenotypes while controlling spatial influence | x/y coordinates and at least three numeric features | Cluster assignments, centers, diagnostics, and spatial cluster map | Built-in spatially weighted k-means adaptation |
| Spatial shallow-neural prediction | Nonlinear prediction from features and optional spatial coordinates | Numeric outcome and at least three numeric predictors | Held-out predictions and performance metrics | Executable single-hidden-layer `nnet` workflow |
| Wide-image scalar regression screen | Screening aligned pixel, voxel, or image-derived features against one scalar outcome | One row per subject, numeric outcome, and at least five aligned image features | Ridge predictions, test error, and ranked coefficients | Fast screening model; not presented as the full SV-NN or ST-CAR Bayesian algorithm |
| Latent image-to-image regression | Screening associations between aligned predictor and outcome images | One row per subject with at least three `input__*` and three `output__*` columns | Low-rank latent factors, reconstructed outcomes, and performance summaries | Fast latent-factor screen; not presented as the full SBLF posterior engine |
| Mediation with sensitivity diagnostics | Estimating a model-based indirect effect | Numeric exposure, mediator, and outcome | Direct, indirect, and total effects with bootstrap uncertainty | Built-in model-based workflow with explicit causal assumptions |
| Tumor contour and shape PCA | Summarizing variation among ordered two-dimensional tumor boundaries | Contour ID and ordered x/y landmark coordinates | Resampled/aligned contours, shape scores, variance summaries, and mean shape | R implementation of the public shape workflow |

Every path returns numerical/statistical output rather than an automated scientific narrative. Network edges are conditional associations unless the selected method and study design establish something stronger. For interpretation, method selection, or a study collaboration, contact [Sagnik Bhadury](mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration).

The application also includes synthetic demonstration data, automatic column suggestions, structural validation, compatibility gating, documented controls, and downloadable result bundles. The complete v0.3.0 source distribution includes local package sources for all three advanced engines; the live deployment includes all three.

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

The [Function-by-Function Pipeline Reference](docs/PIPELINE_REFERENCE.md) documents every executable workflow, control, input contract, output object, assumption, common failure, citation, and programmatic example. It is also rendered as a dedicated tab in the application.

For wide-image workflows, use one row per subject. Scalar-on-image screening treats the selected numeric columns as aligned image features. Image-to-image screening requires selected predictor columns named `input__*` and outcome-image columns named `output__*`. Registration, mask, resolution, and feature order must already be harmonized.

The spatial AI/ML workflows are deployable adaptations informed by the public Microsoft AI/ML curricula and deep-learning course repositories; they do not claim that a curriculum repository is a statistical package. The image-regression screens cite the related SV-NN, ST-CAR, and SBLF publications and state explicitly when the interactive implementation is a faster screening model rather than the full Bayesian posterior engine.

No uploaded values or analysis results are sent to a third-party AI interpretation service. Automated result explanation remains outside the application by policy.

## Install from the Zenodo release

The Zenodo/GitHub source ZIP is the complete distribution of record. It
contains the Shiny application, installable R package API, documentation,
tests, five vignettes, four compact simulated datasets, and the curated public
source bundle under `vendor/`.

### 1. Prerequisites

- R 4.5.x or a compatible current R release;
- a working compiler toolchain for source packages: Rtools on Windows, Xcode
  Command Line Tools on macOS, or standard build tools plus GSL development
  libraries on Linux;
- Git is optional when installing from Zenodo because the exact engine sources
  are already included;
- sufficient memory and runtime for Bayesian, Stan, image, and 3D workflows.

On Windows, install the Rtools generation matching the first two digits of R
(for example, R 4.5.x requires Rtools 4.5). The complete installer detects a
standard `C:\\rtools45` installation and adds its compiler directories to the
active R process.

Download the newest ZIP from the [Zenodo concept
record](https://doi.org/10.5281/zenodo.21763606), extract it, and make the
extracted directory the working directory in R or a terminal. Do not run the
installer from inside the ZIP preview.

### 2. Install the complete executable Workbench

From a terminal in the extracted release directory:

```text
Rscript install-complete.R
```

The installer obtains missing CRAN dependencies and installs the bundled local
copies of `ISPAT`, `GPGHS`, and `ISPAT3D`; it does not clone those method
repositories. It then installs the `SpatialMethodsWorkbench` package and its
vignettes. Use `--skip-vignettes` when LaTeX/Pandoc is unavailable:

```text
Rscript install-complete.R --skip-vignettes
```

The release includes source snapshots for the related SV-NN, STCAR, and SBLF
workflows. To compile the related R packages as well, where supported:

```text
Rscript install-complete.R --include-related
```

SBLF's upstream package is Linux-only. SV-NN is a Python research-code
snapshot and has its own Python requirements under
`vendor/reference-workflows/SV-NN`; it is not silently used by the built-in
ridge screen. The interactive image workflows remain explicitly labelled as
screening implementations rather than the full SV-NN, STCAR, or SBLF
posteriors.

### 3. Verify the installation

```text
Rscript verify-installation.R
Rscript tests/testthat.R
```

The first command checks the Workbench plus all three direct advanced engines
and prints the 12 registered paths. The second runs the automated package and
Shiny-server tests.

### 4. Launch the local application

In R, with the extracted release as the working directory:

```r
shiny::runApp(".", host = "127.0.0.1", port = 3838)
```

Open `http://127.0.0.1:3838`. Binding to `127.0.0.1` keeps the service on the
local machine. Do not upload protected data to an unsecured public service.

### What “complete distribution” means

Users do not need to download the ISPAT, GP–GHS, or ISPAT-3D repositories
separately. License-safe public source connected to the AI/ML, neural,
SV-NN, STCAR, and SBLF provenance is also archived under `vendor/`; exact
revisions and retained scope are recorded in
[`vendor/manifest.json`](vendor/manifest.json). Large curriculum translations,
media, public demonstration datasets, Git history, and platform-specific
compiled objects are intentionally omitted because they do not execute a
Workbench analysis.

Like normal source-distributed R software, the installer may download missing
CRAN dependencies. A source ZIP cannot contain a portable compiled R library
for every operating system. The Docker build below provides the most isolated
installation path and installs the bundled direct engines inside the image.

## Minimal built-in installation

For only the fast built-in workflows, without the advanced Bayesian engines:

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "jsonlite", "markdown", "zip"))
shiny::runApp(".")
```

Or from a terminal with this repository as the working directory:

```text
Rscript run-local.R
```

Open `http://127.0.0.1:3838`.

## Install the package API from GitHub and open the vignettes

```r
install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github(
  "sagnikbhadury/SpatialMethodsWorkbench",
  build_vignettes = TRUE,
  dependencies = NA
)

vignette(package = "SpatialMethodsWorkbench")
vignette("app-walkthrough", package = "SpatialMethodsWorkbench")
```

The installed package exports the four deterministic data generators, validation and readiness functions, all 12 workflow runners, the common `run_analysis()` dispatcher, citation helpers, and reproducibility-bundle functions. The Shiny app and package call the same analysis implementation.

## Install public advanced engines directly from GitHub

This alternative is for a GitHub checkout. Zenodo users should use the bundled
installer above.

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

### Method and source attribution

Cite the Workbench **and** the sources applicable to the workflow you use. The table below is a compact map; [REFERENCES.md](REFERENCES.md) contains the complete bibliography and public repository provenance.

| Workbench workflow | Primary method/software sources to cite |
|---|---|
| Spatial exploration and autocorrelation | Moran (1950), [doi:10.1093/biomet/37.1-2.17](https://doi.org/10.1093/biomet/37.1-2.17) |
| Region-aware conditional network | Hoerl & Kennard (1970), [doi:10.1080/00401706.1970.10488634](https://doi.org/10.1080/00401706.1970.10488634), for ridge regularization; report this Workbench implementation as an exploratory precision/partial-correlation screen |
| ISPat 2D Bayesian network | Bhadury et al. (2026), [doi:10.1038/s41598-026-35341-8](https://doi.org/10.1038/s41598-026-35341-8), plus the public [`ISPAT`](https://github.com/sagnikbhadury/ISPAT) package |
| GP–GHS network | Bhadury, Gaskins & Rao (2026), [doi:10.64898/2026.04.01.715977](https://doi.org/10.64898/2026.04.01.715977), plus the public [`GP-GHS`](https://github.com/sagnikbhadury/GP-GHS) package |
| ISPat 3D volumetric network | Public [`ISPAT-3D`](https://github.com/sagnikbhadury/ISPAT-3D) software; the Workbench does not cite or disclose an unpublished manuscript |
| Linear/logistic and shallow-neural prediction | Roberts et al. (2017), [doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881), for structured validation; Venables & Ripley (2002), [doi:10.1007/978-0-387-21706-2](https://doi.org/10.1007/978-0-387-21706-2), for `nnet` |
| Spatial clustering | Chavent et al. (2018), [doi:10.1007/s00180-018-0791-1](https://doi.org/10.1007/s00180-018-0791-1), as related spatial-clustering methodology; the Workbench path is a distinct weighted-coordinate k-means adaptation |
| Wide-image scalar regression screen | Hoerl & Kennard (1970) for the implemented ridge screen; Wu, Wu & Kang (2025), [JMLR 26(116)](https://www.jmlr.org/papers/v26/22-0246.html), and Xu & Kang (2025), [ICLR 2025](https://proceedings.iclr.cc/paper_files/paper/2025/hash/f418594e90047a10f4c158f70d6701cc-Abstract-Conference.html), as related full spatial Bayesian methods |
| Latent image-to-image regression screen | Guo, Kang & Johnson (2022), [doi:10.1111/biom.13420](https://doi.org/10.1111/biom.13420), as the related full SBLF method; the Workbench path is a separate PCA/score-regression screen |
| Mediation | Imai, Keele & Tingley (2010), [doi:10.1037/a0020761](https://doi.org/10.1037/a0020761) |
| Shape PCA | Dryden & Mardia (2016), [doi:10.1002/9781119072492](https://doi.org/10.1002/9781119072492) |

### Public fork provenance

The following public forks informed workflows or related-method documentation. Attribution goes to the original upstream project as well as the fork used during development. A curriculum or related-method fork is **not** represented as an exact statistical engine unless the Workbench directly calls that package.

| Role in the Workbench | Original upstream source | Public fork used during development |
|---|---|---|
| AI workflow design and educational structure | [microsoft/AI-For-Beginners](https://github.com/microsoft/AI-For-Beginners) | [sagnikbhadury/AI-For-Beginners](https://github.com/sagnikbhadury/AI-For-Beginners) |
| ML workflow design and educational structure | [microsoft/ML-For-Beginners](https://github.com/microsoft/ML-For-Beginners) | [sagnikbhadury/ML-For-Beginners](https://github.com/sagnikbhadury/ML-For-Beginners) |
| Shallow-neural implementation background | [mikexcohen/DeepUnderstandingOfDeepLearning](https://github.com/mikexcohen/DeepUnderstandingOfDeepLearning) | [sagnikbhadury/DeepUnderstandingOfDeepLearning](https://github.com/sagnikbhadury/DeepUnderstandingOfDeepLearning) |
| Related scalar-on-image SV-NN method and software | [benwu233/SV-NN](https://github.com/benwu233/SV-NN) | [sagnikbhadury/SV-NN](https://github.com/sagnikbhadury/SV-NN) |
| Related soft-thresholded CAR image-regression method and software | [yuliangxu/STCAR](https://github.com/yuliangxu/STCAR) | [sagnikbhadury/STCAR](https://github.com/sagnikbhadury/STCAR) |
| Related spatial Bayesian latent-factor method and software | [umich-biostatistics/SBLF](https://github.com/umich-biostatistics/SBLF) | [sagnikbhadury/SBLF](https://github.com/sagnikbhadury/SBLF) |

Only the public sources above are claimed as informing the shipped Workbench workflows. Other repositories are not implied dependencies merely because they exist in the account.

Citation acknowledgement is a strong research-norm and provenance mechanism, not a technical guarantee about a later publication's bibliography. A tagged release can be archived with Zenodo to add a persistent DOI.

- Cite the version-specific DOI shown on the downloaded Zenodo record.
- Use the [concept DOI 10.5281/zenodo.21763606](https://doi.org/10.5281/zenodo.21763606) when referring to the Workbench across all versions.

## License and software publication

The Workbench source code is published under the [MIT License](LICENSE.md). The MIT License allows reuse, modification, redistribution, sublicensing, and commercial use, provided the copyright and permission notice are retained. The software is supplied without warranty.

The archived software releases are citable software publications under the [Zenodo concept DOI 10.5281/zenodo.21763606](https://doi.org/10.5281/zenodo.21763606). The MIT License governs reuse of this repository's code; it does **not** replace scholarly attribution, method-specific citations, data-use obligations, or licenses attached to external packages. Bundled third-party source remains under the license retained with that source. The application requires citation acknowledgement and writes the applicable references into every result bundle, but no open-source license can technically guarantee what a later manuscript includes.

## Design evidence

The product and literature rationale, scope boundaries, and deployment recommendation are documented in [docs/LITERATURE_REVIEW.md](docs/LITERATURE_REVIEW.md).

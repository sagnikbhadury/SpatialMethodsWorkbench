# Spatial Methods Workbench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21763606.svg)](https://doi.org/10.5281/zenodo.21763606)
[![Version](https://img.shields.io/badge/release-v0.2.1-176b68)](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/releases/tag/v0.2.1)
[![License: MIT](https://img.shields.io/badge/License-MIT-d46245.svg)](LICENSE)
[![R tests](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/actions/workflows/test.yml/badge.svg)](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/actions/workflows/test.yml)

A functional Shiny application for guided analysis of spatial and structured biomedical data. Users upload a CSV, map its columns, see which analyses are compatible, configure a method, run it, and download a reproducibility bundle containing results, a figure, settings, citations, and the R result object.

**[Launch the live application](https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/)** · **[Read the complete usage and installation guide](docs/USER_GUIDE.md)** · **[Open the function-by-function reference](docs/PIPELINE_REFERENCE.md)** · **[Cite v0.2.1](https://doi.org/10.5281/zenodo.21764196)**

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

The application also includes synthetic demonstration data, automatic column suggestions, structural validation, compatibility gating, documented controls, and downloadable result bundles. The advanced engines are available only when their corresponding public R package is installed; the live v0.2.1 deployment includes all three.

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

- Cite the current `v0.2.1` release using [DOI 10.5281/zenodo.21764196](https://doi.org/10.5281/zenodo.21764196).
- Use the [concept DOI 10.5281/zenodo.21763606](https://doi.org/10.5281/zenodo.21763606) when referring to the workbench across all versions.

## License and software publication

The Workbench source code is published under the [MIT License](LICENSE). The MIT License allows reuse, modification, redistribution, sublicensing, and commercial use, provided the copyright and permission notice are retained. The software is supplied without warranty.

The archived software release is a citable software publication on Zenodo: [v0.2.1, DOI 10.5281/zenodo.21764196](https://doi.org/10.5281/zenodo.21764196). The MIT License governs reuse of this repository's code; it does **not** replace scholarly attribution, method-specific citations, data-use obligations, or licenses attached to external packages. The application requires citation acknowledgement and writes the applicable references into every result bundle, but no open-source license can technically guarantee what a later manuscript includes.

## Design evidence

The product and literature rationale, scope boundaries, and deployment recommendation are documented in [docs/LITERATURE_REVIEW.md](docs/LITERATURE_REVIEW.md).

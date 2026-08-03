# Development session log

## 2026-08-02 — v0.2.1 public toolkit and documentation release

### Completed

- Built and published a functional R Shiny Workbench around public methods only.
- Added 12 guided executable workflows, input readiness checks, citation acknowledgement, numerical outputs, figures, and reproducibility bundles.
- Integrated the public ISPAT, GP-GHS, and ISPAT-3D engines; fixed and released the public ISPAT-3D `%dopar%` namespace import as ISPAT-3D v0.1.1.
- Installed and smoke-tested the advanced R dependencies locally.
- Enforced the interpretation boundary: no automated result explanation and no LLM prompt in exports.
- Added the complete `docs/USER_GUIDE.md` covering restrictions, real-data preparation, installation, Docker, deployment, privacy, citation, and collaboration.
- Added the 669-line `docs/PIPELINE_REFERENCE.md` covering every workflow/function, input contract, controls, calls, output objects, failures, and citations.
- Rendered both manuals inside the live Shiny application.
- Published GitHub release v0.2.1 and Zenodo DOI `10.5281/zenodo.21764196`.
- Corrected the citation regression contract to contain `SpatialMethodsWorkbench`, version 0.2.1, and the exact DOI.
- Verified 31 automated checks locally with no failures, warnings, or skips.
- Prepared the academic website showcase, launch links, installation guide, function-reference link, v0.2.1 release, and DOI.

### Public endpoints

- Live Shiny app: <https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/>
- Source: <https://github.com/sagnikbhadury/SpatialMethodsWorkbench>
- Release: <https://github.com/sagnikbhadury/SpatialMethodsWorkbench/releases/tag/v0.2.1>
- DOI: <https://doi.org/10.5281/zenodo.21764196>

### Resume prompt for a future coding agent

> Open `LOCAL_DEVELOPMENT_HANDOFF.md`, `SESSION_LOG.md`, `docs/USER_GUIDE.md`, and `docs/PIPELINE_REFERENCE.md` completely. Inspect `git status`, recent commits, and the live endpoints before editing. Add only already-public modules and preserve all privacy, disclosure, citation, and no-interpretation boundaries. Run the complete test suite and update this log after material changes.

### Next safe extension point

New public modules should enter through `analysis_registry()`, `run_analysis()`, a standard result object, Shiny controls, citations, tests, and both manuals. Do not integrate a module until its public status, license, truthful computational scope, and input/output contract have been verified.

### Final publication follow-up

- Expanded the public repository README with a structured catalog of all 12 analyses, their required data, principal outputs, and exact-engine versus screening status.
- Added explicit MIT License and software-publication guidance. The MIT terms govern code reuse; they do not replace scholarly citation or external-package licenses.
- Embedded the live Shiny application directly in the academic website page with a full-screen fallback.
- Updated the website catalog to state the input and output of every analysis and linked the latest v0.2.1 installation guide and complete function reference.
- Synchronized this local development clone with public source commit `a750fcd`; the local handoff commit remains intentionally local and is not pushed to the public repository.

## 2026-08-02 — narrated training and documentation package

### Completed

- Created a detailed 25-slide public usage tutorial covering all 12 workflows, data contracts, preprocessing, readiness, controls, output auditing, local installation, advanced engines, privacy, citation, and reproducibility.
- Created a separate 29-slide local interpretation course grounded in the public method literature.
- The interpretation course distinguishes the quantities used by each method: Moran's I and permutation reference; covariance, precision, and partial correlations; ISPat spatial intensity, shared/zone matrices, edge signs, sample proportions, and RV recovery; GP-GHS spatial edge fields and grouped shrinkage; ISPAT3D volumetric trends and matrices; predictive errors and calibration; clustering geometry; image coefficients or latent factors; mediation effects; and shape PCA scores/loadings.
- Preserved the public no-interpretation boundary. The public tutorial explains operation, statistical purpose, controls, and reproducibility but does not assign biological, clinical, or causal meaning to a user's results.
- Exported both courses as narrated 1080p MP4 files with editable PowerPoint sources, SRT captions, chapter timestamps, and thumbnails.
- Verified media durations with Windows Media Player metadata: public tutorial 20:20; local interpretation course 21:56.
- Added `media/youtube/PUBLIC_USAGE_VIDEO_GUIDE.md`, `METHOD_INTERPRETATION_GUIDE.md`, `LITERATURE_CONTEXT.md`, `YOUTUBE_UPLOAD.md`, and rebuild instructions.
- Kept all videos and interpretation training local. No private repository concept, working manuscript, LLM work, SAE work, or protected data was included.

### Video source of truth

Edit `media/youtube/public_usage_slides.json` or `media/youtube/interpretation_slides.json`, then run `media/youtube/build_videos.ps1`. Generated artifacts are under `media/youtube/output/` and narration WAV files under `media/youtube/audio/`.

## 2026-08-02 — complete method and upstream attribution

### Completed

- Added method-specific citations for all 12 Workbench workflows to the runtime registry and reproducibility bundles.
- Added `REFERENCES.md` as the complete public bibliography and provenance record.
- Added explicit upstream-project and development-fork credit for AI for Beginners, ML for Beginners, Deep Understanding of Deep Learning, SV-NN, STCAR, and SBLF.
- Distinguished exact computational engines from related publications, curriculum sources, and screening implementations so the public claims remain accurate.
- Updated the README with analysis-to-method attribution and public-fork provenance tables.
- Updated and published the existing Zenodo v0.2.1 record with 19 references, eight keywords, the University of Michigan affiliation, the expanded software description, and the method-specific citation requirement.
- Verified that Zenodo retained DOI `10.5281/zenodo.21764196`, version `v0.2.1`, the MIT license, open access, repository relation, and the original archived ZIP checksum.
- Verified the public Zenodo page exposes the updated description, required attribution statement, public upstream/fork acknowledgments, and the explicit private/working-manuscript exclusion.
- GitHub Actions passed all 36 checks for public source commit `87cd1cb`.

## 2026-08-02 — reproducible package and app vignettes

### Completed

- Converted the analysis layer from an app-only project layout into an installable R package while preserving the Shiny source application.
- Added five executable R Markdown vignettes covering the app walkthrough, all 12 analysis workflows, result-bundle auditing, and public fork provenance.
- Added deterministic generators for spatial, 3D/volumetric, aligned image-feature, and ordered-contour inputs.
- Added four upload-ready CSV files under `inst/extdata/`: 72 spatial rows, 72 volumetric rows, 48 image subjects, and six contours with 16 points each.
- Used fast vignette settings: 19 permutations, 25 bootstrap replicates, compact neural networks, and small feature matrices. Advanced Bayesian calls remain opt-in and use clearly labeled smoke-test settings.
- Added an in-app **Reproducible vignettes** tab and a repository vignette index with direct CSV downloads and mapping instructions.
- Documented how every public fork is represented, distinguishing direct R engine adapters, executable built-in adaptations, and related full methods cited alongside screening workflows.
- Exported the public dataset generators, validation/readiness utilities, all 12 workflow runners, dispatcher, citation helpers, and reproducibility-bundle functions.
- Added package manuals, a package-safe MIT license stub plus full `LICENSE.md`, package build exclusions, and continuous-integration vignette building.
- Verified 51 source tests with no failures, warnings, or skips.
- Built all five HTML vignettes and completed `R CMD check --no-manual` with `Status: OK`.
- Added a Shiny-specific deployment exclusion file, reducing the hosted bundle from approximately 526 MB to 95 KB while retaining all runtime code and manuals.
- Redeployed the public shinyapps.io instance and verified in a browser that the vignette tab, fork-representation table, no-interpretation boundary, and all four CSV links are visible.

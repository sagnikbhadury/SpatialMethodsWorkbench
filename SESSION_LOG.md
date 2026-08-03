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

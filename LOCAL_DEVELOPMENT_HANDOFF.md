# Local development handoff

This folder is a full Git clone of the public Workbench. Use it as the durable starting point for future modules in Codex or Claude Code.

## Non-negotiable boundaries

1. Use public methods and public repositories only.
2. Do not expose any non-public source, research direction, result, or other confidential material.
3. Never send uploaded research data to an external interpretation service.
4. Do not generate automated scientific interpretation. Display statistical outputs, assumptions, warnings, and provenance only.
5. Preserve the citation acknowledgement and include both Workbench and method-specific citations in every export.
6. Preserve privacy warnings and the prohibition on PHI in the public service.
7. Describe adaptations honestly; do not label an interactive screen as an exact implementation of a more complex published model.

## Architecture

- `app.R`: Shiny user interface, dynamic controls, readiness display, server orchestration, downloads.
- `R/registry.R`: workflow registry, compatibility rules, Workbench citation, method citations.
- `R/methods_core.R`: built-in executable methods.
- `R/methods_advanced.R`: adapters for public advanced packages plus central `run_analysis()` dispatch.
- `R/utils.R`: upload parsing, validation, demo data, safe column handling.
- `R/reporting.R`: manifest and reproducibility ZIP bundle.
- `docs/PIPELINE_REFERENCE.md`: public function-by-function reference.
- `docs/USER_GUIDE.md`: operational and installation guide.
- `tests/testthat/`: behavioral tests.

The central programmatic interface is:

```r
result <- run_analysis(
  id = "workflow_id",
  data = analysis_data,
  mapping = list(x = "x", y = "y", sample = "sample_id"),
  features = c("feature_1", "feature_2"),
  params = list()
)
```

## Add a module

1. Confirm the method, code, license, and citations are already public.
2. Define the exact input contract and distinguish required from optional roles.
3. Implement `run_<module>()` in `R/methods_core.R` or a public-package adapter in `R/methods_advanced.R`.
4. Return the standard result object documented in `docs/PIPELINE_REFERENCE.md`; do not return prose interpretation.
5. Add a registry entry in `analysis_registry()` with ID, label, description, requirements, and engine type.
6. Add a dispatch branch to `run_analysis()`.
7. Add method controls and parameter construction in `app.R`.
8. Add the method's publications to `method_citations()`.
9. Add readiness, numerical-output, export, failure, and citation tests.
10. Document the function, data layout, assumptions, output schema, limits, and example in both manuals.
11. Run the full suite and a real adapter smoke test when an external engine is involved.
12. Bump the version, update `CITATION.cff`, create a GitHub release, obtain the Zenodo DOI, then update the live app and website.

## Local verification

From this directory in PowerShell:

```powershell
& 'C:\Users\bhadury\AppData\Local\Programs\R\R-4.5.2\bin\Rscript.exe' tests/testthat.R
```

Expected for v0.2.1: `FAIL 0`, `WARN 0`, `SKIP 0`, `PASS 31`.

Run locally:

```powershell
& 'C:\Users\bhadury\AppData\Local\Programs\R\R-4.5.2\bin\Rscript.exe' run-local.R
```

The public advanced engines are installed from:

- `sagnikbhadury/ISPAT`
- `sagnikbhadury/GPGHS`
- `sagnikbhadury/ISPAT-3D`

See `docs/USER_GUIDE.md` for the complete dependency and Docker instructions.

## Deployment and release checklist

1. Confirm `git status` and inspect every diff so unrelated/private files cannot enter the release.
2. Run `tests/testthat.R` and smoke-test external adapters.
3. Deploy only the explicit public file list: `app.R`, `DESCRIPTION`, `R/`, `docs/USER_GUIDE.md`, `docs/PIPELINE_REFERENCE.md`, and `www/`.
4. Verify the live HTTP response and the presence of Usage, Pipeline reference, citation gate, and collaboration controls.
5. Create a versioned GitHub release only after tests pass.
6. Confirm the Zenodo record through its public API and propagate the exact DOI to citation metadata, tests, manuals, application, and website.
7. Verify the GitHub Actions run and GitHub Pages build.

Never store shinyapps.io tokens, GitHub tokens, passwords, uploaded datasets, temporary deployment URLs, or other secrets in this repository or in session logs.

# Literature and product review

## Recommendation

Build a guided spatial-methods workbench, not another monolithic spatial-omics platform. The useful gap is a low-friction path from a researcher's own table to a defensible analysis: validate the data, expose only compatible methods, explain assumptions, run a selected method, and export complete provenance and citations.

This direction matches Sagnik Bhadury's public research identity—statistical methods for structure in cancer data—while keeping unpublished manuscripts, private repositories, and private method ideas outside the application. The public interface does not enumerate source repositories or expose implementation provenance beyond the software and scholarly citations needed for responsible reuse.

## Evidence from existing tools

The ecosystem already has strong broad platforms:

- [Giotto](https://doi.org/10.1186/s13059-021-02286-2) supplies an extensive spatial-expression analysis and visualization toolbox.
- [Squidpy](https://doi.org/10.1038/s41592-021-01358-2) provides scalable spatial molecular analysis in Python.
- [Vitessce](https://doi.org/10.1038/s41592-024-02436-x) demonstrates the value of browser-based coordinated views for multimodal and spatial data, with preprocessing handled outside the browser.
- [histoCAT](https://doi.org/10.1038/nmeth.4391), [SPIAT](https://doi.org/10.1038/s41467-023-37822-0), and the [imcRtools workflow](https://doi.org/10.1038/s41596-023-00881-0) cover interactive and reproducible analysis of highly multiplexed tissue imaging.
- [spatialLIBD](https://pubmed.ncbi.nlm.nih.gov/35689177/) shows how a Shiny application can make a specialist spatial workflow accessible without requiring users to program.
- [cyjShiny](https://doi.org/10.1371/journal.pone.0285339) illustrates interactive network exploration inside Shiny.

These tools make a generic viewer or all-purpose toolkit a crowded product category. The distinctive opportunity is method selection, structured guidance, and network inference for spatial cancer data.

## Scientific center of gravity

[ISPat](https://doi.org/10.1038/s41598-026-35341-8) is the natural flagship method: it estimates shared and region-specific conditional-association networks from multiplexed spatial imaging data. The application therefore treats networks as a first-class result while explicitly warning that an edge is not automatically a physical interaction, molecular mechanism, or causal effect.

The 3D pipeline belongs in a separate, higher-complexity tier. Serial-section or volumetric analyses require x/y/z coordinates, volume identifiers, region labels, adequate features, and preprocessing/alignment checks. The interface gates this engine until those inputs and the server-side package are available; it never silently substitutes a simpler method under the same name.

## Product architecture

The resulting application has four layers:

1. **Intake and validation.** Upload a CSV or use synthetic data; map coordinates, regions, identifiers, outcomes, exposures, mediators, and features.
2. **Compatibility guidance.** Show which workflows are ready and exactly what is missing from blocked workflows.
3. **Analysis.** Offer fast exploratory modules and opt-in public advanced engines. Expensive Bayesian or 3D models remain server-side and are gated by installed dependencies.
4. **Responsible export.** Download results, plots, parameters, an R result object, runtime metadata, required citations, and a collaboration contact—without automated result interpretation.

Public AI and machine-learning educational material informs the guided prediction design, but the deployed app does not expose unpublished LLM/SAE research or upload raw study data to a language model. It does not generate explanations or export prompts for explaining results. Researchers are invited to contact Sagnik Bhadury for interpretation, methodological guidance, or collaboration.

## Deployment conclusion

GitHub should be the public source of record, issue tracker, citation entry point, and release history. GitHub Pages cannot execute an R Shiny process. The runnable service should therefore use shinyapps.io for a lightweight initial release or Posit Connect/a container host for authentication, job control, larger uploads, and compiled Bayesian engines. The academic website should link to both the deployed service and its versioned source.

## Scope and disclosure policy

- Only already-public methods and intentionally reimplemented generic workflows may appear.
- No working manuscript, private repository, private benchmark, experimental result, or unpublished method name may be copied into documentation, examples, UI text, logs, or images.
- Synthetic demonstration data are the default. Public deployments must prohibit protected health information and configure host-level retention and logging accordingly.
- Forked educational repositories are design references, not evidence that every tutorial is a production statistical engine. A workflow is exposed only when its input contract, statistical meaning, license, and runtime behavior can be made explicit and tested.

Review date: 2026-08-02.

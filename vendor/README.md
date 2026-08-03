# Curated public-source bundle

This directory is part of the complete Spatial Methods Workbench source
distribution. It prevents users from having to clone the method repositories
separately and records the exact public revisions used for release testing.

## Directly executed engines

The complete package sources in `engines/` are called directly by Workbench
analysis paths when installed:

- `ISPAT` — ISPat 2D Bayesian networks;
- `GP-GHS` — GP–GHS spatially varying networks (R package name `GPGHS`);
- `ISPAT-3D` — volumetric ISPat networks (R package name `ISPAT3D`).

Run `Rscript install-complete.R` from the repository root to install their
declared R dependencies and then install these local copies. The installer
does not clone the three repositories.

## Reference and related-method source

`reference-workflows/` preserves the license-safe executable source most
closely connected to the Workbench's built-in adaptations and screening
workflows:

- selected neural-network examples from AI for Beginners and A Deep
  Understanding of Deep Learning;
- selected regression and k-means examples from ML for Beginners;
- SV-NN Python source without its large MNIST/Fashion-MNIST demonstration
  datasets;
- the STCAR R/C++ package source without platform-specific compiled objects;
- the SBLF R/C/C++ package source without platform-specific compiled objects.

These sources are included for reproducibility and attribution. They are not
silently substituted for a Workbench method. The current interactive
wide-image and image-to-image paths are explicitly labelled fast screening
workflows rather than the full SV-NN, STCAR, or SBLF posterior algorithms.
STCAR and SBLF can be compiled separately with
`Rscript install-complete.R --include-related`, subject to their system and
platform requirements.

Large curriculum translations, media, Git history, and demonstration datasets
are omitted because they do not execute any Workbench analysis. Their precise
public revisions and the scope retained here are recorded in `manifest.json`.

## Licensing and confidentiality

Every snapshot retains its upstream license or license declaration. The
Workbench itself remains MIT licensed; third-party source remains governed by
its own license. STCAR declares `GPL (>= 2)` in `DESCRIPTION`, and `COPYING`
is included with its curated source.

Only public, release-scoped source is included in this distribution.

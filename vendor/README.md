# Curated public source bundle

This directory makes the Spatial Methods Workbench release a complete,
self-contained public distribution. Every entry is a snapshot of an already
public repository or CRAN package, recorded in `manifest.json` with its
repository, exact commit, package version, license, retained scope, and role.
Nothing here is private research code, and nothing here is executed silently in
place of a Workbench method. The snapshots are the revisions used for release
testing.

## Directly executed engines

The complete package sources in `engines/` are called directly by Workbench
analysis paths when installed:

- `ISPAT-3D` — the CRAN package `ISPAT3D` 0.3.1 (exact CRAN source tarball
  contents). Its `ispat3d_fit_2d()` powers the ISPat 2D path and its
  `ispat3d_fit()` powers the ISPat 3D path; both share the Vecchia GP
  adjustment and the shared-plus-zone covariance likelihood. The compiled
  dependency `gpboost` is installed from CRAN by the installer.
- `GP-GHS` — GP–GHS spatially varying networks (R package name `GPGHS`,
  version 1.0.1, pure R). Its `gpghs_graph()` powers the GP–GHS path. The
  snapshot includes the `analysis/` scripts of the accompanying article.

Run `Rscript install-complete.R` from the repository root to install the
engine dependency and then these local copies. The installer does not clone
the two repositories.

The standalone `ISPAT` package (the historical CAVI/SVI implementation with
Stan) was bundled through v0.3.1 and is retired: no Workbench path calls it,
so it is listed under `retired_sources` in `manifest.json` and is not shipped.

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

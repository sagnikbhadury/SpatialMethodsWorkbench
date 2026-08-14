# Ready-to-analyze synthetic data

Start with `advisor_complete_demo.csv`. It contains 128 fully synthetic rows,
uses fixed seed 2026, contains no patient information, and has been executed
successfully through all nine built-in Spatial Methods Workbench analyses.

Recommended mappings for `advisor_complete_demo.csv`:

| Workbench role | Column(s) |
|---|---|
| X, Y, Z | `x`, `y`, `z` |
| Region / zone | `zone` |
| Numeric outcome | `outcome` |
| Binary outcome | `disease` |
| Exposure, mediator | `exposure`, `mediator` |
| Spatial/network features | `epithelial`, `APC`, `CTL`, `Treg`, `T_helper` |
| Wide-image features | `pixel_1` through `pixel_6` |
| Image-to-image features | all `input__*` and `output__*` columns |
| Shape ID | `contour_id` |
| 3D volume ID | `volume_id` |

Use `contour_id` as the ID for shape PCA and `volume_id` for ISPAT-3D.
ISPAT, GP-GHS, and ISPAT-3D also require their corresponding engines to be
installed on the server.

Focused alternatives:

- `spatial_demo.csv`: spatial QC, networks, prediction, clustering, neural, and mediation
- `image_demo.csv`: scalar-image and image-to-image regression
- `contour_demo.csv`: shape PCA
- `spatial_3d_demo.csv`: volumetric spatial analysis

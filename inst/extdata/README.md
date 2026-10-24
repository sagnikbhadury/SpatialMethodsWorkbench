# Validated example data

`advisor_complete_demo.csv` is the recommended shareable example. It contains
128 fully synthetic observations generated with seed 2026. It has no patient or
protected health information. The automated test suite verifies the generator,
and the CSV has been executed successfully through every built-in Workbench
analysis.

Use these column mappings after upload:

| Role | Column |
|---|---|
| X / Y / Z coordinates | `x` / `y` / `z` |
| Region / zone | `zone` |
| Numeric outcome | `outcome` |
| Binary outcome (optional) | `disease` |
| Exposure / mediator | `exposure` / `mediator` |
| Spatial/network markers | `epithelial`, `APC`, `CTL`, `Treg`, `T_helper` |
| Wide-image features | `pixel_1` through `pixel_6` |
| Image-to-image features | all `input__*` and `output__*` columns |
| Shape ID | `contour_id` |
| 3D volume ID | `volume_id` |

For the shape workflow, change the ID mapping to `contour_id`. For the 3D
workflow, use `volume_id`; its two values are treated as section labels for the
ISPAT3D anchor balancing. Intensive ISPat 2D, GP-GHS, and ISPat 3D analyses also
require the `ISPAT3D` (>= 0.3.0) and `GPGHS` (>= 1.0.0) engines to be installed
on the server. Every zone must contain more rows than selected markers.

The specialized `spatial_demo.csv`, `image_demo.csv`, `contour_demo.csv`, and
`spatial_3d_demo.csv` files remain available for focused examples.

\name{demo_spatial_data}
\alias{demo_spatial_data}
\alias{demo_spatial_3d_data}
\alias{demo_image_data}
\alias{demo_contour_data}
\title{Generate small deterministic Workbench demonstration data}
\usage{
demo_spatial_data(n = 180L, seed = 2026L)

demo_spatial_3d_data(n = 72L, volumes = 2L, seed = 2026L)

demo_image_data(n = 48L, features = 8L, seed = 2026L)

demo_contour_data(contours = 6L, points = 16L, seed = 2026L)
}
\arguments{
  \item{n}{Number of simulated observations or subjects. Minimum safe values are enforced for workflows with sample-size checks.}
  \item{seed}{Integer random seed.}
  \item{volumes}{Number of simulated volume identifiers.}
  \item{features}{Number of features in each simulated image-feature block.}
  \item{contours}{Number of simulated ordered contours.}
  \item{points}{Number of ordered points per contour.}
}
\value{A data frame whose columns satisfy the corresponding vignette and Workbench input contracts.}
\description{
Generate compact synthetic spatial, volumetric, aligned-image, or contour data.
These datasets are intended for examples, tests, and software demonstrations,
not for biological or clinical conclusions.
}
\examples{
d <- demo_spatial_data(72, seed = 2026)
img <- demo_image_data(48, features = 8, seed = 2026)
dim(d)
dim(img)
}

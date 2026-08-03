\name{run_analysis}
\alias{analysis_registry}
\alias{coerce_mapping}
\alias{method_citations}
\alias{recommend_methods}
\alias{result_manifest}
\alias{run_analysis}
\alias{run_conditional_network}
\alias{run_gpghs}
\alias{run_image_to_image}
\alias{run_ispat}
\alias{run_ispat3d}
\alias{run_mediation}
\alias{run_neural_prediction}
\alias{run_scalar_image}
\alias{run_shape_pca}
\alias{run_spatial_clustering}
\alias{run_spatial_ml}
\alias{run_spatial_qc}
\alias{validate_dataset}
\alias{workbench_citation}
\alias{write_result_bundle}
\title{Run and audit Spatial Methods Workbench analyses}
\usage{
analysis_registry()

coerce_mapping(data, mapping)

validate_dataset(data, mapping = list())

recommend_methods(data, mapping, features)

run_analysis(id, data, mapping, features, params = list())

run_spatial_qc(data, mapping, features, params)

run_conditional_network(data, mapping, features, params)

run_ispat(data, mapping, features, params)

run_gpghs(data, mapping, features, params)

run_ispat3d(data, mapping, features, params)

run_spatial_ml(data, mapping, features, params)

run_spatial_clustering(data, mapping, features, params)

run_neural_prediction(data, mapping, features, params)

run_scalar_image(data, mapping, features, params)

run_image_to_image(data, mapping, features, params)

run_mediation(data, mapping, features, params)

run_shape_pca(data, mapping, features, params)

result_manifest(result)

write_result_bundle(result, destination)

workbench_citation()

method_citations(id)
}
\arguments{
  \item{data}{A data frame containing the mapped tabular input.}
  \item{mapping}{A named list mapping roles such as \code{x}, \code{y}, \code{z}, \code{zone}, \code{id}, \code{outcome}, \code{exposure}, and \code{mediator} to column names.}
  \item{features}{Character vector of selected numeric feature columns.}
  \item{id}{Registry identifier returned by \code{names(analysis_registry())}.}
  \item{params}{Named list of workflow controls.}
  \item{result}{A successful result returned by \code{run_analysis()}.}
  \item{destination}{Path of the ZIP archive to create.}
}
\value{
\code{analysis_registry()} and \code{recommend_methods()} return named lists;
\code{validate_dataset()} returns a validation summary; \code{run_analysis()}
returns a standard result list; \code{result_manifest()} returns a JSON-ready
list; \code{write_result_bundle()} writes a ZIP archive; citation functions
return character vectors.
}
\description{
Public programmatic interface shared with the Shiny application. Results contain
statistical output and provenance but do not provide automated scientific
interpretation.
}
\examples{
d <- demo_spatial_data(60, seed = 2026)
mapping <- list(x = "x", y = "y")
result <- run_analysis(
  "spatial_qc", d, mapping, "APC",
  list(feature = "APC", neighbors = 5, permutations = 19, seed = 2026)
)
result$table
}

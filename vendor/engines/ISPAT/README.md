# ISPAT: Informed Spatially aware Patterns for Multiplexed Immunofluorescence Data

An R package implementing ISPAT, a Bayesian framework that identifies shared
and region-specific cell-type interaction patterns in multiplexed
immunofluorescence (mIF) tissue imaging data, with multi-modal integration for
Poisson and Gaussian spatial mixed-effects models.

Methods and application: Bhadury, S. et al., [ISPAT: Informed Spatially Aware
Patterns for Multiplexed Immunofluorescence
data](https://www.nature.com/articles/s41598-026-35341-8), *Nature Scientific
Reports*.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/ISPAT")
```

ISPAT depends on [rstan](https://mc-stan.org/users/interfaces/rstan); make
sure a working C++ toolchain and rstan are installed first (see the
[RStan getting-started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)).

## Usage

```r
library(ISPAT)

mydata <- readRDS(file = paste0(savepath, "data/",
                   "Data - Analysis ready Marginally obtained Spatial intensities ",
                   "of Cells on Tissue Slide_", diseasetype, "_", mypatID, ".rds"),
                   refhook = NULL)
rawdata <- mydata %>% tidyr::pivot_wider(names_from = Cells, values_from = Intensities)
S <- rawdata[, c("X", "Y", "area_category")] %>% as.matrix()
Y <- rawdata %>% dplyr::select("APC", "CTL", "Epi", "THelper", "Treg") %>% t() %>% as.matrix()

saved <- ISPAT(Y, S, ncores = numcores, RefPrior = diag(1, 5, 5), use_ref = FALSE,
               Kernel = "Matern", sGLM_method = "MLE", VB_MSFA = TRUE, MSFA_method = "CAVI")
```

`ISPAT()` returns a shared cell-type interaction network common across all
tissue regions/clusters, plus a region-specific network for each cluster.

## Package layout

* `R/ISPAT.R` — main model-fitting function, `ISPAT()`.
* `R/msfa.R` — the multi-study factor analysis step, `cavi_msfa()` and `svi_msfa()`.
* `R/ext-loop.R` — internal per-marker spatial GLM worker.
* `inst/stan/` — bundled Stan programs for the per-marker spatial fit
  (Gaussian and Poisson spatial mixed models, RBF and Matern kernels).

See `?ISPAT` for full documentation of the main function.

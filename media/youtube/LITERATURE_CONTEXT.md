# Literature context for the tutorial videos

This note records the public literature used to motivate the two tutorial videos. It does not add methods to the Workbench and does not imply that a built-in screening workflow is the full published Bayesian model.

## Spatial tissue and spatial-omics inputs

- Feng et al. (2023), *Spatial analysis with SPIAT and spaSim to characterize and simulate tissue microenvironments*. Nature Communications. <https://doi.org/10.1038/s41467-023-37822-0>. Multiplex platforms including mIHC, MIBI, CODEX, and IMC commonly yield analysis-ready tables of cell coordinates, phenotypes, marker intensities, and morphological features after image processing.
- Palla et al. (2022), *Squidpy: a scalable framework for spatial omics analysis*. Nature Methods. <https://doi.org/10.1038/s41592-021-01358-2>. Spatial molecular data may be spot-, single-cell-, subcellular-, or pixel-resolved and combine coordinates, molecular measurements, neighborhood graphs, and images.

These sources support the Workbench's tabular intake after segmentation, registration, normalization, or sequencing preprocessing. The Workbench does not process raw microscopy, DICOM, whole-slide images, or sequencing reads.

## Method-specific foundations

### Spatial autocorrelation

- Moran (1950), *Notes on Continuous Stochastic Phenomena*. Biometrika. <https://doi.org/10.1093/biomet/37.1-2.17>.

Moran's I compares similarity in a measured feature with a user-defined spatial weight graph. The Workbench uses a symmetrized k-nearest-neighbor graph and a permutation reference distribution. The result depends on the coordinate system, neighbor definition, sampling design, and exchangeability of labels under permutation.

### Gaussian graphical models and partial correlations

- Drton and Maathuis (2017), *Structure Learning in Graphical Modeling*. Annual Review of Statistics and Its Application. Related pedagogic review: <https://arxiv.org/abs/1707.04345>.
- Li, Craig, and Bhadra (2019), *The Graphical Horseshoe Estimator for Inverse Covariance Matrices*. Journal of Computational and Graphical Statistics. Preprint: <https://arxiv.org/abs/1707.06661>.

For multivariate Gaussian data, off-diagonal precision-matrix entries encode conditional relationships. Partial correlations summarize sign and standardized strength after conditioning on the other included variables. They are not automatically physical interactions, signaling events, regulatory mechanisms, or causal effects.

### ISPat 2D

- Bhadury et al. (2026), *Informed spatially aware patterns for multiplexed immunofluorescence data*. Scientific Reports. <https://doi.org/10.1038/s41598-026-35341-8>.

ISPat models spatial cell-type intensity surfaces and constructs shared and region-specific conditional-dependence networks. The paper motivates multiplexed immunofluorescence, tumor microenvironment heterogeneity, biologically defined tissue zones, and optional prior information. It explicitly treats the network as conditional dependence and notes the need for experimental validation before causal or therapeutic claims.

### GP-GHS

- Bhadury, Gaskins, and Rao, *Spatially Varying Graphical Models for Cell-Cell Interaction Networks in Multiplexed Tissue Imaging*. Public full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC13060317/>.
- Public software: <https://github.com/sagnikbhadury/GP-GHS>.

GP-GHS uses spatial basis representations of Gaussian-process fields, grouped horseshoe shrinkage across basis coefficients, and nodewise regressions. Its primary target is evidence for whether an edge varies across space under the fitted model. Edge selection depends on prior, basis, sampling, symmetry, and threshold choices.

### ISPat 3D

- Public software and documentation: <https://github.com/sagnikbhadury/ISPAT-3D>.

The 3D implementation extends spatial adjustment and shared/zone-specific network estimation to registered serial sections or volumetric multiplexed imaging. Valid orientation, common coordinates, z spacing, zone definitions, section registration, volume identifiers, and adequate observations per zone are prerequisites.

### Prediction and validation

- Roberts et al. (2017), *Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure*. Ecography. <https://doi.org/10.1111/ecog.02881>.
- Venables and Ripley (2002), *Modern Applied Statistics with S*. Springer. <https://doi.org/10.1007/978-0-387-21706-2>.

Random train/test splitting can underestimate error when nearby observations or observations from the same patient, slide, site, or time point enter both sets. The Workbench's random held-out split is an exploratory demonstration. Confirmatory prediction requires patient-, slide-, site-, time-, or spatially blocked validation and preferably external validation.

### Spatially weighted clustering

- Chavent et al. (2018), *ClustGeo: an R package for hierarchical clustering with spatial constraints*. Computational Statistics. <https://doi.org/10.1007/s00180-018-0791-1>.

The Workbench does not implement ClustGeo; it uses standardized features plus weighted standardized coordinates in k-means. The literature motivates the general tradeoff between feature homogeneity and spatial coherence. Cluster identity, number, stability, and external biological validation remain separate questions.

### Ridge and scalar-on-image screening

- Hoerl and Kennard (1970), *Ridge Regression: Biased Estimation for Nonorthogonal Problems*. Technometrics. <https://doi.org/10.1080/00401706.1970.10488634>.
- Wu, Wu, and Kang (2025), *Bayesian Scalar-on-Image Regression with a Spatially Varying Single-layer Neural Network Prior*. JMLR. <https://www.jmlr.org/papers/v26/22-0246.html>.
- Xu and Kang (2025), *Bayesian Image Regression with Soft-thresholded Conditional Autoregressive Prior*. ICLR. <https://openreview.net/forum?id=rnL3OafDdw>.

The Workbench uses ridge regression as a computationally fast screen for aligned pixel, voxel, or derived-image columns. It is not the SV-NN or ST-CAR posterior method. Coefficients are conditional on preprocessing, alignment, scaling, selected features, and ridge penalty; spatially adjacent coefficients are not explicitly regularized together in the screen.

### Image-to-image screening

- Guo, Kang, and Johnson (2022), *A Spatial Bayesian Latent Factor Model for Image-on-Image Regression*. Biometrics. <https://doi.org/10.1111/biom.13420>.

The publication motivates using low-dimensional latent factors to connect high-dimensional aligned predictor and outcome images. The Workbench uses separate PCA representations and regression between scores as a fast screen. It is not the full spatial Bayesian latent-factor posterior model.

### Mediation

- Imai, Keele, and Tingley (2010), *A General Approach to Causal Mediation Analysis*. Psychological Methods. <https://doi.org/10.1037/a0020761>.

An indirect-effect estimate is not automatically causal. A causal discussion requires a defensible exposure, mediator, outcome, time order, consistency, positivity, model specification, and identification assumptions concerning exposure–mediator, exposure–outcome, and mediator–outcome confounding. Bootstrap uncertainty does not repair violated identification assumptions.

### Shape analysis

- Dryden and Mardia (2016), *Statistical Shape Analysis, with Applications in R*. Wiley. <https://doi.org/10.1002/9781119072492>.

Shape analysis separates geometric variation from nuisance translation and scale, and often rotation. The Workbench centers, size-normalizes, and resamples ordered contours before PCA, but does not automatically standardize rotation. Comparable boundary ordering, segmentation, anatomical correspondence, and rotation policy must therefore be established before scientific interpretation.

## Communication rule used in both videos

For every workflow, separate five statements:

1. What was observed in the input data.
2. What the fitted statistic or model estimates.
3. What uncertainty or validation was performed.
4. What biological hypothesis the result may motivate.
5. What additional design, replication, or experiment would be required to claim mechanism or causation.

The public usage video stops at statements 1–3. The local interpretation-training video teaches how to discuss statements 4–5 responsibly.

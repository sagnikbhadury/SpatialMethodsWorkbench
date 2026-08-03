# Spatial Methods Workbench: model interpretation guide

**Local consultation and training document.** This guide is for explaining model outputs to collaborators. It is not automated interpretation, clinical guidance, or a substitute for study-specific statistical review. Apply it only after verifying the study design, preprocessing, input contract, and fitted-model diagnostics.

## A five-layer interpretation discipline

For every result, speak in this order:

1. **Design:** observational unit, sampling hierarchy, preprocessing, coordinate system, and scientific question.
2. **Estimand:** the mathematical quantity the workflow estimates.
3. **Uncertainty/validation:** permutation reference, bootstrap interval, posterior rule, held-out error, or stability assessment.
4. **Biological hypothesis:** the tissue organization or mechanism that could be consistent with the result.
5. **Evidence boundary:** alternative explanations and the experiment, cohort, or validation required for a stronger claim.

Never jump directly from an image, edge, coefficient, cluster, or p-value to a mechanism.

## 1. Spatial exploration and Moran's I

### Displayed quantities

- Complete observations used.
- Number of nearest neighbors, `k`.
- Moran's I.
- Two-sided permutation p-value.
- Spatial point map of the selected feature.

### Statistical interpretation

The Workbench constructs a symmetrized k-nearest-neighbor weight matrix. Moran's I measures whether pairs connected by that graph tend to have similar centered feature values. Larger positive values indicate positive global spatial autocorrelation under that graph; negative values indicate neighboring dissimilarity; values near the randomization expectation indicate little global pattern under that definition of neighborhood.

The permutation p-value asks how extreme the observed statistic is relative to randomly reassigning feature values across the retained coordinates. It is only meaningful when that label-exchangeability scheme is defensible.

### Biological framing

For cell or spot data, positive spatial autocorrelation can be consistent with spatial compartments, density gradients, local tissue niches, or unremoved technical fields. Negative autocorrelation can be consistent with alternation or exclusion, but may also arise from segmentation and neighborhood choices.

### Synthetic demonstration

The bundled synthetic data yield Moran's I approximately `0.782` and permutation `p = 0.005` with `k = 8`. Say: “The simulation was intentionally generated with a global spatial pattern, and this configuration detects that pattern.” Do not attach a tumor, immune, or clinical mechanism to this synthetic value.

### Required sensitivity checks

- Repeat scientifically plausible `k` values or distance/radius definitions.
- Check coordinates, duplicate points, boundaries, and sampling density.
- Plot each patient/slide separately when observations are nested.
- Use a permutation scheme that preserves the sampling design.
- Distinguish global autocorrelation from local hotspots.

## 2. Region-aware conditional network

### Displayed quantities

- Ridge-regularized partial-correlation matrices.
- Edge source, target, sign, and weight above the display threshold.
- Number of networks and visible edges.
- Ridge and display thresholds.

### Statistical interpretation

The inverse covariance matrix is converted to standardized partial correlations. An edge summarizes the residual linear association between two selected variables after conditioning on the other selected variables. The ridge term stabilizes inversion; the threshold controls display rather than proving statistical significance.

### Biological framing

A positive partial correlation can be described as coordinated variation conditional on the included panel. A negative value can be described as conditional opposition. Either may motivate hypotheses about shared niches, exclusion, or coordinated states. Neither demonstrates physical contact, ligand–receptor signaling, regulation, or causation.

### Required checks

- Examine transformations, outliers, nonlinearity, and missingness.
- Confirm adequate observations relative to the number of features within each zone.
- Vary ridge and edge thresholds.
- Check whether omitted cell types, batch, patient, or tissue features could induce edges.
- Treat separately fitted zones as exploratory unless formal contrasts are performed.

## 3. ISPat 2D Bayesian networks

### Publication quantities

The ISPat publication interprets spatial intensity estimates, shared covariance/loading patterns, region-specific covariance/loading patterns, partial correlations, thresholded positive/negative edges, and the proportion of samples in which edge signs occur. Simulation recovery is evaluated using the RV coefficient between true and estimated matrices.

### Workbench quantities

- Shared and zone-specific conditional-association matrices.
- Thresholded edge tables and heat maps.
- Full public `ISPAT` fitted object in the R result bundle.

### Statistical interpretation

ISPat first adjusts cell-type intensity surfaces for spatial structure and then separates shared from zone-specific covariance patterns through multi-study factor structure. Partial correlations derived from these patterns describe conditional association after accounting for other included cell types.

### Biological framing

Shared edges summarize organization conserved across the modeled tissue zones. Zone-specific deviations identify conditional organization that differs in a prespecified spatial domain. The ISPat paper used epithelial, APC, CTL, T-helper, and Treg intensity patterns across tumor-density regions in pancreatic tissue. In a new study, biological labels and zones must be justified independently.

### Required checks

- Validate cell segmentation, phenotype classification, intensity estimation, kernel/bandwidth, and zones.
- Check observations per zone and patient/slide replication.
- Review convergence/variational settings and sensitivity to prior information.
- Examine the continuous matrices, not only thresholded graphs.
- Require independent-cohort and functional validation before mechanistic claims.

## 4. GP-GHS spatially varying networks

### Publication quantities

The GP-GHS publication focuses on spatially varying edge fields, grouped shrinkage of all spatial basis coefficients for an edge, posterior/selection rules for edge presence, sensitivity to symmetry and threshold choices, and network-recovery metrics in simulations.

### Workbench quantities

- Selected adjacency matrix and edge list.
- Number of selected edges.
- Full `GPGHS` fit, including spatial basis/posterior objects supplied by the package.

### Statistical interpretation

An edge is selected when the group of coefficients representing its spatial surface survives the configured horseshoe shrinkage and selection rule. The app's adjacency is a summary of model-based edge presence; the detailed spatial surface remains in the fitted object.

### Biological framing

The method is appropriate when a conditional relationship may change continuously across a tissue rather than only between discrete zones. A selected edge means the model found evidence for a spatially structured conditional relationship. It does not by itself establish sign, mechanism, direct contact, or functional signaling.

### Required checks

- Inspect the spatial edge surface, not only adjacency.
- Review basis dimension, smoothness, burn-in, effective samples, chains, symmetry rule, and threshold.
- Evaluate edge stability across images/patients and preprocessing variants.
- Check boundary and sampling-density artifacts.

## 5. ISPat 3D volumetric networks

### Displayed quantities

- Shared and zone-specific matrices for the first displayed volume.
- Thresholded partial-correlation edges.
- Full fitted objects for all volumes in the result bundle.

### Statistical interpretation

The public package estimates marker-specific 3D spatial trends, removes modeled spatial structure, and uses multi-study factor analysis to recover shared and zone-specific covariance patterns. Partial correlations summarize conditional relationships derived from those patterns.

### Biological framing

The 3D path addresses tissue organization that cannot be represented by a single section. Differences across depth or zones can motivate hypotheses about volumetric niches and tissue architecture only if section alignment and z geometry are reliable.

### Required checks

- Verify orientation, registration, section order, section thickness/spacing, missing sections, masks, and units.
- Confirm that apparent z variation is not batch or staining drift.
- Inspect each volume, not only the first display.
- Test zone and registration sensitivity.

## 6. Spatial machine-learning prediction

### Displayed quantities

- Regression: held-out RMSE, MAE, R-squared, observed/predicted values, residuals.
- Binary classification: accuracy, probabilities, observed/predicted classes.
- Fitted linear or logistic model.

### Statistical interpretation

RMSE emphasizes larger errors; MAE is the average absolute error; held-out R-squared summarizes squared predictive correlation in this implementation. Classification probabilities are model-estimated risks/scores, while accuracy depends on the 0.5 cutoff and class prevalence.

### Biological framing

Prediction answers whether measured features and optional coordinates carry reproducible information about an outcome. It does not show that a predictor causes the outcome or that a coordinate feature is a biological mechanism.

### Required checks

- Replace random splitting with patient-, slide-, site-, time-, or spatially blocked validation.
- Report uncertainty and class-specific metrics where relevant.
- Compare with a clinically/scientifically meaningful baseline.
- Use external validation and calibration for translational claims.

## 7. Spatially weighted phenotype clustering

### Displayed quantities

- Cluster assignment for each observation.
- Standardized feature/coordinate centers in the fitted k-means object.
- Cluster count, spatial weight, and between-total sum-of-squares ratio.
- Spatial map of assignments.

### Statistical interpretation

K-means minimizes within-cluster squared distance in a constructed space containing standardized biological features and spatial coordinates multiplied by the chosen spatial weight. Increasing the spatial weight favors geographically coherent clusters. The between-total ratio is an internal compactness/separation summary, not external biological validation.

### Biological framing

Clusters may generate hypotheses about spatial phenotypes or niches. Call them “algorithmic clusters” until marker identity, stability, replication, and external biological annotations support stronger labels.

### Required checks

- Vary cluster count, seeds, feature set, transformations, and spatial weight.
- Measure stability and compare against nonspatial clustering.
- Check whether clusters follow tissue edges, batches, or segmentation artifacts.
- Validate against independent annotations or assays.

## 8. Shallow-neural prediction

### Displayed quantities

- Held-out RMSE, MAE, and R-squared for numeric outcomes.
- Held-out accuracy and probabilities for binary outcomes.
- Hidden-unit count, weight decay, and fitted `nnet` object.

### Statistical interpretation

The single-hidden-layer network estimates nonlinear prediction functions after standardizing predictors. Weight decay regularizes parameters. Performance is interpreted on data not used for fitting; random held-out performance remains exploratory when sampling is spatially or hierarchically dependent.

### Biological framing

Nonlinear predictability may motivate interaction or threshold hypotheses, but network weights are not direct biological effects. Use permutation importance, partial-dependence tools, or carefully designed follow-up models if explanations are required.

### Required checks

- Repeat seeds and tuning choices.
- Compare against transparent linear/logistic baselines.
- Use blocked/external validation and assess calibration.
- Avoid interpreting individual neural weights biologically.

## 9. Wide-image scalar regression screen

### Publication quantities versus Workbench quantities

SV-NN and ST-CAR infer spatially structured effect surfaces, active image regions, posterior inclusion/selection measures, uncertainty, and prediction. The Workbench instead displays standardized ridge coefficients ranked by absolute magnitude and held-out RMSE/MAE.

### Statistical interpretation

A ridge coefficient is the conditional change in predicted scalar outcome per standardized feature unit under the fitted linear screen, with all selected features and the specified penalty. Correlated neighboring pixels/voxels distribute signal across coefficients; rank is not a posterior inclusion probability.

### Biological framing

Coefficient maps or ranked regions can nominate image locations/features for follow-up. They do not establish localized causal effects and can reflect registration, mask, smoothing, scanner, or preprocessing differences.

### Required checks

- Confirm identical registration, mask, resolution, and feature ordering.
- Select penalty without test-set leakage.
- Use subject-level validation and site/scanner blocking.
- Assess coefficient stability under resampling.
- Use the full published spatial models when spatially coherent selection and posterior uncertainty are scientific requirements.

## 10. Latent image-to-image regression screen

### Publication quantities versus Workbench quantities

The SBLF publication models image predictors and outcomes with spatial latent factors and spatially varying regression coefficients, assessing prediction and spatial associations. The Workbench displays the chosen number of PCA factors, per-output held-out RMSE, overall RMSE, and an observed-versus-predicted plot for the first outcome feature.

### Statistical interpretation

Separate PCA bases reduce the predictor and outcome images, and regression connects their subject scores. Reconstruction error measures predictive fidelity in the original outcome-feature space. A latent factor is a variance direction, not automatically a biological pathway.

### Biological framing

The screen can test whether one aligned modality contains information about another—for example, whether one image-derived pattern predicts another. Spatial or mechanistic interpretation requires reconstructing factor loadings, studying stability, and validating across subjects/sites.

### Required checks

- Confirm cross-subject alignment for both modalities.
- Tune latent dimension without test leakage.
- Inspect factor loadings, reconstructed images, residual images, and subject-level errors.
- Compare against null and simpler baselines.
- Use the full SBLF model for spatial posterior inference.

## 11. Mediation with bootstrap uncertainty

### Displayed quantities

- Indirect effect: product of the exposure-to-mediator coefficient and mediator-to-outcome coefficient conditional on exposure.
- Direct effect: exposure coefficient in the outcome model including mediator.
- Total effect: exposure coefficient in the outcome-only model.
- Percentile bootstrap interval for the indirect effect.

### Statistical interpretation

These are model-based linear mediation estimates. A bootstrap interval quantifies resampling variability under the observed data-generating and sampling structure. It does not validate sequential ignorability or repair time-order/confounding problems.

### Biological framing

Use wording such as “the data are consistent with a model in which part of the exposure–outcome association operates through the measured mediator.” Reserve “mediates” as a causal verb for designs and assumptions that support it.

### Required checks

- Establish exposure before mediator before outcome.
- Defend consistency, positivity, and no unmeasured confounding for all relevant relations.
- Include prespecified confounders; the current three-variable app path does not do this.
- Conduct sensitivity analysis and account for clustering/repeated measures.

## 12. Tumor contour and shape PCA

### Displayed quantities

- Resampled, centered, size-normalized contours.
- PC1 and PC2 scores per contour.
- Percent variance explained by each component.
- PCA loadings and mean shape in the fitted object.

### Statistical interpretation

PCA identifies orthogonal directions of greatest variation in the processed landmark-coordinate vectors. A score indicates where one contour lies along a component; a loading describes how landmarks move along that component. Variance explained is geometric variability in this processed sample, not outcome association or biological importance.

### Biological framing

Shape components can nominate recurring morphological modes, but their interpretation requires visualizing mean shape plus and minus component multiples and checking whether they reflect anatomy rather than segmentation, orientation, or landmark-order artifacts.

### Required checks

- Confirm ordered boundary points and consistent orientation.
- Decide whether rotation should be standardized; the Workbench does not do it automatically.
- Assess segmentation/reader repeatability.
- Avoid overinterpreting PCs from fewer than a stable number of contours.
- Test any outcome association in a separate model with appropriate validation.

## How to present a result to collaborators

Use this template:

> We analyzed **[observational unit]** using **[workflow]** to estimate **[quantity]**. Under **[neighbor graph/model/penalty/validation]**, the estimate was **[value and uncertainty]**. This supports the descriptive statement **[statistical statement]**. It is consistent with, but does not establish, **[biological hypothesis]**. The main alternative explanations are **[technical/design/confounding issues]**. A stronger conclusion would require **[blocked/external validation, independent cohort, perturbation, orthogonal assay, or prospective design]**.

## Red-flag phrases to avoid

- “This edge proves the cells interact.”
- “The p-value proves spatial biology.”
- “The cluster is a new cell state.”
- “This coefficient identifies the causal voxel.”
- “The neural network learned the mechanism.”
- “The mediator explains the effect” without causal identification.
- “PC1 is tumor aggressiveness” without an external outcome model.

Prefer: *conditional association*, *model-selected edge*, *algorithmic cluster*, *held-out prediction*, *screening coefficient*, *model-based indirect effect*, *shape variation component*, and *hypothesis requiring validation*.

## References

See `LITERATURE_CONTEXT.md` for the primary publications, software records, and exact relationship between each publication and the Workbench implementation.

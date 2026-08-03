# Spatial Methods Workbench: method and software references

Use this bibliography together with the versioned Spatial Methods Workbench citation. Cite only the sources applicable to the workflow actually used, including the exact external engine when relevant. Screening workflows are explicitly distinguished from related full Bayesian methods.

## Workbench software

Bhadury, S. (2026). *Spatial Methods Workbench* (Version 0.2.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21764196

## Statistical and scientific method references

1. **Global spatial autocorrelation.** Moran, P. A. P. (1950). Notes on continuous stochastic phenomena. *Biometrika, 37*(1–2), 17–23. https://doi.org/10.1093/biomet/37.1-2.17

2. **Ridge regularization used in exploratory precision/network and wide-image screens.** Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. *Technometrics, 12*(1), 55–67. https://doi.org/10.1080/00401706.1970.10488634

3. **ISPat 2D.** Bhadury, S., Peruzzi, M., Acharyya, S., Eliason, J., Pasca Di Magliano, M., Frankel, T. L., Ravikumar, V., Krishnan, S., & Rao, A. (2026). Informed spatially aware patterns for multiplexed immunofluorescence data. *Scientific Reports, 16*(1). https://doi.org/10.1038/s41598-026-35341-8. Software: https://github.com/sagnikbhadury/ISPAT

4. **GP–GHS spatially varying networks.** Bhadury, S., Gaskins, J. T., & Rao, A. (2026). Spatially varying graphical models for cell-cell interaction networks in multiplexed tissue imaging. https://doi.org/10.64898/2026.04.01.715977. Software: https://github.com/sagnikbhadury/GP-GHS

5. **ISPat 3D volumetric engine.** Bhadury, S. (2026). *ISPAT-3D* [Computer software]. https://github.com/sagnikbhadury/ISPAT-3D. No unpublished manuscript is cited or disclosed by the Workbench.

6. **Structured predictive validation.** Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography, 40*(8), 913–929. https://doi.org/10.1111/ecog.02881

7. **Related spatial-clustering methodology.** Chavent, M., Kuentz-Simonet, V., Labenne, A., & Saracco, J. (2018). ClustGeo: An R package for hierarchical clustering with spatial constraints. *Computational Statistics, 33*(4), 1799–1822. https://doi.org/10.1007/s00180-018-0791-1. The Workbench implements a distinct weighted-coordinate k-means adaptation, not ClustGeo.

8. **Single-hidden-layer neural-network implementation.** Venables, W. N., & Ripley, B. D. (2002). *Modern Applied Statistics with S* (4th ed.). Springer. https://doi.org/10.1007/978-0-387-21706-2

9. **Related SV-NN scalar-on-image method.** Wu, B., Wu, K., & Kang, J. (2025). Bayesian scalar-on-image regression with a spatially varying single-layer neural network prior. *Journal of Machine Learning Research, 26*(116), 1–38. https://www.jmlr.org/papers/v26/22-0246.html

10. **Related ST-CAR image-regression method.** Xu, Y., & Kang, J. (2025). Bayesian image regression with soft-thresholded conditional autoregressive prior. *International Conference on Learning Representations 2025*. https://proceedings.iclr.cc/paper_files/paper/2025/hash/f418594e90047a10f4c158f70d6701cc-Abstract-Conference.html

11. **Related SBLF image-on-image method.** Guo, C., Kang, J., & Johnson, T. D. (2022). A spatial Bayesian latent factor model for image-on-image regression. *Biometrics, 78*(1), 72–84. https://doi.org/10.1111/biom.13420

12. **Causal mediation framework.** Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal mediation analysis. *Psychological Methods, 15*(4), 309–334. https://doi.org/10.1037/a0020761

13. **Statistical shape analysis.** Dryden, I. L., & Mardia, K. V. (2016). *Statistical Shape Analysis, with Applications in R* (2nd ed.). Wiley. https://doi.org/10.1002/9781119072492

## Public fork and upstream repository attribution

These public forks informed a shipped workflow or its related-method documentation. Both upstream and fork locations are recorded for provenance.

| Purpose | Original upstream | Development fork | Relationship to shipped Workbench |
|---|---|---|---|
| AI curriculum | https://github.com/microsoft/AI-For-Beginners | https://github.com/sagnikbhadury/AI-For-Beginners | Educational/workflow design influence; not an executable statistical engine |
| ML curriculum | https://github.com/microsoft/ML-For-Beginners | https://github.com/sagnikbhadury/ML-For-Beginners | Educational/workflow design influence; not an executable statistical engine |
| Deep-learning curriculum | https://github.com/mikexcohen/DeepUnderstandingOfDeepLearning | https://github.com/sagnikbhadury/DeepUnderstandingOfDeepLearning | Background for the shallow-neural adaptation; the Workbench executes R `nnet` |
| SV-NN | https://github.com/benwu233/SV-NN | https://github.com/sagnikbhadury/SV-NN | Related full spatial Bayesian method; Workbench implements a separately labeled ridge screen |
| ST-CAR | https://github.com/yuliangxu/STCAR | https://github.com/sagnikbhadury/STCAR | Related full spatial Bayesian method; Workbench implements a separately labeled ridge screen |
| SBLF | https://github.com/umich-biostatistics/SBLF | https://github.com/sagnikbhadury/SBLF | Related full spatial Bayesian method; Workbench implements a separately labeled PCA/score-regression screen |

Repository presence alone is not evidence that a source was used. Forks not listed in this table are not claimed as dependencies or methodological sources of the released Workbench.

## Citation boundary

- Cite the Workbench version and the method(s) actually used.
- Follow the license and citation terms of every exact external engine.
- Do not cite a related full Bayesian method as though the Workbench screening implementation reproduced its posterior model.
- The bibliography does not disclose private repositories, working manuscripts, or unpublished method ideas.

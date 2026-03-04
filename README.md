# pMEM: Predictive Moran's Eigenvector Maps for Spatial Modeling

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/pMEM)](https://CRAN.R-project.org/package=pMEM)
[![R-CMD-check](https://github.com/guenardg/pMEM/workflows/R-CMD-check/badge.svg)](https://github.com/guenardg/pMEM/actions)
[![Methods in Ecology and Evolution](https://img.shields.io/badge/DOI-10.1111/2041--210X.14413-blue)](https://doi.org/10.1111/2041-210X.14413)

**pMEM** implements *Predictive Moran's Eigenvector Maps*, a method for
spatially-explicit prediction of environmental variables using
eigen-decomposition of distance-based spatial weighting matrices.

Guénard, G., & Legendre, P. (2024). Spatially-explicit predictions using spatial
eigenvector maps. *Methods in Ecology and Evolution*.
<https://doi.org/10.1111/2041-210X.14413>

---

## What Does pMEM Do?

pMEM extends classical Moran's Eigenvector Maps (MEM) by:

- **Enabling prediction at unsampled locations** via `predict()` methods
- **Supporting tunable distance weighting functions** (exponential, Gaussian,
  power, linear, spherical, hyperbolic, hole_effect)  
- **Optimizing spatial scale parameters** via cross-validated model selection  
- **ntegrating with regression frameworks**: linear models, GLMs, elastic net  
- **Providing fast C++ backend** via Rcpp for efficient eigenvector selection  
- **Supporting asymmetric (directional) spatial processes** via complex-valued
  distance metrics  

Ideal for ecologists, geographers, and spatial analysts modeling
spatially-autocorrelated variables such as:

- Environmental gradients (depth, temperature, substrate)
- Species distributions and abundances
- Landscape connectivity metrics
- Directional processes (river flow, wind patterns, ocean currents)

---

## Installation

### From CRAN (stable release)

`install.packages("pMEM")`

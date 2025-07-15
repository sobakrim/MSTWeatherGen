# MSTWeatherGen <img src="man/figures/MSTWeatherGen.png" align="right" alt="" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/sobakrim/MSTWeatherGen/actions/workflows/r-package-check.yml/badge.svg)](https://github.com/sobakrim/MSTWeatherGen/actions/workflows/r-package-check.yml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MSTWeatherGen)](https://cran.r-project.org/package=MSTWeatherGen)
<!-- badges: end -->

## Overview

{MSTWeatherGen} Multivariate and Space-Time stochastic Weather Generator R package is designed for the simulation of multivariate spatio-temporal meteorological variables. It provides tools for estimating the model parameters and generating synthetic weather data that can be used for a variety of applications, including climate research, agricultural or hydrological modeling.  

For more details on the paper, see the publication: 
<blockquote>
  <p>
    Obakrim, S. <em>et&nbsp;al.</em> “A multivariate and space-time stochastic weather generator using a latent Gaussian framework,”
    <em> Stochastic Environmental Research and Risk Assessment</em>, 2025.
    <a href="https://doi.org/10.1007/s00477-024-02897-8">doi:10.1007/s00477-024-02897-8</a>
  </p>
</blockquote>


## Installation

You can install `{MSTWeatherGen}` from GitHub as follows:
```r
# install.packages("remotes")
remotes::install_github("sobakrim/MSTWeatherGen")
```

## Getting Started

To learn how to use the `{MSTWeatherGen}` package, please refer to the detailed vignette available [here](https://sobakrim.github.io/MSTWeatherGen/articles/MSTWeatherGen.html).

## Funding
This work was supported by funding from the French National Research Agency (ANR) as part of the BEYOND project (Contract No. 20-PCPA-0002). We acknowledge their support for enabling the development of this package.
## License

The package MSTWeatherGen is under GNU GPL V3.   
See [LICENSE](LICENSE) file.  

## Authors

- Said Obakrim  (Author and Maintainer))
- Jean-François Rey (Author)

## Contributors

- Lionel Benoit 
- Denis Allard



<!-- README.md is generated from README.Rmd. Please edit that file -->

# epistasisGA

<!-- badges: start -->
<!-- badges: end -->

The epistasisGA package implements the GADGETS approach for detecting
gene-gene interactions in case-parent triad or affected/unaffected
sibling studies.

## Installation

The epistasisGA package is now available from Bioconductor (release
version:
<https://www.bioconductor.org/packages/release/bioc/html/epistasisGA.html>,
devel version:
<https://www.bioconductor.org/packages/devel/bioc/html/epistasisGA.html>).
The most current version will be available from the devel version link,
and also from this github page. To install from Bioconductor, see the
links above.

epistasisGA remains available through github. The main functions of the
package rely on C++ code via the Rcpp and RcppArmadillo packages. These
packages require a suitable compiler for installation. The requirements
for different operating systems are described in section 1.3 of the Rcpp
package FAQs:
<https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf>.
Further details for macOS users can be found in the ‘R Administration’
manual, Appendix C.3:
<https://cran.r-project.org/doc/manuals/r-release/R-admin.html#macOS>.
Additionally, to build the package vignette, Pandoc is required. To
check whether you have pandoc installed, see here:
<https://rmarkdown.rstudio.com/docs/reference/pandoc_available.html>.
After ensuring these prerequisites are installed, epistasisGA can be
installed using the `devtools` package with the following commands:

``` r
library(devtools)
devtools::install_github("mnodzenski/epistasisGA", build_vignettes = TRUE, dependencies = TRUE)
```

## Vignette

Please consult the package vignette for example usages of the functions
in epistasisGA. The vignette can be accessed by entering the following
command and clicking on the resulting HTML links:

``` r
browseVignettes("epistasisGA")
```

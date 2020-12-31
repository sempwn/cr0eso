# Covid r0 estimation of outbreaks (cr0eso)

R0 estimation of multiple outbreaks using an SEIR model

*The package name comes from the Welsh word for "welcome"*

## Overview

**cr0eso** fits a Bayesian hierarchical SEIR (susceptible, exposed, infected, recovered) to a collection of outbreak data from different facilities.

## Installation

Before installation, you will need a C++ compiler installed to compile the Stan model. If you are on a Mac or Linux, follow [these instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). If you are on Windows, follow the instructions [here](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows).

The following must return `TRUE` before continuing:

```r
# install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

Then, install the **cr0eso** package with:

```r
# install.packages("remotes")
remotes::install_github("sempwn/cr0eso")
```

See the examples in `?seir_model_fit`.

To install the package with the vignette:

```r
remotes::install_github("sempwn/cr0eso", build_vignettes = TRUE)
```

and run:

```r
browseVignettes("cr0eso")
```

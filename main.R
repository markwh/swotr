# main.R
# Mark Hagemann
# 7/10/2018
# Preparing the package for use

library(devtools)

use_package("dplyr")
use_package("purrr")
use_package("ggplot2")
use_package("tidyr")
use_package("reshape2")

use_package("bamr", type = "suggests")
use_package("rstan", type = "suggests")
use_package("ncdf4", type = "suggests")


use_vignette("swotlists")
use_vignette("test") # Just so I can see the hints in the default file. To be delted.

#' Load the extdata folder with an example netcdf
file.copy("../SWOT/data/ncdata/SacramentoDownstream.nc", to = "inst/extdata")


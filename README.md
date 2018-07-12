# swotr

<img src="https://raw.githubusercontent.com/markwh/swotr/master/logo/logo.png" width=200 alt="Swotr Logo"/>

**swotr** is a package for working with SWOT-like remote-sensed river data (or similar) in R. Broadly, this means data observed at discrete times and locations within a river or river network. Data do not need to come from satellite remote sensing in order to be SWOT-like; the term can include in-situ, model-output, or purely synthetic data. But the underlying motivation is the upcoming SWOT satellite mission, so satellite remote sensing of discharge and hydrologic paramters is the presumed purpose in **swotr**'s conception and development.


**Documentation** including [function reference](https://markwh.github.io/swotr/reference/index.html) and a [tutorial vignette](https://markwh.github.io/swotr/articles/swotlists.html) can be found at the [swotr website](https://markwh.github.io/swotr/index.html).

**Install** as follows:

```
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("markwh/swotr")
```

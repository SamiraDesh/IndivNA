# IndTempNetAna
The IndivNetwork package. An R tool for individual temporal network analysis.


# Installation
This is currently the only available and developmental version of the package. You can install this package from source on Windows. 

```ruby
setwd("path to package folder")
devtools::document()
devtools::install()
```

For successful installation, you need to have `devtools` and `RTools` installed.

# Data
Here is a table that illustrates the required data structure. Y1 through Y10 represent the nodes i.e., the variables that form our network. X1 through X5 represent the characteristics/ covariates of each observation that may affect the associations between nodes i.e. they account for individual variability in the estimated networks. There are two observations (indexed by "id") with values at three time points each. 

<picture>

</picture>



There is no constraint on the number of nodes or covariates. However, only binary nodes and binary or numeric covariates are permissible.  




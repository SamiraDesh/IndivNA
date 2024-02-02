# IndTempNetAna
The IndivNetwork package. An R tool for individual temporal network analysis.


# Installation
This is currently the only available and developmental version of the package. You can install it from source on Windows. 

```ruby
setwd("path to package folder")
devtools::document()
devtools::install()
```

For successful installation, you need to have `devtools` and `RTools` installed.

# Data
Here is a table that illustrates the required data structure. Y1 through Y10 represent the nodes i.e., the variables that form our network. X1 through X5 represent the characteristics/ covariates of each observation that may affect the associations between nodes i.e. they account for individual variability in the estimated networks. There are two observations with values at three time points each. 

<picture>

</picture>



There is no constraint on the number of nodes or covariates. However, only binary nodes and binary or scaled, continuous covariates are permissible.  

# IndivIsing
This function is used to estimate the network. It only supports lag-1 factorization in estimating temporal networks.
In addition to the dataset, the following input parameters need to be specified:
1. `y_index` and `x_index`, the indices of the nodes and individual characteristics, respectively.
2. `networkModel`, 'static' or 'temporal' as appropriate.
3. `timepoint`, column name of the timepoint variable for the construction of temporal networks. If `networkModel='static'`, then need not be specified.
4. `family`, must be 'binomial' since it is the only supported type for nodes.
5. `lowerbound.lambda`, minimum value of tuning parameter lambda (regularization parameter). If one has two networks of different sample sizes but the same number of parameters p, they can be directly compared by setting this value equal to sqrt(log(p)/n), n being the number of observations in the smaller group. 
6. `gamma`, hyperparameter gamma in the extended BIC.
7. `AND`, can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network.

Some important results returned by IndivIsing are: 
1. `estimated_thresholds` and `estimated_coeff_raw`, the estimated intercepts and coefficients in the extended Ising model, respectively.
2. `estimated_bias` The formula of the estimated bias term, including the estimated intercept and the coefficients of the major effects of the individual characteristics. 
3. `estimated_formula` The formula of the estimated edge weights, including the coefficients of the major effects of the nodes and the coefficients of the interaction terms of the nodes and the individual characteristics.

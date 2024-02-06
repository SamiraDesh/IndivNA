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

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/NA_exdata.PNG?raw=true)


There is no constraint on the number of nodes or covariates. However, only binary nodes and binary or scaled, continuous covariates are permissible.  

# IndivIsing
Estimates the static or temporal network. It only supports lag-1 factorization in estimating temporal networks.
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
  
An estimated formula that looks as follows
![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndIsing_example1.PNG)

for the circled formula suggests that the weight of the edge directed from node Y7 to node Y2 consists of the main effect of Y7 on Y2 (0.555) and the effect of the interaction between Y7 and X3 (0.792) on Y2.
# IndivNetwork
Generates the individual static or temporal network for a specific subject. The required inputs are:
1. `IsingResult`, the object returned by the IndivIsing function.
2. `target_id` The ID of the specific subject for whom an individual static or temporal network is to be fitted. 
3. `covar_vec` The vector of individual characteristics for a specific subject of interest needed if `target_id` has not been assigned. The order has to be the same as x_index.
4. `data` and `id`, the original data set and name of the column in it that indicates the ID of each subject.

A subject with all covariate features equal to 1 with the estimated network illustrated above has an estimated adjacent matrix that looks like this:

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivNetwork_example.PNG)

This result is stored as `estimated_network`. Another important output of this function is `node_centrality`, which is the matrix of calculated centrality values. In the temporal network case, the statistics estimated are Betweenness, Closeness, InStrength, OutStrength, OutExpectedInfluence and InExpectedInfluence whereas in the case of a static network, given the non-directionality of edges, only Betweenness, Closeness, Strength and ExpectedInfluence will be estimated.
 
# IndivBootSplitting
 Performs data splitting inference by permutation test on the individual static or temporal network using bootstrap to evaluate its accuracy and stability. Inputs include:
 1. `type`, the type of bootstrap procedure. This can either be `"nonparametric"` for testing the accuracy of the estimated edges and centrality indices or `"casedropping"` for assessing the stability of the calculated centrality indices.
 2. `sample_prob` The vector of sampling proportion for casedropping bootstrap. Default to (0.9,0.8,...0.5).
#' @param nBoots The number of bootstraps. Default to 500.
#' @param nCores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#'   Defaults to 1.
#' @param verbose Logical. Should progress of the function be printed to the console?
#'   Defaults to TRUE.
#' @param library Library location to be used in parallel computing. Default to .libPaths().
#' @param seed The random seed.
 

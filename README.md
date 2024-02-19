# The IndivNetwork package
An R tool for individual static and temporal network analysis.


# Installation
This is currently the only available and developmental version of the package. You can install it from source on Windows. 

```ruby
setwd("path to package folder")
devtools::document()
devtools::install()
```

For successful installation, you need to have `devtools` and `RTools` installed.

# Data
Here is a table that illustrates the required data structure in the temporal network scenario. Y1 through Y10 represent the nodes i.e., the variables that form our network. X1 through X5 represent the characteristics/ covariates of each observation that may affect the associations between nodes i.e. they account for individual variability in the estimated networks. There are two samples with observations at three time points each. In the case of a static network, the structure remains the same, except there is no time column and each sample only has one observation.

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/NA_exdata.PNG?raw=true)


There is no constraint on the number of nodes or covariates. However, only binary nodes and binary or scaled, continuous covariates are permissible.  

# IndivIsing
Estimates the static or temporal network. It only supports lag-1 factorization in estimating temporal networks. The following input parameters need to be specified:
1. `data`, the training portion of the original data set in case of data splitting inference and the entire data set otherwise.
2. `y_index` and `x_index`, the indices of the nodes and individual characteristics, respectively.
3. `networkModel`, 'static' or 'temporal' as appropriate.
4. `timepoint`, column name of the timepoint variable for the construction of temporal networks. If `networkModel='static'`, then need not be specified.
5. `family`, must be 'binomial' since it is the only supported type for nodes.
6. `lowerbound.lambda`, minimum value of tuning parameter lambda (regularization parameter). If one has two networks of different sample sizes but the same number of parameters p, they can be directly compared by setting this value equal to sqrt(log(p)/n), n being the number of observations in the smaller group. 
7. `gamma`, hyperparameter gamma in the extended BIC.
8. `AND`, can be TRUE of FALSE to indicate whether the AND-rule or the OR-rule should be used to define the edges in the network.

Some important results returned by IndivIsing are: 
1. `estimated_thresholds` and `estimated_coeff_raw`, the estimated intercepts and coefficients in the extended Ising model, respectively.
2. `estimated_bias` The formula of the estimated bias term, including the estimated intercept and the coefficients of the major effects of the individual characteristics. 
3. `estimated_formula` The formula of the estimated edge weights, including the coefficients of the major effects of the nodes and the coefficients of the interaction terms of the nodes and the individual characteristics.
  
An estimated formula for the temporal network may look like this

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivIsing_common.PNG)

and the circled formula suggests that the weight of the edge directed from node Y6 to node Y2 is equal to the effect of the interaction between Y6 and X4 (0.096) on Y2. There is no main effect of Y6 on Y2. There also exists an edge that is directed from a node (Y10 and Y6) to itself, which is possible only in the temporal network case since we account for values of the node at a prior timepoint.

# IndivNetwork
Generates the individual static or temporal network for a specific subject. The required inputs are:
1. `IsingResult`, the object returned by the IndivIsing function.
2. `target_id` The ID of the specific subject for whom an individual static or temporal network is to be fitted. 
3. `covar_vec` The vector of individual characteristics for a specific subject of interest needed if `target_id` has not been assigned. The order has to be the same as x_index.
4. `data` and `id`, the original data set and name of the column in it that indicates the ID of each subject.

Considering the estimated temporal network above, at a specific timepoint the results for a subject could look like
Estimated adjacency matrix       |  Estimated network
:-------------------------:|:-------------------------:
![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivNtwrk_temporal.PNG)  |  ![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivNtwrk_temporal2.png)

A static network for another subject may look something like -

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivNtwrk_static.PNG)

with this result translating to this network structure - 

![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivNtwrk_static2.png)


Another important output of this function is `node_centrality`, which is the matrix of calculated centrality values. In the temporal network case, the statistics estimated are Betweenness, Closeness, InStrength, OutStrength, OutExpectedInfluence and InExpectedInfluence whereas in the case of a static network, given the non-directionality of edges, only Betweenness, Closeness, Strength and ExpectedInfluence will be estimated.

# IndivBootSplitting
Evaluates the accuracy and stability of the estimated individual static or temporal network. This is accomplished through bootstrap data splitting inference by permutation test, discussed in our previous work (REF). Inputs include:
 1. `data`, the testing/ inference portion of the original data set.
 3. `type`, the type of bootstrap procedure. This can either be `"nonparametric"` for testing the accuracy of estimated edges and centrality indices or `"casedropping"` for assessing the stability of the calculated centrality indices.
 4. `sample_prob`, the vector of sampling proportion for casedropping bootstrap.
 5. `nBoots` the number of bootstraps.

Two outputs returned by this function that are of note -
1. `IsingResult_new`, the new IndivIsing model fitted on the inference data set.
2. `bootResults`, a list containing the bootstraped networks.


# IndivBootSummary
Summarizes the bootstrap results. Important inputs are:

1. `data`, the inference data set.
2. `IsingBootResult`, the fitted IndivBootSplitting object.
3. `IsingResult`, The fitted IndivIsing object.
4. `sample_covar_mat` The matrix of covariate information for representative subjects to perform casedropping bootstrap.


Key results returned by IndivBootSummary:
1. `test_centrality`, matrix listing the stacked testing results for different centralities between nodes for the first representative subject.
2. `test_centrality_long`, matrix with the long format listing the stacked testing results for different centralities between nodes for the first representative subject.
   
Additional results are returned depending on the type of bootstrap. When running the nonparametric bootstrap, this function also creates -
1. `boot_data`, a matrix containing all the edge weights estimated in each bootstrap replicate.
2. `boot_data_summary`, a matrix containing the summary of all bootstrap results.
3. `edge_data`, a matrix containing all edge weights retained in the original constructed network that were estimated in each bootstrap replicate.
4. `edge_data_summary`, a matrix containing the summary of all edge weights retained in the original constructed network.
5. `edge_test`, a matrix listing the test results of whether the selected edges in the IndivIsing object are significantly different from each other.
6. `edge_quantile`, a matrix listing the range of the differences between two estimated edges out of the selected edges in the IndivIsing object.
7. `boot_centrality`, a matrix summarizing the centrality information for each subject at each bootstrap. 

On the other hand, for the casedropping bootstrap the summary of the correlation results is returned as `cor_results`.

The forest plots below illustrate the summarized edge weights for static and temporal networks - 
Static        |  Temporal
:-------------------------:|:-------------------------:
![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivBootSummary_static.png)  |  ![alt text](https://github.com/SamiraDesh/IndTempNetAna/blob/main/IndivBootSummary_temporal.png)








#' Summarize the estimating results based on replicates.
#' 
#' @param data The data frame containing the binary nodes and the individual characteristics
#'   for network construction. If a static network is to be constructed, there is one observation
#'   for each subject. Otherwise, if a temporal network is to be constructed, multiple observations 
#'   for each subject is required. Only lag-1 factorization is supported in the current version of
#'   estimating temporal networks.
#' @param y_index The indices of the nodes.
#' @param x_index The indices of the individual characteristics. Pre-processing of the 
#'   individual features are required (e.g., scaling of continuous variables, generating 
#'   dummy variables for categorical variables) before performing analysis using this function.
#' @param networkModel Model that should be used for network estimation. Can be either 'static' or 'temporal'.
#' @param timepoint The column name in 'data' that denotes timepoint for the construction 
#'   of temporal network. If `networkModel='static'`, `timepoint` does not need to be specified.Default to NA.
#' @param family Currently only 'binomial' family of data is supported.
#' @param progressbar Logical. Should the pbar be plotted in order to see the progress 
#'   of the estimation procedure? Default to TRUE.
#' @param lowerbound.lambda The minimum value of tuning parameter lambda (regularization parameter). 
#'   Can be used to compare networks that are based on different sample sizes. The lowerbound.lambda 
#'   is based on the number of observations in the smallest group n: sqrt(log(p)/n). 
#'   p is the number of variables, that should be the same in both groups. When both 
#'   networks are estimated with the same lowerbound for lambda (based on the smallest group), 
#'   the two networks can be directly compared. Default to NA. 
#' @param gamma A value of hyperparameter gamma in the extended BIC. Can be anything 
#'   between 0 and 1. Defaults to 0.25.
#' @param AND Logical. Can be TRUE of FALSE to indicate whether the AND-rule or 
#'   the OR-rule should be used to define the edges in the network. Defaults to TRUE
#' @param id The column name in 'data' that denotes the ID of each subject.
#' @param rep The number of replicates of random splitting on the original dataset for network estimation.
#' @param nBoots The number of bootstraps. Default to 500.
#' @param nCores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#'   Defaults to 1.
#' @param sample_prob The vector of sampling proportion for casedropping bootstrap.
#' @param verbose Logical. Should progress of the bootstrap function be printed to the console?
#'   Defaults to TRUE.
#' @param library Library location to be used in parallel computing.
#' @param seed The random seed.
#' @param sample_covar_mat The matrix of covariate information for representative subjects to perform casedropping bootstrap.
#' @return An IndivRepSummary object.
#' @returns `list_estimated_thresholds` The list of estimated thresholds based on replicates.
#' @returns `list_estimated_coeff_raw` The list of estimated raw coefficients based on replicates.
#' @returns `list_boot_data` The list of estimated bootstrapped edge values based on replicates.
#' @returns `list_boot_centrality` The list of estimated bootstrapped centrality values based on replicates.
#' @returns `list_cor_final_results` The list of estimated correlations of centralities based on replicates.
#' @returns `data` The input dataset.
#' @returns `y_index` The input indices of nodes.
#' @returns `x_index` The input indices of individual characteristics.
#' @returns `networkModel` The input networkModel parameter.
#' @returns `timepoint` The input timepoint parameter.
#' @returns `family` The input family parameter.
#' @returns `progressbar` The input progressbar parameter.
#' @returns `lowerbound.lambda` The input lowerbound.lambda parameter.
#' @returns `gamma` The input gamma parameter.
#' @returns `AND` The input AND parameter.
#' @returns `rep` The input rep parameter.
#' @export
IndivRep <- function(data,
                     y_index,
                     x_index,
                     networkModel = c("static", "temporal"),
                     timepoint = NA,
                     family = "binomial",
                     progressbar = TRUE, 
                     lowerbound.lambda = NA, 
                     gamma = 0.25, 
                     AND = TRUE,
                     id,
                     rep = 200,
                     nBoots = 500,
                     nCores = 1,
                     sample_prob = seq(0.9,0.5,-0.1),
                     verbose = TRUE,
                     library = .libPaths(),
                     seed,
                     sample_covar_mat){
  
  ##############################################################################
  ##############################################################################
  t0 <- Sys.time()
  
  ##############################################################################
  ## check conditions ##########################################################
  if(family != "binomial"){
    stop("This procedure currently only supports binary (family='binomial') data.")
  }
  
  if(!(networkModel %in% c("static","temporal"))){
    stop("Model has to be either 'static' or 'temporal'.")
  }else if(networkModel=="temporal"){
    if(is.na(timepoint)){
      stop("Timepoint has to be assigned to construct temporal networks.")
    }else{
      times <- unique(data[,timepoint])  
    }
  }
  
  if(missing(id)){
    stop("The id variable has to be given for bootstrap.")
  }
  
  if(missing(seed)){
    seed <- 1234
  }
  
  ################################################################################
  ################################################################################
  y_names <- colnames(data)[y_index]
  x_names <- colnames(data)[x_index]
  
  list_estimated_thresholds <- vector("list", rep)
  list_estimated_coeff_raw <- vector("list", rep)
  list_boot_data <- vector("list", rep)
  list_boot_centrality <- vector("list", rep)
  list_cor_final_results <- vector("list", rep)
  
  for(i in 1:rep){
    
    set.seed(seed+i)
    train_index <- sample(1:nrow(data), ceiling(nrow(data)/2), replace = FALSE)
    
    data_train <- data[train_index,]
    data_inference <- data[-train_index,]
    
    result_Ising <- IndivIsing(data = data_train, 
                               y_index = y_index, 
                               x_index = x_index, 
                               networkModel = networkModel,
                               timepoint = timepoint, 
                               family = family,
                               progressbar = progressbar, 
                               lowerbound.lambda = lowerbound.lambda, 
                               gamma = gamma, 
                               AND = AND)
    
    list_estimated_thresholds[[i]] <- result_Ising$estimated_thresholds
    list_estimated_coeff_raw[[i]] <- result_Ising$estimated_coeff_raw
    
    boot_nonparametric <- IndivBootSplitting(data=data_inference, 
                                             IsingResult=result_Ising,
                                             id=id, 
                                             type = "nonparametric",
                                             sample_prob = sample_prob,
                                             nBoots = nBoots, 
                                             nCores = nCores, 
                                             verbose = verbose, 
                                             library = library,
                                             seed = seed)
    
    boot_casedropping <- IndivBootSplitting(data=data_inference, 
                                            IsingResult=result_Ising,
                                            id=id, 
                                            type = "casedropping",
                                            sample_prob = sample_prob,
                                            nBoots = nBoots, 
                                            nCores = nCores, 
                                            verbose = verbose, 
                                            library = library,
                                            seed = seed)
    
    boot_nonparametric_summary <- IndivBootSummary(data=data_inference,
                                                   IsingBootResult=boot_nonparametric,
                                                   IsingResult=boot_nonparametric$IsingResult_new,
                                                   nCores = nCores,
                                                   library = library,
                                                   verbose = verbose)
    
    boot_casedropping_summary <- IndivBootSummary(data=data_inference,
                                                  IsingBootResult=boot_casedropping,
                                                  IsingResult=boot_casedropping$IsingResult_new,
                                                  nCores = nCores,
                                                  library = library,
                                                  verbose = verbose,
                                                  sample_covar_mat)
    
    list_boot_data[[i]] <- boot_nonparametric_summary$boot_data
    list_boot_centrality[[i]] <- boot_nonparametric_summary$boot_centrality
    list_cor_final_results[[i]] <- boot_casedropping_summary$cor_final_results
    
  }
  
  ########################################################################################
  ########################################################################################
  return(list("list_estimated_thresholds"=list_estimated_thresholds,
              "list_estimated_coeff_raw"=list_estimated_coeff_raw,
              "list_boot_data"=list_boot_data,
              "list_boot_centrality"=list_boot_centrality,
              "list_cor_final_results"=list_cor_final_results,
              "data"=data,
              "y_index"=y_index,
              "x_index"=x_index,
              "networkModel"=networkModel,
              "timepoint"=timepoint,
              "family"=family,
              "progressbar"=progressbar,
              "lowerbound.lambda"=lowerbound.lambda,
              "gamma"=gamma,
              "AND"=AND,
              "rep"=rep))
}





















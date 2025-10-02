#' Perform permutation test on the individual static or temporal network using bootstrap.
#' 
#' This function performs permutation tests on the constructed individual static or
#' temporal network using bootstrap to evaluate its accuracy and stability. A nonparametric 
#' bootstrap sampled with replacement is performed to assess the accuracy of estimated 
#' network connections. A case-dropping bootstrap sampled without replacement is performed to 
#' assess the stability of the calculated centrality indices. 
#' 
#' @param data The data frame used to fit the IndivIsing object.
#' @param IsingResult The fitted IndivIsing object.
#' @param id The column name in 'data' that denotes the ID of each subject.
#' @param type The type of bootstrap procedure. Can be either "nonparametric" for edge testing or 
#'   "casedropping" for centrality testing.
#' @param sample_prob The vector of sampling proportion for casedropping bootstrap. Default to (0.9,0.8,...0.5).
#' @param nBoots The number of bootstraps. Defaults to 500.
#' @param nCores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#'   Defaults to 1.
#' @param verbose Logical. Should progress of the function be printed to the console?
#'   Defaults to TRUE.
#' @param library Library location to be used in parallel computing. Default to .libPaths().
#' @param seed The random seed.
#' @return An IndivBootnet object.
#' @returns `bootResults` A list containing the bootstrap results. A list of fitted networks for 
#'  "nonparametric" bootstrap, and a list of centrality correlation values for "casedropping" bootstrap.
#' @returns `id` The input id variable.
#' @returns `type` The input type of bootstrap.
#' @returns `sample_prob` The input sampling proportion for casedropping bootstrap.
#' @returns `time` Time spent on running the function.
#' @export 
IndivBootnet <- function(data, 
                         IsingResult,
                         id, 
                         type = c("nonparametric", "casedropping"),
                         sample_prob = seq(0.9,0.5,-0.1),
                         nBoots = 500, 
                         nCores = 1, 
                         verbose = TRUE, 
                         library = .libPaths(),
                         seed) 
{
  ##############################################################################
  ##############################################################################
  t0 <- Sys.time()
  
  ##############################################################################
  ## check conditions ##########################################################
  if(missing(id)){
    stop("The id variable has to be given for bootstrap.")
  }
  
  if(!(type %in% c("nonparametric", "casedropping"))){
    stop("This procedure currently only supports either nonparametric or casedropping bootstrap.")
  }
  
  if(missing(seed)){
    seed <- 1234
  }
  
  ##############################################################################
  id_vec <- unique(data[,id])
  Np <- length(id_vec)
  
  if(nCores == 1){
    set.seed(seed)
    
    bootResults <- vector("list", nBoots)
    
    if(verbose){
      message("Bootstrapping...")
      pb <- txtProgressBar(0, nBoots, style = 3)
    }
    
    for(b in seq_len(nBoots)){
      tryLimit <- 10  #maximum number of errors allowed
      tryCount <- 0  #start number of errors
      
      if(type=="nonparametric"){
        repeat{
          bootSample <- sample(seq_len(Np), Np, replace = TRUE)
          bootData <- do.call(rbind, lapply(1:length(bootSample), function(k) data[which(data[,id]==id_vec[bootSample[k]]),]))
          
          res <- suppressWarnings(try({
            if(IsingResult$networkModel=="temporal"){
              IndivIsing(data = bootData, 
                         y_index = IsingResult$y_index, 
                         x_index = IsingResult$x_index, 
                         networkModel = "temporal",
                         timepoint = IsingResult$timepoint, 
                         family = IsingResult$family, 
                         progressbar = FALSE, 
                         lowerbound.lambda = IsingResult$lowerbound.lambda, 
                         gamma = IsingResult$gamma, 
                         AND = IsingResult$AND)
            }else{
              IndivIsing(data = bootData, 
                         y_index = IsingResult$y_index, 
                         x_index = IsingResult$x_index, 
                         networkModel = "static",
                         timepoint = NA,
                         family = IsingResult$family, 
                         progressbar = FALSE, 
                         lowerbound.lambda = IsingResult$lowerbound.lambda, 
                         gamma = IsingResult$gamma, 
                         AND = IsingResult$AND)
            }
          }))
          
          gc()
          ##########################################################
          if(is(res, "try-error")){
            if(tryCount == tryLimit){
              stop("Maximum number of errors in bootstrap reached.")
            }
            tryCount <- tryCount + 1
          }
          else{
            break
          }
        }
        
        ##########################################################
        bootResults[[b]] <- res
        if(verbose){
          setTxtProgressBar(pb, b)
        }
      }else if(type=="casedropping"){
        repeat{
          
          casedrop_net <- function(i){
            bootSample <- sample(seq_len(Np), round(sample_prob[i]*Np), replace = FALSE)
            bootData <- do.call(rbind, lapply(1:length(bootSample), function(k) data[which(data[,id]==id_vec[bootSample[k]]),]))
            
            res_prob <- suppressWarnings(try({
              if(IsingResult$networkModel=="temporal"){
                IndivIsing(data = bootData, 
                           y_index = IsingResult$y_index, 
                           x_index = IsingResult$x_index, 
                           networkModel = "temporal",
                           timepoint = IsingResult$timepoint, 
                           family = IsingResult$family, 
                           progressbar = FALSE, 
                           lowerbound.lambda = IsingResult$lowerbound.lambda, 
                           gamma = IsingResult$gamma, 
                           AND = IsingResult$AND)
              }else{
                IndivIsing(data = bootData, 
                           y_index = IsingResult$y_index, 
                           x_index = IsingResult$x_index, 
                           networkModel = "static",
                           timepoint = NA,
                           family = IsingResult$family, 
                           progressbar = FALSE, 
                           lowerbound.lambda = IsingResult$lowerbound.lambda, 
                           gamma = IsingResult$gamma, 
                           AND = IsingResult$AND)
              }
            }))
            
            gc()
            
            return(res_prob)
          }
          
          res <- lapply(1:length(sample_prob), function(i) casedrop_net(i))  
          
          ##########################################################
          if(is(res, "try-error")){
            if(tryCount == tryLimit){
              stop("Maximum number of errors in bootstrap reached.")
            }
            tryCount <- tryCount + 1
          }
          else{
            break
          }
        }
        
        ##########################################################
        bootResults[[b]] <- res
        if(verbose){
          setTxtProgressBar(pb, b)
        }
      }
    }
    
    if(verbose){
      close(pb)
    }
  }else{  ##if (nCores>1)
    
    if(verbose){
      message("Bootstrapping...")
    }
    
    # if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
    #     Sys.info()["sysname"] == "Darwin" && gsub("\\..*", "", getRversion()) == "4"){
    #   snow::setDefaultClusterOptions(setup_strategy = "sequential")
    # }
    # 
    # nClust <- nCores
    # cl <- snow::makeSOCKcluster(nClust)
    # 
    # excl <- c("prepFun", "prepArgs", "estFun", "estArgs",
    #           "graphFun", "graphArgs", "intFun", "intArgs")
    # parallel::clusterExport(cl, c(ls()[!ls() %in% c(excl, "cl")]), envir = environment())  #, "glmnet", "colSums", "IndivIsing","IndivNetwork","qgraph","IsingResult","timepoint"
    # 
    # parallel::clusterSetRNGStream(cl, seed)
    
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    
    bootResults <- parallel::mclapply(seq_len(nBoots), function(b) {
      
      .libPaths(library)
      tryLimit <- 10
      tryCount <- 0
      
      if(type=="nonparametric"){
        repeat{
          bootSample <- sample(seq_len(Np), Np, replace = TRUE)
          bootData <- do.call(rbind, lapply(1:length(bootSample), function(k) data[which(data[,id]==id_vec[bootSample[k]]),]))
          
          res <- suppressWarnings(try({
            if(IsingResult$networkModel=="temporal"){
              IndivIsing(data = bootData,
                         y_index = IsingResult$y_index,
                         x_index = IsingResult$x_index,
                         networkModel = "temporal",
                         timepoint = IsingResult$timepoint,
                         family = IsingResult$family,
                         progressbar = FALSE,
                         lowerbound.lambda = IsingResult$lowerbound.lambda,
                         gamma = IsingResult$gamma,
                         AND = IsingResult$AND)
            }else{
              IndivIsing(data = bootData,
                         y_index = IsingResult$y_index,
                         x_index = IsingResult$x_index,
                         networkModel = "static",
                         timepoint = NA,
                         family = IsingResult$family,
                         progressbar = FALSE,
                         lowerbound.lambda = IsingResult$lowerbound.lambda,
                         gamma = IsingResult$gamma,
                         AND = IsingResult$AND)
            }
          }))
          
          gc()
          
          ##########################################################
          if(is(res, "try-error")){
            if(tryCount == tryLimit){
              stop("Maximum number of errors in bootstrap reached.")
            }
            tryCount <- tryCount + 1
          }
          else{
            break
          }
        }
        
        ##########################################################
      }else if(type=="casedropping"){
        
        repeat{
          
          casedrop_net <- function(i){
            bootSample <- sample(seq_len(Np), round(sample_prob[i]*Np), replace = FALSE)
            bootData <- do.call(rbind, lapply(1:length(bootSample), function(k) data[which(data[,id]==id_vec[bootSample[k]]),]))
            
            res_prob <- suppressWarnings(try({
              if(IsingResult$networkModel=="temporal"){
                IndivIsing(data = bootData,
                           y_index = IsingResult$y_index,
                           x_index = IsingResult$x_index,
                           networkModel = "temporal",
                           timepoint = IsingResult$timepoint, 
                           family = IsingResult$family, 
                           progressbar = FALSE, 
                           lowerbound.lambda = IsingResult$lowerbound.lambda, 
                           gamma = IsingResult$gamma, 
                           AND = IsingResult$AND)
              }else{
                IndivIsing(data = bootData,
                           y_index = IsingResult$y_index,
                           x_index = IsingResult$x_index,
                           networkModel = "static",
                           timepoint = NA, 
                           family = IsingResult$family, 
                           progressbar = FALSE, 
                           lowerbound.lambda = IsingResult$lowerbound.lambda, 
                           gamma = IsingResult$gamma, 
                           AND = IsingResult$AND)
              }
            }))
            
            gc()
            
            return(res_prob)
          }
          
          res <- lapply(1:length(sample_prob), function(i) casedrop_net(i))
          
          ##########################################################
          if(is(res, "try-error")){
            if(tryCount == tryLimit){
              stop("Maximum number of errors in bootstrap reached.")
            }
            tryCount <- tryCount + 1
          }
          else{
            break
          }
        }
        
      }
      
      return(res)
    }, cl=cl)
  }
  
  return(list("bootResults" = bootResults,
              "id" = id,
              "type" = type,
              "sample_prob" = sample_prob,
              "time" = Sys.time()-t0))
}
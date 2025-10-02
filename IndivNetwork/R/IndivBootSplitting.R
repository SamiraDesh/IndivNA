#' Perform data splitting inference by permutation test on the individual static or temporal network.
#' 
#' This function performs data splitting inference on the constructed individual static or
#' temporal network using bootstrap to evaluate its accuracy and stability. A nonparametric 
#' bootstrap sampled with replacement on an independent dataset based on the model constructed on the 
#' original dataset is performed to assess the accuracy of estimated network connections. 
#' A case-dropping bootstrap sampled without replacement is performed to 
#' assess the stability of the calculated centrality indices. 
#' 
#' @param data The independent data frame used for data splitting inference.
#' @param IsingResult The fitted IndivIsing object.
#' @param id The column name in 'data' that indicates the ID of each subject.
#' @param type The type of bootstrap procedure. Can be either "nonparametric" for edge testing or 
#'   "casedropping" for centrality testing.
#' @param sample_prob The vector of sampling proportion for casedropping bootstrap. Default to (0.9,0.8,...0.5).
#' @param nBoots The number of bootstraps. Default to 500.
#' @param nCores Number of cores to use in computing results. Set to 1 to not use parallel computing.
#'   Defaults to 1.
#' @param verbose Logical. Should progress of the function be printed to the console?
#'   Defaults to TRUE.
#' @param library Library location to be used in parallel computing. Default to .libPaths().
#' @param seed The random seed.
#' @return An IndivBootSplitting object.
#' @returns `IsingResult_new` The new IndivIsing model fitted on the inference dataset.
#' @returns `bootResults` A list containing the bootstraped networks.
#' @returns `id` The input id variable.
#' @returns `type` The input type of bootstrap.
#' @returns `sample_prob` The input sampling proportion for casedropping bootstrap.
#' @returns `time` Time spent on running the function.
#' @export 
IndivBootSplitting <- function(data, 
                               IsingResult,
                               id, 
                               type = c("nonparametric", "casedropping"),
                               sample_prob = seq(0.9,0.5,-0.1),
                               nBoots = 500, 
                               nCores = 1, 
                               verbose = TRUE, 
                               library = .libPaths(),
                               seed){

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
  y_index <- IsingResult$y_index
  x_index <- IsingResult$x_index
  
  y_names <- IsingResult$nodes
  x_names <- IsingResult$covar
  
  ##############################################################################
  boot_analysis <- function(data_bootstrap, networkModel){
    
    t00 <- Sys.time()

    y_var <- data_bootstrap[,y_index]  
    x_var <- data_bootstrap[,x_index]  
    
    xy_inter <- do.call(cbind, lapply(1:length(y_index), function(k) y_var[,k]*x_var))
    xy_inter_names <- unlist(lapply(1:length(y_index), function(k) paste0(y_names[k], "-", x_names)))
    colnames(xy_inter) <- xy_inter_names
    var_mat <- cbind(y_var, x_var, xy_inter)
    var_mat <- as.matrix(var_mat)
    
    ##############################################################################
    estimated_bias <- IsingResult$estimated_bias
    estimated_bias_index <- IsingResult$estimated_bias_index
    estimated_formula <- IsingResult$estimated_formula
    estimated_formula_index <- IsingResult$estimated_formula_index
    estimated_coeff_raw <- IsingResult$estimated_coeff_raw
    
    refit_thresholds <- numeric(length(y_index))
    refit_coeff_raw <- matrix(0, nrow=length(y_index), ncol=(length(y_index)+length(x_index)+length(y_index)*length(x_index)))
    colnames(refit_coeff_raw) <- colnames(estimated_coeff_raw)
    
    ##############################################################################
    if(networkModel=="static"){
      
      for(i in 1:nrow(estimated_formula)){
        
        data_model <- as.data.frame(var_mat[,c(i,estimated_bias_index[[i]],estimated_formula_index[[i]])])
        colnames(data_model)[1] <- "outcome"
        
        model <- glm(outcome~., data = data_model, family = "binomial")
        
        refit_thresholds[i] <- summary(model)$coefficients[1,1]
        
        if(ncol(data_model)>=2){
          for(l in 2:ncol(data_model)){
            #refit_coeff_raw[i,which(colnames(refit_coeff_raw) %in% colnames(data_model)[l])] <- summary(model)$coefficients[l,1]
            refit_coeff_raw[i,which(colnames(refit_coeff_raw) %in% colnames(data_model)[l])] <- coef(model, complete = TRUE)[l]
          }  
        }
      }
      
    }else if(networkModel=="temporal"){
      timepoint <- IsingResult$timepoint
      times <- unique(data_bootstrap[,timepoint])  
      
      for(i in 1:nrow(estimated_formula)){
        
        data_model <- as.data.frame(cbind(var_mat[which(data_bootstrap[,timepoint] %in% times[2:length(times)]),i],
                                          var_mat[which(data_bootstrap[,timepoint] %in% times[1:(length(times)-1)]), estimated_formula_index[[i]]]))
        colnames(data_model)[1] <- "outcome"
        
        model <- glm(outcome~., data = data_model, family = "binomial")
        
        refit_thresholds[i] <- summary(model)$coefficients[1,1]
        
        if(ncol(data_model)>=2){
          for(l in 2:ncol(data_model)){
            #refit_coeff_raw[i,which(colnames(refit_coeff_raw) %in% colnames(data_model)[l])] <- summary(model)$coefficients[l,1]
            refit_coeff_raw[i,which(colnames(refit_coeff_raw) %in% colnames(data_model)[l])] <- coef(model, complete = TRUE)[l]
          }  
        }
      }
    }
      
    refit_bias <- matrix(NA, ncol=1, nrow=length(y_index))
    refit_formula <- matrix(NA, nrow = length(y_index), ncol = length(y_index))
    
    if(networkModel=="static"){
      for(i in 1:length(y_index)){
        
        ########################
        formula <- format(round(refit_thresholds[i],3), nsmall=3)
        
        x_var_index <- c((length(y_index)+1):(length(y_index)+length(x_index)))
        
        if(length(which(refit_coeff_raw[i,x_var_index]!=0))>0){
          index <- x_var_index[which(refit_coeff_raw[i,x_var_index]!=0)]
          
          for(k in 1:length(index)){
            formula <- paste0(formula,"+",
                              format(round(refit_coeff_raw[i,index[k]],3), nsmall=3),">>",
                              colnames(refit_coeff_raw)[index[k]])  
          }
        }
        refit_bias[i,1] <- formula
        
        if(networkModel=="static"){
          
        }
        ########################
        for(j in i:length(y_index)){
          if(j>i){
            if(IsingResult$AND == FALSE){
              if(refit_coeff_raw[i,j]!=0 | refit_coeff_raw[j,i]!=0){
                formula <- format(round((refit_coeff_raw[i,j]+refit_coeff_raw[j,i])/2,3), nsmall=3)
              }
              else{
                formula <- NULL
              }
              
              covar_index_i <- c((length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
              covar_index_j <- c((length(y_index)+length(x_index)*i+1):(length(y_index)+length(x_index)*(i+1)))
              
              merge_coeff <- as.vector(as.matrix((refit_coeff_raw[i,covar_index_i]+refit_coeff_raw[j,covar_index_j])/2))
              if(length(which(merge_coeff!=0))>0){
                index <- which(merge_coeff!=0)
                for(k in 1:length(index)){
                  formula <- paste0(formula, "+",
                                    format(round(merge_coeff[index[k]],3), nsmall=3),">>",
                                    x_names[index[k]])
                }
              }
            }
            else{
              if(refit_coeff_raw[i,j]!=0 & refit_coeff_raw[j,i]!=0){
                formula <- as.character(round((refit_coeff_raw[i,j]+refit_coeff_raw[j,i])/2,3))
              }
              else{
                formula <- NULL
              }
              
              covar_index_i <- c((length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
              covar_index_j <- c((length(y_index)+length(x_index)*i+1):(length(y_index)+length(x_index)*(i+1)))
              
              check_and <- sapply(1:length(x_index), function(l) ifelse(refit_coeff_raw[i,covar_index_i[l]]!=0 & refit_coeff_raw[j,covar_index_j[l]]!=0, 1, 0))
              merge_coeff <- as.vector(as.matrix((refit_coeff_raw[i,covar_index_i]+refit_coeff_raw[j,covar_index_j])*check_and/2))
              if(length(which(merge_coeff!=0))>0){
                index <- which(merge_coeff!=0)
                for(k in 1:length(index)){
                  formula <- paste0(formula, "+",
                                    format(round(merge_coeff[index[k]],3), nsmall=3),">>",
                                    x_names[index[k]])
                }
              }
            }
            if(!is.null(formula)){
              refit_formula[i,j] <- formula
              refit_formula[j,i] <- formula
            }
            else{
              refit_formula[i,j] <- NA
              refit_formula[j,i] <- NA
            }
          }
        }
      }
    }else if(networkModel=="temporal"){
      for(i in 1:length(y_index)){
        for(j in 1:length(y_index)){
          
          x_var_index <- c(j,(length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
          
          if(length(which(refit_coeff_raw[i,x_var_index]!=0))>0){
            
            index <- x_var_index[which(refit_coeff_raw[i,x_var_index]!=0)]
            
            if(index[1]==j){
              formula <- format(round(refit_coeff_raw[i,index[1]],3), nsmall=3)
            }
            else{
              formula <- paste0(format(round(refit_coeff_raw[i,index[1]],3), nsmall=3),
                                ">>", x_names[index[1]-length(y_index)-length(x_index)*j])
            }
            
            if(length(index)>1){
              for(k in 2:length(index)){
                formula <- paste0(formula,"+",
                                  format(round(refit_coeff_raw[i,index[k]],3), nsmall=3),">>",
                                  x_names[index[k]-length(y_index)-length(x_index)*j])
              }
            }
            
            refit_formula[i,j] <- formula
          }
        }
      }
    }
    
    rownames(refit_bias) <- y_names
    colnames(refit_formula) <- y_names
    rownames(refit_formula) <- y_names
    
    ##############################################################################
    return(list("nodes" = y_names,
                "covar" = x_names,
                "y_index" = y_index,
                "x_index" = x_index,
                "networkModel" = networkModel,
                "timepoint" = IsingResult$timepoint, 
                "AND" = IsingResult$AND,
                "family" = IsingResult$family,
                "progressbar" = IsingResult$progressbar, 
                "lowerbound.lambda" = IsingResult$lowerbound.lambda, 
                "gamma" = IsingResult$gamma, 
                "estimated_thresholds" = refit_thresholds,
                "estimated_coeff_raw" = refit_coeff_raw,
                "estimated_bias"=refit_bias,
                "estimated_bias_index"=estimated_bias_index,
                "estimated_formula" = refit_formula,
                "estimated_formula_index" = estimated_formula_index,
                "time" = Sys.time() - t00))
  }
  
  ##############################################################################
  IsingResult_new <- boot_analysis(data, IsingResult$networkModel)
  
  ##############################################################################
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
              boot_analysis(bootData, "temporal")
            }else{
              boot_analysis(bootData, "static")
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
                boot_analysis(bootData, "temporal")
              }else{
                boot_analysis(bootData, "static")
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
              boot_analysis(bootData, "temporal")
            }else{
              boot_analysis(bootData, "static")
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
                boot_analysis(bootData, "temporal")
              }else{
                boot_analysis(bootData, "static")
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
    }, mc.cores=nCores)
  }

  return(list("IsingResult_new" = IsingResult_new,
              "bootResults" = bootResults,
              "id" = id,
              "type" = type,
              "sample_prob" = sample_prob,
              "time" = Sys.time()-t0))
  
}
  


  
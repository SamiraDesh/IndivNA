#' Build individual-level static or temporal network for binary data
#' 
#' This function takes the indices of the binary nodes and the indices of the
#' individual characteristics to estimate an individual-level network based on 
#' the extended Ising or logistic autogressive model. The network can be either 
#' static based on cross-sectional data or temporal based on longitudinal data.
#' 
#' @param data The dataset containing the binary nodes and the individual characteristics
#'   for network construction. If a static network is to be constructed, there is one observation
#'   for each subject; if a temporal network is to be constructed, multiple observations 
#'   for each subject is required. Only lag-1 factorization is supported in the current version of
#'   estimating temporal networks.
#' @param y_index The indices of the binary nodes.
#' @param x_index The indices of the individual characteristics. Pre-processing of the 
#'   individual features are required (e.g., scaling of continuous variables, generating 
#'   dummy variables for categorical variables).
#' @param networkModel Model for network estimation. Can be either 'static' or 'temporal'.
#' @param timepoint The column name in 'data' that indicates timepoint for the construction 
#'   of temporal network. If `networkModel='static'`, `timepoint` does not need to be specified. Default to NA.
#' @param family Currently only 'binomial' family of data is supportedn for binary nodes.
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
#'   the OR-rule should be used to define the edges in the network. Defaults to TRUE.
#' @return An IndivIsing object. 
#' @returns `nodes` The names of nodes in the fitted network. 
#' @returns `covar` The names of individual characteristics.
#' @returns `y_index` The input index of the nodes.
#' @returns `x_index` The input index of the individual characteristics. 
#' @returns `networkModel` The input model of the network constructed.
#' @returns `timepoint` The input timepoint variable.
#' @returns `AND` The input AND variable.
#' @returns `family` The input family variable.
#' @returns `progressbar` The input progressbar variable.
#' @returns `lowerbound.lambda` The input lowerbound.lambda variable.
#' @returns `gamma` The input gamma variable.
#' @returns `estimated_thresholds` The estimated intercepts in the extended Ising model.
#' @returns `estimated_coeff_raw` The estimated coefficients in the extended Ising model.
#' @returns `estimated_bias` The formula of the estimated bias term, including the estimated 
#'   intercept and the coefficients of the major effects of the individual characteristics. 
#' @returns `estimated_bias_index` The list of the indices (in the estimated_coeff_raw matrix) of the features 
#'   that are involved in the determination of the bias for each node. 
#' @returns `estimated_formula` The formula of the estimated edge weights, including the 
#'   coefficients of the major effects of the nodes and the coefficients of the interaction terms of the
#'   nodes and the individual characteristics.
#' @returns `estimated_formula_index` The list of the indices (in the estimated_coeff_raw matrix) of the nodes and interaction terms 
#' that are involved in the determination of the edges for each node. 
#' @returns `time` Time used for computation.
#' @export 
IndivIsing <- function(data, 
                       y_index, 
                       x_index, 
                       networkModel = c("static", "temporal"),
                       timepoint = NA, ##name of timepoint variable
                       family = "binomial",
                       progressbar = TRUE, 
                       lowerbound.lambda = NA, 
                       gamma = 0.25, 
                       AND = TRUE){
  
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
  
  ##############################################################################
  ## check symptoms ############################################################
  y_names <- colnames(data)[y_index]
  x_names <- colnames(data)[x_index]
  
  y_var <- data[,y_index]
  x_var <- data[,x_index]
  
  checklognet <- function(y){
    res <- c()
    y <- as.factor(y)
    ntab <- table(y)
    minclass <- min(ntab)
    if(minclass <= 1){
      res <- 0
    }else{
      res <- 1
    }
    return(res)
  }
  
  NodesToAnalyze <- apply(y_var, 2, checklognet) != 0
  names(NodesToAnalyze) <- colnames(y_var)
  if(!any(NodesToAnalyze)){ #not at least one of the values is true
    stop("No variance in dataset.")
  } 
  
  if(any(!NodesToAnalyze)){ 
    warning(paste("Nodes with too little variance (not allowed):",
                  paste(colnames(y_var)[!NodesToAnalyze], collapse = ", ")))
  }
  
  ##############################################################################
  ## create interaction terms ##################################################
  xy_inter <- do.call(cbind, lapply(1:length(y_index), function(k) y_var[,k]*x_var))
  xy_inter_names <- unlist(lapply(1:length(y_index), function(k) paste0(y_names[k], "-", x_names)))
  colnames(xy_inter) <- xy_inter_names
  var_mat <- cbind(y_var, x_var, xy_inter)
  var_mat <- as.matrix(var_mat)
  
  ##############################################################################
  ## model fitting #############################################################
  nvar <- length(y_index)
  if(networkModel=="temporal"){
    p <- nvar+length(x_index)+nvar*length(x_index) # intercept not included
  }else{
    p <- nvar-1+length(x_index)+(nvar-1)*length(x_index) # intercept not included
  }
  
  intercepts <- betas <- lambdas <- vector("list", nvar)
  nlambdas <- rep(0, nvar)
  
  for(i in 1:nvar){
    
    if(networkModel=="temporal"){
      model <- glmnet::glmnet(x=var_mat[which(data[,timepoint] %in% times[1:(length(times)-1)]),],
                              y=as.factor(var_mat[which(data[,timepoint] %in% times[2:length(times)]),i]), 
                              family = "binomial", alpha = 1)
    }else{
      model <- glmnet::glmnet(x=var_mat[,-c(i,(nvar+length(x_index)*i+1):(nvar+length(x_index)*(i+1)))],
                              y=as.factor(var_mat[,i]), family = "binomial", alpha = 1)  
    }
    
    intercepts[[i]] <- model$a0
    betas[[i]] <- model$beta
    lambdas[[i]] <- model$lambda
    nlambdas[i] <- length(lambdas[[i]])
  }
  
  if(progressbar == TRUE){
    pb <- txtProgressBar(max = nvar, style = 3)
  }
  
  J <- matrix(0, nrow=max(nlambdas), ncol=nvar)
  
  for(i in 1:nvar){
    J[1:ncol(betas[[i]]), i] <- colSums(as.matrix(betas[[i]]) != 0)
  }
  
  if(networkModel=="temporal"){
    logl_M <- P_M <- array(0, dim = c((nrow(data)/length(times))*(length(times)-1), max(nlambdas), nvar))
    N <- (nrow(data)/length(times))*(length(times)-1)
  }else{
    logl_M <- P_M <- array(0, dim = c(nrow(data), max(nlambdas), nvar))
    N <- nrow(data)  
  }
  
  for(i in 1:nvar){
    betas.ii <- as.matrix(betas[[i]])
    int.ii <- intercepts[[i]]
    y <- matrix(0, nrow = N, ncol = ncol(betas.ii))
    
    ############################################################################
    if(networkModel=="temporal"){
      xi <- var_mat[which(data[,timepoint] %in% times[1:(length(times)-1)]),]
    }else{
      xi <- var_mat[,-c(i,(nvar+length(x_index)*i+1):(nvar+length(x_index)*(i+1)))]  
    }
    
    ############################################################################
    NB <- nrow(betas.ii) ##number of predictors
    
    for(bb in 1:NB){
      y <- y+betas.ii[rep(bb, N), ]*xi[, bb]
    }
    y <- matrix(int.ii, nrow = N, ncol = ncol(y), byrow = TRUE)+y
    
    n_NA <- max(nlambdas) - ncol(y)
    if(n_NA > 0){
      for(vv in 1:n_NA){
        y <- cbind(y, NA)
      }
    }
    
    if(networkModel=="temporal"){
      P_M[, , i] <- exp(y*var_mat[which(data[,timepoint] %in% times[2:length(times)]),i])/(1+exp(y))
    }else{
      P_M[, , i] <- exp(y*var_mat[,i])/(1+exp(y))  
    }
    logl_M[, , i] <- log(P_M[, , i])
    
    if(progressbar == TRUE){
      setTxtProgressBar(pb, i)
    }
  }
  
  logl_Msum <- colSums(logl_M, 1, na.rm = FALSE)
  
  if(progressbar == TRUE){
    close(pb)
  }
  
  sumlogl <- logl_Msum
  sumlogl[sumlogl == 0] = NA
  
  penalty <- J*log(N)+2*gamma*J*log(p)
  EBIC <- -2*sumlogl+penalty
  lambda.mat <- matrix(NA, nrow(EBIC), ncol(EBIC))
  for(i in 1:nvar){
    lambda.mat[, i] <- c(lambdas[[i]], rep(NA, nrow(EBIC)-length(lambdas[[i]])))
  }
  
  if(!is.na(lowerbound.lambda)){
    EBIC <- EBIC/(lambda.mat >= lowerbound.lambda) * 1
  }
  
  ####################################################################################
  ## summarize results ###############################################################
  lambda.opt <- apply(EBIC, 2, which.min)
  lambda.val <- rep(NA, nvar)
  
  estimated_thresholds <- numeric(length(lambda.opt))  ## intercept
  for (i in 1:length(lambda.opt)){
    lambda.val[i] <- lambda.mat[lambda.opt[i], i]
    estimated_thresholds[i] <- intercepts[[i]][lambda.opt[i]]
  }
  
  estimated_coeff_raw <- matrix(NA, nrow = nvar, ncol = ncol(var_mat))
  for(i in 1:nvar){
    if(networkModel=="temporal"){
      estimated_coeff_raw[i,] <- betas[[i]][,lambda.opt[i]]
    }else{
      estimated_coeff_raw[i, seq(1,ncol(var_mat),1)[-c(i,(nvar+length(x_index)*i+1):(nvar+length(x_index)*(i+1)))]] <- betas[[i]][,lambda.opt[i]]  
    }
  }
  estimated_coeff_raw <- as.data.frame(estimated_coeff_raw)
  colnames(estimated_coeff_raw) <- c(y_names, x_names, xy_inter_names)
  
  ####################################################################################
  estimated_bias <- matrix(NA, ncol=1, nrow=nvar)
  estimated_bias_index <- vector("list", nvar)
  
  estimated_formula <- matrix(NA, nrow = nvar, ncol = nvar)
  estimated_formula_index <- vector("list", nvar)
  
  for(i in 1:nvar){
    
    ########################
    formula <- format(round(estimated_thresholds[i],3), nsmall=3)
    
    x_var_index <- c((nvar+1):(nvar+length(x_index)))
    
    if(length(which(estimated_coeff_raw[i,x_var_index]!=0))>0){
      index <- x_var_index[which(estimated_coeff_raw[i,x_var_index]!=0)]
      
      estimated_bias_index[[i]] <- index
      
      for(k in 1:length(index)){
        formula <- paste0(formula,"+",
                          format(round(estimated_coeff_raw[i,index[k]],3), nsmall=3),">>",
                          colnames(estimated_coeff_raw)[index[k]])  
      }
    }
    estimated_bias[i,1] <- formula
    ########################
    
    if(networkModel=="temporal"){
      
      for(j in 1:nvar){
        
        x_var_index <- c(j,(nvar+length(x_index)*j+1):(nvar+length(x_index)*(j+1)))
        
        if(length(which(estimated_coeff_raw[i,x_var_index]!=0))>0){
          
          index <- x_var_index[which(estimated_coeff_raw[i,x_var_index]!=0)]
          
          estimated_formula_index[[i]] <- c(estimated_formula_index[[i]], index)
          
          if(index[1]==j){
            formula <- format(round(estimated_coeff_raw[i,index[1]],3), nsmall=3)
          }else{
            formula <- paste0(format(round(estimated_coeff_raw[i,index[1]],3), nsmall=3),
                              ">>", x_names[index[1]-length(y_index)-length(x_index)*j])
          }
          
          if(length(index)>1){
            for(k in 2:length(index)){
              formula <- paste0(formula,"+",
                                format(round(estimated_coeff_raw[i,index[k]],3), nsmall=3),">>",
                                x_names[index[k]-length(y_index)-length(x_index)*j])
            }
          }
          
          estimated_formula[i,j] <- formula
        }
      }
    }else{ #if networkModel=="static
      for(j in i:nvar){
        if(j>i){
          if(AND == FALSE){
            if(estimated_coeff_raw[i,j]!=0 | estimated_coeff_raw[j,i]!=0){
              formula <- format(round((estimated_coeff_raw[i,j]+estimated_coeff_raw[j,i])/2,3), nsmall=3)
              estimated_formula_index[[i]] <- c(estimated_formula_index[[i]],j)
              estimated_formula_index[[j]] <- c(estimated_formula_index[[j]],i)
            }else{
              formula <- NULL
            }
            
            covar_index_i <- c((nvar+length(x_index)*j+1):(nvar+length(x_index)*(j+1)))
            covar_index_j <- c((nvar+length(x_index)*i+1):(nvar+length(x_index)*(i+1)))
            
            merge_coeff <- as.vector(as.matrix((estimated_coeff_raw[i,covar_index_i]+estimated_coeff_raw[j,covar_index_j])/2))
            if(length(which(merge_coeff!=0))>0){
              index <- which(merge_coeff!=0)
              
              estimated_formula_index[[i]] <- c(estimated_formula_index[[i]],covar_index_i[index])
              estimated_formula_index[[j]] <- c(estimated_formula_index[[j]],covar_index_j[index])
              
              for(k in 1:length(index)){
                formula <- paste0(formula, "+",
                                  format(round(merge_coeff[index[k]],3), nsmall=3),">>",
                                  x_names[index[k]])
              }
            }
          }else{ #if AND==TRUE
            if(estimated_coeff_raw[i,j]!=0 & estimated_coeff_raw[j,i]!=0){
              formula <- as.character(round((estimated_coeff_raw[i,j]+estimated_coeff_raw[j,i])/2,3))
              estimated_formula_index[[i]] <- c(estimated_formula_index[[i]],j)
              estimated_formula_index[[j]] <- c(estimated_formula_index[[j]],i)
            }else{
              formula <- NULL
            }
            
            covar_index_i <- c((nvar+length(x_index)*j+1):(nvar+length(x_index)*(j+1)))
            covar_index_j <- c((nvar+length(x_index)*i+1):(nvar+length(x_index)*(i+1)))
            
            check_and <- sapply(1:length(x_index), function(l) ifelse(estimated_coeff_raw[i,covar_index_i[l]]!=0 & estimated_coeff_raw[j,covar_index_j[l]]!=0, 1, 0))
            merge_coeff <- as.vector(as.matrix((estimated_coeff_raw[i,covar_index_i]+estimated_coeff_raw[j,covar_index_j])*check_and/2))
            if(length(which(merge_coeff!=0))>0){
              index <- which(merge_coeff!=0)
              
              estimated_formula_index[[i]] <- c(estimated_formula_index[[i]],covar_index_i[index])
              estimated_formula_index[[j]] <- c(estimated_formula_index[[j]],covar_index_j[index])
              
              for(k in 1:length(index)){
                formula <- paste0(formula, "+",
                                  format(round(merge_coeff[index[k]],3), nsmall=3),">>",
                                  x_names[index[k]])
              }
            }
          }
          
          if(!is.null(formula)){
            estimated_formula[i,j] <- formula
            estimated_formula[j,i] <- formula
          }else{
            estimated_formula[i,j] <- NA
            estimated_formula[j,i] <- NA
          }
        }
      }
    }
  }
  
  rownames(estimated_bias) <- y_names
  
  estimated_formula <- as.data.frame(estimated_formula)
  colnames(estimated_formula) <- y_names
  rownames(estimated_formula) <- y_names
  
  return(list("nodes" = y_names,
              "covar" = x_names,
              "y_index" = y_index,
              "x_index" = x_index,
              "networkModel" = networkModel,
              "timepoint" = timepoint, 
              "AND" = AND,
              "family" = family,
              "progressbar" = progressbar, 
              "lowerbound.lambda" = lowerbound.lambda, 
              "gamma" = gamma, 
              "estimated_thresholds" = estimated_thresholds,
              "estimated_coeff_raw" = estimated_coeff_raw,
              "estimated_bias"=estimated_bias,
              "estimated_bias_index"=estimated_bias_index,
              "estimated_formula" = estimated_formula,
              "estimated_formula_index" = estimated_formula_index,
              "time" = Sys.time() - t0))
}




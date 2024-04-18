#' Generate individual static or temporal network with binary data for a specific subject.
#' 
#' This function generates the individual static or temporal network with binary nodes 
#' for a specific subject based on the estimated IndivIsing object. 
#' 
#' @param IsingResult The fitted IndivIsing object.
#' @param data The data frame used to fit the IndivIsing object.
#' @param id The column name in 'data' that indicates the ID of each subject.
#' @param target_id The ID of the specific subject that an individual static or temporal network
#'   is to be fitted. 
#' @param covar_vec The vector of individual characteristics for a specific subject to plot. 
#'   The order has to be the same as x_index. Have to be assigned if target_id is not given.
#' @return An IndivNetwork object. 
#' @returns `estimated_network` The estimated adjacent matrix for the specific subject. 
#' @returns `net_centrality` A list containing the centrality statistics of the generated graph.
#' @returns `node_centrality` The matrix of the calculated centrality values. 
#' @returns `time` Time used for computation.
#' @export
IndivNetwork <- function(IsingResult, 
                         data, 
                         id, 
                         target_id,
                         covar_vec){
  #IsingResult <- result #non case dropping
  #target_id <- 2 (demo)
  # kk <- 3; jj <- 2
  #IsingResult = IsingBootResult$bootResults[[kk]][[jj]] #case dropping, kkth bootstrap and jjth sampling probability
  #covar_vec <- as.vector(as.matrix(sample_covar_mat[1,])) <- numeric(5)
  #covar_vec <- as.vector(as.matrix(sample_covar_mat[2,])) <- feature vector with all entries = 1
  #id <- "id"
  ##############################################################################
  ##############################################################################
  t0 <- Sys.time()
  
  if( missing(target_id) & missing(covar_vec)){
    stop("The target id or the covariate vector has to be given. If both are given, the target id will be used.")
  }
  
  if(!missing(target_id)){
    if(missing(data)){
      stop("The dataset has to be given if the individual network is to be generated using target_id.")
    }
    
    if(missing(id)){
      stop("The id variables has to be given if the individual network is to be generated using target_id")
    }
  }
  
  ##############################################################################
  ##############################################################################
  x_index <- IsingResult$x_index
  
  if(!missing(target_id)){
    cov_data <- as.data.frame(data[data[,id]==target_id, x_index])  
  }else if(!missing(covar_vec)){
    cov_data <- as.data.frame(t(covar_vec))
    colnames(cov_data) <- IsingResult$covar
  }
  
  cov_data$intercept <- 1
  cov_data <- cov_data[,c(length(x_index)+1, 1:length(x_index))]  #(x+1)th column is the intercept
  
  estimated_network_target <- as.data.frame(matrix(0, ncol=ncol(IsingResult$estimated_formula), nrow=nrow(IsingResult$estimated_formula)))
  colnames(estimated_network_target) <- colnames(IsingResult$estimated_formula)
  rownames(estimated_network_target) <- colnames(IsingResult$estimated_formula)
  
  if(is.na(IsingResult$timepoint)){ #static
    for(i in 1:(ncol(estimated_network_target)-1)){
      for(j in (i+1):ncol(estimated_network_target)){
        #IsingResult$estimated_formula[i,j] <- "0.642+1.180>>X3"
        #j<-3; i<-1
        if(!is.na(IsingResult$estimated_formula[i,j]) & !(grepl(">>", IsingResult$estimated_formula[i,j], fixed = TRUE))){ 
          #non zero weight exists and if not affected by a feature
          estimated_network_target[i,j] <- as.numeric(IsingResult$estimated_formula[i,j])
          estimated_network_target[j,i] <- as.numeric(IsingResult$estimated_formula[i,j])
          #non zero weight affected by features
        }else if(!is.na(IsingResult$estimated_formula[i,j]) & grepl(">>", IsingResult$estimated_formula[i,j], fixed = TRUE)){
          factors <- strsplit(IsingResult$estimated_formula[i,j], "+", fixed = TRUE)[[1]]
          cov_coeff <- as.data.frame(do.call(rbind, sapply(1:length(factors), function(k) strsplit(factors[k], ">>"))))
          if(cov_coeff[1,1]==cov_coeff[1,2]){
            cov_coeff[1,2] <- "intercept"
          }
          colnames(cov_coeff) <- c("coeff", "cov")
          cov_coeff$coeff <- as.numeric(cov_coeff$coeff)
          cov_coeff$value <- as.vector(as.matrix(cov_data[,which(colnames(cov_data) %in% cov_coeff$cov)]))
          #get the values of those features that influence edges
          estimated_network_target[i,j] <- sum(cov_coeff$coeff*cov_coeff$value) #coefficient estimate x covariate value
          estimated_network_target[j,i] <- sum(cov_coeff$coeff*cov_coeff$value) #to ensure symmetry
        }
      }
    }
  }else{ #if(!is.na(IsingResult$timepoint)) #temporal
    for(i in 1:ncol(estimated_network_target)){
      for(j in 1:ncol(estimated_network_target)){
        if(!is.na(IsingResult$estimated_formula[i,j]) & !(grepl(">>", IsingResult$estimated_formula[i,j], fixed = TRUE))){
          estimated_network_target[i,j] <- as.numeric(IsingResult$estimated_formula[i,j])
        }else if(!is.na(IsingResult$estimated_formula[i,j]) & grepl(">>", IsingResult$estimated_formula[i,j], fixed = TRUE)){
          factors <- strsplit(IsingResult$estimated_formula[i,j], "+", fixed = TRUE)[[1]]
          cov_coeff <- as.data.frame(do.call(rbind, sapply(1:length(factors), function(k) strsplit(factors[k], ">>"))))
          if(cov_coeff[1,1]==cov_coeff[1,2]){
            cov_coeff[1,2] <- "intercept"
          }
          colnames(cov_coeff) <- c("coeff", "cov")
          cov_coeff$coeff <- as.numeric(cov_coeff$coeff)
          cov_coeff$value <- as.vector(as.matrix(cov_data[1,which(colnames(cov_data) %in% cov_coeff$cov)]))

          estimated_network_target[i,j] <- sum(cov_coeff$coeff*cov_coeff$value)
        }
      }
    }
  }
  #estimated network target is a 10x10 grid showcasing the pairwise edge weights 
  #################################################################################
  net_centrality <- qgraph::centrality_auto(as.matrix(estimated_network_target), weighted = TRUE, signed = TRUE)
  #weighted - Logical, set to FALSE to set all edge weights to 1 or -1
  #signed - Logical, set to FALSE to make all edge weights absolute
  node_centrality <- net_centrality$node.centrality
  
  
  return(list("estimated_network"=estimated_network_target,
              
              "net_centrality"=net_centrality,
              
              "node_centrality"=node_centrality,
              
              "time"=Sys.time() - t0))
  
  
}



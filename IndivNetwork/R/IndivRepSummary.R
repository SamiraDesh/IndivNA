#' Summarize the estimating results based on replicates.
#' 
#' @param IndivRepResult An returned IndivRep object.
#' @param cutoff The cutoff of the probability among all the replicatesbeyond which the edge will be 
#' considered to be included in the final estimated network.
#' @return An IndivRepSummary object.
#' @returns `result_rep` An IndivIsing object that contains the network information estimated based on replicates.
#' @returns `boot_data` A matrix containing all the edge weights estimated in each bootstrap replicate.
#'   Only returned for "nonparametric" bootstrap.
#' @returns `boot_data_summary` A matrix containing summary of all bootstrap results. Only returned 
#'   for "nonparametric" bootstrap. Only returned for "nonparametric" bootstrap.
#' @returns `edge_data` A matrix containing all edge weights retained in the original constructed network that 
#'   estimated in each bootstrap replicate.
#' @returns `edge_data_summary` A matrix containing summary of all edge weights retained in the original constructed network.
#' @returns `edge_test` A matrix listing the test results of whether the edges are significantly different
#'   from each other. 
#' @returns `edge_quantile` A matrix listing the range of the differences between two estimated edges.
#' @returns  `boot_centrality` A matrix summarizing the centrality information for each subject at each bootstrap.
#' #' @returns `test_centrality` A matrix listing the stacked testing results for different centralities between 
#'   nodes for the first representative subject.
#' @returns `test_centrality_long` A matrix with the long format listing the stacked testing results for different centralities between 
#'   nodes for the first representative subject.
#' @returns `cor_final_results` A matrix listing the ranges of the correlation of centralities.
#' @export
IndivRepSummary <- function(IndivRepResult,
                            cutoff = 0.5){
  
  ########################################################################################
  ########################################################################################
  list_estimated_thresholds <- IndivRepResult$list_estimated_thresholds
  list_estimated_coeff_raw <- IndivRepResult$list_estimated_coeff_raw
  list_boot_data <- IndivRepResult$list_boot_data
  list_boot_centrality <- IndivRepResult$list_boot_centrality
  list_cor_final_results <- IndivRepResult$list_cor_final_results
  
  data <- IndivRepResult$data
  y_index <- IndivRepResult$y_index
  x_index <- IndivRepResult$x_index
  networkModel <- IndivRepResult$networkModel
  timepoint <- IndivRepResult$timepoint
  family <- IndivRepResult$family
  progressbar <- IndivRepResult$progressbar
  lowerbound.lambda <- IndivRepResult$lowerbound.lambda
  gamma <- IndivRepResult$gamma
  AND <- IndivRepResult$AND
  rep <- IndivRepResult$rep
  
  ########################################################################################
  ########################################################################################
  mat_estimated_thresholds <- do.call(rbind, list_estimated_thresholds)
  estimated_thresholds <- apply(mat_estimated_thresholds, 2, mean)
  
  mat_estimated_coeff_raw <- do.call(rbind, list_estimated_coeff_raw)
  estimated_coeff_raw <- NULL
  for(i in 1:length(y_index)){
    temp_estimated_coeff_raw <- mat_estimated_coeff_raw[seq(i,nrow(mat_estimated_coeff_raw),by=length(y_index)),]
    above_cutoff <- sapply(1:ncol(temp_estimated_coeff_raw), function(j) length(which(!is.na(temp_estimated_coeff_raw[,j]) & temp_estimated_coeff_raw[,j]!=0)))/rep
    
    estimated_coeff_raw <- rbind(estimated_coeff_raw, apply(temp_estimated_coeff_raw, 2, mean)*ifelse(above_cutoff>=cutoff,1,0))
  }
  
  ########################################################################################
  ########################################################################################
  boot_data <- do.call(rbind, list_boot_data)
  boot_data_summary <- boot_data %>% group_by(edge_names) %>% summarise(mean=mean(edge_weights),
                                                                        lb=quantile(edge_weights, probs=0.025, type=6),
                                                                        ub=quantile(edge_weights, probs=0.975, type=6))
  boot_data_summary <- boot_data_summary[order(boot_data_summary$mean),]
  boot_data_summary$edge_names <- factor(boot_data_summary$edge_names, levels = boot_data_summary$edge_names)
  
  boot_centrality <- do.call(rbind, list_boot_centrality)
  
  ###################################
  ## test results for centralities for the first representative subject
  boot_centrality_id <- boot_centrality[boot_centrality$id==1,]
  
  if(IndivRepResult$networkModel=="static"){
    test_centrality_strength <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_strength) <- y_names
    rownames(test_centrality_strength) <- y_names
    
    test_centrality_expected_influence <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_expected_influence) <- y_names
    rownames(test_centrality_expected_influence) <- y_names
    
    test_centrality_betweenness <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_betweenness) <- y_names
    rownames(test_centrality_betweenness) <- y_names
    
    test_centrality_closeness <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_closeness) <- y_names
    rownames(test_centrality_closeness) <- y_names
    
    for(j in 1:(length(y_index)-1)){
      for(k in (j+1):length(y_index)){
        boot_centrality_id_node1 <- boot_centrality_id[seq(j,nrow(boot_centrality_id),length(y_index)),]
        boot_centrality_id_node2 <- boot_centrality_id[seq(k,nrow(boot_centrality_id),length(y_index)),]
        
        diff_vec_strength <- boot_centrality_id_node1$Strength-boot_centrality_id_node2$Strength
        diff_vec_expected_influence <- boot_centrality_id_node1$ExpectedInfluence-boot_centrality_id_node2$ExpectedInfluence
        diff_vec_betweennesse <- boot_centrality_id_node1$Betweenness-boot_centrality_id_node2$Betweenness
        diff_vec_closeness <- boot_centrality_id_node1$Closeness-boot_centrality_id_node2$Closeness
        
        if(as.vector(quantile(diff_vec_strength, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_strength, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_strength[j,k] <- 1
          test_centrality_strength[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_expected_influence, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_expected_influence, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_expected_influence[j,k] <- 1
          test_centrality_expected_influence[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_betweennesse, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_betweennesse, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_betweenness[j,k] <- 1
          test_centrality_betweenness[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_closeness, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_closeness, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_closeness[j,k] <- 1
          test_centrality_closeness[k,j] <- 1
        }
        
      }
    }
    
    test_centrality_strength_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                "Y"=rep(y_names, times=length(y_names)),
                                                "value"=as.vector(as.matrix(test_centrality_strength)))
    test_centrality_strength_long$value <- factor(test_centrality_strength_long$value)
    
    test_centrality_expected_influence_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                          "Y"=rep(y_names, times=length(y_names)),
                                                          "value"=as.vector(as.matrix(test_centrality_expected_influence)))
    test_centrality_expected_influence_long$value <- factor(test_centrality_expected_influence_long$value)
    
    test_centrality_betweenness_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                   "Y"=rep(y_names, times=length(y_names)),
                                                   "value"=as.vector(as.matrix(test_centrality_betweenness)))
    test_centrality_betweenness_long$value <- factor(test_centrality_betweenness_long$value)
    
    test_centrality_closeness_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                 "Y"=rep(y_names, times=length(y_names)),
                                                 "value"=as.vector(as.matrix(test_centrality_closeness)))
    test_centrality_closeness_long$value <- factor(test_centrality_closeness_long$value)
    
    test_centrality <- as.data.frame(rbind(test_centrality_strength,
                                           test_centrality_expected_influence,
                                           test_centrality_betweenness,
                                           test_centrality_closeness))
    test_centrality$centrality <- c(rep("Strength",length(y_index)),rep("ExpectedInfluence",length(y_index)),
                                    rep("Betweenness",length(y_index)),rep("Closeness",length(y_index)))
    
    test_centrality_long <- as.data.frame(rbind(test_centrality_strength_long,
                                                test_centrality_expected_influence_long,
                                                test_centrality_betweenness_long,
                                                test_centrality_closeness_long))
    test_centrality_long$centrality <- c(rep("Strength",length(y_index)*length(y_index)),rep("ExpectedInfluence",length(y_index)*length(y_index)),
                                         rep("Betweenness",length(y_index)*length(y_index)),rep("Closeness",length(y_index)*length(y_index)))
  }else{
    test_centrality_in_strength <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_in_strength) <- y_names
    rownames(test_centrality_in_strength) <- y_names
    
    test_centrality_out_strength <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_out_strength) <- y_names
    rownames(test_centrality_out_strength) <- y_names
    
    test_centrality_in_expected_influence <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_in_expected_influence) <- y_names
    rownames(test_centrality_in_expected_influence) <- y_names
    
    test_centrality_out_expected_influence <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_out_expected_influence) <- y_names
    rownames(test_centrality_out_expected_influence) <- y_names
    
    test_centrality_betweenness <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_betweenness) <- y_names
    rownames(test_centrality_betweenness) <- y_names
    
    test_centrality_closeness <- matrix(0, nrow=length(y_index), ncol=length(y_index))
    colnames(test_centrality_closeness) <- y_names
    rownames(test_centrality_closeness) <- y_names
    
    for(j in 1:(length(y_index)-1)){
      for(k in (j+1):length(y_index)){
        boot_centrality_id_node1 <- boot_centrality_id[seq(j,nrow(boot_centrality_id),length(y_index)),]
        boot_centrality_id_node2 <- boot_centrality_id[seq(k,nrow(boot_centrality_id),length(y_index)),]
        
        diff_vec_in_strength <- boot_centrality_id_node1$InStrength-boot_centrality_id_node2$InStrength
        diff_vec_out_strength <- boot_centrality_id_node1$OutStrength-boot_centrality_id_node2$OutStrength
        diff_vec_in_expected_influence <- boot_centrality_id_node1$InExpectedInfluence-boot_centrality_id_node2$InExpectedInfluence
        diff_vec_out_expected_influence <- boot_centrality_id_node1$OutExpectedInfluence-boot_centrality_id_node2$OutExpectedInfluence
        diff_vec_betweennesse <- boot_centrality_id_node1$Betweenness-boot_centrality_id_node2$Betweenness
        diff_vec_closeness <- boot_centrality_id_node1$Closeness-boot_centrality_id_node2$Closeness
        
        if(as.vector(quantile(diff_vec_in_strength, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_in_strength, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_in_strength[j,k] <- 1
          test_centrality_in_strength[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_out_strength, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_out_strength, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_out_strength[j,k] <- 1
          test_centrality_out_strength[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_in_expected_influence, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_in_expected_influence, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_in_expected_influence[j,k] <- 1
          test_centrality_in_expected_influence[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_out_expected_influence, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_out_expected_influence, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_out_expected_influence[j,k] <- 1
          test_centrality_out_expected_influence[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_betweennesse, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_betweennesse, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_betweenness[j,k] <- 1
          test_centrality_betweenness[k,j] <- 1
        }
        
        if(as.vector(quantile(diff_vec_closeness, 0.025, type = 6, na.rm = TRUE)*quantile(diff_vec_closeness, 0.975, type = 6, na.rm = TRUE))>0){
          test_centrality_closeness[j,k] <- 1
          test_centrality_closeness[k,j] <- 1
        }
        
      }
    }
    
    test_centrality_in_strength_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                   "Y"=rep(y_names, times=length(y_names)),
                                                   "value"=as.vector(as.matrix(test_centrality_in_strength)))
    test_centrality_in_strength_long$value <- factor(test_centrality_in_strength_long$value)
    
    test_centrality_out_strength_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                    "Y"=rep(y_names, times=length(y_names)),
                                                    "value"=as.vector(as.matrix(test_centrality_out_strength)))
    test_centrality_out_strength_long$value <- factor(test_centrality_out_strength_long$value)
    
    test_centrality_in_expected_influence_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                             "Y"=rep(y_names, times=length(y_names)),
                                                             "value"=as.vector(as.matrix(test_centrality_in_expected_influence)))
    test_centrality_in_expected_influence_long$value <- factor(test_centrality_in_expected_influence_long$value)
    
    test_centrality_out_expected_influence_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                              "Y"=rep(y_names, times=length(y_names)),
                                                              "value"=as.vector(as.matrix(test_centrality_out_expected_influence)))
    test_centrality_out_expected_influence_long$value <- factor(test_centrality_out_expected_influence_long$value)
    
    test_centrality_betweenness_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                   "Y"=rep(y_names, times=length(y_names)),
                                                   "value"=as.vector(as.matrix(test_centrality_betweenness)))
    test_centrality_betweenness_long$value <- factor(test_centrality_betweenness_long$value)
    
    test_centrality_closeness_long <- data.frame("X"=rep(y_names, each=length(y_names)),
                                                 "Y"=rep(y_names, times=length(y_names)),
                                                 "value"=as.vector(as.matrix(test_centrality_closeness)))
    test_centrality_closeness_long$value <- factor(test_centrality_closeness_long$value)
    
    test_centrality <- as.data.frame(rbind(test_centrality_in_strength,
                                           test_centrality_out_strength,
                                           test_centrality_in_expected_influence,
                                           test_centrality_out_expected_influence,
                                           test_centrality_betweenness,
                                           test_centrality_closeness))
    test_centrality$centrality <- c(rep("InStrength",length(y_index)),rep("OutStrength",length(y_index)),
                                    rep("InExpectedInfluence",length(y_index)),rep("OutExpectedInfluence",length(y_index)),
                                    rep("Betweenness",length(y_index)),rep("Closeness",length(y_index)))
    
    test_centrality_long <- as.data.frame(rbind(test_centrality_in_strength_long,
                                                test_centrality_out_strength_long,
                                                test_centrality_in_expected_influence_long,
                                                test_centrality_out_expected_influence_long,
                                                test_centrality_betweenness_long,
                                                test_centrality_closeness_long))
    test_centrality_long$centrality <- c(rep("InStrength",length(y_index)*length(y_index)),rep("OutStrength",length(y_index)*length(y_index)),
                                         rep("InExpectedInfluence",length(y_index)*length(y_index)),rep("OutExpectedInfluence",length(y_index)*length(y_index)),
                                         rep("Betweenness",length(y_index)*length(y_index)),rep("Closeness",length(y_index)*length(y_index)))
  }
  
  ###############################
  cor_final_results <- do.call(rbind, list_cor_final_results)
  cor_final_results <- cor_final_results %>% group_by(centrality, sample_prob) %>% summarise(mean=mean(mean, na.rm=TRUE),
                                                                                             sd=(sqrt(sum(sd^2)))/rep)
  
  ################################################################################
  ################################################################################
  estimated_bias <- matrix(NA, ncol=1, nrow=length(y_index))
  estimated_formula <- matrix(NA, nrow = length(y_index), ncol = length(y_index))
  
  estimated_bias_index <- vector("list", length(y_index))
  estimated_formula_index <- vector("list", length(y_index))
  
  for(i in 1:length(y_index)){
    
    ########################
    formula <- format(round(estimated_thresholds[i],3), nsmall=3)
    
    x_var_index <- c((length(y_index)+1):(length(y_index)+length(x_index)))
    
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
      
      for(j in 1:length(y_index)){
        
        x_var_index <- c(j,(length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
        
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
      for(j in i:length(y_index)){
        if(j>i){
          if(AND == FALSE){
            if(estimated_coeff_raw[i,j]!=0 | estimated_coeff_raw[j,i]!=0){
              formula <- format(round((estimated_coeff_raw[i,j]+estimated_coeff_raw[j,i])/2,3), nsmall=3)
              estimated_formula_index[[i]] <- c(estimated_formula_index[[i]],j)
              estimated_formula_index[[j]] <- c(estimated_formula_index[[j]],i)
            }else{
              formula <- NULL
            }
            
            covar_index_i <- c((length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
            covar_index_j <- c((length(y_index)+length(x_index)*i+1):(length(y_index)+length(x_index)*(i+1)))
            
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
            
            covar_index_i <- c((length(y_index)+length(x_index)*j+1):(length(y_index)+length(x_index)*(j+1)))
            covar_index_j <- c((length(y_index)+length(x_index)*i+1):(length(y_index)+length(x_index)*(i+1)))
            
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
  
  result_rep <- list("nodes" = y_names,
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
                     "estimated_formula_index" = estimated_formula_index)
  
  ##################################################################################################
  ##################################################################################################
  select_edges <- NULL
  if(is.na(timepoint)){
    for(i in 1:(length(y_index)-1)){
      for(j in (i+1):length(y_index)){
        
        if(!is.na(estimated_formula[i,j])){
          str_vec <- strsplit(estimated_formula[i,j], "+", fixed = TRUE)[[1]]
          
          for(k in 1:length(str_vec)){
            if(!grepl(">>", str_vec[k], fixed = TRUE) & str_vec[k]!=""){
              select_edges <- c(select_edges, paste0(colnames(estimated_formula)[i], "-", 
                                                     colnames(estimated_formula)[j]))
            }
            else if(grepl(">>", str_vec[k], fixed = TRUE) & str_vec[k]!=""){
              select_edges <- c(select_edges, paste0(colnames(estimated_formula)[i], "-", 
                                                     colnames(estimated_formula)[j], ":",
                                                     sub(".*>>", "", str_vec[k])))
            }
          }
        }
      }
    }
  }else{
    for(i in 1:length(y_index)){
      for(j in 1:length(y_index)){
        
        if(!is.na(estimated_formula[i,j])){
          str_vec <- strsplit(estimated_formula[i,j], "+", fixed = TRUE)[[1]]
          
          for(k in 1:length(str_vec)){
            if(!grepl(">>", str_vec[k], fixed = TRUE) & str_vec[k]!=""){
              select_edges <- c(select_edges, paste0(colnames(estimated_formula)[i], "-", 
                                                     colnames(estimated_formula)[j]))
            }
            else if(grepl(">>", str_vec[k], fixed = TRUE) & str_vec[k]!=""){
              select_edges <- c(select_edges, paste0(colnames(estimated_formula)[i], "-", 
                                                     colnames(estimated_formula)[j], ":",
                                                     sub(".*>>", "", str_vec[k])))
            }
          } 
        }
      }
    }
  }
  
  edge_data <- boot_data[which(boot_data$edge_names %in% select_edges),]
  edge_data$edge_names <- factor(edge_data$edge_names, levels = unique(edge_data$edge_names))
  edge_data_summary <- boot_data_summary[which(boot_data_summary$edge_names %in% select_edges),]
  edge_data_summary$edge_names <- factor(edge_data_summary$edge_names, levels = unique(edge_data_summary$edge_names))
  
  edge_test <- as.data.frame(matrix(NA, ncol=length(select_edges), nrow=length(select_edges)))
  edge_quantile <- as.data.frame(matrix(NA, ncol=length(select_edges), nrow=length(select_edges)))
  colnames(edge_test) <- unique(edge_data$edge_names)
  rownames(edge_test) <- unique(edge_data$edge_names)
  colnames(edge_quantile) <- unique(edge_data$edge_names)
  rownames(edge_quantile) <- unique(edge_data$edge_names)
  for(i in 1:(length(select_edges)-1)){
    for(j in (i+1):length(select_edges)){
      
      edge_diff <- boot_data$edge_weights[which(boot_data$edge_names==colnames(edge_test)[i])]-
        boot_data$edge_weights[which(boot_data$edge_names==colnames(edge_test)[j])]
      
      edge_diff_lb <- quantile(edge_diff, probs = 0.025, type=6)
      edge_diff_ub <- quantile(edge_diff, probs = 0.975, type=6)
      
      if(edge_diff_lb*edge_diff_ub>0){
        edge_test[i,j] <- edge_test[j,i] <- 1
      }else{
        edge_test[i,j] <- edge_test[j,i] <- 0
      }
      
      edge_quantile[i,j] <- edge_quantile[j,i] <- paste0("(", format(round(edge_diff_lb,2), nsmall=2),",",
                                                         format(round(edge_diff_ub,2), nsmall=2),")")
    }
  }
  
  ##################################################################################################
  ##################################################################################################
  return(list("result_rep"=result_rep,
              "boot_data"=boot_data,
              "boot_data_summary"=boot_data_summary,
              "edge_data"=edge_data,
              "edge_data_summary"=edge_data_summary,
              "edge_test"=edge_test,
              "edge_quantile"=edge_quantile,
              "boot_centrality"=boot_centrality,
              "test_centrality"=test_centrality,
              "test_centrality_long"=test_centrality_long,
              "cor_final_results"=cor_final_results))
}








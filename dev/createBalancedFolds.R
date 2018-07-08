createBalancedFolds <- function(ID, phenotype, gender, nFolds){
  minNo <- min(table(phenotype, gender))
  minNoRd <- floor(minNo/nFolds)*nFolds
  sizeOfFolds_perGroup <- minNoRd/nFolds
  sizeOfFolds <- minNoRd/nFolds*4
  
  
  group0f <- ID[phenotype == FALSE & gender == 1]
  group1f <- ID[phenotype == TRUE & gender == 1]
  group0m <- ID[phenotype == FALSE & gender == 2]
  group1m <- ID[phenotype == TRUE & gender == 2]
  
  group0f_sampled <- sample(group0f, minNoRd, replace = FALSE)
  group1f_sampled <- sample(group1f, minNoRd, replace = FALSE)
  group0m_sampled <- sample(group0m, minNoRd, replace = FALSE)
  group1m_sampled <- sample(group1m, minNoRd, replace = FALSE)
  
  group0f_sampled_indeces <- which(ID %in% group0f_sampled)
  group1f_sampled_indeces <- which(ID %in% group1f_sampled)
  group0m_sampled_indeces <- which(ID %in% group0m_sampled)
  group1m_sampled_indeces <- which(ID %in% group1m_sampled)
  
  folds <- NULL
  folds_0f <- createFoldsNonStrat(group0f_sampled_indeces, nFolds)
  folds_1f <- createFoldsNonStrat(group1f_sampled_indeces, nFolds)
  folds_0m <- createFoldsNonStrat(group0m_sampled_indeces, nFolds)
  folds_1m <- createFoldsNonStrat(group1m_sampled_indeces, nFolds)
  
  for(i in 1:nFolds){
    folds[[i]] <- c(group0f_sampled_indeces[folds_0f[[i]]], 
                    group1f_sampled_indeces[folds_1f[[i]]],
                    group0m_sampled_indeces[folds_0m[[i]]],
                    group1m_sampled_indeces[folds_1m[[i]]])
                   
  }
  return(folds)
}
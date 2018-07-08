createFoldsNonStrat <- function(indexset, nFolds){
  folds <- cvFolds(length(indexset), K=nFolds, type="random")
  df <- data.frame(fold = folds$which, index = folds$subsets)
  non_strat <- lapply(split(df, df$fold), FUN=function(x) x$index)
  return(non_strat)
}
# training data to fool around with
cur_dat = feature_clusters_cv[[kk]][[ii]][pred_ind]
cur_lab = feature_clusters_cv[[kk]][[ii]]$HCPG

# list of functions for rfe
rfRFE <-  list(summary = defaultSummary,
               fit = function(x, y, first, last, ...){
                 library(randomForest)
                 randomForest(x, y, importance = first, ...)
               },
               pred = function(object, x)  predict(object, x),
               rank = function(object, x, y) {
                 vimp <- varImp(object)
                 vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
                 vimp$var <- rownames(vimp)                  
                 vimp
               },
               selectSize = pickSizeBest,
               selectVar = pickVars)

# list of functions for rfe: ridge regression
cur_fit = function(x,y,first,last,cur_family,type_measure,cur_fold_id,useModelFrame,doGrouped) {
  library(glmnetUtils)
  cur_form = as.formula(paste('y ~', paste(names(x),collapse='+')))
  cur_cvmoda = cva.glmnet(cur_form,family=cur_family,
                                       data=x,type.measure=type_measure,
                                       alpha = c(0),foldid = cur_fold_id,use.model.frame=useModelFrame,
                                       grouped=doGrouped)
  return(cur_cvmoda)
}

cur_pred = function(object,x) {
  cur_full_mods_fs     = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
  predict(cur_full_mods_fs$winning_model, newx = x, s = "lambda.min")
}


rrRFE <-  list(summary = defaultSummary,
               
               
               fit = cur_fit,
               
               
               
               
               pred = function(object, x)  predict(object, x),
               rank = function(object, x, y) {
                 vimp <- varImp(object)
                 vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
                 vimp$var <- rownames(vimp)                  
                 vimp
               },
               selectSize = pickSizeBest,
               selectVar = pickVars)

# correlation: kill highly correlated vars
library(mlbench)
library(caret)
library(parallel)
correlationMatrix <- cor(cur_dat)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.90)
# print indexes of highly correlated attributes
cur_dat = cur_dat[-highlyCorrelated]

# recursive feature selection
normalization <- preProcess(cur_dat)
x <- predict(normalization, cur_dat)
x <- as.data.frame(x)

set.seed(10)
subsets <- c(1:5, 10, 20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900)

set.seed(12213)
index <- createFolds(cur_lab, k = 10, returnTrain = T)

cur_form = as.formula(paste('cur_lab ~', paste(names(cur_dat),collapse='+')))
           
cur_cvmoda = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                     data=x,type.measure=type_measure,
                                     alpha = c(0),foldid = cur_fold_id,use.model.frame=useModelFrame,
                                     grouped=doGrouped)
cur_full_mods_fs     = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
cur_surv             = colnames(cur_full_mods_fs$coef)[abs( cur_full_mods_fs$coef) >0]


ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   index = index,
                   verbose = TRUE,
                   allowParallel = F)

lmProfile <- rfe(x, rnorm(length(cur_dat[,1])),
                 sizes = subsets,
                 rfeControl = ctrl)

## caret example ==============================================================
n <- 100
p <- 40
sigma <- 1
set.seed(1)
sim <- mlbench.friedman1(n, sd = sigma)
colnames(sim$x) <- c(paste("real", 1:5, sep = ""),
                     paste("bogus", 1:5, sep = ""))
bogus <- matrix(rnorm(n * p), nrow = n)
colnames(bogus) <- paste("bogus", 5+(1:ncol(bogus)), sep = "")
x <- cbind(sim$x, bogus)
y <- sim$y

normalization <- preProcess(x)
x <- predict(normalization, x)
x <- as.data.frame(x)
subsets <- c(1:5, 10, 15, 20, 25)

set.seed(10)

ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE,
                   allowParallel = T)

cl    = makeCluster(7)
registerDoSNOW(cl)

lmProfile <- rfe(x, y,
                 sizes = subsets,
                 rfeControl = ctrl)

lmProfile

## old ========================================================================
hc = list()
for(ll in 1:length(index)){
  crrltn = cor(x[index[[ll]],])     
  hc[[ll]] = findCorrelation(crrltn, cutoff = .80, names = T, verbose = T)
}


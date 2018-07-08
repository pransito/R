# training data to fool around with
cur_dat = feature_clusters_cv[[kk]][[ii]][pred_ind]
cur_lab = feature_clusters_cv[[kk]][[ii]]$HCPG

# list of functions for rfe: ridge regression
# fit function
cur_fit = function(x,y,first,last,cur_family,type_measure,useModelFrame,doGrouped) {
  library(glmnetUtils)
  cur_form = as.formula(paste('y ~', paste(names(x),collapse='+')))
  cur_cvmoda = cva.glmnet(cur_form,family=cur_family,
                                       data=x,type.measure=type_measure,
                                       alpha = c(0),use.model.frame=useModelFrame,
                                       grouped=doGrouped)
  return(cur_cvmoda)
}

# prediction function
cur_pred = function(object,x) {
  cur_full_mods_fs = agk.glmnetUtils.cvalpha.getwinner(object)
  cur_preds        = predict(cur_full_mods_fs$winning_model, newx = as.matrix(x), s = "lambda.min")
  cur_preds        = as.factor(ifelse(cur_preds < 0, 'HC', 'PG'))
  return(cur_preds)
}

# rank of predictors function
cur_rank = function(object,x,y) {
  cur_full_mods_fs = agk.glmnetUtils.cvalpha.getwinner(object)
  var              = colnames(cur_full_mods_fs$coef)
  dor              = order(abs(cur_full_mods_fs$coef),decreasing = T)
  vimp             = data.frame(var[dor],as.numeric(cur_full_mods_fs$coef)[dor])
  return(vimp)
}

# summary function
cur_summary = function (cur_op_data) {
  out = c(Accuracy = mean(cur_op_data$pred == cur_op_data$obs))
} 


rrRFE <-  list(summary = defaultSummary,
               fit = cur_fit,
               pred = cur_pred,
               rank = cur_rank,
               selectSize = pickSizeBest,
               selectVar = pickVars,
               summary = cur_summary)

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

# this still to be changed!!
set.seed(12213)
index <- createFolds(cur_lab, k = 10, returnTrain = T)

ctrl <- rfeControl(functions = rrRFE,
                   method    = "repeatedcv",
                   repeats   = 1,
                   #indexOut  = cur_fold_id,
                   verbose   = TRUE,
                   allowParallel = F)

# rfe(x,cur_lab,
#                  sizes = subsets,
#                  rfeControl = ctrl,
#                  useModelFrame = useModelFrame,
#                  doGrouped = doGrouped,
#                  cur_family = cur_family,
#                  type_measure = 'class',
#                  metric = 'Accuracy')

rfeIter(x,cur_lab,x,cur_lab,sizes=subsets,rfeControl = ctrl,useModelFrame = useModelFrame,
        doGrouped = doGrouped,
        cur_family = cur_family,
        type_measure = 'class')

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

## other example pickSizeBest
## Example where we would like to maximize
example2 <- data.frame(Rsquared = c(0.4, 0.6, 0.94, 0.95, 0.95, 0.95, 0.95),
                       Variables = 1:7)
## Percent Loss in performance (positive)
example2$PctLoss <- (max(example2$Rsquared) - example2$Rsquared)/max(example2$Rsquared)*100

xyplot(Rsquared ~ Variables, data= example2)
xyplot(PctLoss ~ Variables, data= example2)

absoluteBest2 <- pickSizeBest(example2, metric = "Rsquared", maximize = TRUE)
within5Pct2 <- pickSizeTolerance(example2, metric = "Rsquared", maximize = TRUE)

cat("numerically optimal:",
    example2$Rsquared[absoluteBest2],
    "R^2 in position",
    absoluteBest2, "\n")
cat("Accepting a 1.5 pct loss:",
    example2$Rsquared[within5Pct2],
    "R^2 in position",
    within5Pct2, "\n")

## .Funcs example
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



## old ========================================================================
hc = list()
for(ll in 1:length(index)){
  crrltn = cor(x[index[[ll]],])     
  hc[[ll]] = findCorrelation(crrltn, cutoff = .80, names = T, verbose = T)
}


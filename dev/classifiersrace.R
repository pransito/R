


set.seed(1234)

library(pracma)



## functions ==================================================================
agk.glmnetUtils.cvalpha.getwinner = function(cur_cvmoda) {
  # extract winner from cvalpha (glmnetUtils) mod
  cur_res = list()
  # get winning alpha
  agk.glmnetUtils.cvalpha.getwinner.getalpha = function(cur_cvmoda) {
    return(min(cur_cvmoda$cvm))
  }
  cur_mins               = unlist(lapply(cur_cvmoda$modlist,FUN=agk.glmnetUtils.cvalpha.getwinner.getalpha))
  # preferring higher alphas (more sparse models)
  ind_alpha              = which((cur_mins) == cur_mins)[length(which((cur_mins) == cur_mins))]
  cur_res$winning_alpha  = cur_cvmoda$alpha[ind_alpha]
  
  # now get the best lambda
  cur_res$cvme           = min(cur_cvmoda$modlist[[ind_alpha]]$cvm)
  cur_res$coef           = t(as.matrix(coef(cur_cvmoda$modlist[[ind_alpha]],s='lambda.min')))
  
  # get the winning model (winning alpha cvmod)
  cur_res$winning_model  = cur_cvmoda$modlist[[ind_alpha]]
  # return
  return(cur_res)
}

## Parameters =================================================================
n_noise_feats = 100
n = 1000

## Prep =======================================================================
# build a array to store accuracy of every iteration of the loop
# and build a mean of all accuracys -> compare

############################### for loop start

x1 <- rnorm(n)
x2 <- rnorm(n)

# construct a data frame and scale the normal distributed data points
df <- cbind.data.frame(x1, x2)
df = scale(df)

# noise params
noise_params = randn(n,n_noise_feats)

# for making the groups equally separable
y  = 3 + 5.0*x1 -0.3*x2

# seperate the data into two groups. 
y  = as.numeric(1/(1+exp(-y))>0.5)

yn = y

# take 10% indices and replace the values with not normal distributed data (1's and 0's)
yn[sample(1:length(y))[1:(length(y)*0.1)]] = randi(2,1,length(y)*0.1)-1

# count the numbers of 1's and 0's

table(yn)

df <- data.frame( y,yn, x1, x2) 
df$y=as.factor(df$y)
df$yn=as.factor(df$yn)

# add noise params to data frame
df = data.frame(df,noise_params)

## GLM =======================================================================

# generalized linear model - trying to fit yn into x1+x2
mdl <- glm( yn ~ x1+x2 , data = df , family=binomial)

# get the name of the colums where the data and noisy data is stored and formulate it as a model-string
noise_form = as.formula(paste('yn ~ ', paste(names(df)[grep('y',names(df),invert=T)],collapse = '+')))

# second glm # why still binomial? only yn is binomial, x1,x2,X1,...,X100 are normal distr. data points
mdl_noise <- glm( noise_form, data = df , family=binomial)

# compare both models (glm, glm_noise) (sum up all different classifications)

sum((predict(mdl,type = 'response')>=0.5) == y)/n
sum((predict(mdl_noise,type = 'response')>=0.5) == y)/n

# square the mean of the comparisons between x and y
cost_fun = function(x,y) {
  y = as.numeric(1/(1+exp(-y))>=0.5)
  x = as.numeric(as.character(x)) # get rid of factors?
  return(mean((y==x)^2))
}

# cvTools 
library(cvTools)
# K is number of groups our data is split into
# R number of replications for K-fold cv
# seed is a value to do better randomization
# y is the "true" vector (of classification from our data)
# yn has 10% of non-normal distrubuted values inside -> errors!

# # accuracy for the simple model without noise parameters
# mdl_cv_score       = cvFit(mdl, formula = mdl$formula, data = df, cost = cost_fun, 
#                            K = 10, R = 10, seed = 1234,y=df$y)

# accuracy for the expanded model with noise parameters
mdl_noise_cv_score = cvFit(mdl_noise, formula = mdl_noise$formula, data = df, cost = cost_fun, 
                           K = 10, R = 10, seed = 1234,y=df$y)

## ENet =======================================================================


library(glmnet)
library(glmnetUtils)

# # glmnet for the simple model without noise parameters
# mdl_glmnet = cva.glmnet(mdl$formula,data = df,family='binomial',type.measure='class')
# cur_res          = agk.glmnetUtils.cvalpha.getwinner(mdl_glmnet)
# cur_res          = cur_res$cvme
# # accuracy
# glmnet_res = 1-cur_res

# glmnet for the expanded model with noise parameters 
mdl_noise_glmnet = cva.glmnet(noise_form,data = df,family='binomial',type.measure='class')
cur_res          = agk.glmnetUtils.cvalpha.getwinner(mdl_noise_glmnet)
cur_res          = cur_res$cvme
# accuracy
glmnet_noise_res = 1-cur_res


## SVM =======================================================================

#library (MASS)
library(e1071)
library(caret)
library(LiblineaR)

# # svm for the simple model without noise parameters
# mdl_svm    = best.svm(mdl$formula, data = df, family='binomial',type.measure='class')
# summary(mdl_svm)
# pred <- predict(mdl_svm, type = 'response')
# # compare pred with the "truth" y, first get rid of the factors (levels) and build a ratio
# confusionMatrix(pred,y)

# svm for the expanded model with noise parameters 
mdl_noise_svm    = svm(noise_form, data = df,type.measure='class',kernel = 'linear')
summary(mdl_noise_svm)
pred <- predict(mdl_noise_svm, type = 'response')
# compare pred with the "truth" y, first get rid of the factors (levels) and build a ratio
confusionMatrix(pred,y)

# or to do it by hand (in case u don't like packages)
sum((as.numeric(as.character(pred)) == y))/n

# SVM using LIBLINEAR for; better for linear kernel; svm always puts out gamma; and cannot be transformed
# into linear regression equation with respect to features
installed.packages('kernlab')
cur_data   = df[names(df)[grep('y',names(df),invert=T)]]
cur_target = df$yn 
llsvm_noise = LiblineaR(data = cur_data, target = cur_target,type=2)
help(heuristicC)
heuristicC(as.matrix(cur_data))

# SVM using kernlab
dfkl  = df
dfkl$yn = as.character(dfkl$yn)
dfkl$yn = ifelse(dfkl$yn == '1', 'group1','group2')
dfkl$yn = as.factor(dfkl$yn)
klabsvm_noise = ksvm(noise_form,data = df,kernel='vanilladot')
cvCtrl <- trainControl(method = "repeatedcv", repeats = 1,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
set.seed(1)
rpartTune <- train(noise_form, data = dfkl, method = "svmLinear",
                   tuneLength = 30,tuneGrid = data.frame(C = c(0.1,0.5,1,1.5,6)),
                   metric = "Accuracy",
                   trControl = cvCtrl)

# metric could also be 'ROC'

# svm_tune <- tune(method = svm, mdl_noise_svm, data = df, train.x=y, kernel="radial") ???

#print(svm_tune)


#################### end of loop ##########################




# ## --> cross-validation of this model
# 
# slope <- coef(mdl)[2]/(-coef(mdl)[3])
# intercept <- coef(mdl)[1]/(-coef(mdl)[3])
# 
# 
# 
# library(lattice)
# xyplot( x2 ~ x1 , data = df, groups = y,
#         panel=function(...){
#           panel.xyplot(...)
#           panel.abline(intercept , slope)
#           panel.grid(...)
#         })

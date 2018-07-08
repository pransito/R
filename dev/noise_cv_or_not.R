


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
x1 <- rnorm(n)
x2 <- rnorm(n)

df <- cbind.data.frame(x1, x2)
df = scale(df)


# begin loop
# noise params
noise_params = randn(n,n_noise_feats)

#y <- sign(-1 - 2 * x1 + 2 * x2 )

y  = 3 + 5.0*x1 -0.3*x2 
y  = as.numeric(1/(1+exp(-y))>0.5)
# noise in 10% of the cases
yn = y
yn[sample(1:length(y))[1:(length(y)*0.1)]] = randi(2,1,length(y)*0.1)-1

table(yn)

# y <- sign(x1 + x2) + rnorm(length(y))

#y[ y == -1] <- 0

df <- data.frame( y,yn, x1, x2)
df$y=as.factor(df$y)
df$yn=as.factor(df$yn)

# add noise params
df = data.frame(df,noise_params)

mdl <- glm( yn ~ x1+x2 , data = df , family=binomial)

noise_form = as.formula(paste('yn ~ ', paste(names(df)[grep('y',names(df),invert=T)],collapse = '+')))
mdl_noise <- glm( noise_form, data = df , family=binomial)

sum((predict(mdl,type = 'response')>=0.5) == y)/n
sum((predict(mdl_noise,type = 'response')>=0.5) == y)/n

# cv
cost_fun = function(x,y) {
  y = as.numeric(1/(1+exp(-y))>=0.5)
  x = as.numeric(as.character(x))
  return(mean((y==x)^2))
}
library(cvTools)
mdl_cv_score       = cvFit(mdl, formula = noise_form, data = df, cost = cost_fun, 
                           K = 10, R = 10, seed = 1234,y=df$y)
mdl_noise_cv_score = cvFit(mdl_noise, formula = noise_form, data = df, cost = cost_fun, 
                           K = 10, R = 10, seed = 1234,y=df$y)

library(glmnet)
library(glmnetUtils)
# glmnet
mdl_noise_glmnet = cva.glmnet(noise_form,data = df,family='binomial',type.measure='class')
cur_res          = agk.glmnetUtils.cvalpha.getwinner(mdl_noise_glmnet)
cur_res          = cur_res$cvme
glmnet_noise_res = 1-cur_res

# svm
mdl_noise_svm    = best.svm(noise_form,data = df,family='binomial',type.measure='class')
cur_res          = agk.glmnetUtils.cvalpha.getwinner(mdl_noise_glmnet)
cur_res          = cur_res$cvme
1-cur_res

## end of loop




## --> cross-validation of this model

slope <- coef(mdl)[2]/(-coef(mdl)[3])
intercept <- coef(mdl)[1]/(-coef(mdl)[3])



library(lattice)
xyplot( x2 ~ x1 , data = df, groups = y,
        panel=function(...){
          panel.xyplot(...)
          panel.abline(intercept , slope)
          panel.grid(...)
        })

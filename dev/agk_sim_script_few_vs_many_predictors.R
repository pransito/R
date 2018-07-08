# F-Test; p-values few or many predictors
library(pracma)

few_pred = as.data.frame(randn(10,2))
preds    = names(few_pred)[grep('V1',names(few_pred),invert = T)]
crit     = names(few_pred)[grep('V1',names(few_pred))]
form     = as.formula(paste(crit,'~',paste(preds,collapse = '+')))
summary(lm(formula = form,data=few_pred))

many_pred = as.data.frame(randn(148,140))
preds     = names(many_pred)[grep('V1',names(many_pred),invert = T)]
crit      = names(many_pred)[grep('V1',names(many_pred))]
form      = as.formula(paste(crit,'~',paste(preds,collapse = '+')))
summary(lm(formula = form,data=many_pred))

# adding a -hopefully- orthogonal predictor
few_pred = as.data.frame(randn(10,2))
preds    = names(few_pred)[grep('V1',names(few_pred),invert = T)]
crit     = names(few_pred)[grep('V1',names(few_pred))]
form     = as.formula(paste(crit,'~',paste(preds,collapse = '+')))
summary(lm(formula = form,data=few_pred))

new_pred    = zeros(length(few_pred[,1]),1)
new_pred[1] = 1
new_pred    = few_pred$V2
few_pred_plus = cbind(few_pred,new_pred)
preds    = names(few_pred_plus)[grep('V1',names(few_pred_plus),invert = T)]
crit     = names(few_pred_plus)[grep('V1',names(few_pred_plus))]
form     = as.formula(paste(crit,'~',paste(preds,collapse = '+')))
summary(lm(formula = form,data=few_pred_plus))


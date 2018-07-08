# load stuff
# laod results
setwd('C:/Users/Alexander/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006')
load('MRT_PP_on_MRT_predGrp1_rounds_noo_noaddfeat.RData')

# # load results glmnet fmri only 100 rounds
# setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/90')
# load('MRT_predGrp1_rounds_wio_onlyPhys_no_perm.RData')
# CV_res_list_op_90 = CV_res_list_op
# 
# setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/10')
# load('MRT_predGrp1_rounds_wio_onlyPhys_no_perm.RData')
# CV_res_list_op_10 = CV_res_list_op
# 
# CV_res_list_glmnet = c(CV_res_list_op_10,CV_res_list_op_90)

# load results glmnet fmri only 102 rounds
setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/102')
load('MRT_predGrp1_rounds_wio_onlyPhys_no_perm.RData')
CV_res_list_svm    = CV_res_list_op

# load results glmnet fmri only 101 rounds
setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/103')
load('MRT_predGrp1_rounds_wio_onlyPhys_no_perm.RData')
CV_res_list_glmnet = CV_res_list_op


# load libraries
setwd('C:/Users/Alexander/Google Drive/Library/R')
source('agk_library.R')
agk.load.ifnot.install('Matching')

# functions
get.truth = function() {
  sample(c(rep('HC',3),rep('PG',3)))
}

get.truth.4 = function() {
  sample(c(rep('HC',4),rep('PG',4)))
}

CVres = CVnoo_res_list_PP_on_MRT
#CVres = CV_res_list_glmnet
CVres = CV_res_list_svm
#CVres = CV_res_list
#CVres = CVcm_res_list
#CVres = CV_res_list_op

# set runs
runs = length(CVres)
runs = 1000


# under 0
# pooled
all_aucs  = c()
all_aucsl = list()
all_accs  = c()
all_accsl = list()
all_sens  = c()
all_spec  = c()
for (ii in 1:runs) {
  print(ii)
  inner_truths = c()
  inner_resps  = c()
  # 3
  for (jj in 1:8) {
    # get truth
    inner_truths = c(inner_truths,as.character(get.truth()))
    # get response
    inner_resps  = c(inner_resps,as.numeric(randn(1,6)*10))
  }
  # 4
  for (jj in 9:10) {
    # get truth
    inner_truths = c(inner_truths,as.character(get.truth.4()))
    # get response
    inner_resps  = c(inner_resps,as.numeric(randn(1,8)*10))
  }
  # cur_auc
  cur_roc         = roc(inner_truths,inner_resps)
  all_aucs[ii]    = cur_roc$auc
  all_aucsl[[ii]] = cur_roc
  
  # accuracy
  inner_preds     = ifelse(inner_resps<0,'HC','PG')
  all_accsl[[ii]] = inner_truths == inner_preds
  all_accs[ii]    = mean(all_accsl[[ii]])
  
  # sens and spec
  cur_cm = caret::confusionMatrix(table(inner_truths,inner_preds))
  all_sens[ii]    = cur_cm$byClass[1]
  all_spec[ii]    = cur_cm$byClass[2]
}

# get true aucs
cur_fun   = function(x) {return(pROC::auc(x$roc))}
real_aucs = as.numeric(unlist(lapply(CVres,cur_fun)))
print(mean(real_aucs))
agk.density_p.c(sample(real_aucs -all_aucs,size = 100000,replace = T),0)
ks.boot(real_aucs,all_aucs)

real_aucs_glmnet = as.numeric(unlist(lapply(CV_res_list_glmnet,cur_fun)))
print(mean(real_aucs_glmnet))
agk.density_p.c(sample(real_aucs_glmnet -all_aucs[1:length(CV_res_list_glmnet)],size = 100000,replace = T),0)

real_aucs_svm = as.numeric(unlist(lapply(CV_res_list_svm,cur_fun)))
print(mean(real_aucs_svm))
agk.density_p.c(sample(real_aucs_svm -all_aucs[1:length(CV_res_list_svm)],size = 100000,replace = T),0)


# get true accs
cur_acc_fun = function(x) {return(x$acc[1])}
real_accs   = as.numeric(unlist(lapply(CVres,cur_acc_fun))) 
print(mean(real_accs))
agk.density_p.c(real_accs- all_accs,0)
ks.boot(real_accs,all_accs)

# get true sens
cur_acc_fun = function(x) {return(x$acc[2])}
real_sens   = as.numeric(unlist(lapply(CVres,cur_acc_fun))) 
print(mean(real_sens))
agk.density_p.c(real_sens- all_sens,0)
ks.boot(real_sens,all_sens)

# get true spec
cur_acc_fun = function(x) {return(x$acc[3])}
real_spec   = as.numeric(unlist(lapply(CVres,cur_acc_fun))) 
print(mean(real_spec))
agk.density_p.c(real_spec- all_spec,0)
ks.boot(real_spec,all_spec)


# get the true-false vectors
cur_tf_fun = function(x) {return(x$tp[,1]==x$tp[,2])}
real_tf_l  = lapply(CVres,cur_tf_fun) 

## ROC ========================================================================
cur_fun = function(x) {return(x$roc)}
rocl = lapply(CVres,FUN=cur_fun)

# plot the ROC under real/H0
agk.plot.mean.roc(rocl,des_ci = 'quant',add = F)
agk.plot.mean.roc(all_aucsl,des_ci = 'quant',add=T)

# add legend
title('Receiver-Operating-Curve for fMRI Classifier')
legend(1.02, 1.02, c("ROC of classifier with 95% quantiles\nacross CV rounds", "ROC of 0 classifier with 95% quantiles", "hypothetical null"), col = c('blue', 'red', 'black'),
       text.col = "black", lty = c(1, 2, 4),
       merge = TRUE, bg = NULL,box.lty=0)

## McNamar Test (paired chi-square) ===========================================
# https://www.quora.com/How-can-I-compute-the-p-values-for-the-accuracy-of-two-classifiers-using-k-fold-cross-validation 
mnp = c()
for (mm in 1:length(real_tf_l)) {
  cur_mn = mcnemar.test(all_accsl[[mm]],real_tf_l[[mm]])
  mnp[mm] = cur_mn$p.value
}

## 0-hypothesis: permutation ==================================================
# under 0 hypothesis
# https://stats.stackexchange.com/questions/13834/bootstrap-and-randomization-tests-to-compare-paired-data-sets?rq=1
# https://stats.stackexchange.com/questions/38174/why-is-my-bootstrap-function-for-paired-samples-t-test-in-r-not-returning-the-sa
# I think not applicable, because I am discarding again the variance of the distribution under permutation
# or in other words: if only the algorithm was run enough times then any mean difference would be signficant
# it is only a question to estimate the mean under 0 and under Ha exactly enough; and that is possible with only 
# enough reps
x_1  = real_aucs
x_2  = all_aucs
means_under0 = c()
for (bb in 1:1000) {
  print(bb)
  # paired data
  cur_dat = data.frame(x_1,x_2)
  # randomly shuffling within pair assignment of H0 and Ha
  cur_dat = as.data.frame(t(apply(t(cur_dat),MARGIN = 2,FUN=sample, replace=F)))
  cd = cur_dat[,1] - cur_dat[,2]
  cm = mean(cd)
  means_under0[bb] = cm
}

mean_diff = mean(x_1-x_2)
1-agk.density_p.c(means_under0,mean_diff)


## get a bootstrapped large sample of under_0 and real (matched) ==============
# accuracy
x_1 = sample.int(length(real_accs),size = 1000000,replace = T)
x_2 = all_accs[x_1]
x_1 = real_accs[x_1]

agk.density_p.c(x_1-x_2,0)

# auc
x_1 = sample.int(length(real_aucs),size = 1000000,replace = T)
x_2 = all_aucs[x_1]
x_1 = real_aucs[x_1]

agk.density_p.c(x_1-x_2,0)

# diff under 0
x_1_0  = sample(all_aucs,size = 100000,replace = T)
x_2_0  = sample(all_aucs,size = 100000,replace = T)
diff_0 = x_1_0 - x_2_0

# diff under real  
x_1    = sample.int(length(real_aucs),size = 100000,replace = T)
x_2    = all_aucs[x_1]
x_1    = real_aucs[x_1]
diff_r = x_1 - x_2

agk.density_p.c(diff_r-diff_0,0)
1-agk.density_p.c(diff_0,mean(diff_r))
ks.test(diff_0,diff_r)
## get a bootstrapped estimate of the true accuracy/auc =======================
true_mn_accs = c()
true_mn_aucs = c()
acc_mn_0     = c()
auc_mn_0     = c()
for (bb in 1:2000) {
  print(bb)
  # experiment
  bs_real_aucs  = sample(real_aucs,size=length(real_aucs),replace = T)
  bs_real_accs  = sample(real_accs,size=length(real_accs),replace = T)
  true_mn_accs[bb] = mean(bs_real_accs)
  true_mn_aucs[bb] = mean(bs_real_aucs)
  
  # null
  bs_all_aucs  = sample(all_aucs,size=length(all_aucs),replace = T)
  bs_all_accs  = sample(all_accs,size=length(all_accs),replace = T)
  acc_mn_0[bb] = mean(bs_all_accs)
  auc_mn_0[bb] = mean(bs_all_aucs)
}

agk.density_p.c(true_mn_accs-acc_mn_0,0)
agk.density_p.c(true_mn_aucs-auc_mn_0,0)

# get "real" 0 hypothesis distribution
x_1 = c(all_accs,all_accs,real_accs,real_accs)
x_1 = sample(x_1,size = length(x_1),replace = F)
x_2 = c(real_accs,all_accs,real_accs,all_accs)
x_2 = sample(x_2,size = length(x_2),replace = F)

x_1 = c(all_aucs,all_aucs,real_aucs,real_aucs)
x_1 = sample(x_1,size = length(x_1),replace = F)
x_2 = c(real_aucs,all_aucs,real_aucs,all_aucs)
x_2 = sample(x_1,size = length(x_2),replace = F)

# observed difference
mean(real_accs-all_accs)
mean(real_aucs-all_aucs)

# p-test
1-agk.density_p.c(x_1-x_2,mean(real_accs-all_accs))
1-agk.density_p.c(x_1-x_2,mean(real_aucs-all_aucs))

## use a weighted mean model from PDT behav ===================================
setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006')
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')

# get the weights
mod_weights = as.matrix(table(cur_mod_sel_nooCV))
mod_weights = mod_weights/sum(mod_weights)
win_mods    = row.names(mod_weights)

# winning model beta values
mean_PP_mod = list()
for (mm in 1:length(win_mods)) {
  winning_mod = win_mods[mm]
  win_mods_distr_c = list_winning_model_c_nooCV[which(winning_mod == cur_mod_sel_nooCV)]
  win_mods_distr_l = list_winning_model_l_nooCV[which(winning_mod == cur_mod_sel_nooCV)]
  
  # make a data frame of it:
  win_mod_coefs        = as.matrix(win_mods_distr_c[[1]])
  win_mod_coefs        = as.data.frame(t(win_mod_coefs))
  names(win_mod_coefs) = win_mods_distr_l[[1]]
  
  for (ii in 2:length(win_mods_distr_c)) {
    cur_win_mod_coefs        = as.matrix(t(win_mods_distr_c[[ii]]))
    cur_win_mod_coefs        = as.data.frame(cur_win_mod_coefs)
    names(cur_win_mod_coefs) = win_mods_distr_l[[ii]]
    win_mod_coefs = rbind.fill(win_mod_coefs,cur_win_mod_coefs)
  }
  #imp_0 = function(x) {x[is.na(x)] = 0; return(x)}
  #win_mod_coefs = as.data.frame(lapply(win_mod_coefs,FUN=imp_0))
  mean_PP_mod[[mm]] = colMeans(as.matrix(win_mod_coefs))
}

# get the standardization
mod_codes = c(3,9,2,10)
# first load postpilot data
# setwd("C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006")
# mod_codes = c(3,9,2,10)
# for (mm in 1:length(mod_codes)) {
#   pp_b_dat = featmod_coefs[[mod_codes[mm]]]
#   pp_b_dat = pp_b_dat[,grep('pred',names(pp_b_dat))]
#   pp_b_dat = data.frame(pp_b_dat,pred_smoking_ftdt=dat_match$smoking_ftdt)
#   pp_b_dat = scale(pp_b_dat)
#   save(pp_b_dat, file=paste0('POSTPILOT_',win_mods[mm],'_stand.RData'))
# }

# apply the standardization and get decision value
setwd("C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006")
responses = list()
for (mm in 1:length(win_mods)) {
  
  load(paste0('POSTPILOT_', win_mods[mm],'_stand.RData'))
  pp_scale = attributes(pp_b_dat)
  
  mr_b_dat = featmod_coefs[[mod_codes[mm]]]
  mr_b_dat = mr_b_dat[,grep('pred',names(mr_b_dat))]
  mr_b_dat = data.frame(mr_b_dat,pred_smoking_ftdt = dat_match$smoking_ftdt)
  mr_b_dat = scale(mr_b_dat,center = pp_scale$`scaled:center`, scale = pp_scale$`scaled:scale`)
  mr_b_dat = data.frame(ones(length(mr_b_dat[,1]),1),mr_b_dat)
  
  # prediction
  responses[[mm]] = t(as.matrix(mean_PP_mod[[mm]])) %*% t(as.matrix(mr_b_dat))
}

# consensus (weighted sum of decision values)
weighted_responses = mod_weights[1]*responses[[1]]
for (mm in 2:length(win_mods)) {
  weighted_responses = weighted_responses + mod_weights[mm]*responses[[mm]]
}
preds     = ifelse(weighted_responses <= 0, 'HC','PG')

acc = mean(preds == dat_match$HCPG)
roc = pROC::roc(dat_match$HCPG,predictor=weighted_responses)
auc = roc$auc
cm  = confusionMatrix(preds,dat_match$HCPG)
sen = cm$byClass[1]
spe = cm$byClass[2]

# test
1-agk.density_p.c(all_accs,acc)
1-agk.density_p.c(all_aucs,auc)
1-agk.density_p.c(all_sens,sen)
1-agk.density_p.c(all_spec,spe)

# weighted models auc
weighted_mod_mean_auc = auc

## use a consensus of ALL models from PDT behav ===============================
setwd('C:/Users/Alexander/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1009')
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')

mod_names = c("la","lac","laCh","laec")
#mod_codes = c(3,9,2,10) # for featmod_coefs

# predict with each model
responses = list()
for (mm in 1:length(list_winning_model_c_nooCV)) {
  cur_c = list_winning_model_c_nooCV[[mm]]
  cur_l = list_winning_model_l_nooCV[[mm]]
  cur_m = cur_mod_sel_nooCV[mm]
  
  # apply the standardization and get decision value
  setwd("C:/Users/Alexander/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006")
  
  load(paste0('POSTPILOT_', cur_m,'_stand.RData'))
  pp_scale = attributes(pp_b_dat)
  
  mr_b_dat = featmod_coefs[[cur_m]]
  mr_b_dat = mr_b_dat[,grep('pred',names(mr_b_dat))]
  mr_b_dat = data.frame(mr_b_dat,pred_smoking_ftdt = dat_match$smoking_ftdt)
  mr_b_dat = scale(mr_b_dat,center = pp_scale$`scaled:center`, scale = pp_scale$`scaled:scale`)
  mr_b_dat = data.frame(ones(length(mr_b_dat[,1]),1),mr_b_dat)
  
  # prediction
  responses[[mm]] = t(as.matrix(cur_c)) %*% t(as.matrix(mr_b_dat))
}

# consensus (weighted sum of decision values)
weighted_responses = responses[[1]]
for (mm in 2:length(responses)) {
  weighted_responses = weighted_responses + responses[[mm]]
}
preds     = ifelse(weighted_responses <= 0, 'HC','PG')

acc = mean(preds == dat_match$HCPG)
roc = pROC::roc(dat_match$HCPG,predictor=weighted_responses)
auc = roc$auc
cm  = confusionMatrix(preds,dat_match$HCPG)
sen = cm$byClass[1]
spe = cm$byClass[2]

# test
1-agk.density_p.c(all_accs,acc)
1-agk.density_p.c(all_aucs,auc)
1-agk.density_p.c(all_sens,sen)
1-agk.density_p.c(all_spec,spe)

## density plots ==============================================================
# only auc but multiple classifiers
# old mean_auc = 0.6475 (where is this from)
cur_dat_be = data.frame(H_0 = all_aucs,mean_auc = mean(weighted_mod_mean_auc),classifier = 'prev_behav_glmnet')
cur_dat_gl = data.frame(H_0 = all_aucs,mean_auc = mean(real_aucs_glmnet),classifier = 'MRI_glmnet')
cur_dat_sv = data.frame(H_0 = all_aucs,mean_auc = mean(real_aucs_svm),classifier = 'MRI_svm')

cur_dat              = rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
cur_dat              = melt(cur_dat,id.vars = c('classifier'))
cur_dat_H_0          = subset(cur_dat,variable == 'H_0')
cur_dat_H_0$mean_auc = cur_dat$value[cur_dat$variable == 'mean_auc']
cur_dat              = cur_dat_H_0
p = ggplot(cur_dat,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for different classifiers compared to random classifier')
p = p + geom_vline(aes(xintercept = mean_auc),colour = 'green',size= 1.5)
print(p)

## test glmnet and SVM
# glmnet
1-agk.density_p.c(all_aucs,mean(real_aucs_glmnet))
1-agk.density_p.c(all_aucs,mean(real__glmnet))
1-agk.density_p.c(all_aucs,mean(real_aucs_glmnet))
1-agk.density_p.c(all_aucs,mean(real_aucs_glmnet))
1-agk.density_p.c(all_aucs,mean(real_aucs_svm))


## OLD ========================================================================
# use a mean model from PDT behav
setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/1006')
load('POSTPILOT_HCPG_predGrp1_rounds_noo_noaddfeat.RData')

# winning model beta values with CIs
winning_mod = 'laec'
win_mods_distr_c = list_winning_model_c_nooCV[which(winning_mod == cur_mod_sel_nooCV)]
win_mods_distr_l = list_winning_model_l_nooCV[which(winning_mod == cur_mod_sel_nooCV)]

# make a data frame of it:
win_mod_coefs        = as.matrix(win_mods_distr_c[[1]])
win_mod_coefs        = as.data.frame(t(win_mod_coefs))
names(win_mod_coefs) = win_mods_distr_l[[1]]

for (ii in 2:length(win_mods_distr_c)) {
  cur_win_mod_coefs        = as.matrix(t(win_mods_distr_c[[ii]]))
  cur_win_mod_coefs        = as.data.frame(cur_win_mod_coefs)
  names(cur_win_mod_coefs) = win_mods_distr_l[[ii]]
  win_mod_coefs = rbind.fill(win_mod_coefs,cur_win_mod_coefs)
}
#imp_0 = function(x) {x[is.na(x)] = 0; return(x)}
#win_mod_coefs = as.data.frame(lapply(win_mod_coefs,FUN=imp_0))
mean_PP_mod = colMeans(as.matrix(win_mod_coefs))

# apply the standardization
#pp_b_dat = featmod_coefs[[10]]
#pp_b_dat = pp_b_dat[,grep('pred',names(pp_b_dat))]
#pp_b_dat = data.frame(pp_b_dat,pred_smoking_ftdt=dat_match$smoking_ftdt)
#pp_b_dat = scale(pp_b_dat)
#save(pp_b_dat, file='POSTPILOT_laec_stand.RData')

load('POSTPILOT_laec_stand.RData')
pp_scale = attributes(pp_b_dat)

mr_b_dat = featmod_coefs[[10]]
mr_b_dat = mr_b_dat[,grep('pred',names(mr_b_dat))]
mr_b_dat = data.frame(mr_b_dat,pred_smoking_ftdt = dat_match$smoking_ftdt)
mr_b_dat = scale(mr_b_dat,center = pp_scale$`scaled:center`, scale = pp_scale$`scaled:scale`)

mr_b_dat = data.frame(ones(length(mr_b_dat[,1]),1),mr_b_dat)

# prediction
responses = t(as.matrix(mean_PP_mod)) %*% t(as.matrix(mr_b_dat))
preds     = ifelse(responses <= 0, 'HC','PG')

acc = mean(preds == dat_match$HCPG)
roc = pROC::roc(dat_match$HCPG,predictor=responses)
auc = roc$auc
cm  = confusionMatrix(preds,dat_match$HCPG)
sen = cm$byClass[1]
spe = cm$byClass[2]

# test
1-agk.density_p.c(all_accs,acc)
1-agk.density_p.c(all_aucs,auc)
1-agk.density_p.c(all_sens,sen)
1-agk.density_p.c(all_spec,spe)

# under 0
outer_mn_aucs  = c()
all_inner_aucs = list()
for (ii in 1:200) {
  inner_aucs = c()
  for (jj in 1:10) {
    # get truth
    cur_truth = as.factor(get.truth())
    # get response
    cur_resp  = as.numeric(randn(1,6))
    # cur_auc
    cur_roc = roc(cur_truth,cur_resp)
    inner_aucs[jj] = cur_roc$auc
  }
  all_inner_aucs[[ii]] = inner_aucs
  outer_mn_aucs[ii] = mean(inner_aucs)
}
# two sample t-test
under_0 = all_aucs
all_aucs = unlist(lapply(CV_res_list_op,FUN = cur_fun))
tmp = one.boot(sample(c(all_aucs,under_0)),FUN = mean,R=3000)
1-agk.density_p(tmp$t,mean(all_aucs))

# direct paired
agk.density_p(all_aucs-under_0,0)

# unsorted paired
under_0_paired = function(x) {
  return(sample(x,size = 1) - sample(x,size = 1))
}
tmp = one.boot(sample(c(all_aucs,under_0)),FUN = under_0_paired,R=3000)
1-agk.density_p(tmp$t,median(all_aucs-under_0))

# one at a time case
# each cross validation round needs its own permutation cloud
# each gets a p-value
agk.density_p_l = function(x,under_0) {
  return(1-agk.density_p(under_0,x))
}
all_ps = sapply(all_aucs,FUN=agk.density_p_l,under_0 = under_0)
hist(all_ps)
mean(all_ps<0.05)

## old density plots
# auc, sens, spec, acc
cur_dat_ac = data.frame(H_0 = all_accs,H_a = real_accs,measure = 'acc')
cur_dat_au = data.frame(H_0 = all_aucs,H_a = real_aucs,measure = 'auc')
cur_dat_se = data.frame(H_0 = all_sens,H_a = real_sens,measure = 'sen')
cur_dat_sp = data.frame(H_0 = all_spec,H_a = real_spec,measure = 'spe')
cur_dat    = rbind(cur_dat_au,cur_dat_ac,cur_dat_se,cur_dat_sp)
cur_dat = melt(cur_dat)
p = ggplot(cur_dat,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
p = p + facet_grid(measure ~ .)
print(p)

# only auc but multiple classifiers
cur_dat_be = data.frame(H_0 = all_aucs,H_a = real_aucs,classifier = 'prev_behav_glmnet')
cur_dat_gl = data.frame(H_0 = all_accs[1:length(CV_res_list_glmnet)],H_a = real_aucs_glmnet ,classifier = 'MRI_glmnet')
cur_dat_sv = data.frame(H_0 = all_accs[1:length(CV_res_list_svm)],H_a = real_aucs_svm,classifier = 'MRI_svm')

cur_dat    = rbind(cur_dat_be,cur_dat_gl,cur_dat_sv)
cur_dat    = melt(cur_dat)
p = ggplot(cur_dat,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
p = p + facet_grid(classifier ~ .) + ggtitle('AUC densities for different classifiers compared to random classifier')
print(p)

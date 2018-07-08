# ML with feature elimination
# lasso with many different start foldings
# consensus-based selection

# all data: split in training and test

feat_sel = function(cur_dat,cur_lab,cnfolds) {
  
  # more preprocessing: killing correlated variables
  # throw out big correlation
  correlationMatrix = cor(cur_dat)
  highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.90)
  if (!isempty(highlyCorrelated)) {
    cur_dat           = cur_dat[-highlyCorrelated]
  }

  # repetition of CV for feature elimination
  reps = 20
  
  # make formula
  cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))
  
  # attach
  cur_data      = cur_dat
  cur_data$HCPG = cur_lab
  
  # surviviving predictors
  surv_list = list()
  
  # cvme vector
  cvmes = c()
  
  # ct
  ct = 0
  for (ii in 1:reps){
    
    disp(ii)
    
    # get new folds
    cur_flds    = f_K_fold_b2g(cur_data,'HCPG',cnfolds)
    cur_fold_id = agk.get.foldid(cur_flds)
    
    # the alphas
    #alphas = agk.norm.range(log(seq(0.001,1,length.out = length(cur_dat))),0,1)
    alphas = seq(1,0,length.out = 12)
    
    # fit the model
    cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                                 data=cur_data,type.measure=type_measure,
                                                 alpha = alphas,foldid = cur_fold_id,use.model.frame=useModelFrame,
                                                 grouped=doGrouped,standardize = des_stand_glmnet)
    
    # best model
    winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
    
    # get the coef
    cur_surv = colnames(winning_model$coef)[abs(winning_model$coef)>0.0001]
    if (!length(cur_surv) == 1) {
      ct = ct + 1
      surv_list[[ct]] = cur_surv
    }
    
    # get the cvme
    cvmes[ii] = winning_model$cvme
  }
  
  return(list(surv_list = surv_list, cvmes = cvmes))
}
feat_sel.c   = cmpfun(feat_sel)

# consensus rating
agk.consensus = function(surv_list, prop_crit,min_var) {
  # TODO: return threshold used
  # function to return which of the surviving features
  # appear in at least prop_crit times of the 
  # length of surv_list
  # will reduce prop_crit if not at least min_var feat surviving
  all_feat = unlist(surv_list)
  all_feat = all_feat[-which(duplicated(all_feat))]
  
  cnt_fun = function(x,cur_feat) {return(any(cur_feat == x))}
  
  scores = c()
  for (ff in 1:length(all_feat)) {
    scores[ff] = sum(unlist(lapply(surv_list,cnt_fun,cur_feat = all_feat[ff])))
  }
  
  # survival?
  surv_feats = which(scores/length(surv_list) >= prop_crit)
  
  # selection of survival features
  surv_feats = all_feat[surv_feats]
  
  # kill intercept
  surv_feats = surv_feats[-grep('Intercept',surv_feats)]
  
  # check if enough vars survived
  if (length(surv_feats) >= min_var) {
    return(surv_feats)
  } else {
    agk.consensus(surv_list,prop_crit-0.01,min_var)
  }
  
}

agk.get.mean.model = function(all_winning_models) {
  the_coefs = list()
  for (ll in 1:10) {
    # kick out intercept only models
    # TODO this maybe drop!
    cur_coefs = all_winning_models[[ll]]$coef
    if (all(cur_coefs[2:length(cur_coefs)] == 0)) {next}
    # other than that scale the params to make comparable
    # TODO this maybe drop
    the_coefs[[ll]] = agk.scale_range(cur_coefs,-1,1)
  }
  coef_names = colnames(the_coefs[[1]])
  the_coefs = lapply(the_coefs,as.numeric)
  the_coefs_ma = t(as.matrix(the_coefs[[1]]))
  for (ll in 2:length(the_coefs)) {
    the_coefs_ma = rbind(the_coefs_ma,the_coefs[[ll]])
  }
  colnames(the_coefs_ma) = coef_names
  cur_res = (colMeans(the_coefs_ma))
  return(cur_res)
}


# cur_feat_sel = feat_sel.c(cur_dat,cur_lab,cnfolds)
# cur_surv     = agk.consensus(cur_feat_sel$surv_list,prop_crit = 0.8)
# 
# # test the consensus's cvme
# cut_dat  = cur_dat[cur_surv]
# #test_res = feat_sel (cut_dat,cur_lab,cnfolds)

## build complete model =======================================================
# using convar for complete model?
with_convar          = T
# estimate complete model over different tries of inner folds or just one shot 
pred_from_mean_model = T
mean_model_runs      = 10

# get only sel features
cut_dat  = cur_dat[cur_surv]

if (with_convar) {
  # add the convars
  convars  = dat_match[c('VPPG','smoking_ftdt','edu_hollingshead')]
  
  # training convars
  cur_traini      = flds[[kk]]$train
  cur_train_sub   = cur_lev[cur_traini]
  convars         = convars[which(convars$VPPG %in% cur_train_sub),]
  convars_scal_in = scale(convars[c(2,3)])
  convars[c(2,3)] = convars_scal_in
  cut_dat$VPPG    = row.names(cut_dat)
  cut_dat         = merge(cut_dat,convars,by.x='VPPG',by.y='VPPG')
  cut_dat$VPPG    = NULL
}

# make formula
cur_form = as.formula(paste('HCPG ~ ',paste(names(cut_dat),collapse = ' + ')))

# attach
cur_data      = cut_dat
cur_data$HCPG = cur_lab
alphas        = seq(1,0,length.out = 20)

if (pred_from_mean_model == F){
  mean_model_runs = 1
}

all_winning_models = list()
for (ll in 1:mean_model_runs) {
  print(ll)
  cur_flds      = f_K_fold_b2g(cur_data,'HCPG',cnfolds)
  cur_fold_id   = agk.get.foldid(cur_flds)
  
  
  
  
  cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                               data=cur_data,type.measure=type_measure,
                                               alpha = alphas,foldid = cur_fold_id,use.model.frame=useModelFrame,
                                               grouped=doGrouped,standardize = des_stand_glmnet)
  
  winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
  all_winning_models[[ll]] = winning_model
  
  # always check if covariates are not deselected;
  
}

if (pred_from_mean_model) {
  winning_model = agk.get.mean.model(all_winning_models)
} else {
  winning_model = all_winning_models[[1]]
}

# test data
cur_testi             = flds[[kk]]$test
cur_test_sub          = cur_lev[cur_testi]
test_data            = cr_agg_pp[which(cr_agg_pp$subject %in% cur_test_sub),]
row.names(test_data) = test_data$subject
names(test_data)     = paste0('pred_',names(test_data))
test_data            = test_data[cur_surv]
test_data_us         = test_data

# scaling the test data
cur_scaling_info = feature_clusters_cv_msd[[kk]][[2]]
cur_sd           = cur_scaling_info$`scaled:scale`
cur_m            = cur_scaling_info$`scaled:center`

for(gg in 1:length(cur_sd)) {
  # scaling test data
  cur_var              = names(cur_sd)[gg]
  cur_ind              = grep(cur_var,names(test_data))
  if (!isempty(cur_ind)) {
    cur_scaled_test_data = (test_data[[cur_ind]] - cur_m[gg])/cur_sd[gg]
    # inserting the scaled test data
    test_data[[cur_ind]] = cur_scaled_test_data
  }
}

if (with_convar) {
  # adding the convar data
  convars_test = dat_match[c('VPPG','smoking_ftdt','edu_hollingshead')]
  convars_test = convars_test[which(convars_test$VPPG %in% cur_test_sub),]
  
  # scaling the convar test data
  cur_scaling_info = attributes(convars_scal_in)
  cur_sd           = cur_scaling_info$`scaled:scale`
  cur_m            = cur_scaling_info$`scaled:center`
  
  for(gg in 1:length(cur_sd)) {
    # scaling test data
    cur_var              = names(cur_sd)[gg]
    cur_ind              = grep(cur_var,names(convars_test))
    if (!isempty(cur_ind)) {
      cur_scaled_test_data = (convars_test[[cur_ind]] - cur_m[gg])/cur_sd[gg]
      # inserting the scaled test data
      convars_test[[cur_ind]] = cur_scaled_test_data
    }
  }
  
  # merging test pp data and test convar data
  test_data$VPPG       = row.names(test_data)
  test_data            = merge(test_data,convars_test,by='VPPG')
  row.names(test_data) = test_data$VPPG
  test_data$VPPG       = NULL
}

# predicting
test_data_mat = as.matrix(test_data)
if (pred_from_mean_model) {
  cur_mod       = winning_model
  cur_predic    = t(as.matrix(cur_mod)) %*% t(cbind(ones(length(test_data[,1]),1),test_data))
} else {
  cur_mod       = winning_model$winning_model
  cur_predic    = predict.cv.glmnet(cur_mod, s="lambda.min",newx=test_data_mat)
}
cur_response  = as.numeric(cur_predic)

# evaluate
truth = agk.recode.c(row.names(test_data_mat),dat_match$VPPG,dat_match$HCPG)
preds = ifelse(cur_response > 0, 'PG','HC')

# accuracy
print(mean(truth == preds))

# auc
cur_roc = roc(ifelse(truth == 'PG',1,0),cur_response)
print(cur_roc$auc)


# t.test
# cr_agg_pp_t              = cr_agg_pp
# cr_agg_pp_t$HCPG         = agk.recode.c(cr_agg_pp$subject,dat_match$VPPG,dat_match$HCPG)
# cr_agg_pp_t$smoking_ftnd = 
# cur_t                    = function(cr_agg_pp_t,cur_var,cur_grp) {
#   cur_mod = 
# }


# ML with feature elimination
# according to MR
# however: to speed up things: we stop evaluating all models
# as soon as we get a better CV score

# all data: split in training and test

# use training data
# training data to fool around with
ii          = 2
kk          = 4
cur_dat     = feature_clusters_cv[[kk]][[ii]][pred_ind]
cur_lab     = feature_clusters_cv[[kk]][[ii]]$HCPG
cur_fold_id = cur_fold_id_list[[kk]]

# unit test that stratified folds for inner CV
if (unit_test_strat_innCV) {
  for (zz in unique(cur_fold_id)[order(unique(cur_fold_id))]) {
    check_half = table(feature_clusters_cv[[kk]][[ii]]$HCPG[cur_fold_id == zz])
    if(check_half[1] != check_half[2]) {
      stop('Unbalanced test set for inner CV!')
    }
  }
}

# more preprocessing: killing correlated variables
# throw out big correlation
correlationMatrix = cor(cur_dat)
highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.90)
cur_dat           = cur_dat[-highlyCorrelated]

# make formula
cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))

# attach
cur_data      = cur_dat
cur_data$HCPG = cur_lab

# fit the model
cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                             data=cur_data,type.measure=type_measure,
                                             alpha = c(1),foldid = cur_fold_id,use.model.frame=useModelFrame,
                                             grouped=doGrouped,standardize = des_stand_glmnet)

# note down the best fit
winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)

# RFE LOOP
# for all features:
# eliminate feature p and refit the model
# note down the difference in CV to complete fit

# interim for debug
baseline_cvme = winning_model$cvme
cur_x         = cur_dat
cur_y         = cur_lab

# kernel function for core function
kernel_fun = function(ii,cur_dat_cmpl,baseline_cvme,cur_lab) {
  # cut a feature
  cur_dat = cur_dat_cmpl[-ii]
  
  # fit model
  # make formula
  cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))
  
  # attach
  cur_data      = cur_dat
  cur_data$HCPG = cur_lab
  
  # fit the model
  cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                               data=cur_data,type.measure=type_measure,
                                               alpha = c(1),foldid = cur_fold_id,use.model.frame=useModelFrame,
                                               grouped=doGrouped,standardize = des_stand_glmnet)

  
  # note down decrease in cvme
  res = list()
  winning_model            = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
  res$decrease_in_cvme     = baseline_cvme - winning_model$cvme 
  return(res)
}

# compile
kernel_fun.c = cmpfun(kernel_fun)

# core function which fits once for each feature eliminated
rfe_core = function(cur_x,cur_y,baseline_cvme,kernel_fun.c) {
  cur_dat_cmpl     = cur_x
  decrease_in_cvme = c()
  cur_lab          = cur_y
  
  # work in pieces
  total      = length(cur_dat_cmpl)
  jobs_to_do = (detectCores()-1)*2
  pieces_seq = seq(1,total,jobs_to_do)
  if (pieces_seq[length(pieces_seq)] != total) {pieces_seq = c(pieces_seq,total)}
  pp_seq     = 1:(length(pieces_seq)-1)
  pp_seq     = gtools::permute(pp_seq)
  
  for (pp in pp_seq) {
    # in a parallel for loop
    cl    = makeCluster(detectCores()-1)
    registerDoSNOW(cl)
    cur_seq = pieces_seq[pp]:(pieces_seq[pp+1]-1)
    all_res = foreach(ii=cur_seq, .packages=c('glmnetUtils'),.verbose=F,
                      .export = c('useModelFrame','des_stand_glmnet','doGrouped','cur_fold_id',
                                  'type_measure','cur_family','agk.glmnetUtils.cvalpha.getwinner')) %dopar% 
                                  {
                                    kernel_fun.c(ii,cur_dat_cmpl,baseline_cvme,cur_lab)
                                  } 
    stopCluster(cl)
    
    # check if there is a decrease in cvme in this piece
    res              = list()
    decrease_in_cvme = unlist(all_res)
    cur_max          = max(decrease_in_cvme)
    
    if(cur_max <= 0) {
      # no more increase
      res$ind_to_kill       = NA
      res$new_baseline_cvme = baseline_cvme
      next
    } else {
      ind_to_kill           = which(max(decrease_in_cvme) == decrease_in_cvme)[1]
      ind_to_kill           = ind_to_kill+pieces_seq[pp]-1
      res$ind_to_kill       = ind_to_kill
      res$new_baseline_cvme = baseline_cvme - max(decrease_in_cvme)
      return(res)
    }
  }
  return(res)
}

rfe_core.c = cmpfun(rfe_core)

# outer while loop which will recursively call rfe_core.c
elim = T
while (elim) {
  res = rfe_core.c(cur_x = cur_dat,cur_y = cur_lab,baseline_cvme = baseline_cvme,kernel_fun.c)
  if(is.na(res$ind_to_kill)) {
    disp('Not getting better!')
    break
  } else {
    disp('New baseline_cvme is:')
    disp(res$new_baseline_cvme)
    disp('Killing feature:')
    disp(names(cur_dat)[res$ind_to_kill])
  }
  # updating
  cur_dat       = cur_dat[-res$ind_to_kill]
  baseline_cvme = res$new_baseline_cvme
}

# fitting final model
# make formula
cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))

# attach
cur_data      = cur_dat
cur_data$HCPG = cur_lab

# fit the model
cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                             data=cur_data,type.measure=type_measure,
                                             alpha = c(1),foldid = cur_fold_id,use.model.frame=useModelFrame,
                                             grouped=doGrouped,standardize = des_stand_glmnet)

# note down the best fit
winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)

# formula
winning_coef =coef(winning_model$winning_model,s='lambda.min')
sum(as.matrix(winning_coef)>0)
row.names(winning_coef)[(as.matrix(winning_coef)>0)]





# rank all features from a lot increase to no increase to worsening

# elimate the feature the drop of which leads to best decrease


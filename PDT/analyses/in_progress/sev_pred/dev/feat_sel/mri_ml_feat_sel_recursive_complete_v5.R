# ML with feature elimination
# according to MR
# however: to speed up things: we kick out
# variables with very low betas
# (and anneal alpha)
# use different folds on every round
# collect the surviving features

# all data: split in training and test

feat_sel = function(cur_dat,cur_lab,cnfolds) {

  # criterion length of error not getting smaller
  err_length = 15
  
  # more preprocessing: killing correlated variables
  # throw out big correlation
  # correlationMatrix = cor(cur_dat)
  # highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.90)
  # cur_dat           = cur_dat[-highlyCorrelated]
  
  # make formula
  cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))
  
  # attach
  cur_data      = cur_dat
  cur_data$HCPG = cur_lab
  
  # get new folds
  cur_flds    = f_K_fold_b2g(cur_data,'HCPG',cnfolds)
  cur_fold_id = agk.get.foldid(cur_flds)
  
  # the alphas
  alphas = agk.norm.range(log(seq(0.001,1,length.out = length(cur_dat))),0,1)
  #alphas = seq(1,0,length.out = length(cur_dat))
  ct     = 1
  
  # fit the model
  cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                               data=cur_data,type.measure=type_measure,
                                               alpha = alphas[1],foldid = cur_fold_id,use.model.frame=useModelFrame,
                                               grouped=doGrouped,standardize = des_stand_glmnet)
  
  # note down the best fit
  winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
  
  # print the error
  disp('current error')
  print(winning_model$cvme)
  
  # RFE LOOP
  # for all features:
  # eliminate feature p and refit the model
  # note down the difference in CV to complete fit
  
  # interim for debug
  baseline_cvme = winning_model$cvme
  cur_x         = cur_dat
  cur_y         = cur_lab
  
  # elimination
  do_elim = 1
  
  # error did not get smaller
  error_improve_too_small = c()
  
  while(do_elim) {
    # ...
    disp('running elimination')
    # new round
    ct = ct + 1
    
    # kill features (one of the 0 features)
    vars_to_kill = which(abs(winning_model$coef) == 0)
    if (!isempty(vars_to_kill)) {
      vars_to_kill = sample(vars_to_kill,size =round(length(vars_to_kill)/5))
      #disp('Killing variables:')
      #print(colnames(winning_model$coef)[vars_to_kill])
      cur_dat_bcp = cur_dat
      cur_dat     = cur_dat[-vars_to_kill]
      cur_x       = cur_dat
    } else {
      cur_x       = cur_dat
    }
    
    if (length(cur_dat) < 2) {
      disp('No more variables! Putting backed killed variables and stopping.')
      cur_dat = cur_dat_bcp
      cur_x   = cur_dat
      break
    }
    
    # attach
    cur_data      = cur_x
    cur_data$HCPG = cur_lab
    
    # get new folds
    cur_flds    = f_K_fold_b2g(cur_data,'HCPG',cnfolds)
    cur_fold_id = agk.get.foldid(cur_flds)
    
    # make formula
    cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_x),collapse = ' + ')))
    
    cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                                 data=cur_data,type.measure=type_measure,
                                                 alpha = alphas[ct],foldid = cur_fold_id,use.model.frame=useModelFrame,
                                                 grouped=doGrouped,standardize = des_stand_glmnet)
    
    # note down the best fit
    winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
    
    # print the error
    disp('current error')
    print(winning_model$cvme)
    
    # test if the error got smaller
    if (!isempty(vars_to_kill)) {
      delta_error = (baseline_cvme - winning_model$cvme)
      if (delta_error >= 0.001) {
        disp('Success: error got smaller!')
        baseline_cvme = winning_model$cvme
        # print the error
        print(winning_model$cvme)
        # reset "error not getting smaller vector"
        error_improve_too_small = c()
      } else if ((delta_error >= 0) & (delta_error < 0.001)) {
        disp('too small error improvement; I will put back the cut variables')
        # no improvement in error: then we put the cut variables back
        cur_dat = cur_dat_bcp
        error_improve_too_small = c(error_improve_too_small,1)
        if (length(error_improve_too_small) > err_length) {
          disp('Stopping to look for improvement')
          break
        }
      } else {
        disp('error got bigger; I will put back the cut variables, and try again')
        cur_dat = cur_dat_bcp
      }
    } 
  }
  
  # print the final improvement
  disp('Final CV error:')
  print(winning_model$cvme)
  
  # return
  return(winning_model)
}

all_winning_models = list()
for (oo in 1:5) {
  all_winning_models[[oo]] = feat_sel(cur_dat,cur_lab,cnfolds)
}

# extract cvmes
list_fun      = function(x) {return(x$cvme)}
all_cvme      = unlist(lapply(all_winning_models,list_fun))
disp('mean cvme over feat selection')
print(mean(all_cvme))

# extract features
list_fun = function(cur_mod) {
  cur_surv = colnames(cur_mod$coef)[abs(cur_mod$coef)>0]
  if (!length(cur_surv) == 1) {
    return(cur_surv)
  } else {
    return(NULL)
  }
}
all_surv_feat = lapply(all_winning_models,list_fun)

print(agk.consensus(all_surv_feat,0.6))

# all_winning_models = list()
# early_break = F
# for (oo in 1:6) {
#   cur_winning_model = feat_sel(ii,kk)
#   if (cur_winning_model$cvme <= 0.15) {
#     winning_model = cur_winning_model
#     early_break = T
#     break
#   } else {
#     all_winning_models[[oo]] = cur_winning_model
#   }
# }
# 
# # get the overall final model
# if(early_break == F) {
#   list_fun      = function(x) {return(x$cvme)}
#   all_cvme      = unlist(lapply(all_winning_models,list_fun))
#   winning_model = all_winning_models[[which.min(all_cvme)]]
#   print(winning_model$cvme)
# } else {
#   print(winning_model$cvme)
# }


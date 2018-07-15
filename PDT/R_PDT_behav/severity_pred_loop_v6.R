# PREAMBLE ====================================================================
# Version 6.0
# script to run sev pred CV in a loop to see how stable results are
# and through permutation get a null-hypothesis distribution
# also you can run the nooCV version to get the model that is most likely 
# estimated by the given data
# run import_data and select_study before

# author: Alexander Genauck
# email:  alexander.genauck@charite.de
# date:   15.07.2018

# PREPARATION FOR FOREIGN RUN =================================================
rm(list=ls())
root_wd = 'C:/Users/genaucka/Google Drive/Library/R/PDT/R_PDT_behav'
setwd(root_wd)
load('.RData')

# WHAT TO RUN =================================================================
# just the behavioral parameter sets
outer_cv_noaddfeat_noperm = 1 # only behavior
noout_cv_noaddfeat_noperm = 0 # only behavior

# behavior plus peripheral-physiological stuff
outer_cv_wiaddfeat_noperm = 0 # adding physio
noout_cv_wiaddfeat_noperm = 0 # adding physio

# only peripheral-physiological
outer_cv_addfeaton_noperm = 1 # only physio
noout_cv_addfeaton_noperm = 0 # only physio

outer_cv_c_model_noperm   = 0 # control model/null-model for classification

# what to report
do_report                 = 0
do_report_no_added_feat   = 0
do_report_with_added_feat = 0
do_report_feat_only       = 0

# PARAMETERS TO SET: General ==================================================
# number of runs to get the CV results distribution, 1000 recommended
runs                  = 1010 # physio pred: 300 from home; current correct one for behav: 1010
# set some seed to ensure reproducability
des_seed              = 990 # 990 normally (combine 90 runs with 10 runs)
# run the models (for param extraction in exp)
est_models            = 1
# ridge regression binomial;
# measure for loss function; default is "deviance"; check out options
# help("cv.glmnet")
# in case of gaussian (metric), then 'mse' is used
# 'class' is possible; or 'auc': careful auc needs a certain number of test samples
# I programmed to use k = 5 for innerCV in case of auc
type_measure_binomial = "auc"
# should unit tests for balanced CVfolds be used in wioCV script?
# should always be 1
unit_test_strat_outCV = 1
# should unit tests for balanced inner CVfolds be used in wioCV script?
# should always be 1 (see below!)
unit_test_strat_innCV = T
# should scaling be used for innerCV?
# meaning the innerCV training data will be scaled and the scaling params
# will be applied to test data (exp, add feat, convars)
# this is the optional first step of the algorithm (part of the innerCV test)
scale_for_CV          = T
# initial scale: should the whole data set be scaled first
# note that this is a preprocessing step which IS NOT submitted to innerCV testing
# ideally initial_scale and scale_for_CV should either be both on or both off
initial_scale         = F
# standardize in glmnet? coefficients always returned on orig scale
# default is TRUE
des_stand_glmnet      = T
# stratification of inner CV training and(?) test data
# should be T, ideally, but has been not 
strat_innerCV         = T
# estimate the control model (never except below it gets switched on)
c_mod                 = F
# for report use permut (F) or c_mod?
c_mod_report          = T
# all alphas: then no model selection with just ridge but complete elastic net 
# for all models immediately
all_alphas            = F
# message box width
box_width             = 800
# what predictors to control for
if (which_study == 'MRT') {
  #pred_to_control = c('smoking_ftdt','edu_hollingshead')
  pred_to_control = c()
} else {
  pred_to_control = c('smoking_ftdt')
}
#pred_to_control  = c()
# how many times should full model fit repeated with different folds to get
# mean model? glmnet or svm
fullm_reps            = 10
# will be setautomaticall to 'behav_on_MRI' if you wanna apply PDT behav models
# on MRI sample data; ML for glmnet
# can be set to 'SVM' for SVM with radial basis kernel
which_ML               = 'ML' 

# PARAMETERS TO SET: Behavior =================================================
# add nonlinear behav models?
add_nonlinear_behav     = 0
# use PDT behav model and predict MRT?
# here you should only do nooCV estimate and report
PDT_behav_apply_on_MRI  = 0
# repeat the model selection over many kinds of foldings, take modus
ms_reps                 = 10
# ridge the behavioral models after having been fit? (use the ridged versions)
ridge_behav_models      = F
# ridge the behavioral models after having been fit? (fit anew, otherwise load)
ridge_behav_models_anew = F
# what alphas to use when ridging those behavioral models?
ridge_bm_alphas         = c(0,0.5,1)
# how many repetitions to ensure stable params?
ridge_bv_reps           = 20

# PARAMETERS TO SET: Physio ===================================================
# regress out covs (third option; valid option now)
regress_out_covs     = 0
# feature selection:
kill_high_entr       = 0
# only for additional features (physio); 0 for none, 1 for
# lasso regression, 2 caret's rfe with randomForest
# 3 is simply checking if t-test is yielding sig. group result
# 4: experimental feature selection; changes...
des_feat_sel         = 0
# robust t-test in des_feat_sel = 3?; robust does not seem to help
t_robust             = 0
# if feature selection by t-test: p-value cut off?
# use 1.01 if you want no selection
feat_sel_p           = 0.02
# legacy: could be implemented anew as independent of other feat selection
kill_correlated_vars = 0
# if do con should it be held constant (conservative null-hypothesis)
do_con_constant      = 0
# master add cue reactivity predictors to full model: peripheral physiology
# can be set to 0, if feature selection and adding should not be done on this
add_cr_pp_ma         = 1
# master add cue reactivity predictors to full model: ratings
# can be set to 0, if feature selection and adding should not be done on this
add_cr_ra_ma         = 0
# reduce the fMRI features as set in the severity_pred_init_v6 script
reduce_fMRI_data     = 1
# plot and stats should be turned of for phys and rating params
plot_and_stats       = 0

# PROCESS PREPS ===============================================================
# do not make any settings here
pred_grp = 1
if (strat_innerCV == F) {
  unit_test_strat_innCV = F
}

if (initial_scale == T & scale_for_CV == F) {
  stop('initial_scale is on, then should also be scale_for_CV')
}

# SVM
if(which_ML == 'SVM' & fullm_reps > 3) {
  warning('Since we are using SVM will set fullm_reps to 5')
  fullm_reps = 5
}

# do feature selection set it
do_feat_sel = des_feat_sel
add_cr_pp   = add_cr_pp_ma
add_cr_ra   = add_cr_ra_ma

# run the init of (the CV of) severity pred
source('severity_pred_init_v6.R')

# get the original matching subjects group
sub_grp_matching = featmod_coefs_bcp[[1]][c('HCPG')]
if (add_cr_pp_ma == T) {
  stopifnot(all(row.names(featmod_coefs_bcp[[1]]) == row.names(feature_clusters_bcp[[2]])))
}

# prepare permutations
set.seed(des_seed)
myperms    = list()
myperms_ok = FALSE
while (myperms_ok == 0) {
  for (ii in 1:runs) {
    myperms[[ii]] = gtools::permute(1:length(dat_match[,1]))
  }
  myperms_ok = all(duplicated(myperms) == FALSE)
}

# OUTERCV, NO/WITH ADDED FEATURES (PHYSIO) NOPERM =============================
if (outer_cv_noaddfeat_noperm | outer_cv_wiaddfeat_noperm) {
  # doing the looping for outer crossvalidation, no added phys,
  # without permutation
  set.seed(des_seed)
  # CV style: no outer CV
  CV = 'wio'
  # no additional features (physio)
  do_feat_sel = 0
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  if (outer_cv_wiaddfeat_noperm) {
    add_cr_pp   = 1
    add_cr_ra   = 0
  } else {
    add_cr_pp   = 0
    add_cr_ra   = 0
  }
  # run the models (for param extraction in exp)
  est_models  = 1
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # no permutation
  do_permut   = F
  CV_res_list = list()
  cur_title   = "Rounds outer and inner CV"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = 300)
  list_winning_model = list()
  cur_mod_sel_vec    = c()
  for(hh in 1:runs) {
    source('severity_pred_6_wioCV.R')
    CV_res_list[[hh]]        = CV_res
    setWinProgressBar(pb,hh, title=paste(cur_title, round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  if (outer_cv_wiaddfeat_noperm) {
    save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_wiaddfeat_no_perm.RData'),
         list = c('CV_res_list','fm','des_seed'))
  } else {
    save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_noaddfeat_no_perm.RData'),
         list = c('CV_res_list','fm','des_seed'))
  }

  setwd(cur_home)
}

# OUTERCV, NO ADDED FEATURES (PHYSIO) WITHPERM ================================
if (outer_cv_noaddfeat_wiperm) {
  # doing the looping for outer crossvalidation;
  # with permutation
  set.seed(des_seed)
  # no additional features (physio)
  do_feat_sel = 0
  # CV style: no outer CV
  CV = 'wio'
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  add_cr_pp   = 0
  add_cr_ra   = 0
  # run the models (for param extraction in exp)
  est_models  = 1
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # with permutation
  do_permut    = T
  CVp_res_list = list()
  cur_title    = "Rounds outer and inner CV with permutations"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = box_width)
  for(hh in 1:runs) {
    source('severity_pred_5_wioCV.R')
    CVp_res_list[[hh]]        = CV_res
    setWinProgressBar(pb,hh, title=paste(cur_title,round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  save(file =  paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_noaddfeat_wi_perm.RData'),
       list = c('CVp_res_list','fm','des_seed'))
  setwd(cur_home)
}

# NOOUTERCV, NO ADDED FEATURES (PHYSIO) NOPERM ================================
if (noout_cv_noaddfeat_noperm | noout_cv_wiaddfeat_noperm) {
  # doing the looping for noOuter CV;
  # used for model parameter estimation
  # no permutation
  set.seed(des_seed)
  # no additional features (physio)
  do_feat_sel = 0
  # CV style: no outer CV
  CV = 'noo'
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  if (noout_cv_noaddfeat_noperm) {
    add_cr_pp   = 0
    add_cr_ra   = 0  
  } else {
    add_cr_pp   = 1
    add_cr_ra   = 0 
  }
 
  # run the models (for param extraction in exp)
  est_models  = 1
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  if (PDT_behav_apply_on_MRI) {
    fm_MRT = fm
  }
  CVnoo_res_list           = list()
  CVnoo_res_list_PP_on_MRT = list()
  cur_title = "Rounds no outer but inner CV"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = box_width)
  list_winning_model_c_nooCV  = list()
  list_winning_model_l_nooCV  = list()
  cur_mod_sel_nooCV         = c()
  
  for(hh in 1:runs) {
    source('severity_pred_6_wioCV.R')
    if (PDT_behav_apply_on_MRI) {
      stop('PDT_behav_apply_on_MRI: again have to implement in wioCV script; also should not be used anymore; we do consensus prediction.')
      CVnoo_res_list_PP_on_MRT[[hh]] = CV_res
    } else {
      CVnoo_res_list[[hh]]           = CV_res
      cur_mod_sel_nooCV[hh]          = cur_mod_sel_n
      
      # getting the complete final model's coefs
      full_mod_tmp = full_mod_cv[[1]]
      if (!is.numeric(full_mod_tmp)) {
        cur_labs = colnames(full_mod_tmp$coef)
        cur_labs = gsub(cur_labs,pattern="pred_",replacement="")
        #cur_labs = cur_labs[-1]
        cur_coef = full_mod_tmp$coef
        #cur_coef = cur_coef[-1]
      } else {
        cur_labs = names(full_mod_tmp)
        cur_coef = full_mod_tmp
      }
      
      # packing
      list_winning_model_c_nooCV[[hh]] = cur_coef
      list_winning_model_l_nooCV[[hh]] = cur_labs
    }
    setWinProgressBar(pb,hh, title=paste(cur_title, round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  if (noout_cv_noaddfeat_noperm) {
    if (PDT_behav_apply_on_MRI) {
      save(file = paste0(which_study,'_PP_on_MRT_predGrp',pred_grp,'_rounds_noo_noaddfeat.RData'),
           list = c('CVnoo_res_list_PP_on_MRT','fm','fm_MRT','des_seed'))
    } else {
      save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_noo_noaddfeat.RData'),
           list = c('list_winning_model_c_nooCV','list_winning_model_l_nooCV',
                    'CVnoo_res_list','cur_mod_sel_nooCV','fm','des_seed'))
    }
  } else {
    save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_noo_wiaddfeat.RData'),
         list = c('list_winning_model_c_nooCV','list_winning_model_l_nooCV',
                  'CVnoo_res_list','cur_mod_sel_nooCV','fm','des_seed'))
  }
  setwd(cur_home)
}

# OUTERCV, ADDED FEATURES ONLY (PHYSIO) NOPERM ================================
if (outer_cv_addfeaton_noperm) {
  # doing the looping for outer crossvalidation, ONLY phys,
  # without permutation 
  set.seed(des_seed)
  # only for additional features (physio)
  do_feat_sel      = des_feat_sel
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  add_cr_pp        = add_cr_pp_ma
  add_cr_ra        = add_cr_ra_ma   
  # run the models (for param extraction in exp)
  est_models       = 1
  # add RT (only for behav!)
  add_rt           = 0
  # run the init of (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # only physio
  use_behav_params = F
  # do no permutation
  do_permut      = F
  CV_res_list_op = list()
  cur_title      = "Rounds outer and inner CV only physio"
  pb             = winProgressBar(title = cur_title, min = 0,
                                  max = runs, width = box_width)
  for(hh in 1:runs) {
    source('severity_pred_5_wioCV.R')
    CV_res_list_op[[hh]]        = CV_res
    setWinProgressBar(pb,hh, title=paste(cur_title, round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_onlyPhys_no_perm.RData'),
       list = c('CV_res_list_op','fm','des_seed'))
  setwd(cur_home)
}

# OUTERCV, ADDED FEATURES ONLY (PHYSIO) WITHPERM ==============================
if (outer_cv_addfeaton_wiperm) {
  # doing the looping for outer crossvalidation; ONLY physio;
  # with permutation 
  set.seed(des_seed)
  # only for additional features (physio)
  do_feat_sel = des_feat_sel
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding 
  # should not be done on this
  add_cr_pp   = add_cr_pp_ma
  add_cr_ra   = add_cr_ra_ma
  # run the models (for param extraction in exp)
  est_models  = 1
  # add RT (only for behav!)
  add_rt           = 0
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # only physio
  use_behav_params = F
  # do no permutation
  do_permut      = T
  CVp_res_list_op = list()
  cur_title       = "Rounds outer and inner CV with permutations only physio"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = box_width)
  for(hh in 1:runs) {
    source('severity_pred_6_wioCV.R')
    CVp_res_list_op[[hh]] = CV_res
    setWinProgressBar(pb,hh, title=paste(cur_title,round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive = T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_onlyPhysio_wi_perm.RData'),
       list = c('CVp_res_list_op','fm','des_seed'))
  setwd(cur_home)
}

# NOOUTERCV, ADDED FEATURES ONLY (PHYSIO) NOPERM ==============================
if (noout_cv_addfeaton_noperm) {
  # doing the looping for noOuter CV; only physio
  # used for model parameter estimation
  # no permutation
  set.seed(des_seed)
  # only for additional features (physio)
  do_feat_sel      = des_feat_sel
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  add_cr_pp        = add_cr_pp_ma
  add_cr_ra        = add_cr_ra_ma
  # run the models (for param extraction in exp)
  est_models       = 1
  # add RT (only for behav!)
  add_rt           = 0
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # only physio
  use_behav_params = F
  CVnoo_res_list_op              = list()
  list_winning_model_c_nooCV_op  = list()
  list_winning_model_l_nooCV_op  = list()
  cur_mod_sel_nooCV_op           = c()
  cur_title         = "Rounds no outer but inner CV (only physio)"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = box_width)
  for(hh in 1:runs) {
    source('severity_pred_4_nooCV.R')
    CVnoo_res_list_op[[hh]]  = CV_res
    
    # getting the complete final model's coefs
    cur_labs = colnames(full_mod_tmp$coef)
    cur_labs = gsub(cur_labs,pattern="pred_",replacement="")
    cur_labs = cur_labs[-1]
    cur_coef = full_mod_tmp$coef
    cur_coef = cur_coef[-1]
    
    # TODO:
    # inspection of the additional feature selection results
    
    # packing
    list_winning_model_c_nooCV_op[[hh]] = cur_coef
    list_winning_model_l_nooCV_op[[hh]] = cur_labs
    setWinProgressBar(pb,hh, title=paste(cur_title, round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_noo_onlyPhysio.RData'),
       list = c('list_winning_model_c_nooCV_op','list_winning_model_l_nooCV_op',
                'CVnoo_res_list_op','fm','des_seed'))
  setwd(cur_home)
}

# OUTERCV, CONTROL MODEL NOPERM ===============================================
# alternative to permutation
# the control model; intercept only, or only the control variables
if (outer_cv_c_model_noperm) {
  # doing the looping for outer crossvalidation, no added phys,
  # without permutation
  set.seed(des_seed)
  # no additional features (physio)
  do_feat_sel = 0
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  add_cr_pp   = 0
  add_cr_ra   = 0
  # run the models (for param extraction in exp)
  est_models  = 1
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # no permutation
  do_permut   = F
  # CV style
  CV = 'wio'
  # control model?
  c_mod       = T
  CVcm_res_list = list()
  cur_title   = "Rounds outer and inner CV"
  pb = winProgressBar(title = cur_title, min = 0,
                      max = runs, width = 300)
  list_winning_model = list()
  for(hh in 1:runs) {
    source('severity_pred_6_wioCV.R')
    CVcm_res_list[[hh]]        = CV_res
    setWinProgressBar(pb,hh, title=paste(cur_title, round(hh/runs*100),
                                         "% done"))
  }
  close(pb)
  
  # saving
  cur_home = getwd()
  dir.create(file.path(cur_home, paste0('results/',runs)),recursive=T)
  setwd(file.path(cur_home, paste0('results/',runs)))
  save(file = paste0(which_study,'_predGrp',pred_grp,'_rounds_wio_conmod_no_perm.RData'),
       list = c('CVcm_res_list','fm','des_seed'))
  setwd(cur_home)
}

# REPORTING: PREPARATION =====================================================
if (do_report) {
  cur_home = getwd()
  setwd(root_wd)
  setwd('results')
  setwd(as.character(runs))
  if (pred_grp) {
    all_res_files = dir(pattern = paste0(which_study,'_predGrp1'))
  } else {
    all_res_files = dir(pattern = paste0(which_study,'_predGrp0'))
  }
  
  if (do_report_with_added_feat) {
    all_res_files_Ha = all_res_files[grep('wiaddfeat',all_res_files)]
  } else if (do_report_no_added_feat) {
    all_res_files_Ha = all_res_files[grep('noaddfeat',all_res_files)]
  } else {
    stop('Ha not defined for report.')
  }
  
  for (ii in 1:length(all_res_files)) {
    load(all_res_files[ii])
  }
  
  # over-load the correct Ha
  for (ii in 1:length(all_res_files_Ha)) {
    load(all_res_files_Ha[ii])
  }
  
  setwd(cur_home)
  if (c_mod_report) {
    # using control model instead of permutation
    CVp_res_list    = CVcm_res_list
    CVp_res_list_op = CVcm_res_list
  }
}



# REPORTING: NO ADDED PHYS ====================================================
if (do_report_no_added_feat | do_report_with_added_feat) {
  # p-value for the algorithm
  # get the accuracy null-distrib/alt hypothesis distribution
  cor    = c() # correlation
  cor_p  = c() # correlation permuted
  coru   = c() # correlation unpooled (mean of k-fold)
  coru_p = c() # correlation unpooled (mean of k-fold) permuted
  mse    = c() # mean squared error
  mse_p  = c() # mean squared error permuted
  accs_p = c()
  accs   = c()
  sens_p = c()
  sens   = c()
  spec_p = c()
  spec   = c()
  auc_p  = c()
  auc    = c()
  roc_pl = list()
  rocl   = list()
  for (ii in 1:length(CVp_res_list)) {
    accs_p[ii]   = CVp_res_list[[ii]]$acc$accuracy
    accs[ii]     = CV_res_list[[ii]]$acc$accuracy
    sens_p[ii]   = CVp_res_list[[ii]]$acc$cur_sens
    sens[ii]     = CV_res_list[[ii]]$acc$cur_sens
    spec_p[ii]   = CVp_res_list[[ii]]$acc$cur_spec
    spec[ii]     = CV_res_list[[ii]]$acc$cur_spec
    roc_pl[[ii]] = CVp_res_list[[ii]]$roc
    auc_p[ii]    = as.numeric(roc_pl[[ii]]$auc)
    rocl[[ii]]   = CV_res_list[[ii]]$roc
    auc[ii]      = as.numeric(rocl[[ii]]$auc)
    cor[[ii]]    = CV_res_list[[ii]]$cor$estimate
    cor_p[[ii]]  = CVp_res_list[[ii]]$cor$estimate
    coru[[ii]]   = CV_res_list[[ii]]$coru
    coru_p[[ii]] = CVp_res_list[[ii]]$coru
    mse[[ii]]    = CV_res_list[[ii]]$mse
    mse_p[[ii]]  = CVp_res_list[[ii]]$mse
  }
  
  # make new ROC objects
  for (ii in 1:length(rocl)) {
    # non_permuted
    predictor  = rocl[[ii]]$original.predictor
    response   = rocl[[ii]]$original.response
    response   = ifelse(response == 'PG',1,0)
    predictor  = agk.scale_range(predictor,0,1)
    rocl[[ii]] = roc(response,predictor)
    auc[ii]    = auc(rocl[[ii]])
    rocl[[ii]] = smooth( rocl[[ii]],n = 500)
    
    # permuted (or control model)
    predictor    = roc_pl[[ii]]$original.predictor
    response     = roc_pl[[ii]]$original.response
    response     = ifelse(response == 'PG',1,0)
    predictor    = agk.scale_range(predictor,0,1)
    roc_pl[[ii]] = roc(response,predictor)
    auc_p[ii]    = auc(roc_pl[[ii]])
    
    if (var(roc_pl[[ii]]$original.predictor) == 0) {
      cur_obj = list()
      cur_obj$sensitivities = seq(0,1,length.out = 502)
      cur_obj$specificities = seq(0,1,length.out = 502)
      roc_pl[[ii]]          = cur_obj
    } else {
      roc_pl[[ii]] = smooth(roc_pl[[ii]],n = 500)
    }
  }
  
  # print the median of acc, sens, spec
  acc_sens_spec       = data.frame(accs,sens,spec,auc)
  acc_sens_spec_p     = data.frame(accs_p,sens_p,spec_p,auc_p)
  med_ac_se_sp        = lapply(acc_sens_spec,FUN=mean)
  med_ac_se_sp_np_ci  = lapply(acc_sens_spec,FUN=agk.boot.ci,
                               cur_fun=mean,R=3000,lower=0.025,upper=0.975)
  med_ac_se_sp_npp_ci = lapply(acc_sens_spec_p,FUN=agk.boot.ci,
                               cur_fun=mean,R=3000,lower=0.025,upper=0.975)
  disp('Mean accuracy, sensitivity, specificity across CV rounds.')
  print(med_ac_se_sp)
  # acc
  diffs        = accs - accs_p
  cur_p        = agk.density_p.c(diffs,0)
  disp("Probability that permuted and non-permuted CV'd accuracy of algorithms are the same.")
  print(cur_p)
  # sens
  diffs        = sens - sens_p
  cur_p        = agk.density_p.c(diffs,0)
  disp("Probability that permuted and non-permuted CV'd sensitivity of algorithms are the same.")
  print(cur_p)
  # spec
  diffs        = spec - spec_p
  cur_p        = agk.density_p(diffs,0)
  disp("Probability that permuted and non-permuted CV'd specificity of algorithms are the same.")
  print(cur_p)
  # auc
  diffs        = auc - auc_p
  cur_p        = agk.density_p.c(diffs,0)
  disp("Probability that permuted and non-permuted CV'd smoothed ROC-AUC of algorithms are the same.")
  print(cur_p)
  
  # ROC curve (get all the coordinates)
  # not perm
  specificities  = t(as.matrix(rocl[[1]]$specificities))
  sensitivities  = t(as.matrix(rocl[[1]]$sensitivities))
  
  # perm
  specificitiesl = t(as.matrix(roc_pl[[1]]$specificities))
  sensitivitiesl = t(as.matrix(roc_pl[[1]]$sensitivities))
  
  # get all the ROC curves
  for (ii in 2:length(rocl)) {
    if (length(rocl[[ii]]$specificities) != length(specificities[ii-1,])) {
      stop('length not same!')
    }
    
    if (length(roc_pl[[ii]]$specificities) != length(specificitiesl[ii-1,])) {
      stop('length not same!')
    }
    
    specificities = rbind(specificities,t(as.matrix(rocl[[ii]]$specificities)))
    sensitivities = rbind(sensitivities,t(as.matrix(rocl[[ii]]$sensitivities)))
    
    specificitiesl = rbind(specificitiesl,t(as.matrix(roc_pl[[ii]]$specificities)))
    sensitivitiesl = rbind(sensitivitiesl,t(as.matrix(roc_pl[[ii]]$sensitivities)))
  }
  
  # get CI
  specificities  = as.data.frame(specificities)
  sensitivities  = as.data.frame(sensitivities)
  specificitiesl = as.data.frame(specificitiesl)
  sensitivitiesl = as.data.frame(sensitivitiesl)
  #spec_min_specl = specificities - specificitiesl
  #sens_min_sensl = sensitivities - sensitivitiesl
  
  # spec_mci  = lapply(specificities,FUN = agk.boot.ci.c,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  # sens_mci  = lapply(sensitivities,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  # specl_mci = lapply(specificitiesl,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  # sensl_mci = lapply(sensitivitiesl,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  
  # not ci of mean but mean and percentiles over CV rounds
  agk.mean.quantile.c = cmpfun(agk.mean.quantile)
  spec_mci  = lapply(specificities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  sens_mci  = lapply(sensitivities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  specl_mci = lapply(specificitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  sensl_mci = lapply(sensitivitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  
  # # subtract the permuted version
  # spemspel_mci = lapply(spec_min_specl,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  # senmsenl_mci = lapply(sens_min_sensl,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  
  spec_mci     = as.data.frame(matrix(unlist(spec_mci),ncol = 3,byrow = T))
  sens_mci     = as.data.frame(matrix(unlist(sens_mci),ncol = 3,byrow = T))
  specl_mci    = as.data.frame(matrix(unlist(specl_mci),ncol = 3,byrow = T))
  sensl_mci    = as.data.frame(matrix(unlist(sensl_mci),ncol = 3,byrow = T))
  # spemspel_mci = as.data.frame(matrix(unlist(spemspel_mci),ncol = 3,byrow = T))
  # senmsenl_mci = as.data.frame(matrix(unlist(spemspel_mci),ncol = 3,byrow = T))
  
  # plot
  # real
  plot(spec_mci$V1[order(spec_mci$V1)],sens_mci$V1[order(spec_mci$V1)],xlim = c(1,0),type='l',lty=1,
       xlab = 'specificity', ylab = 'sensitivity',col='blue',lwd=4)
  lines(spec_mci$V2[order(spec_mci$V1)],sens_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
  lines(spec_mci$V3[order(spec_mci$V1)],sens_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
  
  #perm
  lines(specl_mci$V1[order(specl_mci$V1)],sensl_mci$V1[order(specl_mci$V1)],xlim = c(1,0),type='l',lty=2,lwd=4,col='red')
  lines(specl_mci$V2[order(spec_mci$V1)],sensl_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
  lines(specl_mci$V3[order(spec_mci$V1)],sensl_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
  
  abline(a=1,b=-1,lty=4,lwd=2)
  title('Receiver-Operating-Curve for Behavioral Classifier')
  legend(1.02, 1.02, c("ROC of classifier with 95% bounds", "ROC of control classifier with 95% bounds", "hypothetical null"), col = c('blue', 'red', 'black'),
         text.col = "black", lty = c(1, 2, 4),
         merge = TRUE, bg = "gray90")
  
  # winning model frequencies
  disp('These models have been chosen:')
  cur_tab = table(cur_mod_sel_nooCV)
  print(cur_tab)
  
  # barplot the winning models
  # prep data frame
  cur_mod_sel_nooCV_freq  = as.data.frame(table(cur_mod_sel_nooCV))
  
  # get complexity score
  complexity = data.frame(names(fm),unlist(lapply(fm,length)))
  names(complexity) = c('modname','complexity')
  
  all_models_num          = as.character(1:length(fm))
  all_models_str          = names(fm)
  # add models that have 0 freq
  cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV = as.character(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV)
  for (ii in 1:length(all_models_str)) {
    if (!all_models_str[ii] %in% cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV) {
      cur_mod_sel_nooCV_freq = rbind(cur_mod_sel_nooCV_freq,c(all_models_str[ii],0))
    }
  }
  
  # transform to numeric
  cur_mod_sel_nooCV_freq$Freq = as.numeric(cur_mod_sel_nooCV_freq$Freq)
  
  # get complexity score
  cur_mod_sel_nooCV_freq$complexity = as.numeric(agk.recode.c(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV,complexity$modname,complexity$complexity))
  
  # sort
  cur_mod_sel_nooCV_freq = cur_mod_sel_nooCV_freq[order(cur_mod_sel_nooCV_freq$complexity),]
  # get the model names
  names_models_orig               = names(fm)
  cur_mod_sel_nooCV_freq$modnames = cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV
  #cur_mod_sel_nooCV_freq$modnames = names_models_orig[as.numeric(cur_mod_sel_nooCV_freq$cur_mod_sel_nooCV)]
  #cur_mod_sel_nooCV_freq$modnames = agk.recode.c(cur_mod_sel_nooCV_freq$modnames,names_models_orig,names_models)
  #cur_mod_sel_nooCV_freq$modnames = factor(cur_mod_sel_nooCV_freq$modnames, levels = cur_mod_sel_nooCV_freq$modnames)
  cur_mod_sel_nooCV_freq$Freq     = as.numeric(cur_mod_sel_nooCV_freq$Freq)
  
  # throwing out the non-linear one
  if (!isempty(grep('lanle',cur_mod_sel_nooCV_freq$modnames))) {
    warning('I am throwing out the non-linear mods in report.')
    cur_mod_sel_nooCV_freq = cur_mod_sel_nooCV_freq[-grep('lanle',cur_mod_sel_nooCV_freq$modnames),]
  }
  
  # # recoding
  #cur_mod_sel_nooCV_freq$modnames = agk.recode.c(cur_mod_sel_nooCV_freq$modnames, c('ac','acc'),c('a','ac'))
  
  # fixing the order by complexity, which we prepped above
  cur_mod_sel_nooCV_freq$modnames = factor(cur_mod_sel_nooCV_freq$modnames, levels = cur_mod_sel_nooCV_freq$modnames)
  
  p = ggplot(data = cur_mod_sel_nooCV_freq, aes(modnames,Freq))
  p = p+geom_bar(stat="identity") + ylab("Frequency model selected") + xlab('models ordered by model complexity')
  p = p + coord_cartesian(ylim=c(0,700))
  p = p + ggtitle(paste0("Frequency of model selection in ",
                         sum(cur_mod_sel_nooCV_freq$Freq)," cross validation rounds"))
  p = p + theme(text = element_text(size=20),
                axis.text.x = element_text(angle=90, hjust=1)) 
  p = p + coord_flip()
  print(p)
  
  # get the most often chosen model
  if (length(which(max(cur_tab) == cur_tab)) > 1) {
    stop('At least two models won equally often.')
  }
  winning_mod = names(cur_tab)[which(max(cur_tab) == cur_tab)]
  
  # winning model beta values with CIs
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
  imp_0 = function(x) {x[is.na(x)] = 0; return(x)}
  win_mod_coefs = as.data.frame(lapply(win_mod_coefs,FUN=imp_0))
  
  # now we get the mean, upper and lower bootstrapped CI
  # ci_res = lapply(win_mod_coefs,FUN=agk.boot.ci,cur_fun=mean,
  #                 lower=0.025,upper=0.975,R=1000)
  
  ci_res = lapply(win_mod_coefs,FUN=agk.mean.quantile.c,
                  lower=0.025,upper=0.975)
  ci_res = as.data.frame(t(as.data.frame(ci_res)))
  names(ci_res) = c('mean','lower','upper')
  
  # FIRST INTERCEPT IS THE GROUP PRED INTERCEPT, SECOND IS FROM THE BEHAV MODEL
  ci_res$coef                                = row.names(ci_res)
  ci_res$coef[ci_res$coef == 'Intercept']    = 'Int_behav_model'
  ci_res$coef[ci_res$coef == 'X.Intercept.'] = 'Int_classifier'
  

  p = ggplot(data = ci_res, aes(coef,mean))
  p = p+geom_bar(stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                        width = 0) + ylab("mean (95% CI over CV rounds)\n")
  
  p <- p + ggtitle("Estimated classification weights of winning model with CIs")
  p = p + theme(text = element_text(size=20),
                axis.text.x = element_text(angle=90, hjust=1)) 
  print(p)
}

# REPORTING PHYS ONLY =========================================================
if (do_report_feat_only) {
  set.seed(des_seed)
  # only for additional features (physio)
  do_feat_sel = des_feat_sel
  # add cue reactivity predictors to full model: peripheral physiology
  # can be set to 0, if feature selection and adding should not be done on this
  add_cr_pp   = add_cr_pp_ma
  add_cr_ra   = add_cr_ra_ma
  # run the models (for param extraction in exp)
  est_models  = 1
  # run the init (the CV of) severity pred
  source('severity_pred_init_v6.R')
  # assign the right CVp list
  CV_res_list = CV_res_list_op
  
  # p-value for the algorithm
  # get the accuracy null-distrib/alt hypothesis distribution
  cor    = c() # correlation
  cor_p  = c() # correlation permuted
  coru   = c() # correlation unpooled (mean of k-fold)
  coru_p = c() # correlation unpooled (mean of k-fold) permuted
  mse    = c() # mean squared error
  mse_p  = c() # mean squared error permuted
  accs_p = c()
  accs   = c()
  sens_p = c()
  sens   = c()
  spec_p = c()
  spec   = c()
  auc_p  = c()
  auc    = c()
  roc_pl = list()
  rocl   = list()
  for (ii in 1:length(CVp_res_list)) {
    accs_p[ii]   = CVp_res_list[[ii]]$acc$accuracy
    accs[ii]     = CV_res_list[[ii]]$acc$accuracy
    sens_p[ii]   = CVp_res_list[[ii]]$acc$cur_sens
    sens[ii]     = CV_res_list[[ii]]$acc$cur_sens
    spec_p[ii]   = CVp_res_list[[ii]]$acc$cur_spec
    spec[ii]     = CV_res_list[[ii]]$acc$cur_spec
    roc_pl[[ii]] = CVp_res_list[[ii]]$roc
    auc_p[ii]    = as.numeric(roc_pl[[ii]]$auc)
    rocl[[ii]]   = CV_res_list[[ii]]$roc
    auc[ii]      = as.numeric(rocl[[ii]]$auc)
    cor[[ii]]    = CV_res_list[[ii]]$cor$estimate
    cor_p[[ii]]  = CVp_res_list[[ii]]$cor$estimate
    coru[[ii]]   = CV_res_list[[ii]]$coru
    coru_p[[ii]] = CVp_res_list[[ii]]$coru
    mse[[ii]]    = CV_res_list[[ii]]$mse
    mse_p[[ii]]  = CVp_res_list[[ii]]$mse
  }
  
  # make new ROC objects
  for (ii in 1:length(rocl)) {
    # non_permuted
    predictor  = rocl[[ii]]$original.predictor
    response   = rocl[[ii]]$original.response
    response   = ifelse(response == 'PG',1,0)
    predictor  = agk.scale_range(predictor,0,1)
    rocl[[ii]] = roc(response,predictor)
    auc[ii]    = auc(rocl[[ii]])
    rocl[[ii]] = smooth( rocl[[ii]],n = 500)
    
    # permuted (or control model)
    predictor    = roc_pl[[ii]]$original.predictor
    response     = roc_pl[[ii]]$original.response
    response     = ifelse(response == 'PG',1,0)
    predictor    = agk.scale_range(predictor,0,1)
    roc_pl[[ii]] = roc(response,predictor)
    auc_p[ii]    = auc(roc_pl[[ii]])
    
    if (var(roc_pl[[ii]]$original.predictor) == 0) {
      cur_obj = list()
      cur_obj$sensitivities = seq(0,1,length.out = 502)
      cur_obj$specificities = seq(0,1,length.out = 502)
      roc_pl[[ii]]          = cur_obj
    } else {
      roc_pl[[ii]] = smooth(roc_pl[[ii]],n = 500)
    }
  }
  
  # print the median of acc, sens, spec
  acc_sens_spec       = data.frame(accs,sens,spec,auc)
  acc_sens_spec_p     = data.frame(accs_p,sens_p,spec_p,auc_p)
  med_ac_se_sp        = lapply(acc_sens_spec,FUN=mean)
  med_ac_se_sp_np_ci  = lapply(acc_sens_spec,FUN=agk.boot.ci,
                               cur_fun=mean,R=3000,lower=0.025,upper=0.975)
  med_ac_se_sp_npp_ci = lapply(acc_sens_spec_p,FUN=agk.boot.ci,
                               cur_fun=mean,R=3000,lower=0.025,upper=0.975)
  disp('Mean accuracy, sensitivity, specificity across CV rounds.')
  print(med_ac_se_sp)
  # acc
  diffs        = accs - accs_p
  cur_p        = agk.density_p(diffs,0)
  disp("Probability that permuted and non-permuted CV'd accuracy of algorithms are the same.")
  print(cur_p)
  # sens
  diffs        = sens - sens_p
  cur_p        = agk.density_p(diffs,0)
  disp("Probability that permuted and non-permuted CV'd sensitivity of algorithms are the same.")
  print(cur_p)
  # spec
  diffs        = spec - spec_p
  cur_p        = agk.density_p(diffs,0)
  disp("Probability that permuted and non-permuted CV'd specificity of algorithms are the same.")
  print(cur_p)
  # auc
  diffs        = auc - auc_p
  cur_p        = agk.density_p(diffs,0)
  disp("Probability that permuted and non-permuted CV'd smoothed ROC-AUC of algorithms are the same.")
  print(cur_p)
  
  # ROC curve (get all the coordinates)
  # not perm
  specificities  = t(as.matrix(rocl[[1]]$specificities))
  sensitivities  = t(as.matrix(rocl[[1]]$sensitivities))
  
  # perm
  specificitiesl = t(as.matrix(roc_pl[[1]]$specificities))
  sensitivitiesl = t(as.matrix(roc_pl[[1]]$sensitivities))
  
  # get all the ROC curves
  for (ii in 2:length(rocl)) {
    if (length(rocl[[ii]]$specificities) != length(specificities[ii-1,])) {
      stop('length not same!')
    }
    
    if (length(roc_pl[[ii]]$specificities) != length(specificitiesl[ii-1,])) {
      stop('length not same!')
    }
    
    specificities = rbind(specificities,t(as.matrix(rocl[[ii]]$specificities)))
    sensitivities = rbind(sensitivities,t(as.matrix(rocl[[ii]]$sensitivities)))
    
    specificitiesl = rbind(specificitiesl,t(as.matrix(roc_pl[[ii]]$specificities)))
    sensitivitiesl = rbind(sensitivitiesl,t(as.matrix(roc_pl[[ii]]$sensitivities)))
  }
  
  
  # get CI
  specificities  = as.data.frame(specificities)
  sensitivities  = as.data.frame(sensitivities)
  specificitiesl = as.data.frame(specificitiesl)
  sensitivitiesl = as.data.frame(sensitivitiesl)
  
  # ...using just quantiles
  spec_mci  = lapply(specificities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  sens_mci  = lapply(sensitivities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  specl_mci = lapply(specificitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  sensl_mci = lapply(sensitivitiesl,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  
  spec_mci   = as.data.frame(matrix(unlist(spec_mci),ncol = 3,byrow = T))
  sens_mci   = as.data.frame(matrix(unlist(sens_mci),ncol = 3,byrow = T))
  specl_mci  = as.data.frame(matrix(unlist(specl_mci),ncol = 3,byrow = T))
  sensl_mci  = as.data.frame(matrix(unlist(sensl_mci),ncol = 3,byrow = T))
  
  # plot
  # real
  plot(spec_mci$V1[order(spec_mci$V1)],sens_mci$V1[order(spec_mci$V1)],xlim = c(1,0),type='l',lty=1,
       xlab = 'specificity', ylab = 'sensitivity',col='blue',lwd=4)
  lines(spec_mci$V2[order(spec_mci$V1)],sens_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
  lines(spec_mci$V3[order(spec_mci$V1)],sens_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
  
  #perm
  lines(specl_mci$V1[order(specl_mci$V1)],sensl_mci$V1[order(specl_mci$V1)],xlim = c(1,0),type='l',lty=2,lwd=4,col='red')
  lines(specl_mci$V2[order(spec_mci$V1)],sensl_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
  lines(specl_mci$V3[order(spec_mci$V1)],sensl_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
  
  abline(a=1,b=-1,lty=4,lwd=2)
  title('Receiver-Operating-Curve for Peripheral Physiological Classifier')
  legend(0.68, 0.2, c("ROC of classifier with 95% perc.", "ROC of control classifier with 95% perc.", "hypothetical null"), col = c('blue', 'red', 'black'),
         text.col = "black", lty = c(1, 2, 4),
         merge = TRUE, bg = "gray90")
  
  
  
  
  # winning model beta values with CIs
  win_mods_distr_c = list_winning_model_c_nooCV_op
  win_mods_distr_l = list_winning_model_l_nooCV_op
  
  
  # make a data frame of it:
  # prep which RT set was added
  win_mod_coefs        = t(as.matrix(win_mods_distr_c[[1]]))
  win_mod_coefs        = as.data.frame(win_mod_coefs)
  names(win_mod_coefs) = win_mods_distr_l[[1]]
  
  for (ii in 2:length(win_mods_distr_c)) {
    cur_win_mod_coefs        = t(as.matrix(win_mods_distr_c[[ii]]))
    cur_win_mod_coefs        = as.data.frame(cur_win_mod_coefs)
    names(cur_win_mod_coefs) = win_mods_distr_l[[ii]]
    win_mod_coefs = rbind.fill(win_mod_coefs,cur_win_mod_coefs)
  }
  imp_0 = function(x) {x[is.na(x)] = 0; return(x)}
  win_mod_coefs = as.data.frame(lapply(win_mod_coefs,FUN=imp_0))
  
  
  # make a data frame of it:
  win_mod_coefs = t(as.matrix(win_mods_distr_c[[1]]))
  for (ii in 2:length(win_mods_distr_c)) {
    win_mod_coefs = rbind(win_mod_coefs,win_mods_distr_c[[ii]])
  }
  win_mod_coefs        = data.frame(win_mod_coefs)
  names(win_mod_coefs) = win_mods_distr_l[[1]]
  # now we get the mean, upper and lower bootstrapped CI
  ci_res = lapply(win_mod_coefs,FUN=agk.mean.quantile.c,
                  lower=0.025,upper=0.975)
  ci_res = as.data.frame(t(as.data.frame(ci_res)))
  names(ci_res) = c('mean','lower','upper')
  ci_res$coef   = row.names(ci_res)
  
  # take out vars that do not work
  ci_res = ci_res[ci_res$mean != 0,]
  
  p = ggplot(data = ci_res, aes(coef,mean))
  p = p+geom_bar(stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, color=cbbPalette[4],
                        width = 0) + ylab("mean (95% CI over CV rounds)\n")
  
  p <- p + ggtitle("Estimated parameters of best model (with ratings only) with CIs \n")
  p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}
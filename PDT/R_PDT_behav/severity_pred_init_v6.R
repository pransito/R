## PREAMBLE ===================================================================
# Version 6.0
# trying to predict severity of the whole group (HC+PG)
# using all LA models with or without cat modulation
# in this Version: selection of which model to extract exp params is
# part of training (inner CV will be used)
# this script initializes the data; gets the behav model parameters

# initialization script to set parameters and get data
# before, please run: import_data_pdt.R, select_study.R
# this script will be called by severity_pred_loop_v2.R

# here focus on lmlist; glmer too calculation intense
# and ecologically invalid; cannot do glmer on single subject
# ML: only glmnet
# dev: new in 6.0: formula is now always "Y ~ ."
# no more "pred_"

# CV
# is done on for model selection, feature selection

# one model full and then CV
# no more full and baseline
# every model is a model in its own right
# this is prep for feature selection + ML + prediction(CV)

# author: Alexander Genauck
# email:  alexander.genauck@charite.de
# date:   09.06.2018
warning('REDO the non-linear bit after review; https://de.wikipedia.org/wiki/Verlustaversion')

# PARAMETERS TO SET: General ==================================================
# in glmnet final model which alphas to be tested?
# CAREFUL: MODEL SELECTION PUT ON RIDGE (alpha=0) BY HAND; ALPHAS IGNORED!
# TO NOT UNSELECT ALL
alphas          = c(0,0.05,0.1,0.2,0.4,0.8,0.95,1)

# desired_lambdas
# should always be NULL
# however there seems to be a bug:
# https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R
# so need to use predefined lambdas
des_lambdas      = exp(seq(log(0.001), log(10), length.out=100))
# use model matrix (glmnetUtils)
useModelFrame    = T
doGrouped        = F
# how many inner CV folds to be used for picking lambda
# k = 10 recommended, 'LOOCV' is leave-one-out CV
#cnfolds         = "LOOCV"
cnfolds          = 10
if (type_measure_binomial == 'auc') {
  cnfolds        = 5
}
# k of the outer CV
# if left empty then LOOCV will be used
outer_k          = c(10)
# do permutations? (for checking that CV procedure is not biased)
# not for significance testing of the CV score: cause there is no certain CV estimate
do_permut        = F

# PARAMETERS TO SET: Behav ====================================================
# in case of lmlist should residual sum of squares be pooled?
do_pool          = F
# use behav params; if FALSE then only phys/rating; normally TRUE
use_behav_params = T

# PARAMETERS TO SET: Physio ===================================================
# which aggregate fun to be used for cue reactivity/physio/MRI variables
agg_fun          = mean.rmna

## PROCESS PREPS ==============================================================
# estimate models backup
est_models_bcp = est_models

# getting data and subsetting if desired
data_pdt = data_pdt_bcp_study_selected
dat_match = dat_match_bcp_study_selected

# prep if peripheral physiology or ratings should be added
if (add_cr_pp  == 1 & which_study != 'MRT' || add_cr_ra  == 1) {
  # excluding subjects because of missing in pp or ra
  # prepping data frames for ra and pp
  cur_path = getwd()
  setwd(root_wd)
  source("get_phys_and_rating_params_and_plots.R")
  setwd(root_wd)
} else if (add_cr_pp  == 1 & which_study == 'MRT') {
  # excluding subjects because of missing in pp
  # prepping data frames for pp
  setwd(root_wd)
  source("get_phys_and_rating_params_MRI.R")
  setwd(root_wd)
}

# cleaning cr_agg_pp
if (add_cr_pp_ma) {
  cr_agg_pp_uncleaned = cr_agg_pp
  if (regress_out_covs) {
    if (!exists('cr_agg_pp_cleaned')) {
      disp('Cleaning cr_agg_pp from influence of unmatched covariates of no interest...')
      # test if subjects aligned
      stopifnot(length(dat_match$VPPG) == length(cr_agg_pp$subject))
      stopifnot(all(dat_match$VPPG == cr_agg_pp$subject))
      if (which_study == 'MRT') {
        # cleaning
        preds = dat_match[c('smoking_ftdt','edu_years','edu_years_voca')]
        cr_agg_pp_cleaned = as.data.frame(lapply(cr_agg_pp,FUN=agk.regress.out.c,rob=F,inf_crit='BIC',preds=preds,ro_info = 0))
        cr_agg_pp         = cr_agg_pp_cleaned
      } else if (which_study == 'POSTPILOT_HCPG') {
        # cleaning
        preds = dat_match[c('smoking_ftdt')]
        cr_agg_pp_cleaned = as.data.frame(lapply(cr_agg_pp,FUN=agk.regress.out.c,rob=F,inf_crit='BIC',preds=preds,ro_info = 0))
        cr_agg_pp         = cr_agg_pp_cleaned
      }
      disp('Done.')
    } else {
      cr_agg_pp = cr_agg_pp_cleaned
    }
  }
}

# reduce the fMRI data
if (reduce_fMRI_data & add_cr_pp_ma & which_study == 'MRT') {
  cr_agg_pp_r = cr_agg_pp[c(names(cr_agg_pp)[grep('PicgamonxaccX',names(cr_agg_pp))],'subject')]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__.*DRN_8',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__.*AIns',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__.*PIns',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('_BA_',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('_ACgG',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__grp01.*_.OrG',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__grp01.*_MFC',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__grp01.*_MSFG',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_.*_MFC',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_.*_MSFG',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_.*_full_midbrain',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_._Acc.*_._.OrG',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_._Amy.*_._Amy',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp_r = cr_agg_pp_r[grep('SS__PPI_._Acc.*_._Acc',names(cr_agg_pp_r),invert = T)]
  cr_agg_pp   = cr_agg_pp_r
}

# set deviance measure
cur_family = "binomial"
type_measure = type_measure_binomial

# cross-valid results
CV_res = list()

# cons
if (do_data_inv == 0) {
  contrasts(data_pdt$cat) = contrasts(data_pdt$cat)
}

# control variables
if (length(pred_to_control) == 1) {
  # exactly one thing to control for
  # then in ridging need to add a random regressor
  do_con = 1
} else if (length(pred_to_control) > 1) {
  # more than one to control for
  do_con = 2
} else {
  # will no con for anything
  do_con = 0
}

## CHECK AVAILABLE DATA =======================================================
# checking data which is there
# data pdt
disp("In current data_pdt I have")
all_subs = unique(data_pdt$subject)
grp_var  = agk.recode.c(all_subs,data_pdt$subject,data_pdt$HCPG)
subgrp   = data.frame(subject = all_subs,group = grp_var)
print(table(subgrp$group))

# dat_match
disp(paste("Of the", length(all_subs), "subjects, I have dat_match info on",
           sum(dat_match$VPPG %in% all_subs), "subjects"))

# dat_match KFG, SOGS, BIG check
cur_dat_match = dat_match[dat_match$VPPG %in% all_subs,]
cur_dat_match = cur_dat_match[c("KFG","SOGS",'BIG')]
disp("Are there any NA's in KFG or SOGS or BIG?")
print(table(is.na(cur_dat_match)))

## MODEL PARAMS PER SUB =======================================================
# report to user
disp('')
disp('Getting the behavioral model parameters per subject.')

# function for charpentier model
# maybe not use anymore cause can be done in the general get lambda function
transform_coefs = function(lml) {
  if (!is.data.frame(lml)) {
    lml      = coef(lml)
  }
  lml$mu   = lml$gain
  lml$LA   = -(lml$loss/lml$gain)
  lml$gain = NULL
  lml$loss = NULL
  return(lml)
}


# for charpentier per category
data_pdt_neu       = data_pdt[data_pdt$cat == "neutral",]
data_pdt_gam       = data_pdt[data_pdt$cat == "gambling",]
data_pdt_neg       = data_pdt[data_pdt$cat == "negative",]
data_pdt_pos       = data_pdt[data_pdt$cat == "positive",]

# estimate models for params extraction
if (acc_num == 1) {
  needed_family = 'gaussian'
} else {
  needed_family = 'binomial'
}

if (est_models == 1) {
  
  # legacy; but needs to stay for now
  clmList = lmList
  
  # lmlist
  a         = clmList(accept_reject ~ 1 | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  ac        = clmList(accept_reject ~ cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  laec      = clmList(accept_reject ~ gain+loss+ed_abs+cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  lac       = clmList(accept_reject ~ gain+loss+cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  laeci     = clmList(accept_reject ~ (gain+loss+ed_abs)*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  laci      = clmList(accept_reject ~ (gain+loss)*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  lae       = clmList(accept_reject ~ gain+loss+ed_abs | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  la        = clmList(accept_reject ~ gain+loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  # ratio lmlist (no mu, please! drop it! it's enough!)
  lar       = clmList(accept_reject ~ ratio | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  larc      = clmList(accept_reject ~ ratio + cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  larci     = clmList(accept_reject ~ ratio*cat | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  
  # charpentier lmlist (same as la but without intercept, then mu and lambda can be analytically computed)
  laCh      = clmList(accept_reject ~ 0 + gain + loss | subject,data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
  # Charpentier per category
  laChci    = clmList(accept_reject ~ 0 + gain + loss + gain.catgambling + gain.catnegative + gain.catpositive +       
                        loss.catgambling + loss.catnegative + loss.catpositive | subject,
                      data = data_pdt,family = needed_family,na.action = NULL,pool = do_pool)
}

# packing (in order of complexity)
fm = list()
fm$a       = a
fm$ac      = ac
fm$lar     = lar
fm$laCh    = laCh
fm$la      = la
fm$lae     = lae
fm$larc    = larc
fm$lac     = lac
fm$laec    = laec
fm$larci   = larci
fm$laci    = laci
fm$laeci   = laeci
fm$laChci  = laChci

# backup the fm
fm_bcp = fm

## RIDGE THE EXP BEHAV MODELS =================================================
if (ridge_behav_models == T) {
  if (ridge_behav_models_anew == T) {
    # warning
    message('update glmnet to perhaps avoid having to provide fixed lambdas directly')
    # what alphas to use?
    cur_alpha = ridge_bm_alphas
    
    disp('Ridging the extracted exp. parameters for all models within every subject')
    # body of function to parallelize
    cur.mod.ridging = function(kk) {
      # decide whether this model is ridgeable
      cur_coefs = coef(fm[[kk]])
      if (!isempty(grep('Intercept',names(cur_coefs)))) {min_length = 3} else {min_length = 2} 
      if (length(names(cur_coefs)) < min_length) {
        fm[[kk]] = coef(fm[[kk]])
        next
      }
      
      disp(paste('Ridging model...',names(fm)[kk]))
      cur_mod    = lapply(fm[[kk]],FUN=agk.glmnet.lm.cvalpha,fam_arg='binomial',
                          type_measure='auc',lambda=des_lambdas,nfolds=5,
                          alphas = cur_alpha,verbosity = 1)
      coef_listr = cur_mod[[1]]$coef
      
      for (ii in 2:length(cur_mod)) {
        coef_listr = rbind(coef_listr,cur_mod[[ii]]$coef)
      }
      row.names(coef_listr) = NULL
      coef_listr = data.frame(coef_listr)
      colnames(coef_listr) = names(coef(fm[[kk]]))
      
      # packing back
      return(coef_listr)
    }
    
    total = length(fm)
    cl    = parallel::makeCluster((detectCores()-1))
    registerDoSNOW(cl)
    fm_ridge = foreach(kk=1:total, .packages=c('glmnetUtils','cvTools'),.verbose=T,.export = c()) %dopar% {
      cur.mod.ridging(kk)
    } 
    stopCluster(cl)
    
    # packing back, saving
    fm = fm_ridge
    save(file = 'ridged_behav_models.RData', list = c('fm'))
  } else {
    disp('loading the ridged behav params')
    load('ridged_behav_models.RData')
  }
}

# getting all the subs
all_subs = names(fm[[1]])

# Compute and ammend lambdas (as extra models)
# state which models should get LA appended
LA_append = c("la","lae","lac","laec","laci","laeci","laCh","laChci") 
fm_LA = list()
for (ii in which(names(fm) %in% LA_append)) {
  fm_LA[[ii]]      = agk.get.lambda(fm[[ii]],'gain','loss','(:|\\.)')
  names(fm_LA)[ii] = names(fm)[ii]
}

# drop models which haven't got LA appended
fm_LA = fm_LA[!unlist(lapply(fm_LA,is.null))]

# bind to none_LA models
names(fm_LA) = paste0(names(fm_LA),'_LA')
fm = c(fm,fm_LA)
cur_is_null = as.logical(unlist(lapply(fm,FUN=is.null))) == FALSE
fm          = fm[cur_is_null]

# coerce all fm to data frame (extract the coefs)
for (ii in 1:length(fm)) {
  df = fm[[ii]]
  if (is.data.frame(df)) {
    df = df
  } else {
    df = coef(df)
  }
  if (!is.data.frame(df)) {
    stop('Cur df is not a data frame nor can I extract coefficients data frame from df object.')
  }
  fm[[ii]] = df
}

# put back the row.names
for (ff in 1:length(fm)) {
  row.names(fm[[ff]]) = row.names(fm[[1]])
}

# Count and order by complexity
disp('Ordering models by complexity.')
complexity = unlist(lapply(fm,length))
fm         = fm[order(complexity)]

# names of models
names_models = names(fm)

# Scaling
if (initial_scale) {
  disp('Scaling the model parameters in fm.')
  fm = lapply(fm,FUN=agk.scale.ifpossible)
}

# getting the coefficient data frames per model: packing into featmod_coefs
disp('Packing fm into featmod_coefs')
featmod_coefs = list()

for (ii in 1:length(fm)) {
  # lmlist coef extraction of complete model
  if (is.data.frame(fm[[1]])) {
    tmp = fm[[ii]]
  } else {
    tmp = coef(fm[[ii]])
  }
  
  # feature: exp params
  if (names(tmp)[1] == '(Intercept)') {
    names(tmp)[1] = 'Intercept'
  }
  tmp_exp       = tmp
  
  # packing
  featmod_coefs[[ii]] = tmp_exp
}

# give a name
names(featmod_coefs) = names(fm)

# prep the coefs tables
for (ii in 1:length(featmod_coefs)) {
  tmp                 = featmod_coefs[[ii]]
  names(tmp)          = gsub(pattern = ":",replacement = "_of_",names(tmp))
  #names(tmp)          = paste0("pred_",names(tmp))
  #tmp$subject         = row.names(tmp)
  featmod_coefs[[ii]] = tmp
}

# getting the dependent variables ((PCA on) KFG and SOGS, BIG, HCPG)
for (ii in 1:length(featmod_coefs)) {
  tmp      = featmod_coefs[[ii]]
  #tmp$KFG  = as.numeric(agk.recode.c(tmp$subject,dat_match$VPPG,dat_match$KFG))
  tmp$HCPG = as.factor(agk.recode.c(row.names(tmp),dat_match$VPPG,dat_match$HCPG))
  featmod_coefs[[ii]] = tmp
}

# make a back up of the featmod_coefs
featmod_coefs_bcp = featmod_coefs

# # make the formulas for model selection
# featmod_sel_forms = list()
# for (ii in 1:length(featmod_coefs)) {
#   tmp        = featmod_coefs[[ii]]
#   pred_vars  = names(tmp)[grep("pred_",names(tmp))]
#   form_cmpl  = as.formula(paste('HCPG','~', paste(pred_vars, collapse=" + ")))
#   featmod_sel_forms[[ii]] = form_cmpl
# }

# getting other features ready for analysis


if (add_cr_pp  == 1 & add_cr_ra  == 0) {
  # ADDING OTHER FEATURES (PHYSIO)
  # getting the phys feature cluster
  tmp_phys            = cr_agg_pp

  # scaling (not done in with outer CV)
  if (initial_scale) {
    disp('Scaling additional features overall (NOT RECOMMENDED!).')
    tmp_phys    = as.data.frame(scale(tmp_phys))
  } else {
    disp('Not scaling additional features overall.')
  }
  
} else if ((add_cr_pp  == 1 || add_cr_ra  == 1) & which_study == 'MRT') {
  # ADDING OTHER FEATURES (E.G. PHYSIO)
  # getting the rat and phys feature cluster
  tmp         = tmp_exp
  exp_params  = names(tmp)
  tmp$subject = row.names(tmp)
  tmp_rat     = merge(tmp,cr_agg_ra,by = "subject",all.x = T)
  tmp_rat     = tmp_rat[!names(tmp_rat) %in% exp_params]
  tmp_rat     = tmp_rat
  rat_params  = names(tmp_rat)
  rat_params  = rat_params[!rat_params %in% c('subject')]
  tmp_phys    = merge(tmp_rat,cr_agg_pp,by = "subject",all.x = T)
  tmp_phys    = tmp_phys[!names(tmp_phys) %in% rat_params]
  tmp_phys    = tmp_phys
  
  # add row.names
  row.names(tmp_rat)  = tmp_rat$subject
  row.names(tmp_phys) = tmp_phys$subject
  tmp_rat$subject     = NULL
  tmp_phys$subject    = NULL
  
  # scaling (not done in with outer CV)
  if (initial_scale) {
    tmp_rat             = as.data.frame(scale(tmp_rat))
    tmp_phys            = as.data.frame(scale(tmp_phys))
  }
  
}


# packing selected model coefs plus other feature clusters
feature_clusters      = list()
feature_clusters[[1]] = NULL
ct                    = 2
if (add_cr_ra) {
  feature_clusters[[ct]] = tmp_rat
  ct                     = ct + 1
}
if (add_cr_pp) {
  feature_clusters[[ct]] = tmp_phys
  ct                     = ct + 1
}

# prep the coef table
for (ii in 1:length(feature_clusters)) {
  if (ii > 1) {
    # ii == 1 is already done
    tmp = feature_clusters[[ii]]
    names(tmp)    = gsub(pattern = ":",replacement = "_of_",names(tmp))
    #names(tmp)    = paste0("pred_",names(tmp))
    #tmp$subject   = row.names(tmp)
    feature_clusters[[ii]] = tmp
  }
}

# getting the dependent variable (PCA on KFG and SOGS and BIG)
for (ii in 1:length(feature_clusters)) {
  if (ii > 1) {
    tmp = feature_clusters[[ii]]
    #tmp$KFG  = as.numeric(agk.recode.c(tmp$subject,dat_match$VPPG,dat_match$KFG))
    tmp$HCPG = as.factor(agk.recode.c(row.names(tmp),dat_match$VPPG,as.character(dat_match$HCPG)))
    feature_clusters[[ii]] = tmp
  }
}

# scale the control variables
if (initial_scale & !isempty(pred_to_control)) {
  dat_match[[pred_to_control]] = scale(dat_match[[pred_to_control]])
}

# back up the unpermuted feature_clusters list
feature_clusters_bcp = feature_clusters

# # make the formulas for feature selection
# feature_sel_forms = list()
# for (ii in 1:length(feature_clusters)) {
#   if (ii > 1) {
#     # ii == 1 already done
#     tmp        = feature_clusters[[ii]]
#     pred_vars  = names(tmp)[grep('pred_',names(tmp))]
#     pred_str   = paste(pred_vars[1])
#     for (jj in 2:length(pred_vars)) {
#       pred_str = paste0(pred_str,paste(" +",pred_vars[jj]))
#     }
#     form_cmpl  = as.formula(paste("HCPG ~",pred_str))
#     feature_sel_forms[[ii]] = form_cmpl
#   }
# }
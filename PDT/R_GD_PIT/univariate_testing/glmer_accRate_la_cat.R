## PREAMBLE ===================================================================
# script to describe acceptance rate
# script to fit glmer models to answer questions on acceptance rate,
# loss aversion, effects of category and group
# model comparison to get significance of overall effect of e.g. group or category
# fixed effects of paramters and CIs will be based on bootstrapping non-parametrically
# the fixed effects
# permutation test to get p-value for parameters
# TODO: let permutation function also return an accurarcy or aic or so measure

# PREPARATION FOR FOREIGN RUN =================================================
rm(list=ls())
root_wd = 'C:/Users/genaucka/Google Drive/Library/R/PDT/R_GD_PIT'
setwd(root_wd)
load('.RData')

# LIBRARIES ===================================================================
agk.load.ifnot.install("psych")
agk.load.ifnot.install("pracma")
agk.load.ifnot.install("pls")
agk.load.ifnot.install("Hmisc")
agk.load.ifnot.install("lme4")
agk.load.ifnot.install("reshape2")
agk.load.ifnot.install("R.matlab")
agk.load.ifnot.install("gtools")
agk.load.ifnot.install("plyr")
agk.load.ifnot.install("ggplot2")
agk.load.ifnot.install("rgl")
agk.load.ifnot.install("gridExtra")
agk.load.ifnot.install("boot")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install("corrplot")
agk.load.ifnot.install("glmnet")
agk.load.ifnot.install("glmnetUtils")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install("parallel")
agk.load.ifnot.install("foreach")
agk.load.ifnot.install("doSNOW")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install("GPArotation")
agk.load.ifnot.install("nnet")
agk.load.ifnot.install("msm")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install('readxl')
agk.load.ifnot.install("rJava")
agk.load.ifnot.install("xlsx")
agk.load.ifnot.install('e1071')
agk.load.ifnot.install('ptw')
agk.load.ifnot.install('lmPerm')
agk.load.ifnot.install('pROC')
agk.load.ifnot.install('cvTools')
agk.load.ifnot.install('matlib')
agk.load.ifnot.install('robust')
agk.load.ifnot.install('e1071')
agk.load.ifnot.install('compiler')

## SETTINGS ===================================================================
# control object for the glmer fitting
cur_control = glmerControl(check.conv.grad="ignore",
                           check.conv.singular="ignore",
                           check.conv.hess="ignore",
                           optCtrl=list(optimizer = "nloptwrap",maxfun=250))
# do the boot or just load prebootstrapped results?
doBoot    = 0
# wd for saving the results of the bootstraps
bootResWd = paste(root_wd,'effects_under_0_la',sep = '/')
# how many bootstraps?
cur_num   = 500
# how many cpus to use?
cur_cpus  = 8
# put in the original fe instead of mean of the bootstrapped fe? 
put_in_original_fe = 1
# fit glmer models
doFitGlmer = 1

## PROCESS SETTINGS ===========================================================
message_fixef_setting =  function(put_in_original_fe) {
  if (put_in_original_fe) {
    message('using original fixef values for plotting and as estimates of fixef')
  } else {
    message('using mean of bootstrap instead of fixef values as means and as estimates of fixef')
  }
}
message_fixef_setting(put_in_original_fe)


## glmer models acceptance rate ===============================================
if (doFitGlmer) {
  moda_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial')
  moda_01  = glmer(accept_reject ~ cat + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_02  = glmer(accept_reject ~ cat*HCPG + (cat|subject) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  moda_01b = glmer(accept_reject ~ HCPG + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
}

## glmer models la ============================================================
if (doFitGlmer) {
  # with | cat
  modla_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|stim) + (1|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_01  = glmer(accept_reject ~ gain + loss + ed_abs + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs|stim) + (gain + loss + ed_abs|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_0g  = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs|stim)  + (gain + loss + ed_abs|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_c0  = glmer(accept_reject ~ (gain + loss + ed_abs) + cat + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs|stim) + (gain + loss + ed_abs|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cg  = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + cat*HCPG + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs|stim)  + (gain + loss + ed_abs|cat),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_cgi = glmer(accept_reject ~ (gain + loss + ed_abs)*cat*HCPG + ((gain + loss + ed_abs)*cat|subject) + (gain + loss + ed_abs|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  modla_ci  = glmer(accept_reject ~ (gain + loss + ed_abs)*cat + ((gain + loss + ed_abs)*cat|subject) + (gain + loss + ed_abs|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
  
  # reformulate the cgi model to make permutation test against cg model possible
  cur_mm        = model.matrix(accept_reject ~ (gain + loss + ed_abs)*cat*HCPG,data = data_pdt)
  cur_mm        = as.data.frame(cur_mm)
  cur_mm        = cur_mm[grep(':cat',names(cur_mm))]
  names(cur_mm) = gsub(':','.',names(cur_mm),fixed=T)
  data_pdt_rf   = data.frame(data_pdt,cur_mm)
  cur_ref_names = names(cur_mm)
  cur_ref_names = cur_ref_names[grep('.HCPGPG',cur_ref_names,fixed=T,invert = T)]
  cgi_form_fe = paste('accept_reject ~ (gain + loss + ed_abs)*HCPG + cat*HCPG + ', paste(names(cur_mm),collapse = ' + '))
  cgi_form_re = paste('(gain + loss + ed_abs + cat +',paste(cur_ref_names,collapse = ' + '), ' | subject',')','+ (gain + loss + ed_abs|stim)')
  cgi_form    = paste(cgi_form_fe,cgi_form_re,sep = ' + ')
  modla_cgirf = glmer(cgi_form,data = data_pdt_rf,family = 'binomial',nAGQ = 0,control=cur_control)
}

## check model fit per subject
cur_dp         = modla_cg@frame
cur_dp$pred_00 = as.numeric(as.numeric(predict(modla_00) >= 0) == cur_dp$accept_reject)
cur_dp$pred_cg = as.numeric(as.numeric(predict(modla_cg) >= 0) == cur_dp$accept_reject)
dens_df        = aggregate(cbind(pred_00,pred_cg) ~ subject + HCPG, data = cur_dp, FUN = mean)

## acceptance rate and under different cue conditions #########################
# acceptance rate graph, descriptives (CIs over subjects; better SD?)
mod_acc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject,data_pdt$cat), FUN=mean.rmna)
names(mod_acc) = c('subject','category','mean_acceptance')
mod_acc$Group  = agk.recode.c(mod_acc$subject,dat_match$VPPG,dat_match$HCPG)
mod_acc        = aggregate(mod_acc$mean_acceptance,by=list(mod_acc$Group,mod_acc$cat),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
mod_acc        = data.frame(mod_acc[[1]], mod_acc[[2]],mod_acc[[3]])
names(mod_acc) = c('Group','category','mean_acceptance','ci_0025','ci_0975')
mod_acc$Group  = agk.recode.c(mod_acc$Group,'PG','GD')

mRat  = ggplot(mod_acc, aes(category, mean_acceptance,fill=Group))
mRat  = mRat + labs(x='category', y=paste('Mean of acceptance (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean acceptance across categories")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# acceptance rate only between group
mod_accnc        = aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject), FUN=mean.rmna)
names(mod_accnc) = c('subject','mean_acceptance')
mod_accnc$Group  = agk.recode.c(mod_accnc$subject,dat_match$VPPG,dat_match$HCPG)
mod_accnc        = aggregate(mod_accnc$mean_acceptance,by=list(mod_accnc$Group),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
mod_accnc        = data.frame(mod_accnc[[1]],mod_accnc[[2]])
names(mod_accnc) = c('Group','mean_acceptance','ci_0025','ci_0975')
mod_accnc$Group  = agk.recode.c(mod_accnc$Group,'PG','GD')

# stats
anova(moda_00,moda_01,moda_02)

# stats without cat (simple acceptance rate difference between groups)
anova(moda_00,moda_01b)

## la overall #################################################################
# stats glmer
anova(modla_00,modla_01,modla_0g,modla_cg,modla_cgi)
anova(modla_00nc,modla_01nc,modla_0gnc,modla_cg)

# bootstrap
setwd(bootResWd)
if (doBoot == 1) {
  # bootstrap p-value modla_0g (permutation)
  effects_under_0_0g = agk.boot.p.mermod(mermod = modla_0g,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,permvars = c('HCPG'),type='perm')
  save(file= 'effects_under_0_0g_perm_1000.RData',list=c('effects_under_0_0g'))
  
  # bootstrap cfint modla_0g (np boot)
  boot_cfint_0g = agk.boot.cfint.mermod(mermod = modla_0g,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')
  save(file = 'boot_cfint_0g_1000.RData',list=c('boot_cfint_0g'))
}

# graph fixed effects la model
# alternative graph using cfint from fixed effects bootstrap
cie    = new.env()
load('boot_cfint_0g_1000_wc.RData',envir = cie)
rcfint = cie$boot_cfint_0g[[1]]
for (cc in 2:length(cie$boot_cfint_0g)) {
  rcfint = rbind(rcfint,cie$boot_cfint_0g[[cc]])
}

# add the original (TODO: should be part of bootstrapping function)
rcfint = rbind(rcfint,fixef(modla_0g))

# df of bootstrap data
rcfint           = data.frame(rcfint)
names(rcfint)[1] = c('Intercept')
rcfint_HC        = rcfint[c('Intercept','gain','loss')]
rcfint_PG        = rcfint[c('HCPGPG','gain.HCPGPG','loss.HCPGPG')] + rcfint_HC
names(rcfint_PG) = c('Intercept','gain','loss')
rcfint_HC$group  = 'HC'
rcfint_PG$group  = 'PG'
rcfint           = rbind(rcfint_HC,rcfint_PG)
rcfint$la        = -rcfint$loss/rcfint$gain

rcfinta         = aggregate(.~group,data=rcfint,FUN=agk.mean.quantile,lower=0.025,upper=0.975)
rcfintdf        = data.frame(rcfinta[[1]],rcfinta[[2]])
names(rcfintdf) = c('Group','mean','ci_0025','ci_0975')
rcfintdf$var    = names(rcfinta)[2] 
for (rr in 3:length(rcfinta)) {
  cur_df        = data.frame(rcfinta[[1]],rcfinta[[rr]])
  names(cur_df) = c('Group','mean','ci_0025','ci_0975')
  cur_df$var    = names(rcfinta)[rr] 
  rcfintdf      = rbind(rcfintdf,cur_df)
}
la_overall = rcfintdf

if (put_in_original_fe) {
  # put in the original fe
  obs_fixef      = get_la_fixef_pdt(modla_0g)
  cur_vars       = unique(la_overall$var)
  cur_fe_HC      = obs_fixef[c(1,3,5,7)]
  cur_fe_PG      = obs_fixef[c(2,4,6,8)]
  for (gg in 1:length(cur_vars)) {
    la_overall$mean[la_overall$Group == 'HC' & la_overall$var == cur_vars[gg]] = cur_fe_HC[gg]
    la_overall$mean[la_overall$Group == 'PG' & la_overall$var == cur_vars[gg]] = cur_fe_PG[gg]
  }
}

# new graph (plot)
la_overall$var = factor(la_overall$var, levels = c('la','gain','loss','Intercept'))
mRat           = ggplot(la_overall, aes(Group, mean,fill = var))
mRat           = mRat + labs(x='Group', y=paste('Mean of LA (',0.95*100,'% CI, bootstrapped)'))
mRat           = mRat + ggtitle("Fixed effects for loss aversion parameters per group")
mRat           = mRat + geom_bar(position="dodge", stat="identity")
dodge          = position_dodge(width=0.9)
mRat           = mRat + geom_bar(position=dodge, stat="identity")
mRat           = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# get the p-value
setwd(bootResWd)
load('effects_under_0_0g_perm_1000_wc.RData')
cur_vars    <- c("bg_HC","bg_PG","bl_HC", "bl_PG", "x_la_HCgrPG", "x_in_HCgrPG",'mse')
extr_params = function(x) {
  res             = list()
  res$in_HC       = x['(Intercept)']
  res$in_PG       = x['(Intercept)'] + x['HCPGPG']
  res$bg_HC       = x['gain']
  res$bl_HC       = x['loss']
  res$bg_PG       = x['gain'] + x['gain:HCPGPG']
  res$bl_PG       = x['loss'] + x['loss:HCPGPG']
  res$la_HC       = -x['loss']/x['gain']
  res$la_PG       = -res$bl_PG/res$bg_PG
  res$x_la_HCgrPG = res$la_HC - res$la_PG 
  res$x_in_HCgrPG = res$in_HC - res$in_PG
  res$mse         = x['mse']
  n_res           = names(res)
  res             = t(as.matrix(as.numeric(res)))
  names(res)      = n_res
  return(res)
}
effects_under_0 = lapply(effects_under_0_0g,FUN = extr_params)

for (ii in 1:length(cur_vars)) {
  cur_fun     <- function(l) {return(l[cur_vars[ii]])}
  eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0,FUN = cur_fun))",sep="")))
}
x_bg_HCgrPG  = bg_HC-bg_PG 
x_bl_HCgrPG  = bl_PG-bl_HC 
nulls        = list(x_la_HCgrPG,x_in_HCgrPG,x_bg_HCgrPG,x_bl_HCgrPG)
of           = la_overall
obs_fixef    = c()
ofvars       = c('la','Intercept','gain','loss')
for (oo in 1:length(ofvars)) {
  obs_fixef[oo] = of$mean[of$var == ofvars[oo] & of$Group == 'HC'] - of$mean[of$var == ofvars[oo] & of$Group == 'PG']
}
dirs         = c(-1,1,-1,1)

for (nn in 1:length(nulls)) {
  disp(paste('Permutation test group for:',ofvars[nn]))
  print(agk.density_p(nulls[[nn]],obs_fixef[nn]*dirs[nn]))
}

# overall permuation test
disp(paste('Permutation test group for:','mse'))
mse_a = mean((1/(1+exp(-predict(modla_0g)))-as.numeric(as.character(modla_0g@frame$accept_reject)))^2)
print(agk.density_p(mse,mse_a))


## la model plus effect of category ===========================================
# TODO: the whole thing needs to work as function for modla_0g modla_cg modla_cgi
if (doBoot == 1) {
  # bootstrap p-value modla_cg (permutation)
  effects_under_0_cg = agk.boot.p.mermod(mermod = modla_cg,mermod0 = modla_01,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,permvars = c('cat'),type='perm')
  save(file = 'effects_under_0_cg_perm_1000_wc.RData',list=c('effects_under_0_cg'))
  
  # bootstrap cfint modla_cg (np boot)
  boot_cfint_cg = agk.boot.cfint.mermod(mermod = modla_cg,num_cpus = cur_cpus,num = cur_num,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')
  save(file = 'boot_cfint_cg_1000_wc.RData',list=c('boot_cfint_cg'))
}

# loading boot cfint
all_saved_boots   = dir(pattern = 'boot_cfint_cg_1000_wc.RData')
all_boot_cfint_cg = c()
for (ss in 1:length(all_saved_boots)) {
  cur_e = new.env()
  load(all_saved_boots[ss],envir = cur_e)
  all_boot_cfint_cg = c(all_boot_cfint_cg,cur_e$boot_cfint_cg)
}

# coef computation
coefs = t(as.matrix(all_boot_cfint_cg[[1]]))
for (tt in 2:length(all_boot_cfint_cg)) {
  coefs = rbind(coefs,all_boot_cfint_cg[[tt]])
}

# add the original (TODO: should be part of bootstrapping function)
coefs = rbind(coefs,fixef(modla_cg))

# putting bootstrapped coefs in df
coefs         = data.frame(coefs)
coefs_HC      = coefs[c('X.Intercept.','gain','loss','catgambling', 'catnegative','catpositive')]
coefs_PG      = coefs[c('HCPGPG','gain.HCPGPG','loss.HCPGPG', 'HCPGPG.catgambling', 'HCPGPG.catnegative','HCPGPG.catpositive')]
coefs_PG      = coefs_PG + coefs_HC
coefs_HC$HCPG = 'HC'
coefs_PG$HCPG = 'PG'
# adjusting names
names(coefs_HC)[1] = 'Intercept'
names(coefs_PG)    = names(coefs_HC)
# getting la coef
coefs         = rbind(coefs_HC,coefs_PG)
coefs$la      = -coefs$loss/coefs$gain
var_names     = names(coefs)[!names(coefs) %in% 'HCPG']

# get the cfint (TODO THIS IS WRONG; notanother boot!)
cfint           = aggregate(.~HCPG,data = coefs,FUN=agk.mean.quantile,lower=0.025,upper=0.975)
cfint_df        = data.frame(cfint[[1]],cfint[[2]])
cfint_df$var    = var_names[1]
names(cfint_df) = c('Group','mean','ci_0025','ci_0975','var') 
for (cc in 3:length(cfint)) {
  cur_df        = data.frame(cfint[[1]],cfint[[cc]])
  cur_df$var    = var_names[cc-1]
  names(cur_df) = c('Group','mean','ci_0025','ci_0975','var') 
  cfint_df      = rbind(cfint_df,cur_df)
}

if (put_in_original_fe) {
  # replace means by real estimated fixed effects
  cur_fe              = fixef(modla_cg)
  cur_fe_HC           = cur_fe[c('(Intercept)','gain','loss','catgambling', 'catnegative','catpositive')]
  cur_fe_PG           = cur_fe[c('HCPGPG','gain:HCPGPG','loss:HCPGPG', 'HCPGPG:catgambling', 'HCPGPG:catnegative','HCPGPG:catpositive')]
  cur_fe_PG           = cur_fe_HC + cur_fe_PG
  cur_fe_PG$HCPG      = 'PG'
  cur_fe_HC$HCPG      = 'HC'
  names(cur_fe_HC)[1] = 'Intercept'
  names(cur_fe_PG)    = names(cur_fe_HC)
  cur_fe_df           = rbind(data.frame(cur_fe_HC),data.frame(cur_fe_PG))
  cur_fe_df$la        = -cur_fe_df$loss/cur_fe_df$gain
  cur_fe_df           = melt(cur_fe_df)
  cfint_df$mean       = cur_fe_df$value
}

# graph fixed effects la model with category values and bootstr cfints
cfint_df$var = factor(cfint_df$var,
                      levels = c('la','gain','loss','catgambling','catpositive','catnegative','Intercept'))
mRat  = ggplot(cfint_df, aes(Group, mean,fill = var))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean parameters per group")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# get the p-value
setwd(bootResWd)
load('effects_under_0_cg_perm_1000_wc.RData')
cur_vars    <- c("bg_HC","bg_PG","bl_HC", "bl_PG", "x_la_HCgrPG","catgam_HC","catgam_PG","catneg_HC","catneg_PG",
                 "catpos_HC","catpos_PG","x_in_HCgrPG","x_bg_HCgrPG","x_bl_HCgrPG","x_cg_HCgrPG","x_cp_HCgrPG","x_cn_HCgrPG",'mse')
extr_params = function(x) {
  res             = list()
  res$in_HC       = x['(Intercept)']
  res$in_PG       = x['(Intercept)'] + x['HCPGPG']
  res$bg_HC       = x['gain']
  res$bl_HC       = x['loss']
  res$bg_PG       = x['gain'] + x['gain:HCPGPG']
  res$bl_PG       = x['loss'] + x['loss:HCPGPG']
  res$la_HC       = -x['loss']/x['gain']
  res$la_PG       = -res$bl_PG/res$bg_PG
  res$cg_HC       = x['catgambling']
  res$cn_HC       = x['catnegative']
  res$cp_HC       = x['catpositive']
  res$cg_PG       = x['catgambling'] + x['HCPGPG:catgambling']
  res$cn_PG       = x['catnegative'] + x['HCPGPG:catnegative']
  res$cp_PG       = x['catpositive'] + x['HCPGPG:catpositive']
  res$mse         = x['mse']
  
  # group comparisons
  res$x_in_HCgrPG = res$in_HC - res$in_PG 
  res$x_bg_HCgrPG = res$bg_HC - res$bg_PG 
  res$x_bl_HCgrPG = res$bl_HC - res$bl_PG 
  res$x_la_HCgrPG = res$la_HC - res$la_PG 
  res$x_cg_HCgrPG = res$cg_HC - res$cg_PG 
  res$x_cn_HCgrPG = res$cn_HC - res$cg_PG 
  res$x_cp_HCgrPG = res$cp_HC - res$cp_PG 
  
  n_res           = names(res)
  res             = t(as.matrix(as.numeric(res)))
  names(res)      = n_res
  return(res)
}
effects_under_0 = lapply(effects_under_0_cg,FUN = extr_params)

for (ii in 1:length(cur_vars)) {
  cur_fun     <- function(l) {return(l[cur_vars[ii]])}
  eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0,FUN = cur_fun))",sep="")))
}

nulls        = list(x_cg_HCgrPG,x_cn_HCgrPG,x_cp_HCgrPG)
of           = cfint_df
obs_fixef    = c()
ofvars       = c('catgambling','catnegative','catpositive')
for (oo in 1:length(ofvars)) {
  obs_fixef[oo] = of$mean[of$var == ofvars[oo] & of$Group == 'HC'] - of$mean[of$var == ofvars[oo] & of$Group == 'PG']
}
dirs         = c(1,-1,1)

for (nn in 1:length(nulls)) {
  disp(paste('Permutation test group for:',ofvars[nn]))
  print(agk.density_p(nulls[[nn]]*dirs[nn],obs_fixef[nn]))
}

# overall permuation test
disp(paste('Permutation test group for:','mse'))
mse_a = mean((1/(1+exp(-predict(modla_cg)))-as.numeric(as.character(modla_cg@frame$accept_reject)))^2)
print(agk.density_p(mse,mse_a))

## la model times effect of category ==========================================
if (doBoot == 1) {
  # get from job scheduler
}

# loading boot cfint
all_saved_boots    = dir(pattern = 'boot_cfint_cgi_1000_wc.RData')
all_boot_cfint_cgi = c()
for (ss in 1:length(all_saved_boots)) {
  cur_e = new.env()
  load(all_saved_boots[ss],envir = cur_e)
  all_boot_cfint_cgi = c(all_boot_cfint_cgi,cur_e$boot_cfint_cgi)
}

# coef computation
coefs = t(as.matrix(all_boot_cfint_cgi[[1]]))
for (tt in 2:length(all_boot_cfint_cgi)) {
  coefs = rbind(coefs,all_boot_cfint_cgi[[tt]])
}

# add the original (TODO: should be part of bootstrapping function)
coefs = rbind(coefs,fixef(modla_cgi))

# putting bootstrapped coefs in df
coefs           = data.frame(coefs)
coefs           = agk.clean.intercept.name(coefs)
cur_ns          = names(coefs)
cur_ns_HC       = cur_ns[grep('HCPGPG',cur_ns,invert = T)]
cur_ns_PG       = cur_ns[grep('HCPGPG',cur_ns,invert = F)]
coefs_HC        = coefs[cur_ns_HC]
coefs_PG        = coefs[cur_ns_PG]
names(coefs_PG) = gsub('.HCPGPG','',names(coefs_PG))
coefs_PG        = coefs_PG + coefs_HC
if (!isempty(grep('gain.',names(coefs_HC)))) {
  # HC
  cur_pdf = coefs_HC[grep('gain.',names(coefs_HC))]
  coefs_HC[grep('gain.',names(coefs_HC))] = cur_pdf + as.data.frame(rep(coefs_HC['gain'],length(cur_pdf)))
  cur_pdf = coefs_HC[grep('loss.',names(coefs_HC))]
  coefs_HC[grep('loss.',names(coefs_HC))] = cur_pdf + as.data.frame(rep(coefs_HC['loss'],length(cur_pdf)))
  # PG
  cur_pdf = coefs_PG[grep('gain.',names(coefs_PG))]
  coefs_PG[grep('gain.',names(coefs_PG))] = cur_pdf + as.data.frame(rep(coefs_PG['gain'],length(cur_pdf)))
  cur_pdf = coefs_HC[grep('loss.',names(coefs_PG))]
  coefs_PG[grep('loss.',names(coefs_PG))] = cur_pdf + as.data.frame(rep(coefs_PG['loss'],length(cur_pdf)))
}


coefs_HC$HCPG = 'HC'
coefs_PG$HCPG = 'PG'
# adjusting names
names(coefs_HC)[1] = 'Intercept'
names(coefs_PG)    = names(coefs_HC)
# putting HC and PG together
coefs         = rbind(coefs_HC,coefs_PG)
# getting la coef
cur_la        = -coefs[grep('los.',names(coefs))]/coefs[grep('gai.',names(coefs))]
names(cur_la) = gsub('loss','la',names(cur_la))
coefs         = data.frame(coefs,cur_la)
# adjusting names
var_names     = names(coefs)[!names(coefs) %in% 'HCPG']


# get the cfint
cfint           = aggregate(.~HCPG,data = coefs,FUN=agk.mean.quantile,lower=0.025,upper=0.975)
cfint_df        = data.frame(cfint[[1]],cfint[[2]])
cfint_df$var    = var_names[1]
names(cfint_df) = c('Group','mean','ci_0025','ci_0975','var') 
for (cc in 3:length(cfint)) {
  cur_df        = data.frame(cfint[[1]],cfint[[cc]])
  cur_df$var    = var_names[cc-1]
  names(cur_df) = c('Group','mean','ci_0025','ci_0975','var') 
  cfint_df      = rbind(cfint_df,cur_df)
}

if (put_in_original_fe) {
  # replace means by real estimated fixed effects
  cur_fe              = fixef(modla_cgi)
  cur_ns              = names(cur_fe)
  cur_ns_HC           = cur_ns[grep('HCPGPG',cur_ns,invert = T)]
  cur_ns_PG           = cur_ns[grep('HCPGPG',cur_ns,invert = F)]
  cur_fe_HC           = cur_fe[cur_ns_HC]
  cur_fe_PG           = cur_fe[cur_ns_PG]
  names(cur_fe_PG)    = gsub(':HCPGPG','',names(cur_fe_PG))
  cur_fe_PG           = cur_fe_HC + cur_fe_PG
  if (!isempty(grep('gain.',names(cur_fe_HC)))) {
    cur_fe_HC[grep('gain.',names(cur_fe_HC))] = cur_fe_HC[grep('gain.',names(cur_fe_HC))] + cur_fe_HC['gain']
    cur_fe_HC[grep('loss.',names(cur_fe_HC))] = cur_fe_HC[grep('loss.',names(cur_fe_HC))] + cur_fe_HC['loss']
    cur_fe_PG[grep('gain.',names(cur_fe_PG))] = cur_fe_PG[grep('gain.',names(cur_fe_PG))] + cur_fe_PG['gain']
    cur_fe_PG[grep('loss.',names(cur_fe_PG))] = cur_fe_PG[grep('loss.',names(cur_fe_PG))] + cur_fe_PG['loss']
  }
  cur_fe_PG$HCPG      = 'PG'
  cur_fe_HC$HCPG      = 'HC'
  names(cur_fe_HC)[1] = 'Intercept'
  names(cur_fe_PG)    = names(cur_fe_HC)
  cur_fe_df           = rbind(data.frame(cur_fe_HC),data.frame(cur_fe_PG))
  cur_la              = -cur_fe_df[grep('los.',names(cur_fe_df))]/cur_fe_df[grep('gai.',names(cur_fe_df))]
  names(cur_la)       = gsub('loss','la',names(cur_la))
  cur_fe_df           = data.frame(cur_fe_df,as.data.frame(cur_la))
  cur_fe_df           = melt(cur_fe_df)
  cfint_df$mean       = cur_fe_df$value
}

# graph fixed effects la model with category values and bootstr cfints
cfint_df$var = factor(cfint_df$var,
                      levels = c('la',unique(cfint_df$var)[grep('la.',unique(cfint_df$var))],
                                 'gain','loss','catgambling','catpositive','catnegative','Intercept',
                                 unique(cfint_df$var)[grep('gain.',unique(cfint_df$var))],
                                 unique(cfint_df$var)[grep('loss.',unique(cfint_df$var))]))
# cut out vars of no interest
cfint_df = cfint_df[grep('^cat',cfint_df$var,invert=T),]
cfint_df = cfint_df[grep('Intercept',cfint_df$var,invert=T),]
# plotting
mRat  = ggplot(cfint_df, aes(Group, mean,fill = var))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean parameters per group")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# get the p-value
setwd(bootResWd)
load('effects_under_0_cgi_perm_1000_wc.RData')
cur_vars    <- c("x_la_HCgrPG","x_in_HCgrPG","x_bg_HCgrPG","x_bl_HCgrPG","x_cg_HCgrPG","x_cp_HCgrPG","x_cn_HCgrPG",
                 'x_bgcg_HCgrPG','x_blcg_HCgrPG','x_lacg_HCgrPG',
                 'x_bgcn_HCgrPG','x_blcn_HCgrPG','x_lacn_HCgrPG',
                 'x_bgcp_HCgrPG','x_blcp_HCgrPG','x_lacp_HCgrPG',
                 'mse')
extr_params = function(x) {
  names(x) = gsub(':','.',names(x),fixed=T)
  res             = list()
  for (rr in 1:length(x)) {
    res[[names(x)[rr]]] = x[rr]
  }
  
  # la per category
  bg_HC       = x['gain']
  bl_HC       = x['loss']
  bg_PG       = x['gain'] + x['gain.HCPGPG']
  bl_PG       = x['loss'] + x['loss.HCPGPG']
  res$la_HC       = -x['loss']/x['gain']
  res$la_PG       = -bl_PG/bg_PG - res$la_HC
  
  lacg_HC     = -(x['loss'] + x['loss.catgambling'])/(x['gain'] + x['gain.catgambling'])
  lacn_HC     = -(x['loss'] + x['loss.catnegative'])/(x['gain'] + x['gain.catnegative'])
  lacp_HC     = -(x['loss'] + x['loss.catpositive'])/(x['gain'] + x['gain.catpositive'])
  
  lacg_PG     = -(x['loss'] + x['loss.catgambling'] + x['loss.catgambling.HCPGPG'])/(x['gain'] + x['gain.catgambling'] + x['gain.catgambling.HCPGPG'])
  lacn_PG     = -(x['loss'] + x['loss.catnegative'] + x['loss.catnegative.HCPGPG'])/(x['gain'] + x['gain.catnegative'] + x['gain.catnegative.HCPGPG'])
  lacp_PG     = -(x['loss'] + x['loss.catpositive'] + x['loss.catpositive.HCPGPG'])/(x['gain'] + x['gain.catpositive'] + x['gain.catpositive.HCPGPG'])
  
  # baseline correct la
  res$lacg_HC     = lacg_HC     - res$la_HC
  res$lacn_HC     = lacn_HC     - res$la_HC
  res$lacp_HC     = lacp_HC     - res$la_HC
  res$lacg_PG     = lacg_PG     - res$la_HC - res$la_PG - res$lacg_HC 
  res$lacn_PG     = lacn_PG     - res$la_HC - res$la_PG - res$lacn_HC
  res$lacp_PG     = lacp_PG     - res$la_HC - res$la_PG - res$lacp_HC
  
  # mse
  if (!isempty(grep('mse',names(x)))) {
    res$mse         = x['mse']
  } else {
    res$mse         = NA
  }
  
  n_res           = names(res)
  res             = t(as.matrix(as.numeric(res)))
  names(res)      = n_res
  return(res)
}
effects_under_0 = lapply(effects_under_0_cgi,FUN = extr_params)

cur_vars     = names(effects_under_0[[1]])
cur_vars     = cur_vars[grep('mse',cur_vars,invert = T)]
cur_vars     = cur_vars[grep('(Intercept)',cur_vars,invert = T)]

for (ii in 1:length(cur_vars)) {
  cur_fun     <- function(l) {return(l[cur_vars[ii]])}
  eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0,FUN = cur_fun))",sep="")))
}

nulls     = cur_vars
obs_fixef = extr_params(fun_extr(modla_cgi))
obs_fixef = obs_fixef[grep('mse',names(obs_fixef),invert = T)]
obs_fixef = obs_fixef[grep('(Intercept)',names(obs_fixef),invert = T)]

for (nn in 1:length(names(obs_fixef))) {
  disp(paste('Permutation test group for:',names(obs_fixef)[nn]))
  cur_fe = names(obs_fixef)[nn]
  if (cur_fe == 'catnegative.HCPGPG') {cur_fe = 'HCPGPG.catnegative'}
  if (cur_fe == 'catpositive.HCPGPG') {cur_fe = 'HCPGPG.catpositive'}
  if (cur_fe == 'catgambling.HCPGPG') {cur_fe = 'HCPGPG.catgambling'}
  nulls_ind = which(cur_fe == nulls)
  print(agk.density_p(eval(parse(text=nulls[nulls_ind])),obs_fixef[nn],type = 'two.sided'))
}

# overall permuation test
disp(paste('Permutation test group for:','mse'))
mse_a = mean((1/(1+exp(-predict(modla_cgi)))-as.numeric(as.character(modla_cgi@frame$accept_reject)))^2)
print(agk.density_p(mse,mse_a,type = 'smaller'))

# just plot the the la per category and group lamod_cgi ======================
cfe    = fixef(modla_cgi)
res_HC = c()
res_PG = c()

# HC
bl_HC     = cfe['loss']
bl_cg_HC  = bl_HC + cfe['loss:catgambling']
bl_cn_HC  = bl_HC + cfe['loss:catnegative']
bl_cp_HC  = bl_HC + cfe['loss:catpositive']

bg_HC     = cfe['gain']
bg_cg_HC  = bg_HC + cfe['gain:catgambling']
bg_cn_HC  = bg_HC + cfe['gain:catnegative']
bg_cp_HC  = bg_HC + cfe['gain:catpositive']

res_HC[1] = -bl_HC/bg_HC
res_HC[2] = -bl_cg_HC/bg_cg_HC
res_HC[3] = -bl_cn_HC/bg_cn_HC
res_HC[4] = -bl_cp_HC/bg_cp_HC

# PG
bl_PG     = bl_HC + cfe['loss:HCPGPG']
bl_cg_PG  = bl_HC + cfe['loss:catgambling'] + cfe['loss:catgambling:HCPGPG']
bl_cn_PG  = bl_HC + cfe['loss:catnegative'] + cfe['loss:catnegative:HCPGPG']
bl_cp_PG  = bl_HC + cfe['loss:catpositive'] + cfe['loss:catpositive:HCPGPG']

bg_PG     = bg_HC + cfe['gain:HCPGPG']
bg_cg_PG  = bg_HC + cfe['gain:catgambling'] + cfe['gain:catgambling:HCPGPG']
bg_cn_PG  = bg_HC + cfe['gain:catnegative'] + cfe['gain:catnegative:HCPGPG']
bg_cp_PG  = bg_HC + cfe['gain:catpositive'] + cfe['gain:catpositive:HCPGPG']

res_PG[1] = -bl_PG/bg_PG
res_PG[2] = -bl_cg_PG/bg_cg_PG
res_PG[3] = -bl_cn_PG/bg_cn_PG
res_PG[4] = -bl_cp_PG/bg_cp_PG

la_cgi = data.frame(HC = res_HC,PG =res_PG, cat = c('overall','gambling','negative','positive'))
la_cgi = melt(la_cgi)
names(la_cgi)[2] = 'Group'

mRat  = ggplot(la_cgi, aes(Group, value,fill=cat))
mRat  = mRat + labs(x='Group', y=paste("Loss aversion computed from fixed effects of gain and loss"))
mRat  = mRat + ggtitle("Loss aversion computed from fixed effects of gain and loss")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
print(mRat)


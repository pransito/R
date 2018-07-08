p_fun = function(x) {
  return(x$p.value)
}

res = lapply(cr_agg_pp[c(2:length(cr_agg_pp_cleaned))],FUN = ttest.group,cur_group=dat_match$HCPG,cur_alt = 'two.sided')
res = unlist(lapply(res,p_fun))
res = res[res<0.01]

# modeling behav
res = lapply(featmod_coefs[[9]][c(1:6)],FUN = ttest.group,cur_group=dat_match$HCPG,cur_alt = 'two.sided')
res = unlist(lapply(res,p_fun))
res = res[res<0.05]

##  checking the cr_agg_pp ====================================================
cur_fun   = function(x) {
  if(!is.list(x)) {
    return(NA)
  }
  return(x$p.value)
}

ttest.group.con = function(x,cur_group,cur_con) {
  if (!is.numeric(x)){return(NA)}
  if (is.null(cur_con)) {ttest.group(x,cur_group)}
  # function to compute a ttest group, with cur_con as covariate
  cur_m0      = lm(x ~ cur_con)
  cur_mt      = lm(x ~ cur_con + cur_group)
  summar      = anova(cur_m0,cur_mt)
  p.value     = summar$`Pr(>F)`[2]
  res         = list()
  res$p.value = p.value
  return(res)
}

tt_res_cl = lapply(cr_agg_pp_cleaned,FUN = ttest.group,cur_group = dat_match$HCPG)
tt_res_cl = unlist(lapply(tt_res_cl,cur_fun))
tt_res_cl = tt_res_cl[tt_res_cl < 0.02]

tt_res_ucl = lapply(cr_agg_pp_uncleaned,FUN = ttest.group,cur_group = dat_match$HCPG)
tt_res_ucl = unlist(lapply(tt_res_ucl,cur_fun))
tt_res_ucl = tt_res_ucl[tt_res_ucl < 0.02]

## checking behav t-tests =====================================================
# laec
behav             = featmod_coefs$laec
behav$model       = 'laec'
tt_res_bv         = lapply(behav,FUN = ttest.group.con,cur_group = dat_match$HCPG, cur_con = dat_match$smoking_ftdt)
tt_res_bv         = unlist(lapply(tt_res_bv,cur_fun))
tt_res_bv         = tt_res_bv[tt_res_bv < 0.05]

# la
behav_la          = featmod_coefs$la
behav_la$model    = 'la'

# lac
behav_lac         = featmod_coefs$lac
behav_lac$model   = 'lac'
tt_res_bv_lac     = lapply(behav_lac,FUN = ttest.group.con,cur_group = dat_match$HCPG, cur_con = dat_match$smoking_ftdt)
tt_res_bv_lac     = unlist(lapply(tt_res_bv_lac,cur_fun))
tt_res_bv_lac     = tt_res_bv_lac[tt_res_bv_lac < 0.05]

# laci
behav_laci        = featmod_coefs$laci
behav_laci$model  = 'laci'


# laCh
behav_laCh            = featmod_coefs$laCh
behav_laCh$model      = 'laCh'
tt_res_bv_laCh        = lapply(behav_lac,FUN = ttest.group.con,cur_group = dat_match$HCPG, cur_con = dat_match$smoking_ftdt)
tt_res_bv_laCh        = unlist(lapply(tt_res_bv_laCh,cur_fun))
tt_res_bv_laCh        = tt_res_bv_laCh[tt_res_bv_laCh < 0.05]

# combine
behav         = rbind.fill(behav_la,behav_lac,behav_laci)

# la
behav$pred_gain_of_catgambling = (behav$pred_gain + behav$pred_gain_of_catgambling)
behav$pred_gain_of_catnegative = (behav$pred_gain + behav$pred_gain_of_catnegative)
behav$pred_gain_of_catpositive = (behav$pred_gain + behav$pred_gain_of_catpositive)

behav$pred_loss_of_catgambling = (behav$pred_loss + behav$pred_loss_of_catgambling)
behav$pred_loss_of_catnegative = (behav$pred_loss + behav$pred_loss_of_catnegative)
behav$pred_loss_of_catpositive = (behav$pred_loss + behav$pred_loss_of_catpositive)

behav$la      = -behav$pred_loss/behav$pred_gain
behav$la_of_catgambling = -(behav$pred_loss_of_catgambling )/(behav$pred_gain_of_catgambling )
behav$la_of_catnegative = -(behav$pred_loss_of_catnegative )/(behav$pred_gain_of_catnegative )
behav$la_of_catpositive = -(behav$pred_loss_of_catpositive )/(behav$pred_gain_of_catpositive )

# melt for plotting
behav_m          = melt(behav,id.vars = c('subject','HCPG','model'))
var_cis          = aggregate(behav_m,by=list(behav_m$HCPG,behav_m$variable,behav_m$model), FUN=agk.boot.ci,R=10000,lower=0.025,upper=0.975,cur_fun=median)
var_cis$subject  = NULL
var_cis$HCPG     = NULL
var_cis$variable = NULL
var_cis$model    = NULL

var_cis_df        = data.frame(var_cis[c(1,2,3)],var_cis$value)
names(var_cis_df) = c('HCPG','behav_param','model','mean','ci_0025','ci_0975')
var_cis_df        = var_cis_df[grep('KFG',var_cis_df$behav_param,invert = T),]
names(var_cis_df)[names(var_cis_df) == 'HCPG'] = 'Group'
var_cis_df$Group = agk.recode.c(var_cis_df$Group,'PG','GD')

# get rid of prefix 'pred_'
needed_levels          = levels(var_cis_df$behav_param)
needed_levels          = gsub('pred_','',needed_levels)
var_cis_df$behav_param = gsub('pred_','',var_cis_df$behav_param)
var_cis_df$behav_param = factor(var_cis_df$behav_param,levels = needed_levels)

mRat  = ggplot(var_cis_df, aes(behav_param, mean,fill=Group))
mRat  = mRat + labs(x='behavioral parameter', y=paste('Median (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Behav params")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) 
mRatg = mRat + facet_grid(model ~ .)
#mRatg = mRat + theme(strip.text.y = element_text(angle = 180))
mRatg = mRat+ facet_grid(model ~ .) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(mRatg)

## acceptance rate and under different cue conditions #########################
# acceptance rate graph, descriptives
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
cur_control = glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=250))
moda_00     = glmer(accept_reject ~ 1 + (1|subject) + (1|cat) + (1|stim),data = data_pdt,family = 'binomial')
moda_01     = glmer(accept_reject ~ cat + (cat|subject) + (1|cat) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
moda_02     = glmer(accept_reject ~ cat*HCPG + (cat|subject) + (1|cat) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)

anova(moda_00,moda_01,moda_02)

# stats without cat (simple acceptance rate difference between groups)
moda_01b    = glmer(accept_reject ~ HCPG + (1|subject) + (1|cat) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
anova(moda_00,moda_01b)

## la overall #################################################################
# data
behav_la       = featmod_coefs$la_LA
behav_la       = behav_la[grep('^pred_LA',names(behav_la))]
behav_la$Group = agk.recode.c(row.names(behav_la),dat_match$VPPG,dat_match$HCPG)

# aggregate
la_overall        = aggregate(behav_la$pred_LA,by=list(behav_la$Group),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
la_overall        = data.frame(la_overall[[1]],la_overall[[2]])
names(la_overall) = c('Group','mean','ci_0025','ci_0975')

# graph
mRat  = ggplot(la_overall, aes(Group, mean))
mRat  = mRat + labs(x='category', y=paste('Mean of LA (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean loss aversion overall per group")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# stats
mod_00 = lm(pred_LA ~ 1,behav_la)
mod_01 = lm(pred_LA ~ Group,behav_la)
anova(mod_00,mod_01)

# stats glmer
modla_00  = glmer(accept_reject ~ 1 + (1|subject) + (1|cat) + (1|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
modla_01  = glmer(accept_reject ~ gain + loss + (gain + loss|subject) + (gain + loss|cat) + (gain + loss|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
modla_0g  = glmer(accept_reject ~ (gain + loss)*HCPG + (gain + loss|subject) + (gain + loss|cat) + (gain + loss|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
modla_cg  = glmer(accept_reject ~ (gain + loss)*HCPG + cat*HCPG + (gain + loss + cat|subject) + (gain + loss|cat) + (gain + loss|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
modla_cgi = glmer(accept_reject ~ (gain + loss)*cat*HCPG + ((gain + loss)*cat|subject) + (gain + loss|cat) + (gain + loss|stim),data = data_pdt,family = 'binomial',nAGQ = 0,control=cur_control)
anova(modla_00,modla_01,modla_0g,modla_cg)
get_la_fixef_pdt(modla_0g)
agk.get.compl.coef(modla_00)

# bootstrapping p-value...
effects_under_0 = agk.boot.p.mermod(mermod = modla_0g,mermod0 = modla_01,num_cpus = 7,num = 7,fun_extract = fixef,cur_control = cur_control,permvars = c('HCPG'),type='perm')

# bootstrapping confint
#boot_cfint_0g <-bootMer(modla_01,use.u = F,FUN = fixef,nsim = 6,type = "parametric",verbose = TRUE,parallel = "no")
boot_cfint_0g = agk.boot.cfint.mermod(mermod = modla_0g,num_cpus = 7,num = 7,fun_extract = fixef,cur_control = cur_control,type = 'non-parametric')


# get the p-value
#setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/effects_under_0_la')
setwd('E:/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/effects_under_0_la')
load('effects_under_0_0g_perm_1000.RData')
cur_vars    <- c("bg_HC","bg_PG","bl_HC", "bl_PG", "x_la_HCgrPG")
extr_params = function(x) {
  res             = list()
  res$bg_HC       = x['gain']
  res$bl_HC       = x['loss']
  res$bg_PG       = x['gain'] + x['gain:HCPGPG']
  res$bl_PG       = x['loss'] + x['loss:HCPGPG']
  res$la_HC       = -x['loss']/x['gain']
  res$la_PG       = -res$bl_PG/res$bg_PG
  res$x_la_HCgrPG = res$la_HC - res$la_PG 
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

con_bg_PG_HC <- bg_HC-bg_PG 
con_bl_PG_HC <- bl_PG-bl_HC 
## p-values for HC gr PG
# LA difference
obs_fixef = get_la_fixef_pdt(modla_0g)
obs_fixef = obs_fixef['x_la_HCgrPG']
1-agk.density_p(x_la_HCgrPG,obs_fixef)
# beta gain difference
obs_fixef = get_la_fixef_pdt(modla_0g)
obs_fixef = obs_fixef['bg_HC'] - obs_fixef['bg_PG']
1-agk.density_p(con_bg_PG_HC,obs_fixef)
# beta loss difference
obs_fixef = get_la_fixef_pdt(modla_0g)
obs_fixef = obs_fixef['bl_PG'] - obs_fixef['bl_HC']
1-agk.density_p(con_bl_PG_HC,obs_fixef)

# graph fixed effects la model
use_cfint_boot = 1
cur_coef       = coef(modla_0g)
cur_coef       = cur_coef$subject
cur_coef$HCPG  = agk.recode.c(row.names(cur_coef),dat_match$VPPG,dat_match$HCPG)
cur_coef$gainc = NA
cur_coef$lossc = NA
cur_coef$gainc = cur_coef$gain
cur_coef$lossc = cur_coef$loss
# with group fixed effect
cur_coef$gainc[cur_coef$HCPG == 'PG'] = cur_coef$gainc[cur_coef$HCPG == 'PG']  + cur_coef$'gain:HCPGPG'[cur_coef$HCPG == 'PG']
cur_coef$lossc[cur_coef$HCPG == 'PG'] = cur_coef$lossc[cur_coef$HCPG == 'PG']  + cur_coef$'loss:HCPGPG'[cur_coef$HCPG == 'PG']
cur_coef$la                           = -cur_coef$lossc/cur_coef$gainc
# new graph (prep)
las     = list()
la_vars = c('la','gainc','lossc')
for (ll in 1:length(la_vars)) {
  las[[ll]]        = aggregate(cur_coef[la_vars[ll]],by=list(cur_coef$HCPG),FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
  las[[ll]]        = data.frame(las[[ll]]$Group.1,las[[ll]][[2]])
  names(las[[ll]]) = c('Group','mean','ci_0025','ci_0975')
  las[[ll]]$var    = la_vars[ll]
}
las_df         = rbind(las[[1]],las[[2]],las[[3]])
la_overall     = melt(las_df,id.vars = c('Group','var','ci_0025',  'ci_0975','mean'))
la_overall$var = factor(la_overall$var,levels = la_vars)
obs_fixef      = get_la_fixef_pdt(modla_0g)
# insert original la fixef (cause they get distorted by the articficially created "random effects")
la_overall$mean[la_overall$Group == 'HC' & la_overall$var == 'la'] = obs_fixef['la_HC']
la_overall$mean[la_overall$Group == 'PG' & la_overall$var == 'la'] = obs_fixef['la_PG']

if(use_cfint_boot) {
  # alternative graph using cfint from fixed effects bootstrap
  cie    = new.env()
  load('boot_cfint_0g_1000.RData',envir = cie)
  rcfint = cie$boot_cfint_0g[[1]]
  for (cc in 2:length(cie$boot_cfint_0g)) {
    rcfint = rbind(rcfint,cie$boot_cfint_0g[[cc]])
  }
  
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
  
  rcfinta         = aggregate(.~group,data=rcfint,FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
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
  # # put in the original fe
  # obs_fixef      = get_la_fixef_pdt(modla_0g)
  # cur_vars       = unique(la_overall$var)
  # cur_fe_HC      = obs_fixef[c(1,3,5,7)]
  # cur_fe_PG      = obs_fixef[c(2,4,6,8)]
  # for (gg in 1:length(cur_vars)) {
  #   la_overall$mean[la_overall$Group == 'HC' & la_overall$var == cur_vars[gg]] = cur_fe_HC[gg]
  #   la_overall$mean[la_overall$Group == 'PG' & la_overall$var == cur_vars[gg]] = cur_fe_PG[gg]
  # }
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


## la under different cue conditions ##########################################
# data
behav_laci                 = featmod_coefs$laci_LA
pLA_ind                    = which(names(behav_laci) == 'pred_LA')
names(behav_laci)[pLA_ind] = 'pred_LA_of_neutral'
names(behav_laci)          = gsub('pred_','',names(behav_laci))
behav_laci                 = behav_laci[grep('^LA',names(behav_laci))]
behav_laci[c(2,3,4)]       = behav_laci[c(2,3,4)] + behav_laci[[c(1)]] # adding neutral
behav_laci$Group           = agk.recode.c(row.names(behav_laci),dat_match$VPPG,dat_match$HCPG)
behav_laci$Group           = agk.recode.c(behav_laci$Group,'PG','GD')

# # simple t-tests (careful, only usable if neutral not added)
# tt_res_bv_laci = lapply(behav_laci,FUN = ttest.group.con,cur_group = dat_match$HCPG, cur_con = dat_match$smoking_ftdt)
# tt_res_bv_laci = unlist(lapply(tt_res_bv_laci,cur_fun))
# tt_res_bv_laci = tt_res_bv_laci[tt_res_bv_laci < 0.05]

# aggregate
behav_laci          = melt(behav_laci)
laci_overall        = aggregate(behav_laci$value,by=list(behav_laci$Group,behav_laci$variable),
                              FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
laci_overall        = data.frame(laci_overall[[1]],laci_overall[[2]],laci_overall[[3]])
names(laci_overall) = c('Group','category','mean','ci_0025','ci_0975')

# graph
mRat  = ggplot(laci_overall, aes(category, mean,fill=Group))
mRat  = mRat + labs(x='category', y=paste('Mean of loss aversion (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Mean loss aversion across categories")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = ci_0025, ymax = ci_0975), position=dodge, width=0.25) + theme_bw()
print(mRat)

# stats
mod_00 = lm(value ~ 1,behav_laci)
mod_01 = lm(value ~ Group,behav_laci)
mod_02 = lm(value ~ variable+Group,behav_laci)
mod_03 = lm(value ~ variable*Group,behav_laci)
anova(mod_00,mod_01,mod_02,mod_03)

## la model under different cue conditions ####################################
setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/effects_under_0_la')

# # bootstrapping under 0
# effects_under_0_cg = agk.boot.p.mermod(mermod = modla_cg,mermod0 = modla_01,num_cpus = 7,num = 28,fun_extract = fixef,cur_control = cur_control)
# save(file="effects_under_0_cg_28.RData",list=c("effects_under_0_cg"))
# 
# # bootstrapping cfint
# boot_cfint_cg = agk.boot.cfint.mermod(mermod = modla_cg,num_cpus = 7,num = 990,fun_extract = fixef,cur_control = cur_control)
# save(file = "boot_cfint_cg_990.RData",list=c("boot_cfint_cg"))

# get the p-value
#setwd('C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/effects_under_0_la')
setwd('E:/Google Drive/Library/R/PDT/analyses/in_progress/sev_pred/results/effects_under_0_la')
load('effects_under_0_cg_perm_1000.RData')
cur_vars    <- c("bg_HC","bg_PG","bl_HC", "bl_PG", "x_la_HCgrPG","catgam_HC","catgam_PG","catneg_HC","catneg_PG",
                 "catpos_HC","catpos_PG","x_in_HCgrPG","x_bg_HCgrPG","x_bl_HCgrPG","x_cg_HCgrPG","x_cp_HCgrPG","x_cn_HCgrPG")
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
  res$cg_PG       = x['catgambling'] + x['HCPGPG:catnegative']
  res$cg_PG       = x['catgambling'] + x['HCPGPG:catpositive']
  
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

con_bg_PG_HC <- bg_HC-bg_PG 
con_bl_PG_HC <- bl_PG-bl_HC 
## p-values for HC gr PG
# LA difference
obs_fixef = get_la_fixef_pdt(modla_0g)
obs_fixef = obs_fixef['x_la_HCgrPG']
1-agk.density_p(x_la_HCgrPG,obs_fixef)
# beta gain difference
obs_fixef = get_la_fixef_pdt(modla_0g)
obs_fixef = obs_fixef['bg_HC'] - obs_fixef['bg_PG']
1-agk.density_p(con_bg_PG_HC,obs_fixef)

# loading boot cfint
all_saved_boots   = dir(pattern = 'boot_cfint_cg_1000.RData')
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

# get the cfint
cfint           = aggregate(.~HCPG,data = coefs,FUN=agk.boot.ci,R=2000,lower=0.025,upper=0.975,cur_fun=mean)
cfint_df        = data.frame(cfint[[1]],cfint[[2]])
cfint_df$var    = var_names[1]
names(cfint_df) = c('Group','mean','ci_0025','ci_0975','var') 
for (cc in 3:length(cfint)) {
  cur_df        = data.frame(cfint[[1]],cfint[[cc]])
  cur_df$var    = var_names[cc-1]
  names(cur_df) = c('Group','mean','ci_0025','ci_0975','var') 
  cfint_df      = rbind(cfint_df,cur_df)
}

# # replace means by real estimated fixed effects
# cur_fe              = fixef(modla_cg)
# cur_fe_HC           = cur_fe[c('(Intercept)','gain','loss','catgambling', 'catnegative','catpositive')]
# cur_fe_PG           = cur_fe[c('HCPGPG','gain:HCPGPG','loss:HCPGPG', 'HCPGPG:catgambling', 'HCPGPG:catnegative','HCPGPG:catpositive')]
# cur_fe_PG           = cur_fe_HC + cur_fe_PG
# cur_fe_PG$HCPG      = 'PG' 
# cur_fe_HC$HCPG      = 'HC' 
# names(cur_fe_HC)[1] = 'Intercept'
# names(cur_fe_PG)    = names(cur_fe_HC)
# cur_fe_df           = rbind(data.frame(cur_fe_HC),data.frame(cur_fe_PG))
# cur_fe_df$la        = -cur_fe_df$loss/cur_fe_df$gain
# cur_fe_df           = melt(cur_fe_df)
# cfint_df$mean       = cur_fe_df$value

# graph fixed effects la model with category values and bootstr cfints
message('using mean of bootstrap instead of fixef values as means')
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




## old ########################################################################


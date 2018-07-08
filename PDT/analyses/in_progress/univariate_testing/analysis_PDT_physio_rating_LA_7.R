## PREAMBLE ===================================================================
# general analysis script building on
# import_data.R
# select_study.R

# compares LA between categories of pictures and between groups overall
# and correlates LA and the change in LA with questionnaires

# for physio and rating variables, it uses median split
# it firstly (A) uses every variable one by one, which may be
# a bad case of multiple testing
# secondly (B) it proposes to use all variables and 
# thridly (C) it works like (B) but performs model selection
# LEDALAB: correlating ledalab peak indices within trials
# with loss aversion there are three different studies here:
# MRI, POSTPILOT_HCPG, POSTPILOT_PGM_PGF (gender), PHYSIOPILOT
# (only HC)

## DATA CLEANING ==============================================================

## TODO: THIS NEEDS TO GO TO data_import.R

# Data cleaning with respect to accept/reject behavior in PDT task
# discard trials with too short RTs
# CAREFUL: in MRI RT means something else! lower threshold!; in PP: 0.5 in in MRI: 0.15
data_pdt_bb3 = data_pdt
data_pdt$subject = as.factor(as.character(data_pdt$subject))
data_pdt$rt[data_pdt$rt == 99999] = NA
data_pdt = subset(data_pdt, rt > 1.3) # USING ONLY FAST TRIALS rt_cut_off
trials_per_sub = as.data.frame(as.matrix(xtabs(~subject,data=data_pdt)))
trials_per_sub$subject = row.names(trials_per_sub)
too_fast_subs = trials_per_sub[trials_per_sub$V1 <120,]
too_fast_subs = too_fast_subs$subject

# check for obvious gambles:
# data_pdt = data.la
# data_pdt$accept_reject = data_pdt$accept.reject
# data_pdt$gain_bcp = data_pdt$Gewinn_bb
# data_pdt$loss_bcp = data_pdt$Verlust_bb

all_subs = unique(data_pdt$subject)
test_gam = c()
bad_subs = c()
na_subs  = c()
fails    = c()

for (ii in 1:length(all_subs)) {
  mg_test      = 0 # has max gain test occured?
  ml_test      = 0 # has max loss test occured?
  cur_test_gam = 0
  cur_fail     = 0
  cur_dat = subset(data_pdt, subject == all_subs[ii])
  if (sum(is.na(cur_dat$accept_reject)) > 20) {
    na_subs = c(na_subs,as.character(all_subs[ii]))
    next
  }
  # get the 12 most extreme gambles
  cur_dat$ratio = abs(cur_dat$gain_bcp/cur_dat$loss_bcp)
  all_ratios    = na.omit(cur_dat$ratio)
  min_ratio     = min(all_ratios)
  max_ratio     = max(all_ratios)
  all_ratios    = all_ratios[-c(which(max_ratio==all_ratios),which(min_ratio==all_ratios))]
  min2_ratio    = min(all_ratios)
  max2_ratio    = max(all_ratios)
  all_ratios    = all_ratios[-c(which(max2_ratio==all_ratios),which(min2_ratio==all_ratios))]
  min3_ratio    = min(all_ratios)
  max3_ratio    = max(all_ratios)
  
  for (jj in 20:length(cur_dat[,1])) {
    if (is.na(cur_dat$accept_reject[jj])) {next}
    if (cur_dat$ratio[jj] == max_ratio | cur_dat$ratio[jj] == max2_ratio | cur_dat$ratio[jj] == max3_ratio) {
      cur_test_gam = cur_test_gam + 1
      if (cur_dat$accept_reject[jj] == 0) {cur_fail = cur_fail + 1}
      next
    }
    if (cur_dat$ratio[jj] == min_ratio | cur_dat$ratio[jj] == min2_ratio | cur_dat$ratio[jj] == min2_ratio) {
      cur_test_gam = cur_test_gam + 1
      if (cur_dat$accept_reject[jj] == 1) {cur_fail = cur_fail + 1}
      next
    }
  }
  test_gam[ii] = cur_test_gam
  fails   [ii] = cur_fail
}

# if they have failed in more than 80 percent of the test gambles
fraction_fails = fails/test_gam
for (kk in 1:length(fraction_fails)) {
  if (fraction_fails[kk] >=0.8) {
    bad_subs = c(bad_subs,as.character(all_subs[kk]))
  }
}

# clean df of bad subs, too fast subs, na subs
data_pdt_bb         = data_pdt

# bad, na, and fast out
data_pdt_nobds      = data_pdt[(!data_pdt$subject %in% bad_subs),]
data_pdt_nobds_nona = data_pdt_nobds[(!data_pdt_nobds$subject %in% na_subs),]
data_pdt_nobds_nona = data_pdt_nobds_nona[(!data_pdt_nobds_nona$subject %in% too_fast_subs),]
data_pdt = data_pdt_nobds_nona

# too fast and na are taken out
# data_pdt = data_pdt[(!data_pdt$subject %in% too_fast_subs),]
# data_pdt = data_pdt[(!data_pdt$subject %in% na_subs),]

## WHAT ABOUT THIS, TAKE IT OR LEAVE IT?! ##
# taking out subs with task cor of gain or loss with cat 
tmp_g = lmList(gain ~ cat|subject,data=data_pdt)
tmp_l = lmList(loss ~ cat|subject,data=data_pdt)
tmp   = c(tmp_g,tmp_l)
unr_subs = c()
for(ii in 1:length(tmp)) {
  cur_sum = summary(tmp[[ii]])$coefficients[,4][c(2:4)]
  #cur_sum = anova(tmp[[ii]])[,5][1]
  if (any(cur_sum < 0.05)) {
    unr_subs = c(unr_subs,ii)
  }
}
unr_subs = (unique(names(tmp)[unr_subs]))
data_pdt = data_pdt[(!data_pdt$subject %in% unr_subs),]

data_pdt$HCPG = droplevels(data_pdt$HCPG)

## check model fit
# taking out subs with task cor of gain or loss with cat 
tmp_la_model = lmList(accept_reject ~ gain+loss|subject,data=data_pdt,family="binomial",pool=F,na.action = NULL)
unf_subs = c()
for(ii in 1:length(tmp_la_model)) {
  cur_sum = summary(tmp_la_model[[ii]])$coefficients[,4][c(2:3)]
  #cur_sum = anova(tmp[[ii]])[,5][1]
  if (!any(cur_sum < 0.1)) {
    unf_subs = c(unf_subs,ii)
  }
}
unf_subs = (unique(names(tmp_la_model)[unf_subs]))
data_pdt = data_pdt[(!data_pdt$subject %in% unf_subs),]

## how many subjects now?
tmp = aggregate(data_pdt$subject, by=list(data_pdt$subject,data_pdt$HCPG), FUN=first)
xtabs(~Group.2,data = tmp)

## PRELIMNARY ANALYSES ========================================================

# ACCEPTANCE RATE
accbygroup        <- as.data.frame(aggregate(as.numeric(as.character(data_pdt$accept_reject)), by=list(data_pdt$subject,data_pdt$HCPG),FUN="mean.rmna"))
names(accbygroup) <- c("subject", "HCPG", "acceptance_rate")
kruskal.test(accbygroup$acceptance_rate~accbygroup$HCPG)
describeBy(accbygroup$acceptance_rate,accbygroup$HCPG)

# RATING and PHYSIO by category
# cpmputes mean per category and subject (using lmlist)
# then performs group comparison using f_difftest
# prints only if significant
for (ii in 3:(length(pic_value_vars)-1)) {
  cur_formula = paste(pic_value_vars[ii], " ~ (0+cat_rec) | subject")
  ratbycatv   = lmList(as.formula(cur_formula), data = data_pdt,na.action=NULL,pool = F)
  cratbycatv  = coef(ratbycatv)
  cratbycatv$HCPG = agk.recode.c(row.names(cratbycatv),as.character(data_pdt$subject),as.character(data_pdt$HCPG))
  cratbycatv$HCPG = as.factor(cratbycatv$HCPG)
  cur_cons = names(cratbycatv)
  tmp = list()
  for(kk in 1:(length(cur_cons)-1)) {
    if (which_study == "PhysioPilot"){
      tmp[[kk]] = f.difftest(cratbycatv,cur_cons[kk],"")
    } else {
      tmp[[kk]] = f.difftest(cratbycatv,cur_cons[kk],"HCPG")
    }
  }
  cur_cnames = names(tmp[[1]])
  cur_nrow   = length(tmp)
  cur_ncol   = length(tmp[[1]])
  tmp = as.numeric(unlist(tmp))
  tmp = matrix(tmp,nrow=cur_nrow,byrow = T)
  colnames(tmp) = cur_cnames
  tmp = as.data.frame(tmp)
  disp(" ")
  disp(paste("results for... ",pic_value_vars[ii]))
  we_have_sig = 0
  for(ll in 1:length(tmp$matched)) {
    if (which_study != "PhysioPilot") {
      if(tmp$k_matched[ll] == 0 & tmp$t_matched[ll] == 0) {
        we_have_sig = 1
        break
      }
    } else {
      if(tmp$'t_p-value' <= 0.05) {
        we_have_sig = 1
        break
      }
    }
  }
  if (we_have_sig == 1) {
    print(tmp)
  } else {
    disp("No sig. group differences!")
  }
}

# imageRating5 only in people who answered lie-bet positively
agk.load.ifnot.install("lmerTest")
data_pdt_PG = subset(data_pdt,HCPG=="PG")
tmp = lmer(imageRating1 ~ (0+cat) + (0+cat |subject), data = data_pdt_PG)
tmp_coef = coef(tmp_2)
tmp_coef = tmp_coef$subject
names(tmp_coef) = c("neutral", "gambling", "negative", "positive")
tmp_coef = melt(tmp_coef)
p <- ggplot(data = tmp_coef, aes(x = variable, y = value))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=0.2)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <- p + xlab("picture category") + ylab(expression("rating score (0-100) [mean with 95% boots CI]"))
p <- p + ggtitle('Reactions of problem gamblers to question:\n"How much craving for gambling does this picture induce in you?"')
p <- p + coord_cartesian(ylim = c(0, 63)) 
p <- p +theme_bw()
p

data_pdt_HC = subset(data_pdt,HCPG=="HC")
tmp = lmer(imageRating1 ~ (0+cat) + (0+cat |subject), data = data_pdt_HC)
tmp_coef = coef(tmp)
tmp_coef = tmp_coef$subject
names(tmp_coef) = c("neutral", "gambling", "negative", "positive")
tmp_coef = melt(tmp_coef)
p <- ggplot(data = tmp_coef, aes(x = variable, y = value))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=0.2)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <- p + xlab("picture category") + ylab(expression("rating score (0-100) [mean with 95% boots CI]"))
p <- p + ggtitle('Reactions of problem gamblers to question:\n"How much craving for gambling does this picture induce in you?"')
p <- p + coord_cartesian(ylim = c(0, 63)) 
p <- p +theme_bw()
p

## MODEL FITTING FOR FIRST COHORT =============================================

data_pdt         = data_pdt
data_pdt$subject = droplevels(data_pdt$subject)
data_pdt$rt[data_pdt$rt == 99999] = NA

# bl model, la model
bl_mod   = glmer(accept_reject ~ 1 + (1 | subject), data = data_pdt,family="binomial") # bl
la_mod   = glmer(accept_reject ~ gain+loss+ed_abs + (gain+loss+ed_abs| subject),data = data_pdt,family = "binomial") # la model
anova(bl_mod,la_mod)

# la model vs. lac-model
lac_mod  = glmer(accept_reject ~ (gain+loss+cat) + (gain+loss+cat|stim) + (gain+loss+cat|subject),data = data_pdt,family = "binomial",
                 nAGQ = 0,control=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                                                       check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=500)))

lacg_mod = glmer(accept_reject ~ (gain+loss+ed_abs+cat)*HCPG + (gain+loss+ed_abs+cat|stim) + (gain+loss+ed_abs+cat|subject),data = data_pdt,family = "binomial",
                 nAGQ = 0,control=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                                               check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=500)))

laec_mod = glmer(accept_reject ~ (gain+loss+ed_abs)+cat*HCPG+(gain+loss+ed_abs+cat| subject),data = data_pdt,family = "binomial",
                 nAGQ = 0,control=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                                               check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=500)))

anova(la_mod,lac_mod)
summary(lac_modc)

# la model vs. lac-model
laci_mod = glmer(accept_reject ~ (gain+loss+ed_abs)*cat+((gain+loss+ed_abs)*cat| subject),data = data_pdt,family = "binomial",
                         nAGQ = 0,control=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                                                       check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=700))) # la model
anova(lac_mod,laci_mod)
summary(laci_mod)

# RT?

## COMPLETE MODELING (B) (glmnet) =============================================
## PLUS MODEL SELECTION

f_HCPG_cm  <- function(x) {
  tmp <- kruskal.test(x,coef_mod$HCPG)
  med_PG = median(x[coef_mod$HCPG == "PG"],na.rm = T)
  med_HC = median(x[coef_mod$HCPG == "HC"],na.rm = T)
  return(as.matrix(c(tmp$p.value,med_PG-med_HC,med_PG,med_HC)))}

# Latency is a "bad" variable, in more than half of the cases it is NaN
la_lmlist_b <- lmList(accept_reject ~ (gain+loss) | subject, data = data_pdt,family="binomial",na.action=NULL,pool = F)
la_lmlist_c <- lmList(accept_reject ~ (gain+loss)*(cat_rec) | subject, data = data_pdt,family="binomial",na.action=NULL,pool = F)
la_lmlist_p <- lmList(accept_reject ~ (gain+loss)*(corr_auc+zygo+SCR_stim+AmpSum_stim) | subject, data = data_pdt,family="binomial",na.action=NULL,pool = F)
la_lmlist_r <- lmList(accept_reject ~ (gain+loss)*(valence+
                                                     arousal+dominance+imageRating1s+imageRating2s+imageRating3s+imageRating4s) | subject,
                      data = data_pdt,family="binomial",na.action=NULL,pool = F) # CAREFUL: PG people have imageRating5s also ("does this pic make u question gambling")

# pruning
la_lmlist_bp <- lapply(la_lmlist_b,FUN=agk.glmnet.lm.cvalpha,fam_arg="binomial",type_measure='class',lambda=(30*(1/c(1:100)))) # the internal lambda generator does not work here
la_lmlist_cp <- lapply(la_lmlist_c,FUN=agk.glmnet.lm.cvalpha,fam_arg="binomial",type_measure='class',lambda=NULL)
la_lmlist_pp <- lapply(la_lmlist_p,FUN=agk.glmnet.lm.cvalpha,fam_arg="binomial",type_measure='class',lambda=NULL)
la_lmlist_rp <- lapply(la_lmlist_r,FUN=agk.glmnet.lm.cvalpha,fam_arg="binomial",type_measure='class',lambda=NULL)

pruned_models = list(la_lmlist_bp=la_lmlist_bp,la_lmlist_cp=la_lmlist_cp,la_lmlist_pp=la_lmlist_pp,la_lmlist_rp=la_lmlist_rp)
if (which_study == "MRT") {
  pruned_models = list(la_lmlist_bp=la_lmlist_bp,la_lmlist_cp=la_lmlist_cp,la_lmlist_rp=la_lmlist_rp)
}
cvme_vecs     = list()
group_vecs    = list()

for (pp in 1:length(pruned_models)) {
  
  # get the coefs completely from a glment pruned list
  cur_lmlist_pglmnet <- pruned_models[[pp]]
  take_in = c()
  for (ii in 1:length(cur_lmlist_pglmnet)) {
    if (is.null(cur_lmlist_pglmnet[[ii]]$coef)) {
      take_in[ii] = 0
    } else {
      take_in[ii] = 1
    }
  }
  cur_lmlist_pglmnet = cur_lmlist_pglmnet[as.logical(take_in)]
  
  cur_lmlist_pglmnet_bcp = cur_lmlist_pglmnet
  cur_lmlist_pglmnet     = list()
  cur_lmlist_pglmnet_cvme = c()
  for(jj in 1:length(cur_lmlist_pglmnet_bcp)) {
    cur_lmlist_pglmnet[[jj]] = cur_lmlist_pglmnet_bcp[[jj]]$coef
    cur_lmlist_pglmnet_cvme[jj] = cur_lmlist_pglmnet_bcp[[jj]]$cvme
  }
  
  coef_mod           = as.data.frame(matrix(unlist(cur_lmlist_pglmnet),nrow=length(cur_lmlist_pglmnet),byrow=T))
  names(coef_mod)    = colnames(cur_lmlist_pglmnet[[1]])
  
  # make LA columns
  cur_loss_cols      = coef_mod[names(coef_mod)[grep(names(coef_mod),pattern="^loss")]]
  cur_gain_cols      = coef_mod[names(coef_mod)[grep(names(coef_mod),pattern="^gain")]]
  
  cur_LA_cols        = cur_loss_cols*(-1)/cur_gain_cols
  cur_names          = names(cur_LA_cols)
  tmp                = sapply(cur_names,FUN=strsplit,split="loss",fixed = T)
  cur_names[1]       = paste0("LA",tmp[[1]][1])
  if (length(cur_names) > 1) {
    for (kk in 2:length(cur_names)) {
      cur_names[kk]    = paste0("LA",tmp[[kk]][2])
    }
  }
  
  names(cur_LA_cols) = cur_names
  coef_mod           = cbind(coef_mod,cur_LA_cols)
  
  # t-tests
  coef_mod$subject   = names(cur_lmlist_pglmnet_bcp)
  coef_mod$HCPG      = agk.recode.c(coef_mod$subject,as.character(data_pdt$subject),as.character(data_pdt$HCPG))
  coef_mod$HCPG[as.logical(unlist(lapply(lapply(as.list(coef_mod$HCPG),grep,pattern="VPPG"),length)))]=NA
  coef_mod$HCPG      = as.factor(coef_mod$HCPG)
  coef_modres        = sapply((coef_mod[1:(length(coef_mod)-2)]),FUN=f_HCPG_cm) ## INCLUDE TEST FOR SIGNIFICANCE
  
  group_vecs[[pp]]   = coef_mod$HCPG
  
#   # correlation
#   # prepare frames for correlation
#   severity                   = severity_bcp
#   coef_mod                   = coef_mod[which(coef_mod$subject %in% severity$subject),]
#   severity                   = severity[which(coef_mod$subject %in% severity$subject),]
#   severity                   = severity[which(severity$subject %in% coef_mod$subject),]
#   coef_mod                   = coef_mod[which(severity$subject %in% coef_mod$subject),]
#   coef_mod_all               = coef_mod
#   coef_mod                   = subset(coef_mod, HCPG=="PG")
#   severity$HCPG              = as.factor(agk.recode.c(as.character(severity$subject),as.character(dat_match$VPPG),as.character(dat_match$HCPG)))
#   severity_all               = severity
#   severity                   = subset(severity,HCPG == "PG")
#   
#   severity$subject           = NULL
#   coef_mod$subject           = NULL
#   coef_mod$HCPG              = NULL
#   severity$HCPG              = NULL
#   
#   # correlate
#   M = corr.test(coef_mod,severity,method = "spearman",adjust = "none")
#   M_r = as.matrix(M$r)
#   M_p = as.matrix(M$p)
#   M_r[M_p >= 0.01] = 0
#   corrplot(M_r,method = "circle",diag = F)
#   title(cur_var)
  
  # record the cvme of this model
  cvme_vecs[[pp]] = cur_lmlist_pglmnet_cvme
}

## WHICH IS THE BEST LA MODULATION MODEL IN EACH GROUP? =======================
cur_names = names(pruned_models)
list_names= list()
for(ii in 1:length(cur_names)) {
  list_names[[ii]] = rep(cur_names[ii],length(cvme_vecs[[ii]]))
}
models          = unlist(list_names)
cvmes           = unlist(cvme_vecs)
group           = unlist(group_vecs)
modcvmes        = data.frame(models,cvmes,group)
modcvmes$models = as.factor(modcvmes$models)
# throw out subjects who have bad model fit
# modcvmes = modcvmes[modcvmes$model_fail == FALSE,]
summary(lm(cvmes ~ models,subset(modcvmes,group=="HC")))
summary(lm(cvmes ~ models,subset(modcvmes,group=="PG")))
summary(lm(cvmes ~ models))

## LMLIST MODELING ============================================================
# effect of category on acceptance rate 
mod00 = lmList(accept_reject ~ 1 | subject, data=data_pdt, family = "binomial", na.action = NULL, pool = F)
mod01 = lmList(accept_reject ~ cat | subject, data=data_pdt, family = "binomial", na.action = NULL, pool = F)

coef_acc=coef(mod01)
coef_acc$HCPG = agk.recode(row.names(coef_acc),as.character(dat_match$VPPG),as.character(dat_match$HCPG))
## NEEDS TO GO CHECK WITH MILAN
coef_acc$HCPG[row.names(coef_acc) == "VPPG0100"] = "HC"
table(coef_acc$HCPG)
cur_mod = lm(catgambling ~ HCPG, data=coef_acc)
kruskal.test(coef_acc$catgambling, as.factor(coef_acc$HCPG))
shapiro.test(resid(tmp))
agk.normality.check(tmp)

# bootstrap p-value
cur_fun = function(d,indices) {
  cur_d = d[indices,]
  fit9              = lm(catgambling ~ 1, data=cur_d) # density under assumption of null-hypothesis
  coef(fit9)
}
tmp = boot(coef_acc,cur_fun,R=5000)
obs_fixef  <- as.numeric(cur_mod$coefficients[1] + cur_mod$coefficients[2])
nulldens   <- density(tmp$t[,1],bw = c("ucv"),n = 1000)
nulldens$y <- nulldens$y/sum(nulldens$y)
nulldens$y <- cumsum(nulldens$y)
f <- approxfun(nulldens, rule=2)
1-f(obs_fixef)

# bootstrap p-value using HC only as null hypothesis distribution
coef_acc_HC = subset(coef_acc, HCPG=="HC")
tmp = boot(coef_acc_HC,cur_fun,R=5000)
obs_fixef  <- as.numeric(cur_mod$coefficients[1] + cur_mod$coefficients[2])
nulldens   <- density(tmp$t[,1],bw = c("ucv"),n = 1000)
nulldens$y <- nulldens$y/sum(nulldens$y)
nulldens$y <- cumsum(nulldens$y)
f <- approxfun(nulldens, rule=2)
1-f(obs_fixef)

# bootstrap p-value checking probability that HC mean comes from the unified distribution
tmp = boot(coef_acc,cur_fun,R=5000)
obs_fixef  <- as.numeric(cur_mod$coefficients[1])
nulldens   <- density(tmp$t[,1],bw = c("ucv"),n = 1000)
nulldens$y <- nulldens$y/sum(nulldens$y)
nulldens$y <- cumsum(nulldens$y)
f <- approxfun(nulldens, rule=2)
f(obs_fixef)


## COMPLETE MODELING (C) (glmer) ==============================================
# now starting the modeling
## check the gamble matrix: has there been some undue oversampling of middle gambles?

# bl model, la model
cur_control = glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=300))
laec_m      = glmer(accept_reject ~ (gain+loss+ed_abs+cat)*HCPG+(gain+loss+ed_abs+cat|subject)+(gain+loss+ed_abs|stim),data = data_pdt,family = "binomial",nAGQ = 0,control=cur_control) # la model
summary(laec_m)

## RT =========================================================================
cur_control = lmerControl(check.conv.grad="ignore",check.conv.singular="ignore",check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=300))
rtec_m      = lmer(rt ~ (gain+loss+ed_abs+cat)*HCPG+(gain+loss+ed_abs+cat|subject)+(gain+loss+ed_abs|stim),data = data_pdt,nAGQ = 0,control=cur_control) # la model
rtec_m      = lmer(rt ~ (gain+loss+ed_abs+cat*accept_reject)*HCPG+(gain+loss+ed_abs+cat*accept_reject|subject)+(gain+loss+ed_abs+accept_reject|stim),data = data_pdt,nAGQ = 0,control=cur_control) # la model
agk.lme.summary(rtec_m)

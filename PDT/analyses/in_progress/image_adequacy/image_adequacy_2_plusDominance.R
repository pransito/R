## PREAMBLE ===================================================================
# checks all hypotheses of image adequacy of 
# ratings and physio

# run import_data.R before
# run select study_before.R
# TODO: add (cat|stim) random effect

## PREPARATIONS ===============================================================
# what segment to do?
do_ratings   = 1
do_physio    = 0
do_bootstrap = 0
do_leda      = 0

cur_control = lmerControl(check.conv.grad='ignore',check.conv.singular='ignore',
                          check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=200))

# what summary stat?
# TODO: do all
sum_stat     = 'max'

# which contrasts to use?
# standard (default of R) or special (see section to set special)
standard_contrast = 1

# intall/load lme4
# install.packages('lme4')
# install.packages('pracma')
# install.packages('multcomp')
library(lme4)
library(pracma)
library(multcomp)

# contrasts
cur_con           = contr.treatment(length(levels(data_pdt$cat)),base = which(levels(data_pdt$cat) == 'neutral'))
rownames(cur_con) = levels(data_pdt$cat)
colnames(cur_con) = levels(data_pdt$cat)[FALSE == (levels(data_pdt$cat) == 'neutral')]
contrasts(data_pdt$cat) = cur_con

# HC and PG selection
data_pdt_HC = subset(data_pdt, HCPG == 'HC')
data_pdt_PG = subset(data_pdt, HCPG == 'PG')

## IMAGE ADEQUACY RATINGS =====================================================
if (do_ratings == 1) {
  # check the sample
  all_subs         = levels(data_pdt$subject)
  rat_incl         = c()
  
  # check whether ratings are there, if not: exclude
  # TODO: use lmlist for this! needs to be done for all rating scales!
  
  # describe the sampl
  help(hist)
  cur_n         = length(subs_taken)
  
  # split data
  data_pdt_HC = subset(data_pdt, HCPG=="HC")
  data_pdt_PG = subset(data_pdt, HCPG=="PG")
  
  # test ratings hypotheses
  rav0              = lmer(valence ~ 1 + (1| subject) + (1|stim),data = data_pdt,REML = F)
  ravc              = lmer(valence ~ (0+cat) + (0+cat| subject) + (0+cat|stim),data = data_pdt, REML = F, nAGQ = 0,control=cur_control)
  lmList(valence ~ 0+cat| subject,data = data_pdt)
  ravcg             = lmer(valence ~ (cat*HCPG) + (cat|subject) + (cat|stim),data = data_pdt,REML = F,nAGQ = 0,control=cur_control)
  anova(rav0,ravc,ravcg)
  anova(rav0,ravc)
  agk.lme.summary(rav0,type='norm')
  agk.lme.summary(ravc,type='norm')
  agk.lme.summary(ravcg,type='norm')  
  agk.normality.check(ravcg)
  
  raa0              = lmer(arousal ~ 1 + (1|subject) + (1|stim),data = data_pdt,REML = F)
  raac              = lmer(arousal ~ (cat) + (cat| subject) + (cat|stim),data = data_pdt,REML = F)
  raacg             = lmer(arousal ~ (cat)*HCPG + (cat| subject) + (cat|stim),data = data_pdt,REML = F)
  anova(raa0,raac,raacg)
  agk.lme.summary(raa0,type='norm') 
  agk.lme.summary(raac,type='norm') 
  agk.lme.summary(raacg,type='norm')  
  agk.normality.check(raacg)
  
  
  
  # additional contrasts needed
  cons = glht(raac, linfct = mcp(cat = "Tukey"))
  summary(cons,test=adjusted('none'))
  
  
  #Effects of images on dominance
  rad0              = lmer(dominance ~ 1 + (1| subject) + (1|stim),data = data_pdt,REML = F)
  radc              = lmer(dominance ~ (cat) + (cat| subject) + (cat|stim),data = data_pdt, REML = F, nAGQ = 0,control=cur_control)
  radc_PG          = lmer(dominance ~ cat + (cat| subject) + (cat|stim), data = data_pdt_PG,REML=F,control=cur_control)       
  radc_HC          = lmer(dominance ~ cat + (cat| subject) + (cat|stim), data = data_pdt_HC,REML=F,control=cur_control)
  radcg             = lmer(dominance ~ (cat*HCPG) + (cat|subject) + (cat|stim),data = data_pdt,REML = F,nAGQ = 0,control=cur_control)
  anova(rad0,radc,radcg)
  agk.lme.summary(rad0,type='norm')
  agk.lme.summary(radc,type='norm')
  agk.lme.summary(radc_PG,type='norm')
  agk.lme.summary(radc_HC,type='norm')
  agk.lme.summary(radcg,type='norm') 
  

  
  agk.normality.check(radcg)
  
  raa0              = lmer(arousal ~ 1 + (1|subject) + (1|stim),data = data_pdt,REML = F)
  raac              = lmer(arousal ~ (cat) + (cat| subject) + (cat|stim),data = data_pdt,REML = F)
  raacg             = lmer(arousal ~ (cat)*HCPG + (cat| subject) + (cat|stim),data = data_pdt,REML = F)
  anova(raa0,raac,raacg)
  agk.lme.summary(raa0,type='norm') 
  agk.lme.summary(raac,type='norm') 
  agk.lme.summary(raacg,type='norm')  
  agk.normality.check(raacg)
  
  # Effect of gambling category on rating "Representative for gambling"
  if (standard_contrast == 1) {
    #cur_con = contrasts(as.factor(as.character(data_pdt$cat)))
  } else {
    cur_con           = matrix(c(1,-1,0,1,0,-1),nrow = 3,byrow = F)
    colnames(cur_con) = c("gam>neg","gam>pos") 
    rownames(cur_con) = c("gam","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  }
  
  rarg0          = lmer(imageRating2s ~ 1 + (1| subject) + (1|stim),data = data_pdt,REML=F) # representative of gamble
  rargc          = lmer(imageRating2s ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt,REML=F) # representative of gamble
  rargcg         = lmer(imageRating2s ~ (cat*HCPG) + (cat| subject) + (cat| stim),data = data_pdt,REML=F) # representative of gamble
  anova(rarg0,rargc,rargcg)
  agk.lme.summary(rarg0,type='norm')  
  agk.lme.summary(rargc,type='norm')  
  agk.lme.summary(rargcg,type='norm')  
  


  # Effect of gambling category on rating "craving for gambling" [better to do within PG group]
  if (standard_contrast == 1) {
    #cur_con = contrasts(as.factor(as.character(data_pdt$cat)))
  } else {
    cur_con           = matrix(c(-1,1,0,0,0,1,-1,0,0,1,0,-1),nrow = 4,byrow = F)
    colnames(cur_con) = c("gam>neu","gam>neg","gam>pos") 
    rownames(cur_con) = c("neu","gam","neg","pos") 
    contrasts(data_pdt$cat) = cur_con
  }
  
  rarc0             = lmer(imageRating1s ~ 1 + (1| subject) + (1|stim),data = data_pdt,REML=F,control=cur_control)              # craving
  rarcc             = lmer(imageRating1s ~ cat + (cat| subject) + (cat|stim), data = data_pdt,REML=F,control=cur_control)       # craving
  rarcc_PG          = lmer(imageRating1s ~ cat + (cat| subject) + (cat|stim), data = data_pdt_PG,REML=F,control=cur_control)       # craving
  rarcc_HC          = lmer(imageRating1s ~ cat + (cat| subject) + (cat|stim), data = data_pdt_HC,REML=F,control=cur_control)       # craving
  rarccg            = lmer(imageRating1s ~ (cat*HCPG) + (cat| subject) + (cat|stim),data = data_pdt,REML=F,control=cur_control) # craving
  anova(rarc0,rarcc,rarccg)
  agk.lme.summary(rarc0,type='norm')  
  agk.lme.summary(rarcc,type='norm')  
  agk.lme.summary(rarccg,type='norm') 
  agk.lme.summary(rarcc_PG,type='norm') 
  
  agk.(rarccg,type='norm') 
  agk.normality.check(rarccg)
  
  # additional contrasts needed
  cons = glht(rarcc_PG, linfct = mcp(cat = "Tukey"))
  summary(cons,test=adjusted('none'))
  
  cons2 = glht(rarccg, linfct = mcp(cat = "Tukey"))
  summary(cons2,test=adjusted('none'))
  
  cons3 = glht(rarcc, linfct = mcp(cat = "Tukey"))
  summary(cons3,test=adjusted('none'))
  
  cons4 = glht(rarcc_HC, linfct = mcp(cat = "Tukey"))
  summary(cons4,test=adjusted('none'))
  
 

  
  # Effect of positive category on rating "representative for positive c."
  if (standard_contrast == 1) {
    #cur_con = contrasts(as.factor(as.character(data_pdt$cat)))
  } else {
    cur_con           = matrix(c(-1,0,0,1,0,-1,0,1,0,0,-1,1),nrow = 4,byrow = F)
    colnames(cur_con) = c("pos>neu","pos>gam","pos>neg") 
    rownames(cur_con) = c("neu","gam","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  }
  
  if (which_study == "sanity" & standard_contrast == 0) {
    cur_con           = matrix(c(-1,0,1,0,-1,1),nrow = 3,byrow = F)
    colnames(cur_con) = c("pos>neu","pos>neg") 
    rownames(cur_con) = c("neu","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  } else if (which_study == "Prestudy" & standard_contrast == 0) {
    cur_con           = matrix(c(-1,0,1,0,-1,1),nrow = 3,byrow = F)
    colnames(cur_con) = c("pos>gam","pos>neg") 
    rownames(cur_con) = c("gam","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  }
  
  rarp0             = lmer(imageRating4s ~ (1) + (1| subject) + + (1| stim),data = data_pdt,REML=F) # representative of gamble
  rarpc             = lmer(imageRating4s ~ (cat) + (cat| subject) + + (cat| stim),data = data_pdt,REML=F) # representative of gamble
  rarpcg            = lmer(imageRating4s ~ (cat)*HCPG + (cat| subject) + + (cat| stim),data = data_pdt,REML=F) # representative of gamble
  anova(rarp0 ,rarpc,rarpcg)
  agk.lme.summary(rarp0,type='norm') 
  agk.lme.summary(rarpc,type='norm') 
  agk.lme.summary(rarpcg,type='norm') 
  agk.normality.check(rarpcg)
  
  # Effect of neg. category on rating "representative for negative c."
  if (standard_contrast == 1) {
    #cur_con = contrasts(as.factor(as.character(data_pdt$cat)))
  } else {
    cur_con           = matrix(c(-1,0,1,0,0,-1,1,0,0,0,1,-1),nrow = 4,byrow = F)
    colnames(cur_con) = c("neg>neu","neg>gam","neg>pos") 
    rownames(cur_con) = c("neu","gam","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  }
  
  if (which_study == "sanity" & standard_contrast == 0) {
    cur_con           = matrix(c(-1,1,0,0,1,-1),nrow = 3,byrow = F)
    colnames(cur_con) = c("neg>neu","neg>pos") 
    rownames(cur_con) = c("neu","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  } else if (which_study == "Prestudy" & standard_contrast == 0) {
    cur_con           = matrix(c(-1,1,0,0,1,-1),nrow = 3,byrow = F)
    colnames(cur_con) = c("neg>gam","neg>pos") 
    rownames(cur_con) = c("gam","neg","pos")
    contrasts(data_pdt$cat) = cur_con
  }
  
  rarn0             = lmer(imageRating3s ~ (1) + (1| subject) + (1|stim),data = data_pdt,REML=F) # representative of neg
  rarnc             = lmer(imageRating3s ~ (cat) + (cat| subject) + (cat|stim),data = data_pdt,REML=F) # representative of neg
  rarnc_PG          = lmer(imageRating3s ~ cat + (cat| subject) + (cat|stim), data = data_pdt_PG,REML=F,control=cur_control)       #representative of neg 
  rarnc_HC          = lmer(imageRating3s ~ cat + (cat| subject) + (cat|stim), data = data_pdt_HC,REML=F,control=cur_control)       #representative of neg 
  rarncg            = lmer(imageRating3s ~ (cat)*HCPG + (cat| subject) + (cat|stim),data = data_pdt,REML=F) # representative of neg
  anova(rarn0,rarnc,rarncg)
  agk.lme.summary(rarnc,type='norm') 
  agk.lme.summary(rarncg,type='norm') 
  agk.normality.check(rarpcg)
  
  # additional contrasts needed
  cons = glht(rarnc, linfct = mcp(cat = "Tukey"))
  summary(cons,test=adjusted('none'))
  
  cons2 = glht(rarnc_PG, linfct = mcp(cat = "Tukey"))
  summary(cons2,test=adjusted('none'))
  
  summary(cons3,test=adjusted('none'))  cons3 = glht(rarnc_HC, linfct = mcp(cat = "Tukey"))
  
  # how much does this image make you question your gambling behavior
  ramq_lml_PG       = lmList(imageRating5s ~ cat | subject,data = data_pdt_PG)
  ramq0             = lmer(imageRating5s ~ (1) + (1| subject) + (1|stim),data = data_pdt_PG,REML=F)
  ramqc             = lmer(imageRating5s ~ (cat) + (cat| subject) + (cat|stim),data = data_pdt_PG,REML=F) 
  ramqcg            = lmer(imageRating5s ~ (cat)*HCPG + (cat| subject) + (cat|stim),data = data_pdt_PG,REML=F) 
  anova(ramq0,ramqc)
  agk.lme.summary(ramqc,type='norm') 

  agk.normality.check(ramqcg)
  
  # additional contrasts needed
  cons = glht(ramqc, linfct = mcp(cat = "Tukey"))
  summary(cons,test=adjusted('none'))


  
  
  
  ## break ##
  
  # Q: is the neutral image set AW IAPS (cat==6) different than 0
  # should not be, and this across all groups
  # careful: only sanity had the neu_AW before behav study
  if (which_study == "sanity" | which_study == "POSTPILOT_HC" | which_study == "POSTPILOT_PG" |
      which_study == "POSTPILOT_HCPG") {
    data_pdt_iapsAW     = data_pdt[data_pdt$cat=='neutral',]
    data_pdt_iapsAW$cat = as.character(data_pdt_iapsAW$cat)
    
    # make the ideal neutral ratings
    all_subs = unique(data_pdt_iapsAW$subject)
    res_list = list()
    for (ii in 1:length(all_subs)) {
      cur_dat        = data_pdt_iapsAW[data_pdt_iapsAW$subject == all_subs[ii],]
      cur_dat$cat    = as.character(cur_dat$cat)
      ideal_val      = rnorm(length(cur_dat[,1]),mean = 0,sd = sd(cur_dat$valence,na.rm = T))
      ideal_dom      = rnorm(length(cur_dat[,1]),mean = 0,sd = sd(cur_dat$dominance,na.rm = T))
      ideal_aro      = rnorm(length(cur_dat[,1]),mean = -0.8,sd = sd(cur_dat$arousal,na.rm = T)) # set ideal mean to what?
      new_df         = cur_dat
      new_df[c(1:length(new_df))] = NA
      new_df$subject = all_subs[ii]
      new_df$cat     = "99"
      new_df$valence = ideal_val
      new_df$dominance = ideal_dom
      new_df$arousal = ideal_aro
      cur_dat        = rbind(cur_dat,new_df)
      res_list[[ii]] = cur_dat
    }
    now_df = res_list[[1]]
    for (ii in 2:length(res_list)) {
      now_df = rbind(now_df,res_list[[ii]])
    }
    
    data_pdt_iapsAW = now_df
    data_pdt_iapsAW$cat = factor(data_pdt_iapsAW$cat,levels=c("neutral","99"),labels=c("neu_AW","neu_ideal")) 
    
    # model fitting
    ravnu = lmer(valence ~ cat + (cat|subject),data=data_pdt_iapsAW)
    radnu = lmer(dominance ~ cat + (cat|subject),data=data_pdt_iapsAW)
    raanu = lmer(arousal ~ cat + (cat|subject),data=data_pdt_iapsAW)
    agk.lme.summary(ravnu,type='norm')
    agk.lme.summary(radnu,type='norm')
    agk.lme.summary(raanu,type='norm')
  }
  
  # Q: is the neutral image set AW IAPS (cat==6) overall less affect soliciting then the positive and negative images
  data_pdt$affLength=sqrt(data_pdt$arousal^2+data_pdt$dominance^2+data_pdt$valence^2)
  data_pdt_aff = data_pdt
  if (which_study == "sanity") {
    data_pdt_aff$cat = factor(as.character(data_pdt_aff$cat),levels=c("6","2","3"),labels=c("neu_AW","neg","pos"))
  } else if (which_study == "POSTPILOT_HC" | which_study == "POSTPILOT_PG") {
    data_pdt_aff$cat = factor(as.character(data_pdt_aff$cat),levels=c("6","1","2","3"),labels=c("neu_AW","gam","neg","pos"))
  }
  
  data_pdt_aff$affLength_log = log(data_pdt_aff$affLength)
  fit7 = lmer(affLength ~ cat + (cat|subject),data=data_pdt_aff)
  summary(fit7)
  agk.normality.check(fit7)
  
  data_pdt$affLength_phys=sqrt(data_pdt$arousal^2+data_pdt$dominance^2+data_pdt$valence^2)
  data_pdt_aff = data_pdt
  if (which_study == "sanity") {
    data_pdt_aff$cat = factor(as.character(data_pdt_aff$cat),levels=c("6","2","3"),labels=c("neu_AW","neg","pos"))
  } else if (which_study == "POSTPILOT_HC" | which_study == "POSTPILOT_PG") {
    data_pdt_aff$cat = factor(as.character(data_pdt_aff$cat),levels=c("6","1","2","3"),labels=c("neu_AW","gam","neg","pos"))
  }
  
  data_pdt_aff$affLength_log = log(data_pdt_aff$affLength)
  fit7 = lmer(affLength ~ cat + (cat|subject),data=data_pdt_aff)
  summary(fit7)
  agk.normality.check(fit7)
}

## IMAGE ADEQUACY PHYSIO ======================================================
# CONTINUE HERE # 
if (do_physio) {
  # TODO: check with lmlist if there are shady subs!
  # TODO: implement with lmlist also and give option
  
  # which sum stat?
  if (sum_stat == 'mean') {
    data_pdt$corr = data_pdt$corr_auc
    data_pdt$zygo = data_pdt$zygo_auc
    data_pdt$eda  = data_pdt$eda_auc
  } else if (sum_stat == 'max') {
    data_pdt$corr = data_pdt$corr_max
    data_pdt$zygo = data_pdt$zygo_max
    data_pdt$eda  = data_pdt$eda_max
  } else if (sum_stat == 'med') {
    data_pdt$corr = data_pdt$corr_median
    data_pdt$zygo = data_pdt$zygo_median
    data_pdt$eda  = data_pdt$eda_median
  } else {
    stop('No adequat sum_stat for physio selected.')
  }
  
  # contrasts
  cur_con                 = contr.treatment(length(levels(data_pdt$cat)),base = which(levels(data_pdt$cat) == 'neutral'))
  rownames(cur_con)       = levels(data_pdt$cat)
  colnames(cur_con)       = levels(data_pdt$cat)[FALSE == (levels(data_pdt$cat) == 'neutral')]
  contrasts(data_pdt$cat) = cur_con
  
  # select group again
  # HC and PG selection
  data_pdt_HC = subset(data_pdt, HCPG == 'HC')
  data_pdt_PG = subset(data_pdt, HCPG == 'PG')
  
  ## PHYSIO ADEQUACY
  # Effect of category on corrugator
  # lmlist test
  fit8clml          = lmList(corr ~ 1 | subject,data = data_pdt)
  fit9clml          = lmList(corr ~ cat | subject,data = data_pdt)
  
  # lmer
  fit8c             = lmer(corr ~ 1 + (1| subject) + (1|stim),data = data_pdt)
  fit9c             = lmer(corr ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt)
  if (!which_study == 'sanity') {
    fit9c_hc          = lmer(corr ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_HC) 
    fit9c_pg          = lmer(corr ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_PG)
    fit9cg            = lmer(corr ~ (cat*HCPG) + (cat| subject),data = data_pdt,REML=F)
    print(anova(fit9c,fit9cg))
  }
  
  if (do_bootstrap == 1) {
    cur_fun = function(d,indices) {
      cur_d = d[indices,]
      fit9              = lmer(corr ~ (cat) + (cat| subject),data=cur_d)
      fixef(fit9)
    }
    
    cl<-makeCluster(5)
    clusterEvalQ(cl, library(lme4))
    tmp = boot(data_pdt,cur_fun,R=1000,parallel = "snow",cl = cl)
    stopCluster(cl)
  }
  
  print(anova(fit8c,fit9c))
  print(anova(fit9c,fit9cg))
  summary(fit9c)
  agk.normality.check(fit9c)
  
  # Effect of category on zygomaticus
  fit8z              = lmer(zygo ~ 1 + (1| subject) + (1| stim),data = data_pdt)
  fit9z              = lmer(zygo ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt)
  if (!which_study == 'sanity') {
    fit9z_hc         = lmer(zygo ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_HC)
    fit9z_pg         = lmer(zygo ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_PG)
    fit9zg           = lmer(zygo ~ (cat*HCPG) + (cat| subject)+ (cat| stim),data = data_pdt)
    print(anova(fit9z,fit9zg))
  }
  print(anova(fit8z,fit9z))
  summary(fit9z)
  agk.normality.check(fit9z)

  # Effect of category on eda
  fit8e              = lmer(eda ~ 1 + (1| subject) +  (1| stim),data = data_pdt)
  fit9e              = lmer(eda ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt)
  if (!which_study == 'sanity') {
    fit9e_hc         = lmer(eda ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_HC)
    fit9e_pg         = lmer(eda ~ (cat) + (cat| subject) + (cat| stim),data = data_pdt_PG)
    fit9eg           = lmer(eda ~ (cat*HCPG) + (cat| subject) + (cat| stim),data = data_pdt)
    print(anova(fit9e,fit9eg))
    summary(fit9eg)
  }
  print(anova(fit8e,fit9e))
  summary(fit9e)
  agk.normality.check(fit9e)
  coefs <- data.frame(coef(summary(fit9eg)))
  # use normal distribution to approximate p-value
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
}

if (do_leda == 1) {
  # Effect of category on LEDA indices
  fit10a              = lmer(AmpSum_stim ~ 1 + (1| subject),data = data_pdt)
  fit11a              = lmer(AmpSum_stim ~ (cat) + (cat| subject),data = data_pdt)
  anova(fit10a,fit11a)
  summary(fit11a)
  agk.normality.check(fit11a)
  #ggplot(data=data_pdt, aes(x=cat, y=AmpSum_stim, color=cat))+ geom_boxplot() + facet_wrap(~subject, ncol=7) 
  
  
  fit10b              = lmer(AmpSum_gam ~ 1 + (1| subject),data = data_pdt)
  fit11b              = lmer(AmpSum_gam ~ (cat) + (cat| subject),data = data_pdt)
  anova(fit10b,fit11b)
  summary(fit11b)
  agk.normality.check(fit11b)
  #ggplot(data=data_pdt, aes(x=cat, y=AmpSum_gam, color=cat)) + geom_boxplot() + facet_wrap(~subject, ncol=7) 
  
  
  fit12a              = lmer(SCR_stim ~ 1 + (1| subject),data = data_pdt)
  fit13a              = lmer(SCR_stim ~ (cat) + (cat| subject),data = data_pdt)
  anova(fit12a,fit13a)
  summary(fit13a)
  agk.normality.check(fit13a)
  #ggplot(data=data_pdt, aes(x=cat, y=SCR_stim, color=cat)) + geom_boxplot() + facet_wrap(~subject, ncol=7)
  
  
  fit12b              = lmer(SCR_gam ~ 1 + (1| subject),data = data_pdt)
  fit13b              = lmer(SCR_gam ~ (cat) + (cat| subject),data = data_pdt)
  anova(fit12b,fit13b)
  summary(fit13b)
  coeficients(fit13b)
  
  agk.normality.check(fit13b)
}



## CAREFUL GET_LOG function SEEMS NOT TO WORK PROPERLY ##
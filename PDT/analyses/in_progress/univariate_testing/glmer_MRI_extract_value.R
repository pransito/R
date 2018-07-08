## running glmer tests ========================================================
library(nloptr)
nlopt <- function(par, fn, lower, upper, control) {
  .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                            opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                        maxeval = 300, xtol_abs = 1e-6, ftol_abs = 1e-6))
  list(par = res$solution,
       fval = res$objective,
       conv = if (res$status > 0) 0 else res$status,
       message = res$message
  )
}

cur_control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)

cur_control_glmer=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                         check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=300))

# select and order data_pdt
data_pdt_bcp_glmer = data_pdt
des_vars_exp       = c("subject","trial","stim","gain","loss","rt","cat","choice","side","st_dur","age","sex","gain_bcp","loss_bcp","gainxloss","gainxloss_bcp","ed_abs","ed",                     
                       "EV","ratio","diff","RiskMin","RiskAG","RiskMar","SkewMin","EV_bcp","accept_reject","HCPG")
data_pdt           = data_pdt[des_vars_exp]
data_pdt           = agk.mult.order(data_pdt)

# previous stimulus (category)
data_pdt$cat_prev    = NA
data_pdt$cat_prev[1] = 'neutral'
for (cc in 2:length(data_pdt$cat)) {
  data_pdt$cat_prev[cc] = as.character(data_pdt$cat[cc-1]) 
}
data_pdt$cat_prev = factor(data_pdt$cat_prev,levels = levels(data_pdt$cat))

# killing missings
data_pdt$choice[data_pdt$choice == 5] = NA

# with cat and group
lac_l   = lmer(choice ~ (gain + loss + cat)*HCPG + (gain + loss + cat|subject) + (gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control,REML = F)
laec_l  = lmer(choice ~ (gain + loss + ed_abs + cat)*HCPG + (gain + loss + ed_abs + cat|subject) + (gain + loss + ed_abs|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
laecp_l = lmer(choice ~ (gain + loss + ed_abs + cat + cat_prev)*HCPG + (gain + loss + ed_abs + cat + cat_prev|subject) + (gain + loss + ed_abs + cat_prev|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
laecr_l = lmer(choice ~ (gain + loss + rt + cat)*HCPG + (gain + loss + rt + cat|subject) + (gain + loss + rt|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
laCh_l  = lmer(choice ~ (0 + gain + loss)*HCPG + (0 + gain + loss|subject) + (0 + gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
la_l    = lmer(choice ~ (gain + loss)*HCPG + (gain + loss|subject) + (gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
lat_l   = lmer(choice ~ (gain + loss+trial)*HCPG + (gain + loss+trial|subject) + (gain + loss+trial|stim),nAGQ=0,data = data_pdt,control=cur_control,REML=F)
anova(lac_l,laec_l,laCh_l,la_l)
anova(laec_l,laec)

# with cat and group
lac   = glmer(accept_reject ~ (gain + loss + cat)*HCPG + (gain + loss + cat|subject) + (gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
laec  = glmer(accept_reject ~ (gain + loss + ed_abs + cat)*HCPG + (gain + loss + ed_abs + cat|subject) + (gain + loss + ed_abs|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
laecp = glmer(accept_reject ~ (gain + loss + ed_abs + cat + cat_prev)*HCPG + (gain + loss + ed_abs + cat + cat_prev|subject) + (gain + loss + ed_abs + cat_prev|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
laecr = glmer(accept_reject ~ (gain + loss + rt + cat)*HCPG + (gain + loss + rt + cat|subject) + (gain + loss + rt|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
laCh  = glmer(accept_reject ~ (0 + gain + loss)*HCPG + (0 + gain + loss|subject) + (0 + gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
la    = glmer(accept_reject ~ (gain + loss)*HCPG + (gain + loss|subject) + (gain + loss|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
lat   = glmer(accept_reject ~ (gain + loss+trial)*HCPG + (gain + loss+trial|subject) + (gain + loss+trial|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
anova(lac,laec,laCh,la)
anova(la,lat)

# laec has best BIC
# get the value per trial and subject; we leave stimulus random effects, but extract only over subjects
laec_nc   = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + (gain + loss + ed_abs|subject) + (gain + loss + ed_abs|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
anova(laec,laec_nc)
anova(laec,laecp)
anova(laec,laecr)
anova(lac,laec)


# get the complete coef
compl_coef = agk.get.compl.coef(laec_nc,'HCPG')

# compute value per trial and subject
des_vars     = c('accept_reject','gain','loss','ed_abs')
des_vars_exp = c('accept_reject','gain','loss','ed_abs','subject','trial') 
preds        = c('(Intercept)','gain','loss','ed_abs')
all_subs     = unique(data_pdt$subject)
fits         = c()
using_glmer  = F
val_vecs     = c()
group_vec    = c()
mfit_vec     = c()

for (vv in 1:length(all_subs)) {
  # get data
  cur_dat = subset(data_pdt, subject == all_subs[vv])
  cur_dat = cur_dat[des_vars]
  print(length(cur_dat[,1]))
  
  # get the group
  group_vec[vv] = agk.recode.c(all_subs[vv],data_pdt$subject,data_pdt$HCPG)
  
  # get the none-glmer model
  cur_nm  = glm(accept_reject ~  gain + loss, cur_dat,family='binomial')
  
  # get a model matrix where NAs are not omitted
  cur_dat_nona = cur_dat
  cur_dat_nona$accept_reject[is.na(cur_dat_nona$accept_reject)] = 1
  cur_mm       = model.matrix(glm(accept_reject ~  gain + loss + ed_abs, cur_dat_nona,family='binomial'))
  
  if (using_glmer) {
    # get model
    cur_mo  = compl_coef[compl_coef$subject == all_subs[vv],]
    cur_mo  = cur_mo[preds]
    
    # compute value
    cur_va  = as.matrix(cur_mo) %*% t(cur_mm)
  } else {
    cur_va = predict(cur_nm,newdata = cur_dat_nona)
    
  }
  
  # record fit
  predic   = as.numeric(cur_va > 0)
  fits[vv] = (mean.rmna(predic  == cur_dat$accept_reject))
  
  # print mean and sd
  print(mean(cur_va))
  print(sd(cur_va))
  
  # record values
  val_vecs = c(val_vecs, cur_va)
}

print(mean(fits))

# differences in model fit between groups
summary(lm(fits~group_vec))


# export
data_pdt$sub_value = val_vecs
setwd('C:/Users/genaucka/Google Drive/Promotion/VPPG/VPPG_Exchange/Experimente/PDT/Daten/pilot')
write.table(x = data_pdt, file = 'data_pdt.csv',sep = ';',row.names = F)
setwd('C:/Users/genaucka/Google Drive/Library/MATLAB/PDT/MRI/ss_models/classical')
write.table(x = data_pdt, file = 'data_pdt.csv',sep = ';',row.names = F)

## SECOND ROUND FIT: VALUE*CATEGORY
# you can also kill certain variables in formular and then use offset to pass an offset variable (one value per observation)
# https://stats.stackexchange.com/questions/82411/how-to-fix-one-coefficient-and-fit-others-using-regression 
lav    = glmer(accept_reject ~ sub_value     + (sub_value|subject) + (sub_value|stim) ,nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
lavci  = glmer(accept_reject ~ sub_value*cat + (sub_value*cat|subject) + (sub_value*cat|stim) ,nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
lavcig = glmer(accept_reject ~ sub_value*cat*HCPG + (sub_value*cat|subject) +(sub_value*cat|stim),nAGQ=0,data = data_pdt,control=cur_control_glmer, family='binomial')
anova(lav,lavci,lavcig)
summary(lavcig)

data_pdt_PG = subset(data_pdt,HCPG == 'PG')
lav_PG    = glmer(accept_reject ~ sub_value     + (sub_value|subject) + (sub_value|stim) ,nAGQ=0,data = data_pdt_PG,control=cur_control_glmer, family='binomial')
lavci_PG  = glmer(accept_reject ~ sub_value*cat + (sub_value*cat|subject) + (sub_value*cat|stim) ,nAGQ=0,data = data_pdt_PG,control=cur_control_glmer, family='binomial')
anova(lav_PG,lavci_PG)
summary(lavci_PG)

data_pdt_HC = subset(data_pdt,HCPG == 'HC')
lav_HC   = glmer(accept_reject ~ sub_value     + (sub_value|subject) + (sub_value|stim) ,nAGQ=0,data = data_pdt_HC,control=cur_control_glmer, family='binomial')
lavci_HC = glmer(accept_reject ~ sub_value*cat + (sub_value*cat|subject) + (sub_value*cat|stim) ,nAGQ=0,data = data_pdt_HC,control=cur_control_glmer, family='binomial')
anova(lav_HC,lavci_HC)
summary(lavci_HC)


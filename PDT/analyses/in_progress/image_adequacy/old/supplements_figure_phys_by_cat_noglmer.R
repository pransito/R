## Plots means and CIs of physio variables (stim)==============================

# ALIGNED
# means and CIs
des_vars = c('valence','arousal','dominance',"imageRating1s","imageRating2s","imageRating3s","imageRating4s")
cur_map = c('Valence', 'Arousal','Dominance', "Elicits_Craving",
            "Representative_for_Gambles","Representative_for_Negative","Representative_for_Positive")

cfint_list = list()
for (ii in 1:length(des_vars)) {
  
  cur_form = as.formula(paste(des_vars[ii],'~ subject + cat'))
  cfint = aggregate(cur_form,data = data_pdt,FUN=mean)
  # unit test no missings
  if (any(is.na(cfint[[3]]))) {
    stop('Missings!')
  }
  cfint$HCPG = agk.recode.c(cfint$subject,dat_match$VPPG,dat_match$HCPG)
  cur_var = names(cfint)[3]
  cur_form = as.formula(paste(cur_var,'~ HCPG + cat'))
  cfint$HCPG = agk.recode(cfint$HCPG,c('PG'),c('GD'))
  cfint$HCPG = factor(cfint$HCPG, levels = c('HC','GD'))
  cfint = aggregate(cur_form,data = cfint,FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  
  cfint = as.data.frame(cbind(cfint[c('HCPG','cat')],cfint[[3]]))
  names(cfint)[c(3:5)] = c('mean','lower','upper')
  cfint$variable = agk.recode(des_vars[ii],des_vars,cur_map)
  cfint_list[[ii]] = cfint
}
cfint = cfint_list[[1]]
for(ii in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[ii]])
}

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Ratings of different categories of stimuli")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(cat ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

## get info which gambles they do by SOGS
des_vars = c('SOGS_1[a-z]')
des_vars = names(dat_match)[grep(des_vars,names(dat_match))]
cur_map = c('cards', 'bet_animals','bet_sports', "dice","casinos",
            "lotteries","bingo",'stock_market','slot_machines','played_sports')

cfint_list = list()
for (ii in 1:length(des_vars)) {
  

  
  x1      = dat_match[dat_match$HCPG == 'HC',des_vars[ii]]
  x2      = dat_match[dat_match$HCPG == 'PG',des_vars[ii]]
  
  # unit test no missings
  if (any(is.na(x1)) | any(is.na(x2))) {
    stop('Missings!')
  }
  
  cfint   = t.test(x1,x2,'less')
  cur_res = list()
  cur_res$meanHC = mean(x1)
  cur_res$meanPG = mean(x2)
  cur_res$p      = round(cfint$p.value,digits = 3)
  
  cfint_list[[cur_map[ii]]] = cur_res
}

cfint = as.data.frame(cfint_list[[1]])
cfint$gamble = names(cfint_list)[1]
for(ii in 2:length(cfint_list)) {
  cur_res        = as.data.frame(cfint_list[[ii]])
  cur_res$gamble = names(cfint_list)[ii]
  cfint = rbind(cfint,cur_res)
}
cfint = cfint[c(4,1,2,3)]

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Ratings of different categories of stimuli")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(cat ~ .)
mRatg = mRat + theme_bw()
print(mRatg)



## Plots means and CIs of physio variables (stim)==============================

# ALIGNED
# means and CIs
des_vars = c('zygo_auc', 'zygo_max','zygo_median','zygo_pmv',
             'corr_auc','corr_max','corr_median','corr_pmv')
cfint_list = list()
for (ii in 1:length(des_vars)) {
  
  cur_form = as.formula(paste(des_vars[ii],'~ subject + cat'))
  cfint = aggregate(cur_form,data = data_pdt,FUN=mean)
  cfint$HCPG = agk.recode.c(cfint$subject,dat_match$VPPG,dat_match$HCPG)
  cur_form = as.formula(paste(des_vars[ii],'~ HCPG + cat'))
  cfint = aggregate(cur_form,data = cfint,FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  
  cfint = as.data.frame(cbind(cfint[c('HCPG','cat')],cfint[[3]]))
  names(cfint)[c(3:5)] = c('mean','lower','upper')
  cfint$variable = des_vars[ii]
  cfint_list[[ii]] = cfint
}
cfint = cfint_list[[1]]
for(ii in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[ii]])
}

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Physio reactions to different categories of stimuli")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(cat ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

## Plots means and CIs of physio variables (GAMBLE Phase) =====================
cr_pp   = c("corr","zygo")
if (physio_sum_fun == 'all') {
  # use all the available summary stats
  all_phys_names = c()
  for (ii in 1:length(cr_pp)) {
    all_phys_names = c(all_phys_names,names(data_pdt)[grep(cr_pp[ii],names(data_pdt))])
  }
  cr_pp = all_phys_names
}
# here take out the phys names for gamble phase
all_phys_names_gamble = all_phys_names[grep('_gamble',all_phys_names)]
cr_pp                 = all_phys_names[grep('_gamble',all_phys_names,invert = T)]

cur_cr_gam            = data_pdt[all_phys_names_gamble]
cur_cr_la             = data_pdt[c('cat','gain','loss','accept_reject','subject')]
cur_cr_stim           = data_pdt[cr_pp]
cur_cr_gamStim        = cur_cr_gam - cur_cr_stim
names(cur_cr_gamStim) = paste0(names(cur_cr_gamStim),'Stim') 
all_vars_of_int       = names(cur_cr_gamStim)
cur_cr_gamStim        = cbind(cur_cr_gamStim,cur_cr_la)
all_coefs_df          = list()

imp_missing = function(x) {x[is.na(x)] = mean.rmna(x); return(x)}
for (gg in 1:length(all_vars_of_int)) {
  cur_var            = all_vars_of_int[gg]
  cur_form           = as.formula(paste(cur_var,'~','accept_reject*cat+gain+loss|subject'))
  cur_coefs          = coef(lmList(cur_form, data = cur_cr_gamStim,pool = F))
  cur_subs           = row.names(cur_coefs)
  cur_coefs          = as.data.frame(lapply(cur_coefs,FUN=imp_missing))
  names(cur_coefs)   = paste(cur_var,names(cur_coefs),sep='X')
  cur_coefs$subject  = cur_subs
  all_coefs_df[[gg]] = cur_coefs
}

# merging
physGam_mf = all_coefs_df[[1]]
for (gg in 2:length(all_coefs_df)) {
  physGam_mf = merge(physGam_mf,all_coefs_df[[gg]],by.x = c('subject'),by.y = ('subject'))
}

# all the vars
physGam_vars = names(physGam_mf)
physGam_vars = physGam_vars[grep('subject',physGam_vars,invert = T)]

physGam_mf$HCPG = agk.recode.c(physGam_mf$subject,dat_match$VPPG,dat_match$HCPG)

cfint_list = list()
for (gg in 1:length(physGam_vars)) {
  cur_form = as.formula(paste(physGam_vars[gg],'~ HCPG'))
  cfint = aggregate(cur_form,data = physGam_mf,FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  
  cfint = as.data.frame(cbind(cfint[c('HCPG')],cfint[[2]]))
  names(cfint)[c(2:4)] = c('mean','lower','upper')
  cfint$variable = physGam_vars[gg]
  cfint_list[[gg]] = cfint
}
cfint = cfint_list[[1]]
for(ii in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[ii]])
}

# cutting out some variables
cfint = cfint[grep('reject1.cat',cfint$variable),]

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Physio reaction shifts in gamble phase: interactions")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
#mRat  = mRat + facet_grid(cat ~ .)
mRatg = mRat + theme_bw() + theme(legend.position="bottom")
print(mRatg)


## plotting the pspml extracts (pic phase) ====================================
setwd(path_dat)
pspml               = read.table('pspml_FIR_out.txt',header=T,sep=',')
f                   = function(x) {return(all(is.nan(x)))}
pspml               = pspml[,as.logical(unlist(lapply(pspml, FUN=f))) == FALSE]
pspml$Constant1     = NULL
# vars_opt            = names(pspml)[grep('Gam_opt',names(pspml))]
# vars_opt_cat        = vars_opt[grep('acc',vars_opt)]
# vars_opt_cat_n      = vars_opt[grep('neu',vars_opt_cat)]
# vars_opt_cat        = vars_opt[grep('neu',vars_opt,invert = T)]
# pspml[vars_opt_cat] = pspml[vars_opt_cat] - rep(pspml[vars_opt_cat_n],length(vars_opt_cat))
pspml_pic           = pspml[grep('Gam_Opt',names(pspml),invert = T,ignore.case = T)]
# pruning
pspml_pic$acceptance_recon = NULL

pspml_pic           = melt(pspml_pic,id.vars = 'sub')

# get category variable
fCat = function(x) {
  tmp = strsplit(as.character(x),split = '_')
  return (tmp[[1]][2])
}
pspml_pic$category = sapply(pspml_pic$variable,fCat)
pspml_pic$category = factor(pspml_pic$category,levels=c('neu','gam','neg','pos'))

# get basis function variable
fBf = function(x) {
  tmp = strsplit(as.character(x),split = '_')
  return (tmp[[1]][3])
}
pspml_pic$basisFunction = sapply(pspml_pic$variable,fBf)

# get VPPG
pspml_pic$VPPG = agk.recode.c(pspml_pic$sub,tnl$PhysioVP,tnl$VPPG)
pspml_pic$HCPG = agk.recode.c(pspml_pic$VPPG,dat_match$VPPG,dat_match$HCPG)

# align with dat_match
pspml_pic = pspml_pic[pspml_pic$VPPG %in% dat_match$VPPG,]

# get cfint
coefs           = coef(lmList(value ~ category | sub,data=pspml_pic, pool = F))
names(coefs)[1] = 'categoryneu' 
physPic_vars    = names(coefs)
coefs$subject   = row.names(coefs)
coefs$VPPG      = agk.recode.c(coefs$subject,tnl$PhysioVP,tnl$VPPG)
coefs$HCPG      = agk.recode.c(coefs$VPPG,dat_match$VPPG,dat_match$HCPG)
cfint_list      = list()
for (gg in 1:length(physPic_vars)) {
  cur_form = as.formula(paste(physPic_vars[gg],'~ HCPG'))
  cfint = aggregate(cur_form,data = coefs,FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  
  cfint = as.data.frame(cbind(cfint[c('HCPG')],cfint[[2]]))
  names(cfint)[c(2:4)] = c('mean','lower','upper')
  cfint$variable = physPic_vars[gg]
  cfint_list[[gg]] = cfint
}
cfint = cfint_list[[1]]
for(ii in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[ii]])
}
names(cfint)[5] = 'category'

# plotting (faceting by basis function)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=category))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("EDA pspml extracts for different categories of stimuli (PICphase)")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
#mRat  = mRat + facet_grid(basisFunction ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

# stats
pspml_pic$category = factor(pspml_pic$category,levels=c('neu','gam','pos','neg'))
mod_cat            = lm(value ~ category, data=pspml_pic)
mod_cat_g          = lm(value ~ category*HCPG, data=pspml_pic)
summary(mod_cat)
summary(mod_cat_g)

## plotting the pspml extracts (GAM.OPT phase) ================================
# ALIGNED
setwd(path_dat)
pspml               = read.table('pspml_FIR_out.txt',header=T,sep=',')
f                   = function(x) {return(all(is.nan(x)))}
pspml               = pspml[,as.logical(unlist(lapply(pspml, FUN=f))) == FALSE]
pspml$Constant1     = NULL
needed_vars         = c(grep('gam_opt',names(pspml),ignore.case = T),
                        grep('acceptance_',names(pspml),ignore.case = T))
pspml_gam           = pspml[,needed_vars]
pspml_gam$sub       = pspml$sub
pspml_gam           = melt(pspml_gam,id.vars = 'sub')

# get VPPG
pspml_gam$VPPG = agk.recode.c(pspml_gam$sub,tnl$PhysioVP,tnl$VPPG)
pspml_gam$HCPG = agk.recode.c(pspml_gam$VPPG,dat_match$VPPG,dat_match$HCPG)

# align with dat_match
pspml_gam = pspml_gam[pspml_gam$VPPG %in% dat_match$VPPG,]

# get cfint
cfint = aggregate(value ~ HCPG + variable,data = pspml_gam,FUN=agk.boot.ci,
                  cur_fun = mean,lower=0.025, upper=0.975,R=5000)
cfint = as.data.frame(cbind(cfint[c('HCPG','variable')],cfint[[3]]))
names(cfint)[c(3:5)] = c('mean','lower','upper')

# plotting (faceting by basis function)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("EDA pspml extracts in GAM.OPT phase")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRatg = mRat + theme_bw()
print(mRatg)

# stats
pspml_pic_bf2          = subset(pspml_pic,basisFunction == 'Bf2')
pspml_pic_bf2$opt_on   = factor(pspml_pic_bf2$opt_on,levels=c('on','onXEd','onXGain','onXLoss'))
mod_cat                = lm(value ~ opt_on, data=pspml_pic_bf2)
mod_cat_g              = lm(value ~ opt_on*HCPG, data=pspml_pic_bf2)
summary(mod_cat)
summary(mod_cat_g)

## plotting the leda extracts (pic phase) ====================================
# ALIGNED
setwd(path_dat)
leda_dat       = read.table('ledalab_out.csv',sep=',',header=T)
kick_out       = which(names(leda_dat) %in% c('Latency','AmpSum','ISCR') == FALSE)
kick_out       = which(names(leda_dat) %in% c('Latency','AmpSum','ISCR','nSCR','PhasicMax','Tonic') == FALSE)
leda_dat       = leda_dat[,kick_out]
leda_dat       = subset(leda_dat,onset != 'gamble')
leda_dat$trial = NULL
leda_dat       = aggregate(.~subject+onset,data = leda_dat,FUN=mean)
leda_dat       = melt(leda_dat,id.vars = c('subject','onset'))

# subject var
leda_dat$VPPG    = agk.recode.c(leda_dat$subject,tnl$PhysioVP,tnl$VPPG)
leda_dat$subject = NULL

# HCPG var
leda_dat$HCPG = agk.recode.c(leda_dat$VPPG,dat_match$VPPG,dat_match$HCPG)

# align with dat_match
leda_dat = leda_dat[leda_dat$VPPG %in% dat_match$VPPG,]

# get cfint
cfint = aggregate(value ~ HCPG + variable + onset,data = leda_dat,FUN=agk.boot.ci,
                  cur_fun = mean,lower=0.025, upper=0.975,R=5000)
cfint = as.data.frame(cbind(cfint[c('HCPG','variable','onset')],cfint[[4]]))
names(cfint)[c(4:6)] = c('mean','lower','upper')

# plotting (faceting by basis function)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=onset))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("EDA LEDA extracts to different categories of stimuli")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(variable ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

## plotting the leda extracts (GAM phase) ====================================
# ALIGNED
setwd(path_dat)
leda_dat         = read.table('ledalab_out.csv',sep=',',header=T)
kick_out         = which(names(leda_dat) %in% c('Latency','AmpSum','ISCR') == FALSE)
leda_dat         = leda_dat[,kick_out]
leda_dat_pic     = subset(leda_dat,onset != 'gamble')
leda_dat         = subset(leda_dat,onset == 'gamble')

# create a difference score
cur_des_vars           = !names(leda_dat) %in% c('subject','trial','onset')
cur_des_vars           = names(leda_dat)[cur_des_vars]
leda_dat[cur_des_vars] = leda_dat[cur_des_vars] - leda_dat_pic[cur_des_vars]

# names of vars
cur_code                  = names(leda_dat) %in% cur_des_vars
names(leda_dat)[cur_code] = paste0('GamOptShift',names(leda_dat)[cur_code])
cur_des_vars              = names(leda_dat)[cur_code]
leda_dat$VPPG             = agk.recode.c(leda_dat$subject,tnl$PhysioVP,tnl$VPPG)
leda_dat$subject          = NULL
dtmp                      = merge(data_pdt,leda_dat,by.x = c('subject','trial'),by.y = c('VPPG','trial'))

# modeling phase GAM.OPT
imp_missing = function(x) {x[is.na(x)] = mean.rmna(x); return(x)}

coefs_dfs = list()
for (gg in 1:length(cur_des_vars)) {
  cur_form  = as.formula(paste(cur_des_vars[gg],'~','accept_reject*cat+gain+loss|subject'))
  cur_coefs = coef(lmList(cur_form, data = dtmp,pool = F))
  cur_subs  = row.names(cur_coefs)
  cur_coefs = as.data.frame(lapply(cur_coefs,FUN=imp_missing))
  
  # adjust names
  names(cur_coefs) = paste(cur_des_vars[gg],names(cur_coefs),sep = '_')
  
  # add subject variable
  cur_coefs$subject = cur_subs
  
  coefs_dfs[[gg]] = cur_coefs
}

coefs_df = coefs_dfs[[1]]
for (gg in 2:length(cur_des_vars)) {
  coefs_df = merge(coefs_df,coefs_dfs[[gg]],by='subject')
}

coefs_df = melt(coefs_df,id='subject')
coefs_df$HCPG = agk.recode.c(coefs_df$subject,dat_match$VPPG,dat_match$HCPG)

# get cfint
cfint = aggregate(value ~ HCPG + variable,data = coefs_df,FUN=agk.boot.ci,
                  cur_fun = mean,lower=0.025, upper=0.975,R=5000)
cfint = as.data.frame(cbind(cfint[c('HCPG','variable')],cfint[[3]]))
names(cfint)[c(3:5)] = c('mean','lower','upper')

# preselect some variables
cfint = cfint[grep('accept_reject1.cat',cfint$variable),]
cfint = cfint[grep('SCR',cfint$variable),]
cfint = cfint[grep('nSCR',cfint$variable,invert = T),]

# plotting
mRat  = ggplot(cfint, aes(HCPG, mean,fill=variable))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("EDA LEDA extracts gambling phase")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
#mRat  = mRat + facet_grid(variable ~ .)
mRatg = mRat + theme_bw() #+ theme(legend.position="bottom")
print(mRatg)


## plotting the pspmnl extracts (pic phase) ===================================
# ALIGNED
# we got 6 paramters per trial, 3 for event 1 (pic) and 3 for event 2 (gamble)
des_vars = c('FlexibleResponse_1_Amplitude','FlexibleResponse_1_PeakLatency',
             'FlexibleResponse_1_Dispersion','FlexibleResponse_2_Amplitude',
             'FlexibleResponse_2_PeakLatency','FlexibleResponse_2_Dispersion')
setwd(path_dat)
pspmnl_dat = read.table('pspmnl_out_tab.csv',sep=',',header=T)
all_subs   = pspmnl_dat$sub
f = function(x) {return(all(is.na(x)))}
all_dfs = list()
for (ii in 1:length(all_subs)) {
  cur_dat     = subset(pspmnl_dat,sub == all_subs[ii])
  cur_sub     = all_subs[ii]
  cur_dat$sub = NULL
  cur_dat[which(as.logical(unlist(lapply(cur_dat,f))))] = NULL
  cur_df = list()
  for (jj in 1:length(des_vars)) {
    sub_dat = cur_dat[grep(des_vars[jj],names(cur_dat))]
    cur_df[[jj]] = as.numeric(sub_dat)
  }
  names(cur_df) = des_vars
  cur_df        = as.data.frame(cur_df)
  cur_df$trial  = 1:length(cur_df[,1])
  cur_df$sub    = as.character(cur_sub)
  all_dfs[[ii]] = cur_df
}

combined_df = all_dfs[[1]]
for (ii in 1:length(all_dfs)) {
  combined_df = rbind(combined_df,all_dfs[[ii]])
}

# get subject var
combined_df$VPPG = agk.recode.c(combined_df$sub,tnl$PhysioVP,tnl$VPPG)
combined_df$sub  = NULL

# merging to data_pdt
data_pdt = merge(data_pdt,combined_df,by.x = c('subject','trial'),by.y = c('VPPG','trial'))

# means and CIs of PIC phase
des_vars = names(combined_df)[c(1:3)]
cfint_list = list()
for (ii in 1:length(des_vars)) {
  cur_form = as.formula(paste(des_vars[ii],'~ subject + cat'))
  cfint = aggregate(cur_form,data = data_pdt,FUN=mean)
  cfint$HCPG = agk.recode.c(cfint$subject,dat_match$VPPG,dat_match$HCPG)
  cur_form = as.formula(paste(des_vars[ii],'~ HCPG + cat'))
  cfint = aggregate(cur_form,data = cfint,FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  cfint = as.data.frame(cbind(cfint[c('HCPG','cat')],cfint[[3]]))
  names(cfint)[c(3:5)] = c('mean','lower','upper')
  cfint$variable = des_vars[ii]
  cfint_list[[ii]] = cfint
}
cfint = cfint_list[[1]]
for(ii in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[ii]])
}

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=cat))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Physio reactions (PSPMNL) to different categories of stimuli")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(variable ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

## plotting the pspmnl extracts (pic phase) ===================================
imp_missing = function(x) {x[is.na(x)] = mean.rmna(x); return(x)}
data_pdt$FlexibleResponse_12_Amplitude   = data_pdt$FlexibleResponse_2_Amplitude - data_pdt$FlexibleResponse_1_Amplitude
data_pdt$FlexibleResponse_12_PeakLatency = data_pdt$FlexibleResponse_2_PeakLatency - data_pdt$FlexibleResponse_1_PeakLatency 
data_pdt$FlexibleResponse_12_Dispersion  = data_pdt$FlexibleResponse_2_Dispersion - data_pdt$FlexibleResponse_1_Dispersion 

pspmnl_12_list = list()
physioPDT_PSPMNL_12_Ampl = as.data.frame(lapply(coef(lmList(FlexibleResponse_12_Amplitude ~ accept_reject*cat+gain+loss|subject, data = data_pdt,pool = F)),FUN=imp_missing))
physioPDT_PSPMNL_12_Lat  = as.data.frame(lapply(coef(lmList(FlexibleResponse_12_PeakLatency ~accept_reject*cat+gain+loss|subject, data = data_pdt,pool = F)),FUN=imp_missing))
physioPDT_PSPMNL_12_Disp = as.data.frame(lapply(coef(lmList(FlexibleResponse_12_Dispersion ~ accept_reject*cat+gain+loss|subject, data = data_pdt,pool = F)),FUN=imp_missing))
pspmnl_12_list$pspmnlAmpl = physioPDT_PSPMNL_12_Ampl
pspmnl_12_list$pspmnlLat  = physioPDT_PSPMNL_12_Lat
pspmnl_12_list$pspmnlDisp = physioPDT_PSPMNL_12_Disp

# means and CIs of STIM phase (shifts in physio response)
cfint_list = list()
des_vars   = names(pspmnl_12_list)
for (kk in 1:length(pspmnl_12_list)) {
  pspmnl_12_list[[kk]]$subject = unique(data_pdt$subject)
  pspmnl_12_list[[kk]]         = melt(pspmnl_12_list[[kk]])
  pspmnl_12_list[[kk]]$HCPG    = agk.recode.c(pspmnl_12_list[[kk]]$subject,dat_match$VPPG,dat_match$HCPG)
  cur_form = as.formula(paste('value','~ HCPG + variable'))
  cfint = aggregate(cur_form,data = pspmnl_12_list[[kk]],FUN=agk.boot.ci,
                    cur_fun = mean,lower=0.025, upper=0.975,R=5000)
  cfint = as.data.frame(cbind(cfint[c('HCPG','variable')],cfint[[3]]))
  names(cfint)[c(3:5)] = c('mean','lower','upper')
  names(cfint)[2] = 'shift_gamble_phase'
  cfint$variable = des_vars[kk]
  cfint_list[[kk]] = cfint
}
cfint = cfint_list[[1]]
for(kk in 2:length(cfint_list)) {
  cfint = rbind(cfint,cfint_list[[kk]])
}

# plotting (faceting by category)
mRat  = ggplot(cfint, aes(HCPG, mean,fill=shift_gamble_phase))
mRat  = mRat + labs(x='Group', y=paste('Mean (',0.95*100,'% CI, bootstrapped)'))
mRat  = mRat + ggtitle("Physio reaction (PSPMNL) shifts in gamble phase")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(variable ~ .)
mRatg = mRat + theme_bw()
print(mRatg)


# # stats
# mod_cat             = lmer(FlexibleResponse_1_Amplitude ~ cat + (cat|subject) + (cat|stim), data=data_pdt)
# mod_catg            = lmer(FlexibleResponse_1_Amplitude ~ cat*HCPG + (cat|subject) + (cat|stim), data=data_pdt)
# agk.lme.summary(mod_cat,'norm')
# agk.lme.summary(mod_catg,'norm')
# mod_catl            = lmer(FlexibleResponse_1_PeakLatency ~ cat + (cat|subject) + (cat|stim), data=data_pdt)
# mod_catlg           = lmer(FlexibleResponse_1_PeakLatency ~ cat*HCPG + (cat|subject) + (cat|stim), data=data_pdt)
# agk.lme.summary(mod_catl,'norm')
# agk.lme.summary(mod_catlg,'norm')
# mod_catd            = lmer(FlexibleResponse_1_Dispersion ~ cat + (cat|subject) + (cat|stim), data=data_pdt)
# mod_catdg           = lmer(FlexibleResponse_1_Dispersion ~ cat*HCPG + (cat|subject) + (cat|stim), data=data_pdt)
# agk.lme.summary(mod_catd,'norm')
# agk.lme.summary(mod_catdg,'norm')
# 
# 
# 


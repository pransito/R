# here we plot the estimated rating means per picture category

# dependence
# need to run import_data.R first
# then run analysis_PDT_physio_rating_LA_6.R until
# model fit3e_0pg

# the models which give fixed effects (means) and cnfints
modlist=list()
modlist[[1]] = fit3_0
modlist[[2]] = fit3b_0
modlist[[3]] = fit3c_0
modlist[[4]] = fit3d_0
modlist[[5]] = fit3e_0

# the models which give fixed effects (means) and cnfints by group (HC)
modlisthc=list()
modlisthc[[1]] = fit3_0hc
modlisthc[[2]] = fit3b_0hc
modlisthc[[3]] = fit3c_0hc
modlisthc[[4]] = fit3d_0hc
modlisthc[[5]] = fit3e_0hc

# ...(PG)
modlistpg=list()
modlistpg[[1]] = fit3_0pg
modlistpg[[2]] = fit3b_0pg
modlistpg[[3]] = fit3c_0pg
modlistpg[[4]] = fit3d_0pg
modlistpg[[5]] = fit3e_0pg

# insert check how many subjects are there in PG and HC
# check in the actual models used in the rating pictures
# generation
agk.check.whoisthere =  function(cur_mod) {
  # function for glmer model
  # to count how many subjects (ranef)
  # are there
  cur_ranef = ranef(cur_mod)
  return(length(cur_ranef[[1]][,1]))
}
# HCs
lapply(modlisthc,FUN=agk.check.whoisthere)
lapply(modlistpg,FUN=agk.check.whoisthere)

# some params to set
method = "Wald"
level  = 0.95

# a dictionary to translate variable names
cur_src = c("imageRating1s","imageRating2s","imageRating3s","imageRating4s")
cur_map = c("Elicits_Craving","Representative_for_Gambling","Representative_for_Negative","Representative_for_Positive")

# function to get the mean and lower and upper CI for a model
agk.get.meanCI.lmer = function(mod, method, level) {
  parm = names(fixef(mod))
  cnfint = confint(mod,method=method,parm = parm)
  cnfint = cbind(fixef(mod),cnfint)
  cnfint = as.data.frame(cnfint)
  cnfint = cbind(cnfint,rep(all.vars(mod@call[[2]])[1],length(parm)))
  colnames(cnfint) = c("me","lo","up","var")
  return(cnfint)
}

# now extract from all models the means and CIs
allcfint = lapply(modlist,FUN = agk.get.meanCI.lmer,method=method,level=level)
cmpcfint = allcfint[[1]]
for (ii in 2:length(allcfint)) {cmpcfint = rbind(cmpcfint,allcfint[[ii]])}
cmpcfint$var = agk.recode.c(cmpcfint$var,cur_src,cur_map)

# make the category variable pretty
ct_var  = row.names(cmpcfint)
ct_var  = strsplit(ct_var,"cat")
ct_vars = c()
for (ii in 1:length(ct_var)) {ct_vars[ii] = ct_var[[ii]][2]}
ct_vars = gsub('[0-9]+', '', ct_vars)
cmpcfint$cat = ct_vars
row.names(cmpcfint) = NULL
names(cmpcfint)[c(4,5)] = c("Rating","Category")

# no group
mRat   = ggplot(cmpcfint, aes(Category, me,fill=Rating))
mRat   = mRat + labs(x='Image Category', y=paste('Fixed Effect (',level*100,'% CI, bootsrapped)'))
mRat   = mRat + ggtitle("Stimuli validation")
mRat   = mRat + geom_bar(position="dodge", stat="identity")
dodge  = position_dodge(width=0.9)
mRat   = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
mRatng = mRat + theme_bw()
mRatng

# bygroup
# now extract from all models the means and CIs
allcfintHC = lapply(modlisthc,FUN = agk.get.meanCI.lmer,method=method,level=level)
allcfintPG = lapply(modlistpg,FUN = agk.get.meanCI.lmer,method=method,level=level)
for (ii in 1:length(allcfintHC)) {allcfintHC[[ii]]$group = "HC"}
for (ii in 1:length(allcfintPG)) {allcfintPG[[ii]]$group = "PG"}

cmpcfintHC = allcfintHC[[1]]
for (ii in 2:length(allcfintHC)) {cmpcfintHC = rbind(cmpcfintHC,allcfintHC[[ii]])}
cmpcfintPG = allcfintPG[[1]]
for (ii in 2:length(allcfintPG)) {cmpcfintPG = rbind(cmpcfintPG,allcfintPG[[ii]])}
cmpcfint   = rbind(cmpcfintHC,cmpcfintPG)
cmpcfint$var = agk.recode.c(cmpcfint$var,cur_src,cur_map)

# make the category variable pretty
ct_var  = row.names(cmpcfint)
ct_var  = strsplit(ct_var,"cat")
ct_vars = c()
for (ii in 1:length(ct_var)) {ct_vars[ii] = ct_var[[ii]][2]}
ct_vars = gsub('[0-9]+', '', ct_vars)
cmpcfint$cat = ct_vars
row.names(cmpcfint) = NULL
names(cmpcfint)[c(4,6)] = c("Rating","Category")

# plotting
mRat = ggplot(cmpcfint, aes(Category, me,fill=Rating))
mRat = mRat + labs(x='Image Category', y=paste('Fixed Effect (',level*100,'% CI, Walt method)'))
mRat = mRat + ggtitle("Stimuli validation")
mRat = mRat + geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim = c(-1, 3.5)) 
dodge <- position_dodge(width=0.9)
mRat = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
mRat = mRat + facet_grid(group ~ .)
mRat + theme_bw()

# plotting (faceting by category)
mRat  = ggplot(cmpcfint, aes(group, me,fill=Rating))
mRat  = mRat + labs(x='Group', y=paste('Fixed Effect (',level*100,'% CI, Wald method)'))
mRat  = mRat + ggtitle("Stimuli validation")
mRat  = mRat + geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim = c(-1, 3.5)) 
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
mRat  = mRat + facet_grid(Category ~ .)
mRatg = mRat + theme_bw()
mRatg = mRatg + theme(text = element_text(size=20)) 
mRatg
















cur_vars = c('aff.1', 'aff.2', 'aff.3','aff.4','affup.1', 'affup.2', 'affup.3','affup.4', 'afflo.1', 'afflo.2', 'afflo.3','afflo.4',
             'aff.5','affup.5','afflo.5','aff.6','affup.6','afflo.6','aff.7','affup.7','afflo.7',
             'HCPG','cat')
data_pdt_aff1=data_pdt_meci2[cur_vars]
data_pdt_aff=reshape(data_pdt_aff1, varying =  cur_vars[-c(length(cur_vars)-1,length(cur_vars))],  direction = "long" )
data_pdt_aff$var=factor(data_pdt_aff$time, labels = c('Arousal', 'Dominance', 'Valence','Craving','Gambling','Negative','Positive'))
curvar='aff'

mRat= ggplot(data_pdt_aff, aes_string('cat', paste0(curvar,''), fill='HCPG')) +  labs(title=curvar, x='Image Category', y='Mean (90% CI, bootsrapped)') 
mRat = mRat + geom_bar(position="dodge", stat="identity")
mRat = mRat+facet_grid(var ~ .)
dodge <- position_dodge(width=0.9)
mRat = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes_string(ymin = paste0(curvar,'lo'), ymax=paste0(curvar,'up')), position=dodge, width=0.25)
mRat + theme_bw()

# here we plot the estimated rating means per picture category

# dependence
# need to run import_data.R first
# then run select_study.R
# then run image_adequacy.R for ratings and physio until
# model fit3e_0pg

# the models which give fixed effects (means) and cnfints
modlist=list()
modlist[[1]] = fit9c
modlist[[2]] = fit9z
modlist[[3]] = fit9e

# TODO: MAKE ALSO WORK WITH NORMAL MEANS (LMLIST)

if (!which_study == 'sanity') {
  # the models which give fixed effects (means) and cnfints by group (HC)
  modlisthc=list()
  modlisthc[[1]] = fit9c_hc
  modlisthc[[2]] = fit9z_hc
  modlisthc[[3]] = fit9e_hc
  
  # ...(PG)
  modlistpg=list()
  modlistpg[[1]] = fit9c_pg
  modlistpg[[2]] = fit9z_pg
  modlistpg[[3]] = fit9e_pg
}

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
# HCs and PGs
print(lapply(modlist,FUN=agk.check.whoisthere))
if (which_study != 'sanity') {
  print(lapply(modlisthc,FUN=agk.check.whoisthere))
  print(lapply(modlistpg,FUN=agk.check.whoisthere))
}

# some params to set
method = "Wald"
level  = 0.95

# a dictionary to translate variable names
cur_src = c("corr","zygo","eda")
cur_map = c("Corrugator","Zygomaticus","EDA")

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
names(cmpcfint)[c(4,5)] = c("Physio","Category")

# cleaning in case "intercept" went to NA
cmpcfint[is.na(cmpcfint)] = 'neutral'

# now add neutral category to all means so that
# we have actual means per cat
cat_of_int = unique(cmpcfint$Category)[which(!unique(cmpcfint$Category) %in% 'neutral')]
neu_to_add = cmpcfint[cmpcfint$Category == 'neutral',c('me','lo','up')]
for (ii in 1:length(cat_of_int)) {
  cmpcfint[cmpcfint$Category == cat_of_int[ii],c('me','lo','up')] + neu_to_add
}


# no group
mRat   = ggplot(cmpcfint, aes(Category, me,fill=Physio))
mRat   = mRat + labs(x='Image Category', y=paste('Fixed Effect (',level*100,'% CI, bootsrapped)'))
mRat   = mRat + ggtitle("Stimuli validation")
mRat   = mRat + geom_bar(position="dodge", stat="identity")
dodge  = position_dodge(width=0.9)
mRat   = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
mRatng = mRat + theme_bw()
print(mRatng)

if (!which_study == 'sanity') {
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
  names(cmpcfint)[c(4,6)] = c("Physio","Category")
  
  # cleaning in case "intercept" went to NA
  cmpcfint[is.na(cmpcfint)] = 'neutral'
  
  # plotting
  mRat = ggplot(cmpcfint, aes(Category, me,fill=Physio))
  mRat = mRat + labs(x='Image Category', y=paste('Fixed Effect (',level*100,'% CI, Wald method)'))
  mRat = mRat + ggtitle("Stimuli validation")
  mRat = mRat + geom_bar(position="dodge", stat="identity")
  dodge <- position_dodge(width=0.9)
  mRat = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
  mRat = mRat + facet_grid(group ~ .)
  print(mRat + theme_bw())
  
  # plotting (faceting by category)
  mRat  = ggplot(cmpcfint, aes(group, me,fill=Physio))
  mRat  = mRat + labs(x='Group', y=paste('Fixed Effect (',level*100,'% CI, Wald method)'))
  mRat  = mRat + ggtitle("Stimuli validation")
  mRat  = mRat + geom_bar(position="dodge", stat="identity")
  dodge = position_dodge(width=0.9)
  mRat  = mRat + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin = lo, ymax = up), position=dodge, width=0.25)
  mRat  = mRat + facet_grid(Category ~ .)
  mRatg = mRat + theme_bw()
  print(mRatg)
}

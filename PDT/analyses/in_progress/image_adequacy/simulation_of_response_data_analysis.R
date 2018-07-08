# Plots means and CIs of Katha's TEST session
data_pdt$HCPG = as.character(data_pdt$HCPG)
data_pdt$HCPG = as.factor(agk.recode(data_pdt$HCPG,c('VPPG7777','VPPG8888'),c('nosim','sim')))

# means and CIs
des_vars = c('zygo_auc_stim', 'zygo_max_stim','zygo_median_stim',
                   'corr_auc_stim','corr_max_stim','corr_median_stim',
                   'eda_auc_stim','eda_max_stim','eda_median_stim')
cfint_list = list()
for (ii in 1:length(des_vars)) {
  cur_form = as.formula(paste(des_vars[ii],'~ HCPG + cat'))
  cfint = aggregate(cur_form,data = data_pdt,FUN=agk.boot.ci,
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
mRat  = mRat + ggtitle("Simulation of reaction in vivo")
mRat  = mRat + geom_bar(position="dodge", stat="identity")
dodge = position_dodge(width=0.9)
mRat  = mRat + geom_bar(position=dodge, stat="identity")
mRat  = mRat + geom_errorbar(aes(ymin = lower, ymax = upper), position=dodge, width=0.25) 
mRat  = mRat + facet_grid(cat ~ .)
mRatg = mRat + theme_bw()
print(mRatg)

# analysis
for (ii in 1:length(des_vars)) {
  print(des_vars[ii])
  cur_form = as.formula(paste(des_vars[ii],'~cat*HCPG'))
  model    = lm(cur_form, data=data_pdt)
  print(summary(model))
}

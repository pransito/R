# function to perform median difference tests (using group label permutation)
median_rmna = function(x) {
  if(!(is.numeric(x))) {return(NA)}
  if(all(is.na(x))) {return(NA)}
  return(median(x,na.rm = T))
}
var_cis = aggregate(behav_m,by=list(behav_m$HCPG,behav_m$variable,behav_m$model), FUN=median_rmna)

var_cis$subject  = NULL
var_cis$HCPG     = NULL
var_cis$variable = NULL
var_cis$model    = NULL

var_cis_df        = data.frame(var_cis)
names(var_cis_df) = c('HCPG','behav_param','model','median')
var_cis_df        = var_cis_df[grep('KFG',var_cis_df$behav_param,invert = T),]
names(var_cis_df)[names(var_cis_df) == 'HCPG'] = 'Group'
var_cis_df$Group = agk.recode.c(var_cis_df$Group,'PG','GD')

var_cis_df_GD = subset(var_cis_df,Group == 'GD')
var_cis_df_HC = subset(var_cis_df,Group == 'HC')

var_cis_df_di        = var_cis_df_HC
var_cis_df_di        = merge(var_cis_df_di,var_cis_df_GD,by=c('behav_param','model'))
var_cis_df_di$median = var_cis_df_di$median.x - var_cis_df_di$median.y

# difference under 0-hypothesis
diff_collected = list()
for (pp in 1:1000) {
  var_cis = aggregate(behav_m,by=list(behav_m$HCPG,behav_m$variable,behav_m$model), FUN=median_rmna)
  
  var_cis$subject  = NULL
  var_cis$HCPG     = NULL
  var_cis$variable = NULL
  var_cis$model    = NULL
  
  var_cis_df        = data.frame(var_cis)
  names(var_cis_df) = c('HCPG','behav_param','model','median')
  var_cis_df        = var_cis_df[grep('KFG',var_cis_df$behav_param,invert = T),]
  names(var_cis_df)[names(var_cis_df) == 'HCPG'] = 'Group'
  var_cis_df$Group = agk.recode.c(var_cis_df$Group,'PG','GD')
  
  var_cis_df_GD = subset(var_cis_df,Group == 'GD')
  var_cis_df_HC = subset(var_cis_df,Group == 'HC')
  
  # permutation
  for (ff in 1:length(var_cis_df_GD[,1])) {
    if (randn() > 0) {
      cur_line_GD        = var_cis_df_GD[ff,]
      var_cis_df_GD[ff,] = var_cis_df_HC[ff,]
      var_cis_df_HC[ff,] = cur_line_GD
    }
  }
  
  var_cis_df_di        = var_cis_df_HC
  var_cis_df_di        = merge(var_cis_df_di,var_cis_df_GD,by=c('behav_param','model'))
  var_cis_df_di$median = var_cis_df_di$median.x - var_cis_df_di$median.y
  
  if(is.na(var_cis_df_di$median[35])) {stop()}
  diff_collected[[pp]] = var_cis_df_di
}

diffs = t(as.matrix(diff_collected[[1]]$median))
for (pp in 1:length(diff_collected)) {
  diffs = rbind(diffs,t(as.matrix(diff_collected[[pp]]$median)))
}
diffs      = data.frame(diffs)

# now test the hypotheses (group differences)
var_cis_df_di$p_value_two_sided = NA
for (pp in 1:length(var_cis_df_di$p_value_two_sided)) {
  var_cis_df_di$p_value_two_sided[pp] = agk.density_p(diffs[[pp]],var_cis_df_di$median[pp],type = 'two.sided')
}


# make values for stimuli
# just takes all the variables rating and physio and computes the mean across subjects per stimulus
# uses subjects not included in the current study (i.e. need the data_pdt_inv)
# give a data frame
# data frame must have a stimulus variable and multiple contributing to value
# will do it for PG and HC split

# prep
data_pdt_bb2 = data_pdt
smf      = function(x) {return(median(x,na.rm = T))}
gs       = function(x) {return(sum(x,na.rm = T))}
std_fun  = function(x) {
  mn_v       = mean.rmna(x)
  sd_v       = sd(x,na.rm = T)
  x          = get_log((x-mn_v)/sd_v)
  x          = agk.scale_range(x,-1,1)
  return(x)
}

# HC or PG or both?
cur_tab = as.data.frame(table(data_pdt_inv$HCPG))
if (cur_tab[cur_tab$Var1 == "HC",]$Freq == 0) {
  doHC = 0
} else {
  doHC = 1
}
if (cur_tab[cur_tab$Var1 == "PG",]$Freq == 0) {
  doPG = 0
} else {
  doPG = 1
}

if (doHC == 1) {
  # case HC
  data_pdt       = subset(data_pdt_inv,HCPG=="HC")
  pic_value_vars = c("arousal","dominance","valence","imageRating1s",
                     "imageRating2s","imageRating3s","imageRating4s","imageRating5s")
  
  by_stim               = data_pdt$stim
  cur_dat               = data_pdt[pic_value_vars]
  cur_agg               = aggregate(cur_dat,by = list(by_stim),FUN = smf)
  #cur_agg[-1]           = apply(MARGIN = 2,X = cur_agg[-1],FUN = std_fun)
  cur_aggHC             = cur_agg
  names(cur_aggHC)[-1]  = paste0('vHC_',names(cur_aggHC)[-1]) 
  names(cur_aggHC)[1]   = 'stim'
}

if (doPG == 1) {
  # case PG
  data_pdt       = subset(data_pdt_inv,HCPG=="PG")
  pic_value_vars = c("arousal","dominance","valence","imageRating1s",
                     "imageRating2s","imageRating3s","imageRating4s","imageRating5s")
  
  by_stim               = data_pdt$stim
  cur_dat               = data_pdt[pic_value_vars]
  cur_agg               = aggregate(cur_dat,by = list(by_stim),FUN = smf)
  #cur_agg[-1]           = apply(MARGIN = 2,X = cur_agg[-1],FUN = std_fun)
  cur_aggPG             = cur_agg
  names(cur_aggPG)[-1]  = paste0('vPG_',names(cur_aggPG)[-1]) 
  names(cur_aggPG)[1]   = 'stim'
}


# case PG and HC together
data_pdt              = data_pdt_inv
pic_value_vars        = c("arousal","dominance","valence","imageRating1s",
                         "imageRating2s","imageRating3s","imageRating4s","imageRating5s")
by_stima              = data_pdt$stim
cur_data              = data_pdt[pic_value_vars]
cur_agga              = aggregate(cur_data,by = list(by_stima),FUN = smf)
#cur_agga[-1]          = apply(MARGIN = 2,X = cur_agga[-1],FUN = std_fun)
names(cur_agga)[-1]   = paste0('v_oa_',names(cur_agga)[-1]) 
names(cur_agga)[1]    = 'stim'

# merging
data_pdt   = data_pdt_bb2
if(doHC == 1) {
  data_pdt   = merge(data_pdt,cur_aggHC,by.x = c("stim"),by.y = c("stim"),all.x = T)
}
if (doPG == 1) {
  data_pdt   = merge(data_pdt,cur_aggPG,by.x = c("stim"),by.y = c("stim"),all.x = T)
}
data_pdt   = merge(data_pdt,cur_agga,by.x = c("stim"),by.y = c("stim"),all.x = T)

# selecting according to PG or HC
v_value_vars           = paste0('v_',pic_value_vars)
data_pdt[v_value_vars] = NA
ct = 0
ct_2 = 0

for(ii in 1:length(data_pdt[,1])) {
  if(is.na(data_pdt$HCPG[ii])) {
    ct = ct + 1
    disp(paste('NA cannot assign v_variable',ct))
    next
  }
  if(data_pdt$HCPG[ii] == "HC" & doHC == 1) {
    data_pdt[v_value_vars][ii,] = data_pdt[names(cur_aggHC)[-1]][ii,]      
  } else if (data_pdt$HCPG[ii] == "PG" & doPG == 1) {
    data_pdt[v_value_vars][ii,] = data_pdt[names(cur_aggPG)[-1]][ii,]
  } else {
    ct_2 = ct_2 +1
    disp(paste('cannot assign v_variable',ct_2))
    next
  }
}

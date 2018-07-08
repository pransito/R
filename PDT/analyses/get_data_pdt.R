# message to user
print("Getting behavioral choice data")

# path
path_cur_home = getwd()

# prep some variables
subject                        = c()
missing_people                 = c()

# get the Behav PDT folders of the MRI subjects
all_behav_dir_mri = c()
if (get_MRI_behav_data == 1) {
  setwd(path_mrt)
  cur_folders = dir(pattern = "VPPG")
  for (ii in 1:length(cur_folders)) {
    setwd(path_mrt)
    setwd(cur_folders[ii])
    setwd("PDT")
    all_behav_dir_mri[ii] = getwd()
  }
}

# get all the necessary .txt files from the current directory
# add the PG data to it (from second directory)
setwd(path_dat)
all_behav_dir       = dir(pattern="PhysioVP")
all_behav_dir       = paste(path_dat,all_behav_dir,sep="/")

setwd(path_pg)
all_behav_dir_pg    = dir(pattern="PhysioVP")
all_behav_dir_pg    = paste(path_pg,all_behav_dir_pg,sep="/")

setwd(path_postpilot_pg)
all_behav_dir_pp_pg = dir(pattern="PhysioVP")
all_behav_dir_pp_pg = paste(path_postpilot_pg,all_behav_dir_pp_pg,sep="/")

setwd(path_postpilot_hc)
all_behav_dir_pp_hc = dir(pattern="PhysioVP")
all_behav_dir_pp_hc = paste(path_postpilot_hc,all_behav_dir_pp_hc,sep="/")

all_behav_dir = c(all_behav_dir, all_behav_dir_pg, all_behav_dir_pp_hc, all_behav_dir_pp_pg)
if (get_MRI_behav_data == 1) {
  all_behav_dir = c(all_behav_dir,all_behav_dir_mri)
}

all_behav_files = list()
all_demog_files = list()
all_stim_names  = list()
age             = c()
sex             = c()

# setting up progress bar
total = length(all_behav_dir)
pb    = txtProgressBar(min = 0, max = total, style = 3)

# get the behav data
for (ii in 1:total) {
  # cd in directory and get output file
  setwd(all_behav_dir[ii])
  cur_behav_file           = list.files(all.files =T, pattern="*_output_extr.txt")
  cur_demog_file           = list.files(all.files = T, pattern = "*_demograph_extr.txt")
  all_behav_files[[ii]]    = cur_behav_file
  all_demog_files[[ii]]    = cur_demog_file
  
  # read in choice data for this subject
  setwd(all_behav_dir[ii])
  cur_df = read.table(all_behav_files[[ii]], header=T, sep="\t")
  
  if (length(grep('VPPG0115',all_behav_dir[ii])) != 0) {
    ## EXTRA READ IN FOR VPPG0115
    cur_df_ex = list()
    # build its own behavioral data frame
    cur_home = getwd()
    cd('0115a')
    P_1 = readMat('P_0115.mat')
    
    # session one
    gainex = c()
    lossex = c()
    chcex  = c()
    rtex   = c()
    catex  = c()
    sideex = c()
    st_dux = c()
    stimex = c()
    for (ff in 1:101) {
      # cur_gamble
      cur_gam  = P_1[[1]][[45]][[3]][[ff]]
      cur_gain = as.numeric(P_1[[1]][[18]][[1]][[cur_gam[[1]][1]]])
      cur_loss = as.numeric(P_1[[1]][[19]][[1]][[cur_gam[[1]][2]]])

      # fill
      gainex[ff] = cur_gain
      lossex[ff] = cur_loss
      chcex[ff]  = P_1[[1]][[45]][[grep('choice',row.names(P_1[[1]][[45]]))]][ff]
      rtex[ff]   = P_1[[1]][[45]][[grep('rt',row.names(P_1[[1]][[45]]))]][ff]
      catex[ff]  = P_1[[1]][[45]][[grep('cat',row.names(P_1[[1]][[45]]))]][ff]
      sideex[ff] = P_1[[1]][[45]][[grep('side',row.names(P_1[[1]][[45]]))]][ff]
      st_dux[ff] = P_1[[1]][[45]][[grep('dur',row.names(P_1[[1]][[45]]))]][ff]
      stimex[ff] = P_1[[1]][[45]][[grep('^stim$',row.names(P_1[[1]][[45]]))]][[ff]][[1]]
    }
    
    # session two
    cd(cur_home)
    cd('0115b')
    P_1 = readMat('P_0115.mat')
    
    for (ff in 1:101) {
      # cur_gamble
      cur_gam  = P_1[[1]][[45]][[3]][[ff+101]]
      cur_gain = as.numeric(P_1[[1]][[18]][[1]][[cur_gam[[1]][1]]])
      cur_loss = as.numeric(P_1[[1]][[19]][[1]][[cur_gam[[1]][2]]])
      
      # fill
      gainex[ff+101] = cur_gain
      lossex[ff+101] = cur_loss
      chcex[ff+101]  = P_1[[1]][[45]][[grep('choice',row.names(P_1[[1]][[45]]))]][ff+101]
      rtex[ff+101]   = P_1[[1]][[45]][[grep('rt',row.names(P_1[[1]][[45]]))]][ff+101]
      catex[ff+101]  = P_1[[1]][[45]][[grep('cat',row.names(P_1[[1]][[45]]))]][ff+101]
      sideex[ff+101] = P_1[[1]][[45]][[grep('side',row.names(P_1[[1]][[45]]))]][ff+101]
      st_dux[ff+101] = P_1[[1]][[45]][[grep('dur',row.names(P_1[[1]][[45]]))]][ff+101]
      stimex[ff+101] = P_1[[1]][[45]][[grep('^stim$',row.names(P_1[[1]][[45]]))]][[ff+101]][[1]]
    }
    
    cur_df_ex$gain   = gainex
    cur_df_ex$loss   = lossex
    cur_df_ex$choice = chcex
    cur_df_ex$rt     = rtex
    cur_df_ex$cat    = catex
    cur_df_ex$side   = sideex
    cur_df_ex$st_dur = st_dux
    cur_df_ex$stim   = stimex
    cur_df = as.data.frame(cur_df_ex)
    cd(cur_home)
  }
  
  # read in demogr data for this sub
  # account for multiple header lines; last line should be data
  cur_dm = read.table(all_demog_files[[ii]], header=T, sep="\t")
  sex[ii]    = as.character(cur_dm$sex[length(cur_dm$sex)])
  age[ii]    = as.numeric(as.character(cur_dm$age[length(cur_dm$age)]))
  
  ## TAKE OUT
  # # add the stimulus name info
  # cur_df$stim = all_stim_names[[ii]]
  ## TAKE OUT
  
  # add a trial number var
  cur_df$trial = seq(1:length(cur_df[,1]))
  
  # add age and sex
  cur_df$age = age[ii]
  cur_df$sex = sex[ii]
  
  # formatting
  cur_df$gain          = as.numeric(cur_df$gain)
  cur_df$loss          = as.numeric(cur_df$loss)
  cur_df$gain_bcp      = cur_df$gain
  cur_df$loss_bcp      = cur_df$loss
  cur_df$gainxloss     = cur_df$gain_bcp*cur_df$loss_bcp
  cur_df$gainxloss_bcp = cur_df$gainxloss 
  
  # aggregating
  cur_df$gain      = agk.aggregate.data.la.gain.loss.c(cur_df$gain,agg=cur_agg)
  cur_df$loss      = agk.aggregate.data.la.gain.loss.c(cur_df$loss,agg=cur_agg)
  cur_df$gainxloss = agk.aggregate.data.la.gain.loss.c(cur_df$gainxloss,agg=cur_agg)
  
  # euclidean distance
  # loss must be used as abs(loss)
  # result is automatically absolute values
  # note that if aggregation has taken place; ed is performed on aggr values
  # as it should be
  ed_neu = c()
  for (ll in 1:length(cur_df[,1])) {
    ed_neu[ll] = agk_get_ed.c(c(cur_df$gain[ll],abs(cur_df$loss[ll]),0),sp = c(26,13,0),vec = c(2,1,0))
  }
  
  cur_df$ed_abs = ed_neu
  cur_df$ed     = NA
  
  # save abs(loss) if desired
  if (loss_abs == 1) {
    cur_df$loss = abs(as.numeric(cur_df$loss)) # work as Tom et al. (2007)
  }
  
  # Risk
  # die Formeln von Minati und Martino sind gleich; allerdings finde ich,
  # man muesste den Term noch durch die Anzahl an Optionen teilen
  # also durch 2
  if (loss_abs == 1) {
    cur_df$EV      = 0.5*cur_df$gain + 0.5*cur_df$loss*(-1)
    cur_df$ratio   = abs(cur_df$gain/cur_df$loss*(-1))
    cur_df$diff    = cur_df$gain-cur_df$loss*(-1)
  } else {
    cur_df$EV      = 0.5*cur_df$gain + 0.5*cur_df$loss
    cur_df$ratio   = abs(cur_df$gain/cur_df$loss)
    cur_df$diff    = cur_df$gain-cur_df$loss
  }
  
  if (loss_abs == 1) {
    cur_df$RiskMin = 0.5*(cur_df$gain-cur_df$EV)^2+0.5*(cur_df$gain-cur_df$EV)^2 # typo??
    cur_df$RiskAG  = 0.5*(cur_df$gain-cur_df$EV)^2+0.5*(cur_df$loss*(-1)-cur_df$EV)^2
    cur_df$RiskMar = (0.5*cur_df$gain-0.5*cur_df$loss*(-1))^2
  } else {
    cur_df$RiskMin = 0.5*(cur_df$gain-cur_df$EV)^2+0.5*(cur_df$gain-cur_df$EV)^2 # typo??
    cur_df$RiskAG  = 0.5*(cur_df$gain-cur_df$EV)^2+0.5*(cur_df$loss-cur_df$EV)^2
    cur_df$RiskMar = 0.5*cur_df$gain-0.5*cur_df$loss^2
  }
  
  
  # Skew aufbauend auf Formel von Minati
  cur_df$SkewMin  = 0.5*(cur_df$gain-cur_df$EV)^3+0.5*(cur_df$gain-cur_df$EV)^3
  
  # scaling
  if (use_z == 1) {
    cur_df$EV_bcp    = cur_df$EV
    cur_df$gain      = scale(cur_df$gain, center = T, scale = F)
    cur_df$loss      = scale(cur_df$loss, center = T, scale = F)
    cur_df$gainxloss = scale(cur_df$gainxloss, center = T, scale = F)
    cur_df$EV        = scale(cur_df$EV, center = T, scale = F)
    cur_df$diff      = scale(cur_df$diff, center = T, scale = F)
    cur_df$ratio     = scale(cur_df$ratio, center = T, scale = F)
    cur_df$ed_abs    = scale(cur_df$ed_abs, center = T, scale = F)
  }
  
  # the choice variable
  accept_reject = cur_df$choice
  
  # checking for missings
  # if there are too many, that person will be dropped
  if (missing_check == 1) {
    m_test = which(accept_reject==5)
    if (length (m_test) > missing_cutoff*length(cur_df[,1])) {
      missing_people[length(missing_people)+1] = all_behav_files[[ii]]
      next
    }
  }
  
  # create subject variable
  cur_sub = all_behav_dir[[ii]]
  cur_sub = strsplit(cur_sub,"/",fixed = T)
  
  if (sum(cur_sub[[1]] == "MRT") > 0) {
    cur_sub = cur_sub[[1]][(length(cur_sub[[1]])-1)]
  } else {
    cur_sub = cur_sub[[1]][length(cur_sub[[1]])]
  }
  
  subject = rep(cur_sub, length = length(cur_df[,1]))
  
  # it so happens that we have missing data in the categorical response variable,
  # we have to correct for this by setting "5" in accept_reject to "NA"
  # and we recode the ordinal response scale into a dichotomous response scale
  # recode accept_reject "5" into NA
  # recode 1,2 into 1 and 3,4 into 0
  accept_reject = as.factor(agk.recode.c(accept_reject,c("1","2","3","4","5"),c("1","1","0","0",NA)))
  cur_df = data.frame(subject,cur_df, accept_reject)
  
  # get the BOLD extract data
  is_MRI = as.logical(length(grep(all_behav_dir[ii],pattern = "MRT")))
  if(get_MRI_extr == 1 & is_MRI & as.logical(length(grep(dir(),pattern = "bold_scores_")))) {
    rm(list=c("tmp"))
    
    # now we have to merge all the bold extract
    extracts  = dir(pattern = "bold_scores_")
    read_extr = list()
    for (kk in 1:length(extracts)) {
      cur_extract     = extracts[kk]
      read_extr[[kk]] = read.table(cur_extract,header=T,sep="\t")
      cur_name = gsub('[.]txt','',cur_extract)
      names(read_extr[[kk]])[3] = cur_name
    }
    tmp = read_extr[[1]]
    if (length(read_extr) > 1) {
      for (kk in 2:length(read_extr)) {
        tmp = merge(tmp,read_extr[[kk]],by = c("trial","stim"))
      }
    }
    
    cur_df        = merge(cur_df,tmp,by=c("trial","stim"))
    
  } else {
    cur_df$bold_scores_l_precuneus = NA
    cur_df$bold_scores_post_cing   = NA
    cur_df$bold_scores_r_precuneus = NA
    cur_df$bold_scores_sup_temp    = NA
  }
  
  if (!exists ("data_pdt")) {data_pdt = cur_df} else if (exists ("data_pdt")) {
    data_pdt = rbind.fill(data_pdt, cur_df)}
  rm(list = c("cur_df"))
  
  # update progress bar
  #Sys.sleep(0.05)
  setTxtProgressBar(pb, ii)
}

# close progress bar
close(pb)

## data handling
data_pdt$gain          = as.numeric(data_pdt$gain)
data_pdt$loss          = as.numeric(data_pdt$loss)
data_pdt$subject       = as.factor(data_pdt$subject)
data_pdt$diff          = as.numeric(data_pdt$diff)
data_pdt$subject       = as.factor(data_pdt$subject)
data_pdt$cat           = as.factor(data_pdt$cat)
data_pdt$accept_reject = as.factor(data_pdt$accept_reject)
data_pdt$side          = factor(as.character(data_pdt$side), levels=c("0","1"), labels=c("gainleft","gainright"))

## cleaning RT
data_pdt$rt[data_pdt$rt == 99999] = NA

## here we clean the stim variable
cur_stim        = strTrim(as.character(data_pdt$stim))
cur_stim        = unlist(strsplit(cur_stim,split = ".jpg"))
cur_stim        = strTrim(cur_stim)
cur_stim        = gsub(cur_stim,pattern = " \\(1\\)",replacement = "")
data_pdt$stim   = cur_stim

# throw out NA rows of subject
data_pdt_sub_na = data_pdt[is.na(data_pdt$subject),]
if(length(data_pdt_sub_na[,1]) > 0) {
  stop('data_pdt has NAs in subject variable. Check this!')
}
data_pdt        = data_pdt[!is.na(data_pdt$subject),]

## here we clean the subject variable
data_pdt$subject = gsub('(1)','',data_pdt$subject, fixed=T)
data_pdt$subject = trimws(data_pdt$subject)

## here we recode the subject variable to VPPG numbers only
data_pdt$subject_new = agk.recode.c(data_pdt$subject,tnl$PhysioVP,tnl$VPPG)
data_pdt$subject_new = as.factor(data_pdt$subject_new)
data_pdt$subject     = data_pdt$subject_new
data_pdt$subject_new = NULL

# prepping tnl
include_df = tnl[c("VPPG","Einschluss","Bilderrating_gemacht")]

# getting all subs that are included in tnl
expected_VPPG = tnl$VPPG[tnl$Einschluss == 1]
excluded_VPPG = tnl$VPPG[tnl$Einschluss == 0]

# some subs may be exempt
expected_VPPG = expected_VPPG[!expected_VPPG %in% behav_exempt]

# dropping subs according to inclusion variable in tnl
data_pdt           = merge(data_pdt, include_df, by.x=c("subject"),by.y=c("VPPG"),all.x = T,all.y=F)
data_pdt_noex      = data_pdt
data_pdt_nona      = subset(data_pdt, !is.na(Einschluss))
data_pdt_na        = subset(data_pdt, is.na(Einschluss))
data_pdt           = subset(data_pdt_nona, Einschluss == "1")
# reporting
exl_subs_nona = unique(data_pdt_nona$subject[which(!data_pdt_nona$subject %in% data_pdt$subject)])
message = paste('Just excluded these subs according to include/excl variable in TN-Liste:\n',
                paste(exl_subs_nona,collapse = ' '))
cat(message)
# reporting on subjects in data_pdt which are not in TN-Liste
exl_subs_na = unique(data_pdt_na$subject)
if (length(exl_subs_na) != 0) {
  mes_1 = 'There are subs that are in data_pdt which are not in TN-List. Check!\n'
  mes_2 = paste(exl_subs_na,collapse=' ')
  stop(paste(mes_1,mes_2))
}

# checking whether data_pdt has all the expected (according to tnl) VPPG subjects
cur_tfv = expected_VPPG %in% unique(data_pdt$subject)
if (!all(cur_tfv)) {
  cur_missing_subs = expected_VPPG[which(!cur_tfv)]
  # excluding the behav_exempt subs
  cur_missing_subs = cur_missing_subs[!cur_missing_subs %in% behav_exempt]
  mes_1 = 'Some subs in tnl which are included have no behav data.'
  mes_2 = paste(cur_missing_subs,collapse = ' ')
  stop(paste(mes_1,mes_2))
}

# for Pretest we need to establish a cat variable
# do we need this?!?!?
# str_fun = function(x) {return(substr(x,1,1))}
# data_pdt$cat_of_stim = sapply(data_pdt$stim, FUN=str_fun)
# data_pdt$cat_bcp     = data_pdt$cat
# data_pdt$cat         = data_pdt$cat_of_stim

## here we select only the odd or even trials
data_pdt_odd  = subset (data_pdt, mod(as.numeric(row.names(data_pdt)),2)==0)
data_pdt_even = subset (data_pdt, mod(as.numeric(row.names(data_pdt)),2)!=0)

# report on missings
if (length(missing_people) != 0) {
  mes = paste(missing_people,collapse=' ')
  warning(paste('Getting behavioral data I dropped these subs due to too many missings:\n',
                mes,'\nCan befound in "dropped_people_due_to_missings"'))
  dropped_people_due_to_missings          = missing_people
  dropped_people_due_to_missings_VPPGcode = agk.recode.c(missing_people,tnl$PhysioVP,tnl$VPPG)
}

# go home
setwd(path_cur_home)

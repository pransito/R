# add physio data
print("Adding model-free physio data")

# path
path_cur_home = getwd()

# set the path
setwd(path_dat)

# read in the extracted model-free physio data
corr_datf = dir(pattern = '^corr_.*csv')
eda_datf  = dir(pattern = '^eda_.*csv')
zygo_datf = dir(pattern = '^zygo_.*csv')

# get the summary stats used from filename
agk.second = function(x) {return(x[2])}
agk.extr.sumfuns = function(v) {
  
  # getting the current channel (corr, eda, zygo)
  tmp = suppressWarnings(regexp(v,'.*_'))
  chn = tmp$match
  
  # getting the function used
  v_funs = sapply(v,FUN=strsplit,split=chn)
  v_funs = unlist(lapply(v_funs, FUN=agk.second))
  v_funs = unlist(sapply(v_funs,FUN=strsplit,split='.csv'))
  return(as.character(v_funs)) 
}
corr_funs = agk.extr.sumfuns(corr_datf)
eda_funs  = agk.extr.sumfuns(eda_datf)
zygo_funs = agk.extr.sumfuns(zygo_datf)

# get the data from all channels and sum stats
corr_dat = list()
eda_dat  = list()
zygo_dat = list()
for (ii in 1:length(corr_datf)) {
  corr_dat[[ii]] = read.table(corr_datf[ii],sep=",", header=T)
  eda_dat[[ii]]  = read.table(eda_datf[ii],sep=",", header=T)
  zygo_dat[[ii]] = read.table(zygo_datf[ii],sep=",", header=T)
}

# merging within channels
physio_dat = corr_dat[[1]]
eda_datm   = eda_dat[[1]]
zygo_datm  = zygo_dat[[1]]
for (ii in 2:length(corr_datf)) {
  physio_dat = merge(physio_dat,corr_dat[[ii]],by.x=c("sub","stim","trial"),by.y=c("sub","stim","trial"),all = T)
  eda_datm   = merge(eda_datm,eda_dat[[ii]],by.x=c("sub","stim","trial"),by.y=c("sub","stim","trial"),all = T)
  zygo_datm  = merge(zygo_datm,zygo_dat[[ii]],by.x=c("sub","stim","trial"),by.y=c("sub","stim","trial"),all = T)
}

# give proper names
chns = c('corr','eda','zygo')
names(physio_dat)[c(4:length(physio_dat))]  = paste0('corr_',corr_funs)
names(eda_datm)[c(4:length(physio_dat))]    = paste0('eda_',corr_funs)
names(zygo_datm)[c(4:length(physio_dat))]   = paste0('zygo_',corr_funs)

# merging the channels
physio_dat    = merge(physio_dat,eda_datm,by.x=c("sub","stim","trial"),by.y=c("sub","stim","trial"),all = T)
physio_dat    = merge(physio_dat,zygo_datm,by.x=c("sub","stim","trial"),by.y=c("sub","stim","trial"),all = T)

# unit test: all names must be unique
if (any(duplicated(names(physio_dat)))) {
  stop('DUPLICATED NAMES IN PHYSIO_DAT')
}

# stripping white and stripping " (1)"
cur_stim        = strTrim(as.character(physio_dat$stim))
cur_stim        = unlist(strsplit(cur_stim,split = ".jpg"))
cur_stim        = strTrim(cur_stim)
cur_stim        = gsub(cur_stim,pattern = " \\(1\\)",replacement = "")
physio_dat$stim = as.factor(cur_stim[cur_stim!=""])
data_pdt$stim   = gsub(data_pdt$stim,pattern = " \\(1\\)",replacement = "")

# trimming (1) suffixes
physio_dat$sub   = gsub('(1)','',physio_dat$sub, fixed=T)

# trimming white space
physio_dat$sub   = trimws(physio_dat$sub)

# recoding
physio_dat$sub_new = agk.recode.c(physio_dat$sub,tnl$PhysioVP,tnl$VPPG)
physio_dat$sub_new = as.factor(physio_dat$sub_new)
physio_dat$sub     = physio_dat$sub_new
physio_dat$sub_new = NULL

# getting all the lengths before
cur_len        = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(cur_len) = c('subject','len')
subs_with_phys = cur_len

# checking if there are duplicates in physio_dat
physio_dat = physio_dat[!duplicated.data.frame(physio_dat),]

# merging with data_pdt
data_pdt         = merge(data_pdt,physio_dat,by.x=c("subject","stim","trial"),by.y=c("sub","stim","trial"),all = T)
len_after        = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(len_after) = c('subject','len')

# check if length stayed the same
# are all subs still there?
cur_tf = cur_len$subject %in% len_after$subject
if (all(cur_tf)) {
  disp('All subs are still in data_pdt that were there before physio attachment')
} else {
  cur_missing_subs = cur_len$subject[!cur_len$subject %in% len_after$subject]
  mes              = paste(paste(cur_missing_subs,collapse=' '))
  stop(paste('Subs in data_pdt were lost during physio attachment',
             mes))
}

# check if subjects were put on that weren't there before
# this means there are physio data without proper PDT behavioral data
# this must not be
cur_tf = len_after$subject %in% cur_len$subject
if (all(cur_tf)) {
  disp('All subs with physio were attached to their proper behav data.')
} else {
  cur_missing_subs = len_after$subject[!cur_tf]
  # ignoring subs that have been excluded
  cur_missing_subs = cur_missing_subs[!cur_missing_subs %in% excluded_VPPG]
  if (length(cur_missing_subs) != 0) {
    mes              = paste(paste(cur_missing_subs,collapse=' '))
    stop(paste('Some subs with physio did not have behav data to match to! Check this!',
               mes))
  }
}

# excluding those subs with physio data which do not have behav data (excluded subs)
data_pdt         = data_pdt[!data_pdt$subject %in% excluded_VPPG,]
len_after        = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(len_after) = c('subject','len')


# check the lengths before and after
if (all(cur_len$subject == len_after$subject)) {
  if (all(cur_len$len == len_after$len)) {
    disp('Number of trials was not changed due to physio attachment.')
  } else {
    stop('In some subs number of trials was changed due to physio attachment.')
  }
} else {
  stop('Subs before physio attachment and after do not allign.')
}

# add the cohort variable
data_pdt$cohort = agk.recode.c(data_pdt$subject,tnl$VPPG,tnl$Cohort)
data_pdt$cohort = as.factor(data_pdt$cohort)

# get the MATLAB report on good and existent physio
setwd(path_dat)
setwd('..')
subs_good_phys        = read.table('good_physio_data.txt',header=T)
subs_good_phys        = subs_good_phys$subs_with_good_physio_data
subs_good_phys        = agk.recode.c(subs_good_phys,tnl$PhysioVP,tnl$VPPG)

# check if it aligns with what we just imported
data_pdt_physio       = subset(data_pdt, (cohort != 'MRT' & cohort != 'MRT_Pilot'))
dpp_subs              = unique(data_pdt_physio$subject)
# ignoring excl subs
subs_good_phys = subs_good_phys[!subs_good_phys %in% excluded_VPPG]              
if (!all(subs_good_phys %in% dpp_subs)) {
  stop('Not all subs which are on the MATLAB physio report list have been attached. Check this!')
}

# check what physio data is there and good
cur_phys_failed_check = which(!dpp_subs %in% subs_good_phys)
cur_phys_failed_check = dpp_subs[cur_phys_failed_check]
tb_excl               = cur_phys_failed_check

# first check against the phys_exempt list if subjects are known non-physios
cur_phys_failed_check = cur_phys_failed_check[!cur_phys_failed_check %in% phys_exempt]
if (length(cur_phys_failed_check) != 0) {
  mes = paste(cur_phys_failed_check,collapse = ' ')
  stop(paste('There are unexpected subs which have no physio.\n',
       mes))
}

# if so desired then subs without physio data will be excluded
if (physio_excl) {
  mes = paste(tb_excl,collapse = ' ')
  disp(paste0('Dropping these subs because of missing physio data\n',mes))
  data_pdt         = data_pdt[-which(data_pdt$subject %in% tb_excl),]
  data_pdt$subject = droplevels(data_pdt$subject)
} else {
  disp('Dropping missing physio subs is switched off.')
  tb_excl = c()
}

excluded_VPPG = c(excluded_VPPG,tb_excl)

# go home
setwd(path_cur_home)

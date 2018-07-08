## message to user
print("Adding image ratings")

## TODO
# add image properties
# at path_imp

# path
path_cur_home = getwd()


## get all the necessary .txt files from the current directory
setwd(path_rat)
all_behav_files = list.files(all.files =T, pattern="*_organized.csv*")
ratings         = data.frame()

# setting up progress bar
total = length(all_behav_files)
pb    = txtProgressBar(min = 0, max = total, style = 3)

## getting all ratings into one data frame
for (ii in 1:total) {
  
  cur_sub = strsplit(all_behav_files[[ii]],"_")
  cur_sub = cur_sub[[1]][3]
  cur_dat = read.csv(all_behav_files[ii])
  
  # check if data is empty
  if (length(cur_dat[,1]) == 0) {
    warning(paste(cur_sub, "has no ratings data. I skip the sub."))
    next
  }
  
  cur_sub = rep(cur_sub,length(cur_dat[,1]))
  cur_dat$subject = cur_sub
  

  
  # identify missings
  cur_dat[cur_dat == 0] = NA
  
  # throw out SLM
  cur_dat = cur_dat[cur_dat$imageID < 90000,]
  
  # factor imageGroup
  cur_dat$imageGroup = as.factor(cur_dat$imageGroup)
  
  # save raw data
  cur_dat_raw = cur_dat
  
  # scaling
  cur_dat$arousal_bcp    = cur_dat$arousal
  cur_dat$valence_bcp    = cur_dat$valence
  cur_dat$dominance_bcp  = cur_dat$dominance
  cur_dat$arousal        = scale(cur_dat$arousal,center = 50)
  cur_dat$valence        = scale(cur_dat$valence, center = 50)
  cur_dat$dominance      = scale(cur_dat$dominance,center = 50)
  
  # first pretest subject had different questions...
  if (length(cur_dat$imageRating2abc) > 0) {
    cur_dat$imageRating2 = NA
    cur_dat$imageRating3 = NA
    cur_dat$imageRating4 = NA
    for(zz in 1:length(cur_dat[,1])) {
      if (cur_dat$imageGroup[zz] == 1) {cur_dat$imageRating2[zz] = cur_dat$imageRating2abc[zz]}
      if (cur_dat$imageGroup[zz] == 2) {cur_dat$imageRating3[zz] = cur_dat$imageRating2abc[zz]}
      if (cur_dat$imageGroup[zz] == 3) {cur_dat$imageRating4[zz] = cur_dat$imageRating2abc[zz]}
    }
  }
  
  # now the variables will be scaled together (i.e. a compund variance and mean will be used)
  # NAs stay NAs
  vars   = c("imageRating1","imageRating2","imageRating3","imageRating4","imageRating5")
  int_df = agk.varselect(cur_dat,vars)
  
  # batch scaling
  int_dfs = agk.batch.scale(int_df,fixed_mean = 0)
  # attaching on df
  cur_dat = cbind(cur_dat,int_dfs)
  
  # collect the current data
  if (ii == 1) {ratings = cur_dat}
  if (ii > 1) {
    ratings = join(ratings,cur_dat,type="full", by="subject",match="first")
  }
  
  # update progress bar
  Sys.sleep(0.05)
  setTxtProgressBar(pb, ii)
}

# close progress bar
close(pb)

# imageID merger
if (!is.null(ratings$NewImageID)) {
  for (kk in 1:length(ratings[,1])) {
    if (!is.na(ratings$NewImageID[kk])) {ratings$imageID[kk] = ratings$NewImageID[kk]}
  }
}

ratings_bck = ratings
ratings = ratings_bck[!is.na(ratings_bck$imageID),] # Milan: added this line on 20.06.16 because for some reason there are NA's here that block the script

# image ID clean; this way imageID only has the 5 digit imageIDs;
# old images which were never used again are NA
for (kk in 1:length(ratings[,1])) { 
  ratings$imageID_dirty[kk] = ratings$imageID[kk]
  if (nchar(as.character(ratings$imageID[kk])) != 5) {
    ratings$imageID[kk] = NA
  }
}

# factorize subject
ratings$subject = as.factor(ratings$subject)

# recoding to VPPG numbers
ratings$subject_new  = agk.recode.c(as.character(ratings$subject),tnl$PhysioVP,tnl$VPPG)
ratings$subject_new  = as.factor(as.character(ratings$subject_new))
ratings$subject      = ratings$subject_new
ratings$subject_new  = NULL

# prep putting ratings data onto the pdt data
colnames(ratings)[which(colnames(ratings)  == "imageID")] = "stim"
ratings$subject      = as.factor(ratings$subject)

# getting all the number of trials before
cur_len        = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(cur_len) = c('subject','len')
subs_with_phys = cur_len

# try merging
data_pdt_bcpr    = data_pdt
if (tnl_compl_excl == 1) {
  # only ratings data matching to behav data will be attached
  data_pdt         = merge(data_pdt, ratings, by=c("subject", "stim"),all.x = T,all.y=F)
  cur_incl         = tnl$VPPG[which(tnl$Einschluss_compl == 1)]
  data_pdt         = data_pdt[which(data_pdt$subject %in% cur_incl),]
} else {
  # all ratings data will be attached
  # unmatched ratings and also subjects' ratings that do not have ANY behav data
  data_pdt         = merge(data_pdt, ratings, by=c("subject", "stim"),all.x = T,all.y=T)
}
len_after        = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(len_after) = c('subject','len')

# are all subs still there?
cur_tf = cur_len$subject %in% len_after$subject
if (all(cur_tf)) {
  disp('All subs are still in data_pdt that were there before ratings attachment')
} else {
  cur_missing_subs = cur_len$subject[!cur_len$subject %in% len_after$subject]
  mes              = paste(paste(cur_missing_subs,collapse=' '))
  if (tnl_compl_excl == 0) {
    stop(paste('Subs in data_pdt were lost during ratings attachment',
               mes))
  } else {
    disp(paste('These subs in data_pdt were cut during ratings attachment because of tnl_compl_excl',
               mes))
  }
}

# check if subjects were put on that weren't there before
# this means there are ratings without proper PDT behavioral data
cur_tf = len_after$subject %in% cur_len$subject
if (all(cur_tf)) {
  disp('All subs now having ratings attached were there before in data_pdt.')
} else {
  cur_missing_subs = len_after$subject[!cur_tf]
  # check against behavioral exemption
  cur_missing_subs = cur_missing_subs[!cur_missing_subs %in% behav_exempt]
  # check against excluded subs
  if (exists('excluded_VPPG')) {
    cur_missing_subs = cur_missing_subs[!cur_missing_subs %in% excluded_VPPG]
    if (length(cur_missing_subs) != 0) {
      mes              = paste(paste(cur_missing_subs,collapse=' '))
      stop(paste('Some subs with ratings did not have behav data in data_pdt! Check this!',
                 mes))
    } else {
      disp('There are new subs now after ratings attachment, but they are accounted for.')
    }
    
  }
}

# exclude the "excluded_VPPG" subs
# Pretest subs should be unaffected since they are not excluded even if they have no behav data
if (exists('excluded_VPPG')) {
  data_pdt  = data_pdt[!data_pdt$subject %in% excluded_VPPG,]
}
len_after = aggregate(data_pdt$subject, by=list(data_pdt$subject),FUN=length)
names(len_after) = c('subject','len')
len_after = len_after[len_after$subject %in% cur_len$subject,]
if (tnl_compl_excl) {
  cur_len = cur_len[which(cur_len$subject %in% cur_incl),]
}

# TODO: in case tnl_excl_compl is turned off, then ratings attachments can lead to increase in
# trial length; so far no completely coded check if those length increases are okay or due to
# some douple import;
# check the number of trials before and after
if (all(droplevels(cur_len$subject) == droplevels(len_after$subject))) {
  if (all(cur_len$len == len_after$len)) {
    disp('Number of trials was not changed due to ratings attachment.')
  } else {
    View(data.frame(cur_len$subject,cur_len$len,len_after$len)[which(cur_len$len != len_after$len),])
    warning('In some subs number of trials was changed due to ratings attachment. Check in viewer.')
  }
} else {
  stop('Subs before ratings attachment and after do not allign.')
}

# go home
setwd(path_cur_home)


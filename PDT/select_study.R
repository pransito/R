## PREAMBLE ===================================================================
# script that selects from data_pdt
# the study (cohort) you need
# run data_import.R before
# will output data_pdt to be used for further analysis

# get the original data_pdt from data_import.R
data_pdt     = data_pdt_bcp
data_pdt_inv = data_pdt

## PARAMETER SETTINGS =========================================================

#which study to look at (Cohorts)?
#which_study = "MRT"
#which_study = "MRT_and_POSTPILOT" # lumping those together (for KFG prediction e.g.)
which_study = "POSTPILOT_HCPG" # CAREFUL: had different set of neutral pictures (?!?!)
#which_study = "TEST" # when K.Brehm used POSTPILOT and simulated facial expression (8888) or not (7777)
#which_study = "Prestudy" # HC groups before core behav study; for image adequacy (PhysioPilot)
#which_study = "sanity"
#which_study = "POSTPILOT_HC" # this for testing image adequacy within HC of core behav study
#which_study = "POSTPILOT_PG" # this for testing image adequacy within PG of core behav study
#which_study = "MRT_HC" # this for testing image adequacy within HC of core behav study
#which_study = "POSTPILOT_PGxGENDER" # this for testing the effect of gender in pg group
#which_study = "MRT_LB" # subsample MA

# default data_inv
if (which_study == 'MRT') {
  do_data_inv = 0 
} else {
  do_data_inv = 0
}

# desired data_inv
# use of cohort all subjects but that subject of the group to determine value
# imageRating1s is craving
do_data_inv      = 0
data_pdt_inv_var = 'imageRating1s'

## PREPARATIONS ===============================================================
if (physio_sum_fun == 'mean') {
  data_pdt$corr = data_pdt$corr_auc
  data_pdt$eda  = data_pdt$eda_auc
  data_pdt$zygo = data_pdt$zygo_auc
} else if (physio_sum_fun == 'max') {
  data_pdt$corr = data_pdt$corr_max
  data_pdt$eda  = data_pdt$eda_max
  data_pdt$zygo = data_pdt$zygo_max
} else if (physio_sum_fun == 'median') {
  data_pdt$corr = data_pdt$corr_median
  data_pdt$eda  = data_pdt$eda_median
  data_pdt$zygo = data_pdt$zygo_median
} else if (physio_sum_fun == 'all') {
  # do nothing
  # here all of the stuff will be used
} else {
  stop('No proper physio_sum_fun provided.')
}

# ## rt cut off (are we still using this? should not be in here!)
# # must go, if so, to data cleaning
# if(which_study == "MRT" | which_study == "MRT_HC"| which_study == "MRT_LB") {
#   rt_cut_off = 0.15
# } else {
#   rt_cut_off = 1.0
# }

## FUNCTIONS ==================================================================
# why is this here?!?!?
cur_summary_function = function(x) median(x, na.rm=TRUE)
# needs the f.difftest function from the import_data file
f = function(x) {
  tmp <- t.test(x)
  return(as.matrix(c(tmp$p.value,mean(x,na.rm = T))))
  }

## SUBSETTING DATA ============================================================
data_pdt$Cohort[data_pdt$Cohort ==  "PhysioIAPS"] = "sanity"
data_pdt$Cohort[is.na(data_pdt$Cohort)]           = "Pretest"

if (which_study == "Prestudy") {
  data_pdt = subset(data_pdt, Cohort == "Pretest" | Cohort == "PhysioPilot")
} else if (which_study == "POSTPILOT_HCPG" | which_study == "POSTPILOT_HC" | 
           which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER") {
  data_pdt = subset(data_pdt, Cohort == "POSTPILOT" | Cohort == "PGPilot")
} else if (which_study == "MRT" | which_study == "MRT_HC" | which_study == "MRT_PG"| which_study == "MRT_LB") {
  data_pdt = subset(data_pdt,Cohort == "MRT")
} else if (which_study == "PhysioPilot") {
  data_pdt = subset(data_pdt,Cohort == "PhysioPilot")
} else if (which_study == "MRT_LB") {
  includelist=read.table("E:/MATLAB/info_mri_selection.csv")
  dat_match_MRT_only=subset(dat_match, dat_match$Cohort=="MRT")
  dat_match_MRT_only=subset(dat_match_MRT_only, dat_match_MRT_only$Einschluss==1)
  data_pdt= subset(data_pdt,data_pdt$subject %in% dat_match_MRT_only$VPPG)
  data_pdt= subset(data_pdt,data_MRT_only$subject %in% includelist$V1)
} else if (which_study == "sanity") {
  data_pdt = subset(data_pdt,Cohort == "sanity")
} else if (which_study == "MRT_and_POSTPILOT") {
  data_pdt = subset(data_pdt,Cohort == "POSTPILOT" | Cohort == "PGPilot" | Cohort == "MRT")
} else if (which_study == "TEST") {
  data_pdt = subset(data_pdt, Cohort == "TEST")
} else {
  stop('No valid cohort selected with var which_study!')
}

# get a HCPG variable
data_pdt$HCPG   = as.factor(as.character(data_pdt$HCPG))

# only one group, i.e. HC or PG?
if ((length(grep(which_study, pattern = "HC")) != 0) & (length(grep(which_study, pattern = "PG")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "HC")
} else if ((length(grep(which_study, pattern = "PG")) != 0) & (length(grep(which_study, pattern = "HC")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "PG")
}

if ((length(grep(which_study, pattern = "HC")) != 0) & (length(grep(which_study, pattern = "PG")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "HC")
} else if ((length(grep(which_study, pattern = "PG")) != 0) & (length(grep(which_study, pattern = "HC")) == 0)) {
  data_pdt = subset(data_pdt,HCPG == "PG")
}

## DATA_INV ===================================================================
# prepare a data_pdt_inv df
# this df uses all data EXCEPT the data in which_study
if (do_data_inv == 1) {
  all_subs          = unique(data_pdt$subject)
  data_pdt$cat_orig = data_pdt$cat
  data_pdt$cat      = NA
  
  for (ss in 1:length(all_subs)) {
    # get everything but subject == cur_sub data
    cur_dat = subset(data_pdt, subject != all_subs[ss])
    
    # select group
    #cur_dat = subset(cur_dat, HCPG == first(data_pdt$HCPG[data_pdt$subject == all_subs[ss]]))
    
    # aggregate by stimulus
    cur_dat = cur_dat[c('subject','stim',data_pdt_inv_var)]
    cur_dat = cur_dat[-which(duplicated(cur_dat[c(1,2)])),]
    cur_dat = aggregate(cur_dat[data_pdt_inv_var],by=list(cur_dat$subject,cur_dat$stim),FUN = mean)
    cur_dat = aggregate(cur_dat[data_pdt_inv_var],by=list(cur_dat$Group.2),FUN=mean)
    
    # put it back
    cur_compl_dat = subset(data_pdt, subject == all_subs[ss])
    for (cc in 1:length(cur_compl_dat[,1])) {
      cur_compl_dat$cat[cc] = cur_dat[which(cur_dat$Group.1 %in% cur_compl_dat$stim[cc]),data_pdt_inv_var]
    }
    data_pdt[data_pdt$subject == all_subs[ss],] = cur_compl_dat
  }
  data_pdt_inv = data_pdt
} else {
  data_pdt_inv = data_pdt
}

## END DATA_INV

## AGE AND GENDER within the sample ##
# length(unique(data_pdt_finCat$subject))
# length(unique(data_pdt_finCat$subject[!is.na(data_pdt_finCat$imageRating1)]))
# age=c()
# gender=c()
# for (i in 1:length(unique(data_pdt_finCat$subject))){
#   subs=unique(data_pdt_finCat$subject)
#   sub=subs[i]
#   age[i]=data_pdt_finCat[data_pdt_finCat$subject== sub,]$age[1]
#   gender[i]=data_pdt_finCat[data_pdt_finCat$subject== sub,]$sex[1]
# }
# mean(age)
# sd(age)
# table(gender)

# CATEGORY LABELS =============================================================
if (do_data_inv == 0) {
  # Main effect of the final experimental categories: gam, pos, neg, neu_aw
  data_pdt_finCat = data_pdt
  if (which_study == "Prestudy") {
    data_pdt_finCat$cat = agk.recode.c(as.character(data_pdt_finCat$cat),c("1","2","3"),c("1","2","3"))
    data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(1,2,3),
                                 labels = c('gambling','negative', 'positive')) ## CAREFUL: I TOOK OUT ALL NEUTRAL PICTURES HERE!
  } else if (which_study == "sanity") {
    data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(6,2,3,7,8),
                                 labels = c('neutral_IAPS','negative_VPPG', 'positive_VPPG','negative_IAPS','positive_IAPS'))
  } else if (which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER" | 
             which_study == "POSTPILOT_HC" | which_study == "MRT_HC" | 
             which_study == "MRT_PG" | which_study == "POSTPILOT_HCPG" | 
             which_study == "MRT" | which_study == "MRT_and_POSTPILOT" |
             which_study == "TEST") {
    data_pdt_finCat$cat = factor(as.numeric(as.character(data_pdt_finCat$cat)),levels = c(6,1,2,3),
                                 labels = c('neutral','gambling','negative', 'positive'))
  }
}

if(sum(is.na(data_pdt$cat))) {
  stop('There are NAs in the data_pdt$cat variable!')
}

## VARIABLE TRANSFORMATIONS ===================================================
if (do_data_inv == 0) {
  data_pdt_finCat$valence_log       = get.log(data_pdt_finCat$valence)
  data_pdt_finCat$imageRating2s_log = get.log.base(data_pdt_finCat$imageRating2s,10)
  data_pdt_finCat$imageRating1s_log = get.log(data_pdt_finCat$imageRating1s)
  data_pdt_finCat$imageRating4s_log = get.log(data_pdt_finCat$imageRating4s)
  data_pdt_finCat$imageRating3s_log = get.log(data_pdt_finCat$imageRating3s)
}

## GET CAT LABELS AND TRANSFORMATIONS INTO DATA_PDT ===========================
# DO THIS BFORE STARTING ANY ANALYSIS OF BEHAVIORAL PDT TASK DATA
if (do_data_inv == 0) {
  data_pdt         = data_pdt_finCat
}
data_pdt$subject = droplevels(data_pdt$subject)

## UNIT TEST GROUP VAR ========================================================
# test of group variable
if (which_study != 'sanity' & which_study != 'TEST') {
  data_pdt = data_pdt[!(data_pdt$HCPG != "PG" & data_pdt$HCPG != "HC"),]
  if (length(levels(data_pdt$HCPG))>2) {
    stop("WOULD NEED TO DROP SUBS DUE TO NO GROUP INFO!")
  }
}

# SUBSET DAT_MATCH ============================================================
# also select dat_match
dat_match = dat_match_bcp
if (which_study == "POSTPILOT_HCPG" | which_study == "POSTPILOT_HC" | 
    which_study == "POSTPILOT_PG" | which_study == "POSTPILOT_PGxGENDER") {
  dat_match = subset(dat_match, Cohort == "POSTPILOT" | Cohort == "PGPilot")
}
# or better yet, align
dat_match = dat_match[dat_match$VPPG %in% data_pdt$subject,]

## ADD CATEGORY VARIABLES FOR laCh ============================================
data_pdt = data_pdt[!is.na(data_pdt$accept_reject),]
all_subs = unique(data_pdt$subject)
enh_dpdt = list()

for (ss in 1:length(all_subs)) {
  # get data and make model matrix
  cur_dat = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_mm  = model.matrix.lm(accept_reject ~ (gain + loss)*cat - cat,data = cur_dat)
  cur_mm  = data.frame(cur_mm)
  
  # dropping some unneeded columns
  undes = c('X.Intercept.','gain','loss')
  cur_mm  = cur_mm[-which(names(cur_mm) %in% undes)]
  
  # attaching
  cur_dat        = data.frame(cur_dat,cur_mm)
  enh_dpdt[[ss]] = cur_dat 
}

# making new data_pdt
data_pdt = enh_dpdt[[1]]
for (ss in 2:length(enh_dpdt)) {
  data_pdt = rbind(data_pdt,enh_dpdt[[ss]])
}

## SAVE TO WORKSPACE ==========================================================
# saving this result
data_pdt_bcp_study_selected  = data_pdt
dat_match_bcp_study_selected = dat_match

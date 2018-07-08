# make values for stimuli
# just takes all the variables rating and physio and computes the mean across subjects per stimulus
# uses subjects not included in the current study (i.e. need the data_pdt_inv)
# give a data frame
# data frame must have a stimulus variable and multiple contributing to value
# will do it for PG and HC split

# STUDIES
# PhysioPilot:
# first try to use our selection of images on students
# idea was to really check the physio set up in the lab
# PhysioIAPS (sanity):
# idea was to use well-known IAPS pictures (pos, neg, neu)
# together with our pos and our neg images, and new Polish neutral pics (NAPS)
# we wanted to check if our neg, pos images function well and also we wanted to check
# if our physio set up was ok to get known effects on pos, neg IAPS
# PGPilot (core study part 1)
# here we have sected the pictures for the core study: pos, neg, neuIAPS, gambling
# and started with PG subjects
# POSTPILOT (core study part 2)
# same as PGPilot; just more PG subjects and matched HC subjects
# MRT (MRI)
# same as PGPilot and POSTPILOT, but in the MRI; 

## set parameters =============================================================
# repetitions bootstrap (should be 2000, for debug can be lower)
des_R       = 100 
# which group should be analyzed and reported on?
which_group = 'HC'

## get the data ===============================================================
setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/Adlershof/Daten/PDT/pilot')
load('data_pdt_Maja.rda')

# functions
agk.boot.ci = function(x,cur_fun,lower,upper,R,min_data_points = 10) {
  # bootstrap confint of a parameter, gives a vector x
  # c(mean, ci_lo,ci_up,number_of_data_points_going_into_estimation,
  # size of CI in percent)
  
  # first check if enough data points are there
  # function that checks if at least 10 data points are there
  if (sum(!is.na(x)) < min_data_points) {return(c(NA,NA,NA,sum(!is.na(x)),NA))}
  
  tmp     = one.boot(x,FUN=cur_fun,R)
  agg_res = as.numeric(agk.mean.quantile(tmp$t,lower,upper))
  res     = c(agg_res,sum(!is.na(x)),(agg_res[3]-agg_res[2])/100)
  return(res)
}
agk.boot.ci.c = compiler::cmpfun(agk.boot.ci)

agk.load.ifnot.install <- function(package_name){
  # tries to load library; if not possible, will download and install it
  # input is package as string variable
  if(require(package_name,character.only = T,quietly = T)){
    print(paste (package_name,"is loaded correctly"))
  } else {
    print(paste("trying to install", package_name))
    install.packages(pkgs = c(package_name))
    if(require(package_name,character.only = T)){
      print(paste(package_name,"installed and loaded"))
    } else {
      stop(paste("could not install",package_name))
    }
  }
}

#agk.load.ifnot.install('rJava')
#agk.load.ifnot.install("xlsx")
agk.load.ifnot.install("boot")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install('R.matlab')
#agk.load.ifnot.install('dplyr')

# prep
data_pdt_bb2 = data_pdt
smf      = function(x) {return(median(x,na.rm = T))}
gs       = function(x) {return(sum(x,na.rm = T))}

# HC or PG or both?======================================================================================================
# which variables are of interest
pic_value_vars        = c("arousal","arousal_bcp","dominance","dominance_bcp","valence","valence_bcp","imageRating1s","imageRating1",
                          "imageRating2s","imageRating2","imageRating3s","imageRating3","imageRating4s","imageRating4","imageRating5s","imageRating5")

# get a valid HCPG variable
for (hh in 1:length(data_pdt$HCPG)) {
  
  if (any(c('HC','PG') %in% data_pdt$HCPG[hh])) {
    cur_res = data_pdt$HCPG[hh]
  } else {
    cur_res = 'HC'
  }
  data_pdt$HCPG[hh] = cur_res  
}

# do it for HC or PG? PG or everything else (== HC)
data_pdt = subset(data_pdt,HCPG == which_group)

# get rid of all subjects that just haven't rated anything at all
sub_ok = c()
all_subs = unique(data_pdt$subject)
for (ss in 1:length(all_subs)) {
  cur_dat  = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_dat  = cur_dat[pic_value_vars]
  if(all(is.na(cur_dat))) {
    sub_ok[ss] = F
  } else {
    sub_ok[ss] = T
  }
}
data_pdt = data_pdt[data_pdt$subject %in% all_subs[sub_ok],]


# get rid of double stimuli within subjects
all_subs = unique(data_pdt$subject)
cur_dat  = data_pdt[data_pdt$subject == all_subs[1],]
cur_dat  = cur_dat[!duplicated(cur_dat$stim),]
dpdt_cln = cur_dat
for (ss in 2:length(all_subs)) {
  cur_dat  = data_pdt[data_pdt$subject == all_subs[ss],]
  cur_dat  = cur_dat[!duplicated(cur_dat$stim),]
  dpdt_cln = rbind(dpdt_cln,cur_dat) 
}
data_pdt = dpdt_cln

# get rid of stimuli not being rated at all
stim_ok = c()
all_stim = unique(data_pdt$stim)
for (ss in 1:length(all_stim)) {
  cur_dat  = data_pdt[data_pdt$stim == all_stim[ss],]
  cur_dat  = cur_dat[pic_value_vars]
  if(all(is.na(cur_dat))) {
    stim_ok[ss] = F
  } else {
    stim_ok[ss] = T
  }
}
data_pdt = data_pdt[data_pdt$stim %in% all_stim[stim_ok],]

by_stima              = data_pdt$stim
cur_data              = data_pdt[pic_value_vars]
cur_agga              = aggregate(cur_data,by = list(by_stima),FUN = agk.boot.ci.c,
                                  cur_fun=smf,lower=0.025,upper=0.975,R=des_R,min_data_points = 9)

# find out for each stimulus in which studies it was used
all_stim   = unique(data_pdt$stim)
study_list = list() 
for (ss in 1:length(all_stim)) {
  cur_dat          = data_pdt[data_pdt$stim == all_stim[ss],]
  cur_stu          = unique(cur_dat$Cohort)
  study_list[[ss]] = cur_stu
}
study_list = lapply(study_list, FUN=agk.recode.c,y=c("PhysioIAPS",'PGPilot','POSTPILOT','PhysioPilot','MRT'),
                    z=c('PI','PP','PO','PH','MR'))
collapse_fun = function(x) {return(paste(x,collapse='_'))}
study_list = lapply(study_list,collapse_fun)
study_list = unlist(study_list)

# give proper names to variables
cmb_cur_agga          = as.data.frame(cur_agga$Group.1)
for (vv in 2:length(cur_agga)) {
  cur_mat = cur_agga[[vv]]
  cur_nm  = names(cur_agga)[vv]
  cur_nms = paste0(cur_nm,c('_mn_oa','_lower_oa','_upper_oa','_N','_size_of_CI_in_percent'))
  cur_mat = as.data.frame(cur_mat)
  names(cur_mat) = cur_nms
  cmb_cur_agga = cbind(cmb_cur_agga,cur_mat)
}
#cur_agga[-1]             = apply(MARGIN = 2,X = cur_agga[-1],FUN = std_fun) # TODO: descide on across stim standardization
cur_agga_oa               = cmb_cur_agga
names(cur_agga_oa)[1]     = 'stim'
cur_agga_oa$which_studies = study_list

#delete the unnecessary columns (delete unstandardized mean,upper,lower,N/leave standardized mean,upper,lower,N and unstardized size of CI in percent)
cur_agga_oa             = subset(cur_agga_oa, select = -c(arousal_size_of_CI_in_percent, arousal_bcp_mn_oa,arousal_bcp_upper_oa,arousal_bcp_lower_oa, 
                                                           arousal_bcp_N,dominance_size_of_CI_in_percent,dominance_bcp_mn_oa,dominance_bcp_upper_oa,
                                                           dominance_bcp_lower_oa, dominance_bcp_N,valence_size_of_CI_in_percent,valence_bcp_mn_oa,
                                                           valence_bcp_upper_oa,valence_bcp_lower_oa,valence_bcp_N,imageRating1s_size_of_CI_in_percent,imageRating1_mn_oa,
                                                           imageRating1_upper_oa,imageRating1_lower_oa,imageRating1_N,imageRating2s_size_of_CI_in_percent,imageRating2_mn_oa,
                                                           imageRating2_upper_oa,imageRating2_lower_oa,imageRating2_N,imageRating3s_size_of_CI_in_percent,imageRating3_mn_oa,
                                                           imageRating3_upper_oa,imageRating3_lower_oa,imageRating3_N,imageRating4s_size_of_CI_in_percent,imageRating4_mn_oa,
                                                           imageRating4_upper_oa,imageRating4_lower_oa,imageRating4_N,imageRating5s_size_of_CI_in_percent,imageRating5_mn_oa,
                                                           imageRating5_upper_oa,imageRating5_lower_oa,imageRating5_N))

# imageProperties===============================================================
#Import ImageProperties(all pictures but slot machine pics and neutral NAPS)
ImageProperties <- read.delim('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/Adlershof/Daten/PDT/pilot/Bilderrating_Maja/ImagePropertie_without5and9.dat')
View(ImageProperties)
colnames(ImageProperties)[colnames(ImageProperties)=='Image'] = 'stim'

#delete .jpg of stim_variable in ImageProperties_dataframe
ImageProperties$stim = trimws(gsub(".jpg","",ImageProperties$stim, fixed=T))

#merge tables imageproperties and ratings
ImageProperties_merge  = merge(cur_agga_oa, ImageProperties,by='stim') 
View(ImageProperties_merge)

# print the result to excel table ==============================================
# ... write.table
# delim (sep = ';') filename = name.csv

R.Version()
#would be an option, because doesnt need java, but requires Rtools ?
agk.load.ifnot.install("openxlsx")
library(openxlsx)
write.xlsx(ImageProperties_merge, 'ImageProperties_merge_HC.xlsx')


# I used this
write.csv2(ImageProperties_merge, file="S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/Adlershof/Daten/PDT/pilot/Bilderrating_Maja/tablescorrect/ImageProperties_merge_PG.csv", row.names=F )



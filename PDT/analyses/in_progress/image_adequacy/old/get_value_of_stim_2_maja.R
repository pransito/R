# make values for stimuli
# just takes all the variables rating and physio and computes the mean across subjects per stimulus
# uses subjects not included in the current study (i.e. need the data_pdt_inv)
# give a data frame
# data frame must have a stimulus variable and multiple contributing to value
# will do it for PG and HC split

# get the data
setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/Adlershof/Daten/PDT/pilot')
load('data_pdt_Maja.rda')

# functions
agk.boot.ci = function(x,cur_fun,lower,upper,R,min_data_points = 10) {
  # bootstrap confint of a parameter, given a vector x
  # c(mean, ci_lo,ci_up,number_of_data_points_going_into_estimation)
  
  # first check if enough data points are there
  # function that checks if at least 10 data points are there
  if (sum(!is.na(x)) < min_data_points) {return(c(NA,NA,NA,sum(!is.na(x))))}
  
  tmp = one.boot(x,FUN=cur_fun,R)
  res = c(as.numeric(agk.mean.quantile(tmp$t,lower,upper)),sum(!is.na(x)))
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



agk.load.ifnot.install('rJava')
agk.load.ifnot.install("xlsx")
agk.load.ifnot.install("boot")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install('R.matlab')

# prep
data_pdt_bb2 = data_pdt
smf      = function(x) {return(median(x,na.rm = T))}
gs       = function(x) {return(sum(x,na.rm = T))}



# HC or PG or both?======================================================================================================
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
                     "imageRating2s","imageRating3s","imageRating4s")
  
  by_stim               = data_pdt$stim
  cur_dat               = data_pdt[pic_value_vars]
  cur_agg               = aggregate(cur_dat,by = list(by_stim),FUN = agk.boot.ci.c,cur_fun=smf,lower=0.025,upper=0.975,R=2000)
  cmb_cur_agg           = as.data.frame(cur_agg$Group.1)
  for (vv in 2:length(cur_agg)) {
    cur_mat = cur_agg[[vv]]
    cur_nm  = names(cur_agg)[vv]
    cur_nms = paste0(cur_nm,c('_mn_HC','_lower_HC','_upper_HC'))
    cur_mat = as.data.frame(cur_mat)
    names(cur_mat) = cur_nms
    cmb_cur_agg = cbind(cmb_cur_agg,cur_mat)
  }
  #cur_agg[-1]           = apply(MARGIN = 2,X = cur_agg[-1],FUN = std_fun) # TODO: descide on across stim standardization
  cur_aggHC             = cmb_cur_agg
  names(cur_aggHC)[1]   = 'stim'
}

if (doPG == 1) {
  # case PG
  data_pdt       = subset(data_pdt_inv,HCPG=="PG")
  pic_value_vars = c("arousal","dominance","valence","imageRating1s",
                     "imageRating2s","imageRating3s","imageRating4s","imageRating5s")
  
  by_stim               = data_pdt$stim
  cur_dat               = data_pdt[pic_value_vars]
  cur_agg               = aggregate(cur_dat,by = list(by_stim),FUN = agk.boot.ci.c,cur_fun=smf,lower=0.025,upper=0.975,R=2000)
  cmb_cur_agg           = as.data.frame(cur_agg$Group.1)
  for (vv in 2:length(cur_agg)) {
    cur_mat = cur_agg[[vv]]
    cur_nm  = names(cur_agg)[vv]
    cur_nms = paste0(cur_nm,c('_mn_PG','_lower_PG','_upper_PG'))
    cur_mat = as.data.frame(cur_mat)
    names(cur_mat) = cur_nms
    cmb_cur_agg = cbind(cmb_cur_agg,cur_mat)
  }
  #cur_agg[-1]           = apply(MARGIN = 2,X = cur_agg[-1],FUN = std_fun) # TODO: descide on across stim standardization
  cur_aggPG             = cmb_cur_agg
  names(cur_aggPG)[1]   = 'stim'
}


# case PG and HC together
data_pdt              = data_pdt_inv

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

pic_value_vars        = c("arousal","dominance","valence","imageRating1s",
                         "imageRating2s","imageRating3s","imageRating4s","imageRating5s")
by_stima              = data_pdt$stim
cur_data              = data_pdt[pic_value_vars]
cur_agga              = aggregate(cur_data,by = list(by_stima),FUN = agk.boot.ci.c,
                                  cur_fun=smf,lower=0.025,upper=0.975,R=100,min_data_points = 9)
cmb_cur_agga          = as.data.frame(cur_agga$Group.1)
for (vv in 2:length(cur_agga)) {
  cur_mat = cur_agga[[vv]]
  cur_nm  = names(cur_agga)[vv]
  cur_nms = paste0(cur_nm,c('_mn_oa','_lower_oa','_upper_oa'))
  cur_mat = as.data.frame(cur_mat)
  names(cur_mat) = cur_nms
  cmb_cur_agga = cbind(cmb_cur_agga,cur_mat)
}
#cur_agga[-1]           = apply(MARGIN = 2,X = cur_agga[-1],FUN = std_fun) # TODO: descide on across stim standardization
cur_agga_oa           = cmb_cur_agga
names(cur_agga_oa)[1]   = 'stim'

# merge PG and HC tables===============================================================================================
cur_agg_HC_PG = merge(cur_aggHC, cur_aggPG, by=c('stim'))
#merge PG_HC + together
cur_agg_HC_PG_oa = merge(cur_agg_HC_PG,cur_agga_oa,by = c('stim'))

# make a category
cur_agg_HC_PG_oa$cat = floor(as.numeric(cur_agg_HC_PG_oa$stim)/10000)
cur_agg_HC_PG_oa$cat = factor(cur_agg_HC_PG_oa$cat,levels=c(6,1,2,3),labels=c('neu','gam','neg','pos'))

# modeling
cur_summary = lm(v_oa_arousal ~ 0+cat, data=cur_agg_HC_PG_oa)

#imageProperties=======================================================================================================
#Import ImageProperties by hand
ImageProperties <- read.delim("T:/R_Aufgaben/PDT_ALEX/ImageProperties.dat")
View(ImageProperties)
colnames(ImageProperties)[colnames(ImageProperties)=='Image'] = 'stim'

#delete .jpg of stim_variable in ImageProperties_dataframe
ImageProperties$stim = gsub(".jpg","",ImageProperties$stim, fixed=T)

#merge tables imageproperties and ratings
ImageProperties_merge  = cbind(cur_agg_HC_PG_oa, ImageProperties) 
View(ImageProperties_merge)

#Prozente für das Konfidenzintervall

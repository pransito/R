## PREAMBLE ===================================================================
# prep behavioral analyses PDT
# 19.09.2018
# Alexander Genauck

# IMPORT OF PDT non-fMRI DATA
# PRETEST,PILOT-HC,SANITY,PILOT-PG,MRT-PG,MRT-HC

## clear workspace
rm(list = ls())

warning('VPPG0115 still has two P structs. Behav data now only from first. Adapt import behav and ss models MRI')

## PARAMETERS =================================================================
## parameters that may be set
# use last exisiting import
import_existing_imp = 1
# import from scratch (choice data, ratings, etc.; takes a bit)
# if 0 will take an older saved version
import_from_scratch = 0
# scaling: currently centered-only (gain, loss)
use_z               = 1 
# accept rject is not binary but 1 through 4 (metric)
acc_num             = 0
# use absolute value of loss as predictor
loss_abs            = 1
# if set to 1 then people with too many missings will be dropped
missing_check       = 1  
# if there are more than x% percent missings in response then drop-out
missing_cutoff      = 0.8
# exclude if physio is missing?
# CAREFUL: ONLY SO FAR WORKS IF DATA IS IMPORTED FROM SCRATCH
# TODO!!!
physio_excl         = 1
# which physio sum stat data to use (mean or max or median or...)
# 'all' if you wanna put all summary stas in the CV (preferrable)
physio_sum_fun      = 'all'
# exclude according to TN-Liste? (should be always 1)
# this is to exclude according to include variable
# ADD Pretest cohort to TN-Liste, or get it from second sheet
# excluding before adding rating; rating can bring in subs that were excluded
# before or are not even in TN list
tnl_excl            = 1
# exclude after adding ratings, (...)  to only have subs which have complete data
# should be 1 if you need complete data, should be 0 if you are interested in ratings
# and other data and not in all the complete behav and others data;
# if 0 then check get_data_ratings_2.R for the TODO at the length check
tnl_compl_excl      = 1
# what kind of pca to extract features? (if no kernel pca then normal pca)
do_kernelpca        = 0
# write the matching tables (should always be 1)
write_match_tables  = 1
# should bootstrapping be used instead of normal t-test when checking for mathcing?
# 1: yes; 2: TODO! permutation test will be used
match_boot          = 2
# get the MRI behav data as well? (only makes sense if the data is accessible)
get_MRI_behav_data  = 1
# estimate LA classical for MRI ss analysis? (using lmlist; quick!)
get_LAcl_for_MRI    = 0
# estimate LA Charpentier for MRI ss analysis? (CAREFUL: TAKES A WHILE)
get_LAch_for_MRI    = 0
# get MRI BOLD extracts (to be done in matlab and moved to google drive)
get_MRI_extr        = 0
# how much aggregate (3 is from 12by12 to 4by4 because 12/3 == 4)
cur_agg             = 3
# matching criterion p-value (for matching experimental groups)
m_crit              = 0.15
# KFG cut off equal and above is PG (16 or 25)
KFG_cutoff          = 16
# for matching: which studies to do it on, and what are the group sizes?
# set desired_n to too high value if you do not want to do the matching for a particular study
which_studies       = c("MRT","POSTPILOT")
desired_n           = list(c(32,32),c(30,30))
# for matching (dom: do matching variables, do matching variables narrowed for elimination of couples)
cur_names_dom          = c('edu_years','edu_years_voca','edu_hollingshead','income_personal','smoking_ftdt','Age','audit','dem_gender','handedness','unemployed')
cur_names_dom_narrowed = c('Age')

## EXCLUSION/EXEMPTION LISTS ==================================================
# subjects that are exempt of physio; have been checked; do not have physio
# and physio is not recoverable, these subs will still be excluded if phys_excl == 1
# however, if a sub has no physio and is not on this list, an error will
# be thrown
phys_exempt  = c("VPPG0401","VPPG0705","VPPG0842",'VPPG0257a','VPPG0036','VPPG0326')
behav_exempt = paste0('Pretest',sprintf("%02d",seq(1,16)))
if (physio_excl) {
  behav_exempt = c(behav_exempt,phys_exempt)
}
# TODO: adding some more to behav_exempt; VPPG0289a is MRT and his data from S: needs to be arranged!
# TODO: VPPG0896 is a new MRI sub; needs to be arranged, behav data exported; etc.
# NOTE: subs without behav data will be disregarded in all further analyses/reports/matching
behav_exempt = c(behav_exempt,'baseNAPS0008', 'baseNAPS0013', 'baseNAPS0014',
                 'baseNAPS0015', 'baseNAPS0016', 'baseNAPS09','baseNAPS99999')

# excl subs by hand (should be always empty!!!)
# only for quick and dirty emergency
subs_excl_hand  = c()

## PATHS ======================================================================
## paths (set according to your system)
user              = paste0(as.character(Sys.info()["login"]),"/")
base_dat          = 'S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/Adlershof/Daten/PDT/'
#base_dat         = 'E:/Google Drive/Promotion/VPPG/VPPG_Exchange/Experimente/PDT/Daten/'
base_dat_GD       = paste0('C:/Users/',user,'Google Drive/Promotion/VPPG/VPPG_Exchange/Experimente/PDT/Daten/')
base              = paste("C:/Users/",user,"Google Drive/Promotion/VPPG/VPPG_Exchange/",sep="")
base_bgg          = paste("C:/Users/",user,"Google Drive/Promotion/VPPG/VPPG_Exchange/",sep="")
base_lib          = paste("C:/Users/",user,"Google Drive/",sep="")

# other working locations
#base              = "E:/Google Drive/Promotion/VPPG/VPPG_Exchange/"
#base_lib          = "E:/Google Drive/"
#base_dat_GD       = 'E:/Google Drive/Promotion/VPPG/VPPG_Exchange/Experimente/PDT/Daten/'


# some other absolute paths to be set with brute force
# path to the matching of names and subject codes
path_mtk          = 'S:/AG/AG-Spielsucht2/Daten/Probanden'
path_bgg          = 'C:/Users/genaucka/Google Drive/Diplom/LA/daten_behav_test_finale_SP_Diplom/Results'

# paths that build on above base paths; DO NOT CHANGE
path_ana          = paste0(base_lib,"Library/R/PDT/analyses/")
path_anp          = paste0(path_ana,"in_progress/sev_pred")
path_dat          = paste0(base_dat,"pilot/")
path_dat_GD       = paste0(base_dat_GD,"pilot/")
path_postpilot_pg = paste0(base_dat,"POSTPILOT/PG/")
path_postpilot_hc = paste0(base_dat,"POSTPILOT/HC/")
path_pg           = paste0(base_dat,"PG/")
path_res          = paste0(base,"Experimente/PDT/analysis/results")
path_mtc          = paste0(base,"BCAN/Probandenlisten/matching")
path_rat          = paste0(base,"Bilderrating/Results_Pretest/Result files/")
path_rrs          = paste0(base,"Bilderrating/analysis/results")
path_que          = paste0(base,"Bilderrating/Results_Pretest/Result files/questionnaires - organize in R")
path_que_pp       = paste0(base,"Bilderrating/Results_Pretest/Result files/import old physio pretest questions")
#path_scr          = paste0(base,'Screening/Screening_Export')
path_scr          = path_que
path_mrt          = paste0(base_dat,'/MRT/')
path_led          = paste0(base,"Experimente/PDT/ledaLab_anal/")
path_lib          = paste0(base_lib,"Library/R")
path_plb          = paste0(base_lib,"Library/R/PDT/analyses/in_progress")
path_mod          = paste0(base_lib,"Library/R/PDT/analyses/in_progress/model_calculation")
path_imp          = paste0(base,'Bilderrating/Bildmaterial/VPPG_stim_04_reresized')
path_mes          = 'C:/Users/Alexander/Google Drive/Library/MATLAB/PDT/MRI/sl/ROIs/from_ext_HD/ROIs/ROIs_ss_model'
path_mep          = 'C:/Users/Alexander/Google Drive/Library/MATLAB/PDT/MRI/sl/ROIs/from_ext_HD/ROIs/ROIs_gPPI_targets'

# if working at home
#warning('path_mep and path_mes are set to working at home')
path_mes          = path_dat_GD
path_mep          = path_dat_GD

## LIBRARIES AND FUNCTIONS ====================================================
setwd(path_lib)
source ('agk_library.R')
setwd(path_ana)
source('perform_matching_tests.R')
#source('get_fixef_functions.R')

## which function to use for matching tests?
if (match_boot==1) {
  warning("Bootstrapping is used for matching. So f.summary is not used and instead determined by f.difftest (case 'boot'); Default there is 'mean'")
  warning("Bootstrapping is used for matching. So p(boot) is enough to determine if matched or not. Non-param. k ignored.")
  warning("Means deterimned and reported by imputing group means for missings.")
} else if (match_boot==2) {
  warning("Permutation lm is used for matching. So f.summary is not used and instead determined by f.difftest (case 'permute'); hence mean is used")
  warning("Permutation lm is used for matching. So p(perm) is enough to determine if matched or not. Non-param. k ignored.")
  warning("Means deterimned and reported by imputing group means for missings.")
} else if (match_boot==0) {
  warning('Simple t.tests used for matching plus non-parametric tests. NOT RECOMMENDED.')
} else {
  stop('match_boot is set to an uninterpretable value.')
}

## LOAD PARTICIPANT'S LIST
## participant's list (to get the VPPG numbers)
tnl           = read_excel(paste0(base, 'BCAN/Probandenlisten/Teilnehmerliste_VPPG_ONLY_USE_THIS_ONE.xlsx'),
                           sheet = 1,col_names = T)
tnl$VPPG      = as.character(tnl$VPPG)
tnl           = tnl[!is.na(tnl$VPPG),]
tnl$VPPG      = trimws(tnl$VPPG)
tnl$PhysioVP  = paste0("PhysioVP",tnl$PhysioVP)
tnl$PhysioVP  = trimws(tnl$PhysioVP)

# prelimn check that there are no duplicates in tnl
if (sum(duplicated(tnl$VPPG))) {
  mes_1 = 'You have these duplicate VPPG numbers in Teilnehmerliste. Fix first!\n'
  mes_2 = paste(tnl$VPPG[duplicated(tnl$VPPG)],collapse= ' ')
  stop(paste(mes_1,mes_2))
}

# prelimn check that inlcude/exclude variable is ok
if (any(is.na(tnl$Einschluss))) {
  stop('You have NAs in tnl$Einschluss. Fix first!')
}
if (any(!(tnl$Einschluss == 1 || tnl$Einschluss == 0))) {
  stop('Other values than 0 and 1 used in tnl$Einschluss. Fix first!')
}

# check who is eligible for POSTPILOT
subs_eligible_POSTPILOT = tnl$VPPG[(tnl$Cohort == 'POSTPILOT' | tnl$Cohort == 'PGPilot') & tnl$Einschluss == 1]

if (import_existing_imp == 0) {
  ## GETTING DATA (TAKES A WHILE)  ============================================
  if (import_from_scratch == 1) {
    ## GETTING DATA FROM SCRATCH ==============================================
    # get all the data in long format
    setwd(path_ana)
    ## IS IT OKAY THAT ALL SUBS HAVE SLIGHTLY DIFFERENT GAIN LOSS AGG STEP?
    source("get_data_pdt.R")
    source("get_physio_aggregates.R")
    
    # INTERIM
    #data_pdt        = (data_pdt[,grep('.x',names(data_pdt),fixed = T,invert=T)])
    #names(data_pdt) = gsub('.y','',names(data_pdt),fixed=T)
    #data_pdt_bcp = data_pdt
    
    source("get_data_ratings_2.R")
    
    # INTERIM
    data_pdt        = (data_pdt[,grep('.x',names(data_pdt),fixed = T,invert=T)])
    names(data_pdt) = gsub('.y','',names(data_pdt),fixed=T)
    data_pdt_bcp = data_pdt
    
    # save the result of import from scratch
    setwd(path_dat)
    save(file = "pure_data_pdt.RData",list = c("data_pdt","ratings"))
  } else {
    setwd(path_dat)
    load("pure_data_pdt.RData")
  }
  ## END OF DATA_PDT GATHERING
  
  ## DATA CLEANING ============================================================
  # dropping subjects due to always accept or reject
  disp("checking if have to exclude due to always accepted or rejected")
  
  ## subjects that ALWAYS ACCEPT or REJECT
  # ignores subs that have only NAs
  often_acc = c()
  often_rej = c()
  all_subs  = as.character(unique(data_pdt$subject))
  for (kk in 1:length(all_subs)) {
    cur_dat = subset(data_pdt, subject == all_subs[kk])
    cur_dat = subset(cur_dat, !is.na(accept_reject))
    if (length(cur_dat[,1]) == 0) {
      next
    }
    cur_acc = sum(cur_dat$accept_reject==1)/length(cur_dat$accept_reject)
    cur_rej = sum(cur_dat$accept_reject==0)/length(cur_dat$accept_reject)
    
    if (cur_acc == 1) {
      often_acc = c(often_acc,all_subs[kk])
    }
    if (cur_rej == 1) {
      often_rej = c(often_rej,all_subs[kk])
    }
  }
  
  ## throw always reject out
  data_pdt      = data_pdt[(!data_pdt$subject %in% often_rej),]
  data_pdt      = data_pdt[(!data_pdt$subject %in% often_acc),]
  if (length(often_rej) != 0 | length(often_acc) != 0) {
    warning(paste("Excluded subject(s) due to always accepted/rejected.",
                  paste0(paste(often_rej,collapse=''),paste(often_acc,collapse=''))))
  }
  
  # create workspace backup
  data_pdt_pure = data_pdt
  
  # end data cleaning and excluding subjects
  # after here YOU ARE NOT ALLOWED ANYMORE TO EXCLUDE
  # SUBS (?!)
  
  ### END DATA CLEANING
  
  ## LEDA IMPORT ==============================================================
  ## TODO: LEDA IMPORT!
  ## ADD LEDALAB
  # data_leda      = read.csv(paste0(path_led,'Data/EDA_ResultsR.csv'))
  # data_leda_cont = read.csv(paste0(path_led,'Data_control/EDA_ResultsR.csv'))
  # data_leda_pp   = read.csv(paste0(path_led,'Postpilot/Data/EDA_ResultsR.csv'))
  # data_leda      = rbind(data_leda,data_leda_cont,data_leda_pp)
  # data_leda$subject_new = agk.recode.c(as.character(data_leda$subject),tnl$PhysioVP,tnl$VPPG)
  # data_leda$subject_new = as.factor(as.character(data_leda$subject_new))
  # data_leda$subject     = data_leda$subject_new
  # data_leda$subject_new = NULL
  
  # data_leda_long      = read.csv(paste0(path_led,'Data_onsetFixStimGam/EDA_ResultsR.csv'))
  # 
  # data_leda_fix = data_leda_long[data_leda_long$onset==0,]
  # colnames(data_leda_fix)[4:10] = paste0(colnames(data_leda_long), '_fix')[4:10]
  # 
  # data_leda_stim = data_leda_long[data_leda_long$onset==1,]
  # colnames(data_leda_stim)[4:10] = paste0(colnames(data_leda_long), '_stim')[4:10]
  # 
  # data_leda_gam = data_leda_long[data_leda_long$onset==2,]
  # colnames(data_leda_gam)[4:10] = paste0(colnames(data_leda_long), '_gam')[4:10]
  # 
  # data_leda_m1=merge(data_leda_fix, data_leda_stim,  by=c('subject','trial'))
  # data_leda=merge(data_leda_m1, data_leda_gam, by=c('subject','trial'))
  # 
  # data_leda$subject_new = agk.recode.c(as.character(data_leda$subject),tnl$PhysioVP,tnl$VPPG)
  # data_leda$subject_new = as.factor(as.character(data_leda$subject_new))
  # data_leda$subject     = data_leda$subject_new
  # data_leda$subject_new = NULL
  # 
  # data_leda$onset=NULL
  # data_leda$onset.x=NULL
  # data_leda$onset.y=NULL
  # 
  # data_pdt_leda  = merge(data_pdt, data_leda, by=c('subject','trial')) 
  # data_pdt       = merge(data_pdt, data_leda, by=c('subject','trial'), all=TRUE)
  
  ## QUESTIONNAIRE DATA =======================================================
  load(paste0(path_que, '/Bilderrating_VPPG_Quest.Rda'))
  
  # Fagerström neu zusammenrechnen #
  # TODO: move to quest import
  FTND_set=data_quest[c("FTND_1","FTND_2","FTND_3","FTND_4","FTND_5","FTND_6")]
  # get number of levels
  num_levels = c()
  for (ff in 1:length(FTND_set)) {num_levels[ff] = (length(levels(FTND_set[[ff]]))-1)}
  FTND_num=(sapply(FTND_set, as.numeric))*-1
  FTND_num[FTND_set=="[NA] nicht beantwortet"] = NA
  for (ff in 1:length(FTND_num[1,])) {FTND_num[,ff] = FTND_num[,ff] + num_levels[ff]}
  FTND_num[is.na(FTND_num)]=0
  data_quest$FTND=rowSums(FTND_num)
  
  # getting the questionnaires var names
  questionnaires_vars_names = c()
  for(ii in 1:ncol(data_quest)) {
    cur_com = comment(eval(parse(text=paste0("data_quest$", attr(data_quest[ii],'name')))))
    if (length(cur_com) == 0) {
      questionnaires_vars_names[ii] = NA
      next 
    } else {
      questionnaires_vars_names[ii] =comment(eval(parse(text=paste0("data_quest$", attr(data_quest[ii],'name')))))
    }
  }
  
  # get GBQ split by subscales
  names_GBQ            <- grep("^GBQ_", names(data_quest), value=TRUE)
  GBQ_persi_ids        = c(4, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21)
  GBQ_illus_ids        = c(1, 2, 3, 5, 7, 8, 9, 19)
  tmp                  <- data_quest[names_GBQ]
  data_quest$GBQ_illus <- apply(tmp[,c(1,2,3,5,7,8,9,19)],FUN = mean.rmna,MARGIN = 1)
  data_quest$GBQ_persi <- apply(tmp[,c(4,6,10,11,12,13,14,15,16,17,18,20,21)],FUN = mean.rmna,MARGIN = 1)
  
  ## QUESTIONNAIRE DATA from the pretests and physio pretests
  data_quest_pp                  = readRDS(paste0(path_que_pp, '/data_quest_physplt.rds'))
  data_quest_pp$data_par.P104_01 = data_quest_pp$PhysioVP
  data_quest_pp                  =  data_quest_pp[rownames(data_quest_pp)!='508',]
  
  data_quest_bckp=data_quest
  data_quest= rbind.fill(data_quest, data_quest_pp)
  
  #  make a VPPG variable
  data_quest$VPPG = data_quest$data_par.P104_01
  data_quest$VPPG = agk.recode.c(data_quest$VPPG,tnl$PhysioVP,tnl$VPPG)
  
  
  ## SCREENING DATA ===========================================================
  setwd(path_scr)
  load("Screening_VPPG_Data-label.Rda")
  
  # getting the screening var names
  screening_vars_names = c()
  for(ii in 1:ncol(data_ano)) {
    cur_com = (comment(eval(parse(text=paste0("data_ano$", attr(data_ano[ii],'name'))))))
    if (length(cur_com) == 0) {
      screening_vars_names[ii] = NA
      next 
    } else {
      screening_vars_names[ii] =(comment(eval(parse(text=paste0("data_ano$", attr(data_ano[ii],'name'))))))
    }
  }
  
  # ... and levels
  screening_vars_levels = list()
  for(ii in 1:ncol(data_ano)) {
    cur_com = attr(data_ano[[ii]],'levels')
    if (length(cur_com) == 0) {
      screening_vars_levels[[ii]] = NA
      next 
    } else {
      screening_vars_levels[[ii]] =cur_com
    }
  }
  
  # categorize subjects into HC and PG based on KFG score
  data_ano$HCPG   = factor(as.numeric(data_ano$KFG>=KFG_cutoff), levels=c(NA, 0,1), labels=c('HC','PG'))
  data_ano$Cohort = agk.recode.c(as.character(data_ano$VPPG),as.character(tnl$VPPG),as.character(tnl$Cohort))
  
  # recode VPPG0051 to VPPG0051a
  # TODO: should be done in screening
  # wrote to Jan Philipp Albrecht
  data_ano$VPPG[data_ano$VPPG == "VPPG_0051"] = "VPPG_0051a"
  
  # cleaning
  # TODO: needs to be done in screening export
  data_ano$VPPG = gsub('_', '', data_ano$VPPG)
  
  ## GET MATCHING TABLES  =====================================================
  # get other variables of interest (for matching: years of edu, income) from screening
  dat_match        = data_ano[c('VPPG','S209','HCPG', 'S201','KFG', 'S313','S314_01', 'S314_02','S320','S321_01', 'S322', 'S323')]
  names(dat_match) = c("VPPG", "gender",'HCPG','handedness',"KFG", "edu_highest_deg","edu_tert_years","edu_tert_months","income_household","income_for_how_many",
                       "debt_overall","debt_gambling")
  if(!all(duplicated(dat_match$VPPG) == FALSE)) {
    stop('duplicates in dat_match$VPPG! after just reading data_ano')
  }
  
  #dat_match$debt_gambling[is.na(dat_match$debt_gambling)] = 0
  dat_match$debt_gambling = factor(dat_match$debt_gambling,levels = c(0:6),
                                   labels = c("0","10000","25000","50000","100000","200000","300000"))
  dat_match$edu_tert_months     = as.numeric(dat_match$edu_tert_months) 
  dat_match$edu_tert_years      = as.numeric(dat_match$edu_tert_years)
  names_match_quest  = names(data_quest[,c("dem_age","dem_edu_bildungsgrad","dem_occupation" ,"dem_houseIncome","dem_housePpl","dem_smoke","AUDIT","FTND")])
  rnames_match_quest = c("Age", "edu_deg_secondary", "occup","income_hh_quest", "income_people_hh_quest","smoking","audit","ftdt") ## CAREFUL: from 23.05.2016 audit is in screening!
  
  if(!all(duplicated(dat_match$VPPG) == FALSE)) {
    stop('duplicates in dat_match$VPPG!')
  }
  
  ## QUESTIONNAIRE DATA WITH MATCH DATA
  dat_match       = merge(dat_match,tnl, by=c("VPPG"))
  dat_match       = merge(dat_match,data_quest,by=c("VPPG"))
  names(dat_match)[which(names(dat_match) %in% names_match_quest)] = rnames_match_quest
  tmp = data.frame(dat_match[grep(pattern="^edu",names(dat_match))])
  
  # add KFG to quest
  #add a line in data_quest for everyone in the teilnehmerliste who also has an entry in data_ano
  # TODO: DUNNO WHAT THIS IS!!!
  # screened_tested_noquest=data_ano$VPPG[data_ano$VPPG %in% tnl$VPPG[!tnl$VPPG %in% data_quest$VPPG]]
  # data_quest[(length(data_quest$VPPG)+1):(length(data_quest$VPPG)+length(screened_tested_noquest)),]$VPPG=screened_tested_noquest
  # DUNNO WHAT THIS IS!!!
  data_quest = merge(data_quest,as.data.frame(data_ano[c("VPPG","KFG")]),by=c("VPPG"),all.x = TRUE)
  
  # Resolve AUDIT: take from screening, if NA use the quest
  # TODO: AUDIT: items: some are numeric some are labels
  data_ano_audit=data_ano[c("VPPG","E101","E102","E103","E104","E105","E106","E107",
                            "E108","E109","E110", "AUDIT")]
  names(data_ano_audit)=c("VPPG",'AUDIT_1','AUDIT_2','AUDIT_3','AUDIT_4','AUDIT_5','AUDIT_6','AUDIT_7','AUDIT_8','AUDIT_9','AUDIT_10','AUDIT')
  data_quest_bckp2=data_quest
  data_quest=merge(data_quest, data_ano_audit, by = 'VPPG', suffixes = c(".que",".scr"))
  data_quest$AUDIT=data_quest$AUDIT.scr
  data_quest$AUDIT[is.na(data_quest$AUDIT.scr)]=data_quest$AUDIT.que[is.na(data_quest$AUDIT.scr)]
  for (i in 1:10) {
    eval(parse(text=paste0('data_quest$AUDIT_',i, '=data_quest$AUDIT_', i, '.scr')))
    eval(parse(text=paste0('data_quest$AUDIT_',i, '=as.character(data_quest$AUDIT_', i,')')))
    eval(parse(text=paste0('data_quest$AUDIT_',i, '[is.na(data_quest$AUDIT_', i,
                           '.scr)]=data_quest$AUDIT_', i, '.que[is.na(data_quest$AUDIT_', i, '.scr)]')))
  }
  warning('Still double-check AUDIT and other quests if correctly tallied!')
  
  if(!all(duplicated(dat_match$VPPG) == FALSE)) {
    stop('duplicates in dat_match$VPPG!')
  }
  
  # get matching ready for analysis
  # education (make a years of education variable)
  # using a dictionary with all known(!) combinations
  setwd(path_mtc)
  cur_dic_tab = xlsx::read.xlsx("dic_edu.xlsx",1)
  dat_match$edu_years        = NA
  dat_match$edu_years_voca   = NA
  dat_match$edu_hollingshead = NA
  cur_dic_tab$dat_match.edu_highest_deg   = as.character(cur_dic_tab$dat_match.edu_highest_deg)
  cur_dic_tab$dat_match.edu_deg_secondary = as.character(cur_dic_tab$dat_match.edu_deg_secondary)
  agk.as.char.NA = function(x){
    if (is.na(x)) {
      return("NA")
    } else {
      as.character(x)
    }
  }
  
  for (ii in 1:length(dat_match$VPPG)) {
    cur_education = c(as.character(dat_match$edu_highest_deg[ii]),as.character(dat_match$edu_deg_secondary[ii]))
    cur_education[grep('[NA]',cur_education,fixed = T)] = 'NA'
    first_match                    = cur_dic_tab[,c(1)] == agk.as.char.NA(cur_education[1])
    secnd_match                    = cur_dic_tab[,c(2)] == agk.as.char.NA(cur_education[2])
    cur_match                      = first_match*secnd_match
    cur_years                      = cur_dic_tab$transl[which(cur_match==1)]
    cur_years_voca                 = cur_dic_tab$transl_vocational[which(cur_match==1)]
    cur_years_holl                 = cur_dic_tab$Hollingshead[which(cur_match==1)]
    dat_match$edu_years[ii]        = as.numeric(as.character(cur_years))
    dat_match$edu_years_voca[ii]   = as.numeric(as.character(cur_years_voca))
    dat_match$edu_hollingshead[ii] = as.numeric(as.character(cur_years_holl))
  }
  
  # income
  tmp = data.frame(dat_match[c(grep(pattern="VPPG",names(dat_match)),
                               grep(pattern="^income",names(dat_match)))])
  cur_source          = attr(tmp$income_household,"levels")
  cur_transl          = c("500","1000","1500","2000","2500","3000","3500","4000","4500","5000","7000")
  tmp$income_personal = as.numeric(agk.recode(as.character(tmp$income_household),cur_source,cur_transl))
  # fixing income_hh_quest
  tmp$income_hh_quest = gsub('k.A.',NA,tmp$income_hh_quest)
  tmp$income_hh_quest = gsub(',00','',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub('.','',tmp$income_hh_quest,fixed = T)
  tmp$income_hh_quest = gsub(',-','',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub('2000-2500','2250',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub('ALG II','700',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub('€','',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub(' ','',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub(',','',tmp$income_hh_quest)
  tmp$income_hh_quest = gsub("~",'',tmp$income_hh_quest,fixed=T)
  tmp$income_hh_quest[tmp$VPPG == 'VPPG0045'] = '770'
  tmp$income_hh_quest[tmp$VPPG == 'VPPG0842'] = '451'
  tmp$income_hh_quest = gsub('NA',NA,tmp$income_hh_quest)
  cur_num_na_before = sum(is.na(tmp$income_hh_quest))
  cur_num_na_after  = sum(is.na(as.numeric(tmp$income_hh_quest)))
  if(cur_num_na_after>cur_num_na_before) {
    stop('NAs were created during income processing. Check tmp$income_hh_quest')
  }
  # replacing income personal with more refined income_hh_quest if not NA
  for (ii in 1:length(tmp$income_hh_quest)) {
    if (!is.na(tmp$income_hh_quest[ii])) {
      tmp$income_personal[ii] = tmp$income_hh_quest[ii]
    }
  }
  
  # replacing people household if NA
  for (ii in 1:length(tmp$income_for_how_many)) {
    if (is.na(tmp$income_for_how_many[ii])) {
      tmp$income_for_how_many[ii] = tmp$income_people_hh_quest[ii]
    }
  }
  # dividing by people in household
  cur_num_na_before = sum(is.na(tmp$income_personal))
  cur_num_na_after  = sum(is.na(as.numeric(tmp$income_personal)))
  if(cur_num_na_after>cur_num_na_before) {
    stop('NAs were created during income processing. Check tmp$income_hh_quest')
  }
  cur_num_na_before = sum(is.na(tmp$income_for_how_many))
  cur_num_na_after  = sum(is.na(as.numeric(tmp$income_for_how_many)))
  if(cur_num_na_after>cur_num_na_before) {
    stop('NAs were created during income processing. Check tmp$income_hh_quest')
  }
  
  # correction of tmp$income_for_how_many: 0 is not allowed; will be 1
  tmp$income_for_how_many = as.numeric(tmp$income_for_how_many)
  tmp$income_for_how_many = ifelse(tmp$income_for_how_many == 0,1,tmp$income_for_how_many)
  
  # dividing by number of people in household
  tmp$income_personal = as.numeric(tmp$income_personal)/as.numeric(tmp$income_for_how_many)
  
  dat_match = merge(dat_match,tmp,by=c('VPPG'))
  
  # smoking
  names(dat_match)[which(names(dat_match) == "ftdt")] = "smoking_ftdt"
  
  # debt
  cur_source = attr(dat_match$debt_overall,"levels")
  cur_transl = c("0","10000","25000","50000","100000","200000","300000","10000","25000","50000","100000","200000","300000","NA")
  dat_match$debt_overall_num = as.numeric(agk.recode(as.character(dat_match$debt_overall),cur_source,cur_transl))
  
  # TODO: should work with levels again; wrote Jan Philipp! 31.05.2016
  #cur_source = attr(dat_match$debt_gambling,"levels")
  #cur_transl = c("0","10000","25000","50000","100000","200000","300000")
  dat_match$debt_gambling_num = as.numeric(as.character(dat_match$debt_gambling))
  #dat_match$debt_gambling_num = as.numeric(agk.recode(as.character(dat_match$debt_gambling),cur_source,cur_transl))
  dat_match$debt_gambling_num = ifelse(is.na(dat_match$debt_gambling_num),0,dat_match$debt_gambling_num)
  dat_match$debt_gambling_num = ifelse(is.na(dat_match$debt_overall_num),NA,dat_match$debt_gambling_num)
  
  # occupation
  dat_match$unemployed = ifelse(dat_match$occup == "arbeitslos","unemployed","not unemployed")
  
  # gender
  cur_levels       = levels(dat_match$gender)[c(1,2)]
  dat_match$gender = ifelse(dat_match$gender == "männlich" | dat_match$gender == "weiblich",dat_match$gender,NA)
  dat_match$gender = as.factor(as.character(dat_match$gender))
  levels(dat_match$gender) = cur_levels
  
  # dem_gender: 'anderes' will be weiblich
  # TODO: take this out here and into quest export
  cur_ind = which(names(dat_match) == 'dem_gender')
  dat_match[[cur_ind]][grep('anderes', dat_match[[cur_ind]])] = 'weiblich'
  dat_match[[cur_ind]] = droplevels(dat_match[[cur_ind]])
  
  # now throw out all subs in dat_match that are not in data_pdt
  # this way dat_match always just worries about data_pdt subjects
  # YOU SHOULD NOT EXCLUDE ANY SUBS LATER AS THIS WILL NOT 
  # FIGURE IN DAT_MATCH THEN
  dat_match = dat_match[dat_match$VPPG %in% data_pdt$subject,]
  
  # cleaning those NA. variables
  cur_vars_drop  = grep('^NA.',names(dat_match))
  if (length(cur_vars_drop)) {
    dat_match      = dat_match[,-cur_vars_drop] 
  }
  
  ## END OF GETTING MATCHING TABLES
  
  ## PERFORM MATCHING TESTS ===================================================
  # performing the matching tests and writing results
  # getting the info on available phys data; based on subs_good_phys variable
  # get the MATLAB report on good and existent physio
  setwd(path_dat)
  setwd('..')
  subs_good_phys        = read.table('good_physio_data.txt',header=T)
  subs_good_phys        = subs_good_phys$subs_with_good_physio_data
  subs_good_phys        = agk.recode.c(subs_good_phys,tnl$PhysioVP,tnl$VPPG)
  
  # handedness "beide" will be "rechts"
  # unemployed as factor
  dat_match$handedness[dat_match$handedness == "beide"] = "rechts"
  dat_match$handedness = droplevels(dat_match$handedness)
  dat_match$unemployed = as.factor(dat_match$unemployed)
  dat_match$smoking[dat_match$smoking == "manchmal"] = "ja"
  dat_match$smoking = droplevels(dat_match$smoking)
  
  # splitting the dat_match
  cur_res      = agk.perform.matching.splitdfs(dat_match)
  dfs          = cur_res$dfs
  cur_groups   = cur_res$cur_groups
  cur_matching = cur_res$cur_matching
  cur_names    = cur_res$cur_names
  cur_gr_levs  = list(c("HC","PG"),c("HC","PG"),c("männlich","weiblich"),c("HC","PG"))
  
  # interpolating using mean per group
  dfs = agk.interpolating.dat_match(dfs,cur_groups,cur_names,cur_gr_levs)
  
  # do matching: find best matching subject for each subject
  matching_res          = agk.domatch(which_studies,desired_n,dfs,cur_groups,cur_names_dom)
  dfs                   = matching_res$dfs
  dropped_subs_matching = matching_res$dropped_HCs_PGs
  
  # reporting who was dropped due to matching
  for (ii in 1:length(which_studies)) {
    cur_text = paste("In study", which_studies[ii],"these subjects were dropped to improve matching:\n",
                     "HC:",paste(dropped_subs_matching[[ii]]$dropped_HC_matching,collapse = " "),"\n",
                     "PG:",paste(dropped_subs_matching[[ii]]$dropped_PG_matching,collapse = " "))
    warning(cur_text)
  }
  
  # align data_pdt and dat_match after do matching
  warning(paste0("dat_match and data_pdt subjects are aligned after dropping subjects due to matching.\n",
                 "But dat_match has no interpolation of missing data as was used for printing demography tables."))
  for (ii in 1:length(which_studies)) {
    data_pdt  = data_pdt[!data_pdt$subject %in% dropped_subs_matching[[ii]]$dropped_HC_matching,]
    data_pdt  = data_pdt[!data_pdt$subject %in% dropped_subs_matching[[ii]]$dropped_PG_matching,]
    dat_match = dat_match[!dat_match$VPPG %in% dropped_subs_matching[[ii]]$dropped_HC_matching,]
    dat_match = dat_match[!dat_match$VPPG %in% dropped_subs_matching[[ii]]$dropped_PG_matching,]
  }
  
  # MRI cohort: drop due to Age, edu years, lefthandedness
  # dat_match_MRT     = subset(dat_match,Cohort == 'MRT')
  # subs_drop_Age_MRT = dat_match_MRT$VPPG[dat_match_MRT$Age <20] # relevant in HC
  # subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$Age >60 & dat_match_MRT$HCPG == 'PG'])
  # subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$Age ==57 & dat_match_MRT$HCPG == 'PG' & dat_match_MRT$handedness == 'links'])
  # subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$Age ==43 & dat_match_MRT$HCPG == 'PG' & dat_match_MRT$handedness == 'links'])
  # subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$edu_years > 20]) # relevant in HC
  
  
  #subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$edu_years_voca == 8 & dat_match_MRT$HCPG == 'HC']) # relevant in HC
  #subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$edu_years_voca == 5 & dat_match_MRT$HCPG == 'HC']) # relevant in HC
  #subs_drop_Age_MRT = c(subs_drop_Age_MRT,dat_match_MRT$VPPG[dat_match_MRT$edu_years_voca == 4 & dat_match_MRT$HCPG == 'HC']) # relevant in HC
  
  # align the dfs MRI
  # dfs[[1]] = subset(dfs[[1]],!VPPG %in% subs_drop_Age_MRT)
  
  # # reporting who was dropped due to by-hand matching
  # cur_grp               = agk.recode(subs_drop_Age_MRT,dat_match$VPPG,as.character(dat_match$HCPG))
  # elim_matching_by_hand = data.frame(subs_drop_Age_MRT,cur_grp,stringsAsFactors = F)
  # cur_text = paste('In MRI study these subs were dropped by hand to improve matching',paste(elim_matching_by_hand,collapse=' '))
  # warning(cur_text)
  # 
  # # align data_pdt and dat_match after do matching
  # warning(paste0("dat_match and data_pdt subjects are aligned after dropping subjects due to matching by hand.\n",
  #                "But dat_match has no interpolation of missing data as was used for printing demography tables."))
  # data_pdt  = data_pdt[!data_pdt$subject %in% subs_drop_Age_MRT,]
  # dat_match = dat_match[!dat_match$VPPG %in% subs_drop_Age_MRT,]
  
  # core perfoming matching tests
  disp('Checking matching and printing tables.')
  match_result_tables = agk.perform.matching.tests(dfs,cur_groups,cur_matching,path_mtc,
                                                   write_match_tables = 1,cur_names)
  
  ## SETTING COHORT AND GROUP VARS FOR DATA_PDT ===============================
  
  # TODO: not so pretty yet!
  # data_pdt$cohort (lower case c) is already needed in physio aggregate
  # we use it here
  # make the Cohort-Variable and HCPG variable according to KFG
  data_pdt$Cohort     = as.character(data_pdt$cohort)
  data_pdt$cohort     = NULL
  data_pdt$HCPG       = agk.recode.c(data_pdt$subject,dat_match$VPPG,dat_match$HCPG)
  # set all the physio pilot people to HC
  data_pdt[data_pdt$Cohort == "PhysioPilot",]$HCPG = "HC" 
  complete_data       = data_pdt
  
  # TODO: where are those Pretest subs?
  # Resolve issues with the Cohort variable  - Milan 6.7.2016
  # - Make Pretest a single cohort
  data_pdt$Cohort[grep("Pretest", data_pdt$Cohort)]="Pretest"
  # - Separate the Pretest01, Pretest02, Pretest03 from the Pretest cohort because those did a different questionnaire - pretest01 who
  # had extra different gambling images, and pretest02 and 03 who gave an online feedback during the ratings, thus the PretestOral
  data_pdt$Cohort[data_pdt$subject=='Pretest01']= 'PretestO1'
  data_pdt$Cohort[data_pdt$subject=='Pretest02' | data_pdt$subject=='Pretest03']= 'PretestOral'
  
  # explanation of variables
  # TODO! FIX! with comments
  # data_pdt_vars = names(data_pdt)
  # data_pdt_vars_explanations = c("subject_ID","trial number within task","stimulus ID","gain value presented",
  #                                "loss value presented","choice by subject 1: yes 2: rather yes 3: rather no 3: no",
  #                                "reaction time", "category of stimulus 1:gam 2:neg 3:pos 4:neutr 5:Polish neutral 6:weinreich neg 7: weinreich pos 8: weinreich neutr",
  #                                "side: where was the gain value shown?", "how long was the stimulus shown?","age of participant", "sex of partic",
  #                                "gain uncentered", "loss uncentered", "gainxloss agg centered and agg", "gainxloss in orig form","euclidean distance of gamble from gamble matrix diagonal","euclidean distance of gamble from gamble matrix diagonal stand.",
  #                                "expected value of gamble","ratio of gamble options","difference between options","Risk according to Minati", "risk according to Alex G",
  #                                "risk according to Martino","Skew of gamble according to Minati","expected value uncentered",
  #                                "choice dichotomized","bold extracted from ROI z-score_l_precun","bold_post_cing","bold_r_precun","bold_sup_temp","corrugator mean activity","zygomaticus mean activity","eda mean activity","number when picture of this trial was shown in post hoc rating",
  #                                "category of picture but not used anymore","craving for gambling question for picture","representative of gambling question of pic","pic representative for negative consequences of gambling",
  #                                "pic repres for positive consequences of abstaining","pic makes you question your gambling habits?","arousal rating of pic SAM","dominance rating of pic SAM","valence ratings of pic SAM",
  #                                "how many seconds did subject spend on rating (not SAM) this picture","how long did they spend on SAM ratings","arousal rating unstand","valence rating unstand","dominance rating unst",
  #                                "craving for gambling question for picture stand.","representative of gambling question of pic stand.","pic representative for negative consequences of gambling stand.",
  #                                "pic repres for positive consequences of abstaining stand.","pic makes you question gambling stand.","an old image ID vars, not used anymore","question that changed depending on cat, neg, pos, gambling question see above; only used with early subjects",
  #                                "old image ID variable, only useful in processing, not for analysis",
  #                                "fixation cross - number of estimated SCR peaks according to LEDA","fixation cross - latency of first SCR peak (LEDA)","fixation cross - sum of amplitudes of all reconvolved SCR with onset in response window",
  #                                "fixation cross - average phasic driver activity within response window","fixation cross - integrated phasic driver activity within response window","fixation cross - maximum phasic response(?)","fixation cross - tonic response",
  #                                "image stimulus onset - number of estimated SCR peaks according to LEDA","image stimulus onset - latency of first SCR peak (LEDA)","image stimulus onset - sum of amplitudes of all reconvolved SCR with onset in response window",
  #                                "image stimulus onset - average phasic driver activity within response window","image stimulus onset - integrated phasic driver activity within response window","image stimulus onset - maximum phasic response(?)","image stimulus onset - tonic response",
  #                                "gamble onset - number of estimated SCR peaks according to LEDA","gamble onset - latency of first SCR peak (LEDA)","gamble onset - sum of amplitudes of all reconvolved SCR with onset in response window",
  #                                "gamble onset - average phasic driver activity within response window","gamble onset - integrated phasic driver activity within response window","gamble onset - maximum phasic response(?)","gamble onset - tonic response",
  #                                "which part of study was this subject part of","HC or PG according to KFG>16")
  # 
  # 
  # variables_explained_pdt = data.frame(data_pdt_vars,data_pdt_vars_explanations)
  
  # get an MRI Sjinfo for SPM
  cur_MRI        = aggregate(data_pdt[c("HCPG","Cohort")],by=list(data_pdt$subject),FUN=first)
  cur_MRI        = subset(cur_MRI, Cohort == "MRT")
  cur_MRI$Cohort = NULL
  names(cur_MRI) = c("subject","group") 
  cur_MRI$group  = ifelse(cur_MRI$group == "PG",1,0)
  # add the covariates info
  dat_match$edu_years_sum = dat_match$edu_years + dat_match$edu_years_voca 
  cur_MRI                 = merge(cur_MRI,dat_match[c('VPPG','edu_years_sum','smoking_ftdt')],by.x = c('subject'),by.y = c('VPPG'))
  write.table(file="info_mri_selection.csv",x = cur_MRI,sep = "\t",row.names = F,quote = F)
  
  setwd(path_dat)
  write.table(file="info_mri_selection.csv",x = cur_MRI,sep = "\t",row.names = F,quote = F)
  
  # save the import
  setwd(path_dat)
  save(file="data_pdt.rda",list = c("data_pdt","dat_match"))
  # save the workspace
  setwd(path_dat)
  save(file="data_pdt_Maja.rda",list = ls())
  
  # get the backup
  data_pdt_bcp = data_pdt
  dat_match_bcp = dat_match
  
  tryCatch({
  # google drive
  setwd(path_dat_GD)
  write.table(file="info_mri_selection.csv",x = cur_MRI,sep = "\t",row.names = F,quote = F)
  setwd(path_dat_GD)
  save(file="data_pdt_Maja.rda",list = ls())
  setwd(path_dat_GD)
  save(file="data_pdt.rda",list = c("data_pdt","dat_match"))
  }, error = function (e) {
    disp('No saving to Google Drive possible.')
  })
  
} else {
  # import from existing import
  # check what paths work
  path_GD = tryCatch({
    setwd(path_dat_GD)
    path_GD = T
  }, error = function(e) {
    path_GD = F
  })
  
  path_S = tryCatch({
    setwd(path_dat)
    path_SD = T
  }, error = function(e) {
    path_SD = F
  })
  
  if (path_GD == F & path_S == F) {
    stop('Cannot access neither Google Drive nor S:')
  }
  
  if (path_S == T) {
    # access to S:?
    disp('Trying to access S:')
    setwd(path_dat)
    load("data_pdt.rda")
  }
  
  if (path_GD == T & path_S == F) {
    disp('Loading from google drive because S path not accessible.')
    setwd(path_dat_GD)
    load("data_pdt.rda")
  }
  
  # get an MRI Sjinfo for SPM
  cur_MRI        = aggregate(data_pdt[c("HCPG","Cohort")],by=list(data_pdt$subject),FUN=first)
  cur_MRI        = subset(cur_MRI, Cohort == "MRT")
  cur_MRI$Cohort = NULL
  names(cur_MRI) = c("subject","group") 
  cur_MRI$group  = ifelse(cur_MRI$group == "PG",1,0)
  # add the covariates info
  dat_match$edu_years_sum = dat_match$edu_years + dat_match$edu_years_voca 
  cur_MRI                 = merge(cur_MRI,dat_match[c('VPPG','edu_years_sum','smoking_ftdt')],by.x = c('subject'),by.y = c('VPPG'))
  write.table(file="info_mri_selection.csv",x = cur_MRI,sep = "\t",row.names = F,quote = F)
  
  if (path_S == F) {
    disp('loaded GD data_pdt_bcp, dat_match_bcp')
  }
  
  if (path_GD == T) {
    # save the import
    setwd(path_dat_GD)
    save(file="data_pdt.rda",list = c("data_pdt","dat_match"))
    # save the workspace
    setwd(path_dat_GD)
    save(file="data_pdt_Maja.rda",list = ls())
  }
  
  if (path_S == T) {
    # save the import
    setwd(path_dat)
    save(file="data_pdt.rda",list = c("data_pdt","dat_match"))
    # save the workspace
    setwd(path_dat)
    save(file="data_pdt_Maja.rda",list = ls())
  }
}

# here we will add the ledalab data
# init
data_pdt$SCR          = NULL
data_pdt$SCR_gamble   = NULL
# get data
setwd(path_dat_GD)
leda_dat       = read.table('ledalab_out.csv',sep=',',header=T)
kick_out       = which(names(leda_dat) %in% c('Latency','AmpSum','ISCR','nSCR','PhasicMax','Tonic') == FALSE)
leda_dat       = leda_dat[,kick_out]
leda_dat$VPPG  = agk.recode.c(leda_dat$subject,tnl$PhysioVP,tnl$VPPG)

leda_dat_st    = subset(leda_dat,onset != 'gamble')
leda_dat_gam   = subset(leda_dat,onset == 'gamble')
leda_dat       = leda_dat_st

# change name to SCR_gamble
names(leda_dat_gam)[names(leda_dat_gam) == 'SCR'] = 'SCR_gamble'

kick_out       = which(names(leda_dat) %in% c('subject','onset') == FALSE)
leda_dat       = leda_dat[,kick_out]
leda_dat_gam   = leda_dat_gam[,kick_out]
data_pdt       = merge(data_pdt,leda_dat,by.x = c('subject','trial'),by.y = c('VPPG','trial'),all.x = T,all.y = F)
data_pdt       = merge(data_pdt,leda_dat_gam,by.x = c('subject','trial'),by.y = c('VPPG','trial'),all.x = T,all.y = F)


# scale per subject SCR
all_subs = unique(data_pdt$subject)
for (ss in 1:length(all_subs)) {
  cur_dat                                        = data_pdt$SCR[data_pdt$subject == all_subs[ss]]
  cur_dat                                        = scale(get.log(cur_dat))
  data_pdt$SCR[data_pdt$subject == all_subs[ss]] = cur_dat
  
  cur_dat                                               = data_pdt$SCR_gamble[data_pdt$subject == all_subs[ss]]
  cur_dat                                               = scale(get.log(cur_dat))
  data_pdt$SCR_gamble[data_pdt$subject == all_subs[ss]] = cur_dat
}

# numeric accept_reject (no good!)
if (acc_num == 1) {
  data_pdt$accept_reject = NA
  data_pdt$accept_reject = data_pdt$choice
  data_pdt$accept_reject[data_pdt$accept_reject == 5] = NA
  data_pdt_bcp = data_pdt
}

# create bcp for select study and other analysis scripts
data_pdt_bcp  = data_pdt
dat_match_bcp = dat_match





# # some descriptives
# data_pdt_agg = aggregate(cbind(age,sex,pilot) ~ subject, data = data_pdt,first)
# describeBy(as.numeric(as.character(data_pdt_agg$age)),group = data_pdt_agg$pilot)
# xtabs(~pilot+sex,data = data_pdt_agg)

## fit behavioral models and export model params for MRI analysis
# path
# if user allows
# source optimization scripts based on data_pdt

# if (get_LAch_for_MRI) {
#   setwd(path_mod)
#   setwd("LA_Charpentier")
#   source("MLE_Charpentier_updated.R")
# }
# 
# if (get_LAch_for_MRI) {
#   setwd(path_mod)
#   setwd("LA_cl")
#   source("MLE_Genauck.R")
# }

# if (get_MRI_extr == 1) {
#   data_pdt_MRI = subset(data_pdt,Cohort== "MRT")
#   bold_overall = aggregate(data_pdt_MRI$bold,by=list(data_pdt_MRI$stim),FUN=median,na.rm=T)
#   names(bold_overall) = c("stim","bold_overall")
#   bold_bygroup = aggregate(data_pdt_MRI$bold,by=list(data_pdt_MRI$stim,data_pdt_MRI$HCPG),FUN=median,na.rm=T)
#   names(bold_bygroup) = c("stim","HCPG","bold_bygroup")
#   data_pdt = merge(data_pdt,bold_overall,by=c("stim"))
#   data_pdt = merge(data_pdt,bold_bygroup,by=c("stim","HCPG"))
#   data_pdt_bcp = data_pdt
#   }


## temp ##
#cur_list = lmList(accept_reject ~ gain + loss |subject, pool =F, na.action = NULL, family = "binomial",data =data_pdt)
#BICtry = function(x) {
#  out <- tryCatch(
#    {
#      BIC(x)
#    },
#    error=function(cond)
#    {
#      NA
#    })
#  return(out)
#}


# for Guillaume
# # BGG study
# setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/OFC_project_Guillaume/data_sent_to_Guillaume_Sescousse/Data/VPPG/Behav')
# cur_dat = xlsx::read.xlsx("VPPG_01_all_behav_AG.xlsx",1)
# desired_vars=c('VPPG','BIS')
# cur_dat_match = dat_match[desired_vars]
# cur_dat=merge(cur_dat,cur_dat_match,by.x="ID",by.y="VPPG",all.x = T,all.y = F)
# xlsx::write.xlsx(cur_dat,file="VPPG_01_all_behav_AG_ed_2016_11_25.xlsx")


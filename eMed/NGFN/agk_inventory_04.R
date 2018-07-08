# read in the inventroy list to get infos on
# NGFN data
library(haven)
library(pracma)
library(stringr)

## paths ======================================================================
# user
user            = paste0(as.character(Sys.info()["login"]),"/")

warning('MANNHEIM: USE ALWAYS ID_NGFN (28...) CODES FOR CONSISTENCY')
warning('MANNHEIM SPSS HAS BOTH ID INFOS')

# inventory
path_inventory  = 'S:/AG/AG-Emotional-Neuroscience/Restricted/NGFN/eMed/Inventarlisten/Inventarliste_Daten_all_Sites_22.12.2016.xls'
path_inv_ponly  = 'S:/AG/AG-Emotional-Neuroscience/Restricted/NGFN/eMed/Inventarlisten/'

# spss masks
path_spss_mnh   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Mannheim/Screen_Base_Neuro_SUMMENWERTE_ALL_FINAL_2013_Multicentre_gelöscht_23.01.2013.sav'
path_spss_bnn   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Bonn/NGFN13_Bonn_alle_Daten.sav'
path_spss_bna   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Übergabe_NGFN_Bonn/Finale Abgabe/Fragebîgen/NGFN13_Bonn_alle_Fragebîgen_ausgewertet.sav'
path_spss_bln   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Berlin/NGFNplus_Datenmaske_Berlin_05.03.13_KCH.sav'

# relapse data
path_relabs_bln = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Berlin/NGFNplus_Datenmaske_Follow-ups_13.01.14_KCH_FINAL.sav'
path_relabs_mnh = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Mannheim/NGFN13_RÜCKFÄLLE_NEU_ALLE_20140521.sav'
path_relabs_bnn = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Bonn/Follow Up/NGFNplus_Follow-Ups_Bonn.sav'

# genetics PsychChip
path_genetics   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik/Genetik-PsychChip_S.Erk_09.03.2017'
path_genetic_oe = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik/Genetik_OmniExpress_S_Erk_2018_Feb_01'

# genetics Zuordnung USI-Code
path_gen_usi_bm = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik/Zuordnung/Berlin_Mannheim/Genetik_USI_Codes_Gesamtliste_Berlin_Bonn_24.03.2013.xlsx'
path_gen_usi_bn = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik/Zuordnung/Bonn/Genetik Zuordnung.xlsx'
base_lib        = paste("C:/Users/",user,"Google Drive/",sep="")
path_lib        = paste0(base_lib,"Library/R")

## settings ===================================================================
# what to predict? relapse or diagnosis
des_pred  = 'diagnosis'
# use rules or just take any where there is at least not "NA" everywhere
# LEAVE T
use_rules = T
sites           = c('Mannheim','Bonn','Berlin')
## packages and functions =====================================================
library(pracma)
library(compiler)
library(plyr)
library(readxl)
library(foreign)
agk.recode <- function(x,y,z) {
  # a recode function
  # function to recode x, given a source vector y
  # and a translated vector z
  x = as.character(x)
  y = as.character(y)
  z = as.character(z)
  for (ii in 1:length(x)) {
    done <- 0
    for (jj in 1:length(y)) {
      # NA in x will be NA
      if(is.na(x[ii])) {
        x[ii] <- NA
        break
      }
      if (x[ii] == y[jj]) {
        x[ii] <- z[jj]
        done <- 1
      }
      if (done == 1) {break}
    }
  }
  return(x)
}
agk.recode.c = cmpfun(agk.recode)

## inventories of three places ================================================
inventories = list()
for (vv in 1:3) {
  invM          = as.data.frame(read_excel(path_inventory,sheet = vv,col_names = T))
  
  # the second line is still part of the header
  newHeader = c()
  for (ii in 1:length(invM[1,])) {
    if(!is.na(invM[1,ii])) {
      newHeader[ii] = paste(names(invM)[ii],invM[1,ii],sep='_')
    } else {
      newHeader[ii] = names(invM)[ii]
    }
  }
  
  # clean the newHeader
  newHeader         = gsub('DICOM__.','DICOM',newHeader)
  invM              = invM[-1,]
  names(invM)       = newHeader
  
  # add site
  invM$site    = sites[vv]
  
  # add scanner
  invM$scanner = paste0(sites[vv],'_scanner') 
  
  # pack
  inventories[[vv]] = invM
}

inv = rbind.fill(inventories)

# clean X__1
for (ii in 1:length(inv$X__1)) {
  if (is.na(inv$X__1[ii])) {
    next
  } else {
    inv$Notes[ii] = paste(inv$Notes[ii],inv$X__1[ii])
    inv$X__1[ii]  = NA
  }
}

# kill all variables that have only NA
cur_fun = function(x) {return(all(is.na(x)))}
inv     = inv[unlist(lapply(inv,cur_fun) == FALSE)]

# trim white space
inv = (as.data.frame(lapply(inv,trimws)))

# make one ID variable
inv$ID = NA
# Berlin
inv_bln    = subset(inv, site == 'Berlin')
inv_bln$ID = as.character(inv_bln$Co.ID)
inv[inv$site == 'Berlin',] = inv_bln

# Mannheim
# NGFN_initials counts! (28...)
inv_mnm    = subset(inv, site == 'Mannheim')
inv_mnm$ID = as.character(inv_mnm$ID_initials_NGFN)
for (ii in 1:length (inv_mnm$ID)) {
  if(is.na(inv_mnm$ID[ii])) {
    # check if a ^28 code is available
    cur_re = regexp(as.character(inv_mnm$ID_initials[ii]),'^28')
    if (!is.null(cur_re$start)) {
      inv_mnm$ID[ii] = as.character(inv_mnm$ID_initials[ii])
    }
  }
}
inv[inv$site == 'Mannheim',] = inv_mnm

# Bonn (using the ADD_P variable)
inv_bnn    = subset(inv, site == 'Bonn')
inv_bnn$ID = as.character(inv_bnn$ID_initials_NGFN)
for (ii in 1:length(inv_bnn$ID)) {
  if (!isempty(grep('_5[0-9][0-9]',inv_bnn$ID[ii]))) {
    inv_bnn$ID[ii] = str_pad(inv_bnn$ID[ii],width = 7,pad = '0')
  } else {
    inv_bnn$ID[ii] = str_pad(inv_bnn$ID[ii],width = 3,pad = '0')
  }
}

inv_bnn$ID = paste0('ADD_P_',inv_bnn$ID)
inv[inv$site == 'Bonn',] = inv_bnn

# add scanner variable
inv$scanner       = ifelse(inv$scanner == 'Berlin_scanner','Berlin_scanner_PTBnew',inv$scanner)
spss_bln          = read_spss(path_spss_bln)

# write a csv of the spss data (Berlin)
cur_home = getwd()
spss_bln_p  = strsplit(path_spss_bln,'/')
spss_bln_p  = paste(spss_bln_p[[1]][1:(length(spss_bln_p[[1]])-1)],collapse='/')
setwd(spss_bln_p)
spss_mtlb = as.data.frame(spss_bln)
spss_mtlb = as.data.frame(lapply(spss_mtlb,trimws))
write.table(spss_mtlb,file = 'spss_bln.csv',sep = '\t',row.names = F)
setwd(cur_home)

# write a csv of the spss data (Bonn)
spss_bnn = read_spss(path_spss_bnn)
cur_home = getwd()
spss_bnn_p  = strsplit(path_spss_bnn,'/')
spss_bnn_p  = paste(spss_bnn_p[[1]][1:(length(spss_bnn_p[[1]])-1)],collapse='/')
setwd(spss_bnn_p)
spss_mtlb = as.data.frame(spss_bnn)
spss_mtlb = as.data.frame(lapply(spss_mtlb,trimws))
write.table(spss_mtlb,file = 'spss_bnn.csv',sep = '\t',row.names = F)
setwd(cur_home)


for (ii in 1:length(inv$Co.ID)) {
  cur_id      = inv$Co.ID[ii]
  if (!is.na(cur_id)) {
  } else {
    next
  }
  cur_ind = grep(paste0('^',cur_id,'$'),spss_bln$ID)
  if (length(cur_ind) != 0) {
    cur_note = spss_bln$notes[cur_ind]
    if (length((grep('TRIO',cur_note)))>0) {
      inv$scanner[ii] = 'Berlin_scanner_PTBold'
    }
  }
}

# fix the Drop-out variable
inv$Drop.Out = agk.recode.c(inv$Drop.Out,c('Einschlus','Einschluss'),c('include','include'))

# recode group variable
inv$group = agk.recode.c(inv$group,c('Kontrollgruppe','Patient','Patienten'),c('controls','patients','patients'))

## REL ABS ====================================================================
warning('TODO: still fix code so that only $ID is used')

# get additional data (Berlin)
relabs_bl          = read.spss(path_relabs_bln)
vars_of_int        = c('ID','t6_Schwerer_Rückfall')
relabs_bl          = relabs_bl[vars_of_int]
relabs_bl          = (as.data.frame(lapply(relabs_bl,trimws)))
relabs_bl$REL      = NA
# code to relapse
relabs_bl$REL[grep('Schwer rückfällig:',relabs_bl$t6_Schwerer_Rückfall,ignore.case = FALSE)] = 'REL'
relabs_bl$REL[grep('Abstinent',relabs_bl$t6_Schwerer_Rückfall)] = 'ABS'
relabs_bl$REL[grep('Rückfällig',relabs_bl$t6_Schwerer_Rückfall)] = 'ABS'

# get additional data (Mannheim)
# relapse (heavy, 6 months after t0)
relabs_mn                = read.spss(path_relabs_mnh)
vars_of_int              = c('ID','ID_ngfn','ID_sfb','sRF_180Tage_t0')
relabs_mn                = relabs_mn[vars_of_int]
relabs_mn                = (as.data.frame(lapply(relabs_mn,trimws)))
relabs_mn$sRF_180Tage_t0 = agk.recode.c(relabs_mn$sRF_180Tage_t0,c('ja','nein'),c('REL','ABS'))

# get additional info (Bonn)
relabs_bn                = read.spss(path_relabs_bnn)
vars_of_int              = c('ADD','MP','t6_Schwerer_Rückfall')
relabs_bn                = relabs_bn[vars_of_int]
relabs_bn                = (as.data.frame(lapply(relabs_bn,trimws)))
relabs_bn$REL            = NA
# code to relapse
relabs_bn$REL[grep('Schwer rückfällig:',relabs_bn$t6_Schwerer_Rückfall,ignore.case = FALSE)] = 'REL'
relabs_bn$REL[grep('Abstinent',relabs_bn$t6_Schwerer_Rückfall)] = 'ABS'
relabs_bn$REL[grep('Rückfällig',relabs_bn$t6_Schwerer_Rückfall)] = 'ABS'

# fill it into REL ABS
inv$relapse = NA
disp('info on relpase before filling in from external source')
inv_pat = subset(inv, group == 'patients')
disp('info on relpase after filling in from external source')
print(table(inv_pat$relapse %in% c('REL','ABS')))

# Berlin
disp('RELAPSE BERLIN')
for (ii in 1:length(inv$relapse)) {
  if (is.na(inv$relapse[ii])) {
    cur_id  = inv$ID[ii]
    if (inv$site[ii] != 'Berlin') {
      next
    }
    if (!is.na(inv$group[ii])) {
      if (inv$group[ii] == 'controls') {
        next
      }
    }
    cur_ind = grep(cur_id,relabs_bl$ID)
    # fill in from external relabs
    if (!isempty(cur_ind)) {
      inv$relapse[ii] = relabs_bl$REL[cur_ind]
    } else {
      disp('no relapse info match')
    }
  }
}

# Mannheim (here going from two ID variables)
disp('RELAPSE MANNHEIM')
for (ii in 1:length(inv$relapse)) {
  if (is.na(inv$relapse[ii])) {
    cur_id  = inv$ID[ii]
    if (is.na(cur_id) | inv$site[ii] != 'Mannheim') {
      next
    }
    if (!is.na(inv$group[ii])) {
      if (inv$group[ii] == 'controls') {
        next
      }
    }

    cur_ind = grep(paste0('^',cur_id,'$'),relabs_mn$ID_ngfn)
    # fill in from external relabs
    if (!isempty(cur_ind)) {
      inv$relapse[ii] = relabs_mn$sRF_180Tage_t0[cur_ind]
    } else {
      disp('try to resort to sfb id_initials')
      cur_id  = inv$ID_initials[ii]
      if (is.na(cur_id)) {
        disp('no relapse info match')
        next
      }
      cur_ind = grep(paste0('^',cur_id,'$'),relabs_mn$ID_sfb)
      # fill in from external relabs
      if (!isempty(cur_ind)) {
        inv$relapse[ii] = relabs_mn$sRF_180Tage_t0[cur_ind]
      } else {
        disp('no relapse info match found')
      }
    }
  }
}

# Bonn
disp('RELAPSE BONN')
relabs_bn$ADD = gsub('ADDP_P_','ADD_P_',relabs_bn$ADD)
for (ii in 1:length(inv$relapse)) {
  if (is.na(inv$relapse[ii])) {
    if (is.na(cur_id) | inv$site[ii] != 'Bonn') {
      next
    }
    if (!is.na(inv$group[ii])) {
      if (inv$group[ii] == 'controls') {
        next
      }
    }
    cur_id  = inv$ID[ii]
    if (is.na(cur_id)) {
      stop('no ID info in inventory')
    }
    cur_ind = grep(cur_id,relabs_bn$ADD)
    # fill in from external relabs
    if (!isempty(cur_ind)) {
      if (length(cur_ind) > 1) {stop('too many matches')}
      inv$relapse[ii] = relabs_bn$REL[cur_ind]
    } else {
      disp('no relapse info match found')
    }
  }
}

inv_pat = subset(inv, group == 'patients')
disp('info on relpase after filling in from external source')
print(table(inv_pat$relapse %in% c('REL','ABS')))

# recode the relapse variable
#inv$relapse = as.character(inv$relapse)
#inv$relapse = ifelse(is.na(inv$relapse),'NA',inv$relapse)

## add gender and age information =============================================
spss_bnn = as.data.frame(lapply(as.data.frame(read.spss(path_spss_bnn)),trimws))
spss_mnm = as.data.frame(lapply(as.data.frame(read.spss(path_spss_mnh)),trimws))

# age
inv$sex = NA
inv$age = NA
for (ii in 1:length(inv$ID)) {
  cur_sub    = inv$ID[ii]
  cur_sub_ed = gsub('ADD_P_0+','',cur_sub)
  cur_sub_ed = gsub('ADD_P_','',cur_sub_ed)
  
  if (sum(spss_mnm$ID_initials_ngfn %in% cur_sub) == 1) {
    # Mannheim NGFN
    disp('found Mannheim NGFN')
    inv$age[ii] = as.character(spss_mnm$age[which(spss_mnm$ID_initials_ngfn %in% cur_sub)])
    inv$sex[ii] = as.character(spss_mnm$sex[which(spss_mnm$ID_initials_ngfn %in% cur_sub)])
  } else if (sum(spss_mnm$ID_initials_sfb %in% cur_sub) == 1) {
    # Mannheim SFB
    disp('found Mannheim SFB')
    inv$age[ii] = as.character(spss_mnm$age[which(spss_mnm$ID_initials_ngfn %in% cur_sub)])
    inv$sex[ii] = as.character(spss_mnm$sex[which(spss_mnm$ID_initials_ngfn %in% cur_sub)])
  } else if (sum(spss_bln$ID %in% cur_sub) == 1) {
    # Berlin
    disp('found Berlin')
    inv$age[ii] = as.character(spss_bln$age[which(spss_bln$ID %in% cur_sub)])
    inv$sex[ii] = as.character(factor(spss_bln$sex,labels = c('weiblich','männlich'))[which(spss_bln$ID %in% cur_sub)])
  } else if (sum(spss_bnn$ADD_P %in% cur_sub_ed) == 1) {
    # Bonn
    disp('found Bonn')
    inv$age[ii] = as.character(spss_bnn$age[which(spss_bnn$ADD_P %in% cur_sub_ed)])
    inv$sex[ii] = as.character(spss_bnn$sex[which(spss_bnn$ADD_P %in% cur_sub_ed)])
  } else {
    disp(paste(ii,'no match found'))
  }
  
}

inv$sex = agk.recode.c(inv$sex,c('männlich','weiblich'),c('m','f'))

## add psychosocial data ======================================================
# Berlin
spss_bln = as.data.frame(lapply(as.data.frame(read.spss(path_spss_bln)),trimws))
psd_bln  = spss_bln[]



## print the complete inventory ===============================================
setwd(path_inv_ponly)
write.table(inv,file = 'ngfn_inventory.csv',sep=';',row.names = F)

## TEST COMPLETENESS OF DATA SET FOR PREDICTION OF GROUP ======================
warning('TODO completeness check seems faulty; the rel abs info that remains is very little')
if (des_pred == 'relapse') {
  inv      = subset(inv, group == 'patients')
}

# complete data MRI tasks and out-of MRI tasks; for prediction HC vs. AD
if (des_pred == 'diagnosis') {
  all_vars = names(inv)
  int_vars = all_vars[c(3,4,9,10:16,17,18,19,20,21,22,23)]
  rul_vars = list(c('controls','patients'),c('include'),c(),c(),c(),c(),
                  c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
} else if (des_pred == 'relapse') {
  all_vars = names(inv)
  int_vars = all_vars[c(3,4,9,10:16,17,18,19,20,21,22,23,29)]
  rul_vars = list(c('controls','patients'),c('include'),c(),c(),c(),c(),
                  c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),
                  c('REL','ABS'))
}

test_rules = function(inv,int_vars,rul_vars) {
  # test all the rules
  sub_ok = c()
  for (ii in 1:length(inv[,1])) {
    ## ii counts rows
    sub_rules_met = c()
    for(jj in 1:length(int_vars)) {
      # jj for columns
      # cur rule
      cur_rule = rul_vars[[jj]]
      
      # current observed value
      cur_val = inv[ii,grep(pattern = int_vars[jj], names(inv))]
      
      # pre-test
      if (is.na(cur_val)) {
        sub_rules_met[jj] = F
        next
      }
      # test the rule
      if (use_rules) {
        if (!isempty(cur_rule)) {
          if (any(cur_rule == cur_val)) {
            sub_rules_met[jj] = T
          } else {
            sub_rules_met[jj] = F
          }
        } else {
          sub_rules_met[jj] = T
        }
      } else {
        sub_rules_met[jj] = T
      }
    }
    if (all(sub_rules_met)) {
      sub_ok[ii] = T
    } else {
      sub_ok[ii] = F
    }
  }
  return(sub_ok)
}

sub_ok = test_rules(inv,int_vars,rul_vars)

## REPORT inventory that is okay for ==========================================
# first let's look how many people we have complete at t1
if (des_pred == 'relapse') {
  disp('')
  disp('Distribution of complete data at t1; NA means no t2 data.')
  print(table(inv$relapse,useNA = 'always'))
}

# complete data MRI tasks and out-of MRI tasks; for prediction diagnosis
# with few cases where fmri data is irregular
disp('')
disp(paste('complete data MRI tasks and out-of MRI tasks; for prediction', des_pred))
if (use_rules) {
  disp('with few cases where fmri data is irregular')
} else {
  disp('with few cases where fmri data is irregular and few where behav is irregular; at least nowhere NA')
}
disp('')
print(table(sub_ok))

# who is in/out?
invIN  = inv[which(sub_ok),]
invOUT = inv[which(sub_ok==FALSE),]

# how many per location?
disp('')
disp('who is in/out?')
disp('how many per location?')
print(table(invIN$site))
print(xtabs(~site+group, data=invIN))
print(xtabs(~site+group, data=invIN))

# for case relapse
disp('')
disp('how relapse-abstain per location')
if (des_pred == 'relapse') {
  print(xtabs(~site+relapse, data=invIN))
  print(xtabs(~site+relapse, data=invIN))
}

# how many dropout per location?
disp('')
disp('how many dropout per location?')
print(table(invOUT$site))

# in case of relapse we want to know who is DO just because no relapse data
int_vars = all_vars[c(3,4,9,10:16,17,18,19,20,21,22,23)]
rul_vars = list(c('controls','patients'),c('include'),c(),c(),c(),c(),
                c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
sub_ok_AD_but_no_relapse_data = test_rules(invOUT,int_vars,rul_vars)
table(sub_ok_AD_but_no_relapse_data)
invOUT_cmpl_but_no_rel_info = invOUT[sub_ok_AD_but_no_relapse_data,]
table(invOUT_cmpl_but_no_rel_info$site)

## GENETICS ===================================================================
# function needed
agk.strfind.gsub = function(pattern,vector) {
  # function to ...
  found_it = c()
  for (ff in 1:length(vector)) {
    if (is.null(strfind(pattern,vector[ff]))) {
      found_it[ff] = F
    } else {
      found_it[ff] = T
    }
  }
  return(which(found_it))
}

setwd(path_genetics)
gen1   = read.table("NGFN-ALK-PsychChip_15048346_B.fam")
gen2   = read.table("2016-201-ILL_NGFN_Psych1-1_N_52.fam")
gen_ps = rbind(gen1,gen2)
setwd(path_genetic_oe)
gen3   = read.table("NGFN-ALC-OExpress12v1-0_J.fam")
gen_ps = rbind(gen_ps,gen3)
stopifnot(all(duplicated(gen_ps$V1) == FALSE))

# some more spss files
spss_bnn = as.data.frame(lapply(as.data.frame(read.spss(path_spss_bnn)),trimws))
spss_mnm = as.data.frame(lapply(as.data.frame(read.spss(path_spss_mnh)),trimws))

# USI mapping files
gen_usi_bm = as.data.frame(read_excel(path_gen_usi_bm))
gen_usi_bm = (as.data.frame(lapply(gen_usi_bm,trimws)))
gen_usi_bn = as.data.frame(read_excel(path_gen_usi_bn))
gen_usi_bn = as.data.frame(lapply(gen_usi_bn,trimws))

# Berlin: using USI
# strip usi from control numbers
gen_usi_bln                   = gen_usi_bm[grep('CHB',gen_usi_bm$USICODE),]
gen_usi_bln$USICODE_stripped  = gsub('DE[0-9][0-9]','DE',gen_usi_bln$USICODE)
gen_usi_bln$USICODE_stripped  = gsub('CHBB1','CHB',gen_usi_bln$USICODE_stripped)
gen_usi_bln$USICODE_stripped  = gsub('CHBI1','CHB',gen_usi_bln$USICODE_stripped)
gen_usi_bln$SCREENINGNUMMER   = as.character(gen_usi_bln$SCREENINGNUMMER)
gen_usi_bln$BESCHRIFTUNG      = as.character(gen_usi_bln$BESCHRIFTUNG)
# check for right form of screening number
for (ss in 1:length(gen_usi_bln[,1])) {
  cur_sn = as.character(gen_usi_bln$SCREENINGNUMMER[ss])
  cur_rs = regexp(s=cur_sn,pat = '[0-9]{3}_[0-9]_[0-9]{3}')
  if(isempty(cur_rs$start)) {
    gen_usi_bln$SCREENINGNUMMER[ss] = gen_usi_bln$BESCHRIFTUNG[ss] 
  }
}
# put subject code on genetic data
gen_bln                     = gen_ps[grep('chb',gen_ps$V1),]
gen_bln$USICODE_stripped    = toupper(gen_bln$V1)
subinfo                     = gen_usi_bln[c('USICODE_stripped','SCREENINGNUMMER')]
gen_bln                     = merge(gen_bln,subinfo, by='USICODE_stripped',all.x=T,all.y=T)
gen_bln                     = merge(gen_bln,spss_bln[c('ID','sex','Diagnose')],by.x = 'SCREENINGNUMMER',by.y='ID',all.x = T,all.y=T)
gen_bln$site                = 'Berlin'
# cleaning names
names(gen_bln)[grep('Diagnose',names(gen_bln))] = 'group'
names(gen_bln)[grep('SCREENINGNUMMER',names(gen_bln))] = 'ID'

# Mannheim
# strip usi from control numbers
gen_usi_mnm                   = gen_usi_bm[grep('ZIM',gen_usi_bm$USICODE),]
gen_usi_mnm$USICODE_stripped  = gsub('DE[0-9][0-9]','DE',gen_usi_mnm$USICODE)
gen_usi_mnm$USICODE_stripped  = gsub('ZIMB1','ZIM',gen_usi_mnm$USICODE_stripped)
gen_usi_mnm$USICODE_stripped  = gsub('ZIMI1','ZIM',gen_usi_mnm$USICODE_stripped)
gen_usi_mnm$USICODE_stripped  = gsub('AAA','A',gen_usi_mnm$USICODE_stripped)
gen_usi_mnm$SCREENINGNUMMER   = as.character(gen_usi_mnm$SCREENINGNUMMER)
gen_usi_mnm$BESCHRIFTUNG      = as.character(gen_usi_mnm$BESCHRIFTUNG)
# sfb to ngfn code
gen_usi_mnm$SCREENINGNUMMER = agk.recode.c(gen_usi_mnm$SCREENINGNUMMER,spss_mnm$ID_initials_sfb,spss_mnm$ID_initials_ngfn)
# fix some of the screeningnummers
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '2801012'] = '2801012_KS'
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '2202044'] = as.character(na.omit(inv$ID_initials_NGFN[inv$ID_initials == '2202044_JA']))
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '2202066'] = as.character(na.omit(inv$ID_initials_NGFN[inv$ID_initials == '2202066_DB']))
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '2802056'] = as.character(na.omit(inv$ID_initials_NGFN[inv$ID_initials_NGFN == '2802056_EL']))
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '22020591_RW'] = as.character(na.omit(inv$ID_initials_NGFN[inv$ID_initials == '2202059_RW']))
gen_usi_mnm$SCREENINGNUMMER[gen_usi_mnm$SCREENINGNUMMER == '2202063_HP'] = as.character(na.omit(inv$ID_initials_NGFN[inv$ID_initials == '2202063_HP']))
# put subject code on genetic data
gen_mnm                      = gen_ps[grep('zim',gen_ps$V1),]
gen_mnm$USICODE_stripped     = toupper(gen_mnm$V1)
gen_mnm$USICODE_stripped     = gsub('DE[0-9][0-9]','DE',gen_mnm$USICODE_stripped)
gen_mnm$USICODE_stripped     = gsub('ZIMD1','ZIM',gen_mnm$USICODE_stripped)
gen_mnm$USICODE_stripped     = gsub('AAA','A',gen_mnm$USICODE_stripped)
subinfo                      = gen_usi_mnm[c('USICODE_stripped','SCREENINGNUMMER')]
gen_mnm                      = merge(gen_mnm,subinfo, by='USICODE_stripped',all.x=T,all.y=T)
gen_mnm                      = merge(gen_mnm,spss_mnm[c('ID_initials_ngfn','sex','group')],by.x = 'SCREENINGNUMMER',by.y='ID_initials_ngfn',all.x = T,all.y=T)
gen_mnm$site                 = 'Mannheim'  
# cleaning names
names(gen_mnm)[grep('SCREENINGNUMMER',names(gen_mnm))] = 'ID'


# Bonn: using USI
# strip usi from control numbers
gen_usi_bnn                   = gen_usi_bn
gen_usi_bnn$USICODE_stripped  = gsub('DE[0-9][0-9]','DE',gen_usi_bnn$Lab_Code)
gen_usi_bnn$USICODE_stripped  = gsub('UKBD1','UKB',gen_usi_bnn$USICODE_stripped)
gen_usi_bnn$ADD_P             = as.character(gen_usi_bnn$ADD_P)
# put subject code on genetic data
gen_bnn                      = gen_ps[grep('ukb',gen_ps$V1),]
gen_bnn$USICODE_stripped     = toupper(gen_bnn$V1)
subinfo                      = gen_usi_bnn[c('USICODE_stripped','ADD_P')]
gen_bnn                      = merge(gen_bnn,subinfo, by='USICODE_stripped',all.x=T,all.y=T)
spss_bnn$ADD_P_stripped      = gsub('_[0-9]{3}','',spss_bnn$ADD_P)
gen_bnn                      = merge(gen_bnn,spss_bnn[c('ADD_P_stripped','sex','group')],by.x = 'ADD_P',by.y='ADD_P_stripped',all.x = T,all.y=T)
gen_bnn$site                 = 'Bonn'
# cleaning names
names(gen_bnn)[grep('ADD_P',names(gen_bnn))] = 'ID'

# complete genetic metainfo
gen_cmpl       = rbind(gen_bln,gen_mnm,gen_bnn)
gen_cmpl$sex   = agk.recode.c(gen_cmpl$sex,c('mÃ¤nnlich'),c('männlich'))
gen_cmpl$sex   = agk.recode.c(gen_cmpl$sex,c('männlich','weiblich','2','1'),c('m','f','m','f'))
gen_cmpl$group = agk.recode.c(gen_cmpl$group,c('gesund','alkoholabhängig','controls','patients','Kontrollgruppe','Patienten','0','102','patient'),
                              c('control','case','control','case','control','case','control','case','case'))
# work with characters
gen_cmpl       = data.frame(lapply(gen_cmpl,FUN = as.character),stringsAsFactors = F)

# print to file
setwd('S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik')
write.table(gen_cmpl,file = 'genetics_info_case_control_sex_site.csv',sep=';',row.names = F)

# explore and report Mannheim's problems
gen_mnm    = subset(gen_cmpl,site=='Mannheim')
cur_fun    = function(x) {return(any(is.na(x)))}
gen_mnm_pr = gen_mnm[apply(gen_mnm,cur_fun,MARGIN = 1),]
gen_mnm_pr[c('V2','V3','V4','V5','V6')] = NULL
names(gen_mnm_pr)[grep('V1',names(gen_mnm_pr))] = 'CHIP_USI'
names(gen_mnm_pr)[1] = 'ID_mapping_table'
names(gen_mnm_pr)[2] = 'USICODE_mapping_table'
write.table(gen_mnm_pr,file = 'genetics_info_Mannheim_problematic_subjects.csv',sep=';',row.names = F)
disp("Number of problematic genetic cases in Mannheim data before Vagelis's information")
print(length(gen_mnm_pr[,1]))

# get infos by Vagelis Zois to fill in info gaps
setwd('S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik/Zuordnung/Berlin_Mannheim')

# gen_mnm_corr gives us correct mappings from Mannheim ids to USIs
gen_mnm_corr = as.data.frame(read_excel('genetics_info_Mannheim_problematic_subjects_VagelisZois.xlsx'))

# correct the usicode_MA (this is the variable by Vagelis Zois)
gen_mnm_corr$USICODE_MA = gsub('DE[0-9][0-9]','DE',gen_mnm_corr$USICODE_MA)
gen_mnm_corr$USICODE_MA = gsub('MI1','M',gen_mnm_corr$USICODE_MA)
gen_mnm_corr$USICODE_MA = gsub('MB1','M',gen_mnm_corr$USICODE_MA)
gen_mnm_corr$USICODE_MA = gsub('MB1','M',gen_mnm_corr$USICODE_MA)
gen_mnm_corr$USICODE_MA = tolower(gen_mnm_corr$USICODE_MA)

# rename variables to align with names of gen_cmpl
source = names(gen_mnm_corr)
target = c('ID','USICODE_stripped','V1',"USICODE_MA","sex","group",'site')
names(gen_mnm_corr) = agk.recode(names(gen_mnm_corr),source,target)

# checks
stopifnot(isempty(grep('^22',gen_mnm$ID)))
stopifnot(isempty(grep('^22',gen_mnm_pr$ID_mapping_table)))

# recode
cur_map         = data.frame(ID_initials = inv$ID_initials, ID_initials_NGFN = inv$ID_initials_NGFN)
cur_map         = na.omit(cur_map)
gen_mnm_corr$ID = agk.recode(gen_mnm_corr$ID,cur_map$ID_initials,cur_map$ID_initials_NGFN)

# prep the CHIP_USI
gen_cmpl$CHIP_USI_stripped = gen_cmpl$V1
gen_cmpl$CHIP_USI_stripped = gsub('de[0-9][0-9]','de',gen_cmpl$CHIP_USI_stripped)
gen_cmpl$CHIP_USI_stripped = gsub('zimd1','zim',gen_cmpl$CHIP_USI_stripped)

# now use gen_mnm_corr to correct the problematic subjects
# CHIP_USI is what we need, because therein is the genetic data
for (cc in 1:length(gen_cmpl$CHIP_USI_stripped)) {
  # problem child row
  cur_chip_usi = gen_cmpl$CHIP_USI_stripped[cc]
  
  # does it have a cur_chip_usi?
  if (is.na(cur_chip_usi)) {next}
  
  # if yes then we check if vagelis sent us info on that person
  cur_ind = which(cur_chip_usi == gen_mnm_corr$USICODE_MA)
  
  # get the info by Vagelis
  if (length(cur_ind) == 1) {
    cur_pr = gen_cmpl[cc,]
    cur_cr = gen_mnm_corr[cur_ind,]
    for (pp in 1:length(cur_pr)) {
      if(is.na(cur_pr[pp])) {
        if (any(names(cur_cr) %in% names(cur_pr)[pp])) {
          cur_pr[pp] = cur_cr[names(cur_pr)[pp]]
        }
      }
    }
    gen_cmpl[cc,] = cur_pr
  } else {
    warning('not found')
  }
}

# recode again
gen_cmpl$group = agk.recode.c(gen_cmpl$group,c('gesund','alkoholabhängig','controls','patients','Kontrollgruppe','Patienten','0','102','patient'),
                              c('control','case','control','case','control','case','control','case','case'))

# now report again
gen_mnm    = subset(gen_cmpl,site=='Mannheim')
cur_fun    = function(x) {return(any(is.na(x)))}
gen_mnm_pr = gen_mnm[apply(gen_mnm,cur_fun,MARGIN = 1),]
gen_mnm_pr[c('V2','V3','V4','V5','V6')] = NULL
names(gen_mnm_pr)[grep('V1',names(gen_mnm_pr))] = 'CHIP_USI'
names(gen_mnm_pr)[1] = 'ID_mapping_table'
names(gen_mnm_pr)[2] = 'USICODE_mapping_table'
write.table(gen_mnm_pr,file = 'genetics_info_Mannheim_problematic_subjects.csv',sep=';',row.names = F)
disp("Number of problematic genetic cases in Mannheim data before Vagelis's information")
print(length(gen_mnm_pr[,1]))

# now write anew
gen_cmpl$CHIP_USI_stripped = NULL
setwd('S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/Genetik')
write.table(gen_cmpl,file = 'genetics_info_case_control_sex_site.csv',sep=';',row.names = F)




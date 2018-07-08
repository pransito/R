# read in the inventroy list to get infos on
# NGFN data

## paths ======================================================================
user            = paste0(as.character(Sys.info()["login"]),"/")
path_inventory  = 'S:/AG/AG-Emotional-Neuroscience/Restricted/NGFN/eMed/Inventarlisten/Inventarliste_Daten_all_Sites_22.12.2016.xls'
sites           = c('Mannheim','Bonn','Berlin')
#path_spss_mnh  = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Mannheim/Screen_Base_Neuro_SUMMENWERTE_ALL_FINAL_2013_Multicentre_gelöscht_23.01.2013.sav'
#path_spss_bnn  = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Bonn/NGFN13_Bonn_alle_Daten.sav'
path_spss_bln   = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Berlin/NGFNplus_Datenmaske_Berlin_05.03.13_KCH.sav'
path_relabs_bln = 'S:/AG/AG-Emotional-Neuroscience/Restricted/NGFN/eMed/TP6_TU Analysen/ABS_REL_Zuordnung_Transfer_TUserver_04.06.2014.xlsx'
path_relabs_mnh = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Mannheim/NGFN13_RÜCKFÄLLE_NEU_ALLE_20140521.sav'
base_lib        = paste("C:/Users/",user,"Google Drive/",sep="")
path_lib        = paste0(base_lib,"Library/R")

## settings ===================================================================
# what to predict? relapse or diagnosis
des_pred  = 'relapse'
# use rules or just take any where there is at least not "NA" everywhere
# LEAVE T
use_rules = T

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
inv$ID = paste0('ID__',inv$ID_initials,'__IDN__',inv$ID_initials_NGFN,'__CoID__',inv$`Co-ID`) 

# add scanner variable
inv$scanner = ifelse(inv$scanner == 'Berlin_scanner','Berlin_scanner_PTBnew',inv$scanner)
spss_bln = read.spss(path_spss_bln)

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
if (des_pred == 'relapse') {
  inv = subset(inv, group == 'patients')
  # get additional data (Berlin)
  relabs           = as.data.frame(read_excel(path_relabs_bln,sheet = 1,col_names = T))
  relabs           = (as.data.frame(lapply(relabs,trimws)))
  cur_fun          = function(x) {return(all(is.na(x)))}
  relabs           = relabs[unlist(lapply(relabs,cur_fun) == FALSE)]
  names(relabs)[1] = 'ID'
  str_fun = function(x) {
    cur_str = strsplit(as.character(x),'_NGFN_')
    cur_str = strsplit(cur_str[[1]][2],'_ALCUE')
    return(cur_str[[1]][1])
  }
  relabs$ID        = sapply(relabs$ID,FUN = str_fun)
  
  # get additional data (Mannheim)
  # relapse (heavy, 6 months after t0)
  relabs_mn                = read.spss(path_relabs_mnh)
  vars_of_int              = c('ID','ID_ngfn','ID_sfb','sRF_180Tage_t0')
  relabs_mn                = relabs_mn[vars_of_int]
  relabs_mn                = (as.data.frame(lapply(relabs_mn,trimws)))
  relabs_mn$sRF_180Tage_t0 = agk.recode.c(relabs_mn$sRF_180Tage_t0,c('ja','nein'),c('REL','ABS'))
  
  # fill it into REL ABS
  inv$relapse = inv$X6..follow.up
  disp('info on relpase before filling in from external source')
  print(table(inv$relapse %in% c('REL','ABS')))
  
  # Berlin
  for (ii in 1:length(inv$relapse)) {
    if (is.na(inv$relapse[ii])) {
      cur_id  = inv$Co.ID[ii]
      if (is.na(cur_id)) {
        next
      }
      cur_ind = grep(cur_id,relabs$ID)
      # fill in from external relabs
      if (!isempty(cur_ind)) {
        inv$relapse[ii] = relabs$X__4[cur_ind]
      }
    }
  }
  
  # Mannheim (here foing from two ID variables)
  for (ii in 1:length(inv$relapse)) {
    if (is.na(inv$relapse[ii])) {
      cur_id  = inv$ID_initials_NGFN[ii]
      if (is.na(cur_id)) {
        next
      }
      cur_ind = grep(paste0('^',cur_id,'$'),relabs_mn$ID_ngfn)
      # fill in from external relabs
      if (!isempty(cur_ind)) {
        inv$relapse[ii] = relabs_mn$sRF_180Tage_t0[cur_ind]
      } else {
        cur_id  = inv$ID_initials[ii]
        if (is.na(cur_id)) {
          next
        }
        cur_ind = grep(paste0('^',cur_id,'$'),relabs_mn$ID_sfb)
        # fill in from external relabs
        if (!isempty(cur_ind)) {
          inv$relapse[ii] = relabs_mn$sRF_180Tage_t0[cur_ind]
        }
      }
    }
  }
  
  disp('info on relpase after filling in from external source')
  print(table(inv$relapse %in% c('REL','ABS')))
  
  # recode the relapse variable
  #inv$relapse = as.character(inv$relapse)
  #inv$relapse = ifelse(is.na(inv$relapse),'NA',inv$relapse)
  
}

## TEST COMPLETENESS OF DATA SET FOR PREDICTION OF GROUP ======================
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

# first let's look how many people we have complete at t1
if (des_pred == 'relapse') {
  disp('Distribution of complete data at t1; NA means no t2 data.')
  print(table(inv$relapse,useNA = 'always'))
}

# REPORT inventory that is okay for:
# complete data MRI tasks and out-of MRI tasks; for prediction diagnosis
# with few cases where fmri data is irregular
disp(paste('complete data MRI tasks and out-of MRI tasks; for prediction', des_pred))
if (use_rules) {
  disp('with few cases where fmri data is irregular')
} else {
  disp('with few cases where fmri data is irregular and few where behav is irregular; at least nowhere NA')
}
print(table(sub_ok))

# who is in/out?
invIN  = inv[which(sub_ok),]
invOUT = inv[which(sub_ok==FALSE),]

# how many per location?
table(invIN$site)
xtabs(~site+group, data=invIN)
xtabs(~site+group+scanner, data=invIN)

# for case relapse
xtabs(~site+relapse, data=invIN)
xtabs(~site+relapse+scanner, data=invIN)

# how many dropout per location?
table(invOUT$site)

# in case of relapse we want to know who is DO just because no relapse data
int_vars = all_vars[c(3,4,9,10:16,17,18,19,20,21,22,23)]
rul_vars = list(c('controls','patients'),c('include'),c(),c(),c(),c(),
                c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
sub_ok_AD_but_no_relapse_data = test_rules(invOUT,int_vars,rul_vars)
table(sub_ok_AD_but_no_relapse_data)
invOUT_cmpl_but_no_rel_info = invOUT[sub_ok_AD_but_no_relapse_data,]
table(invOUT_cmpl_but_no_rel_info$site)




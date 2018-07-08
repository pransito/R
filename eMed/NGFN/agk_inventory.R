# read in the inventroy list to get infos on
# NGFN data

# packages
library(readxl)
library(plyr)
library(foreign)

path_inventory = 'S:/AG/AG-Emotional-Neuroscience/Restricted/NGFN/eMed/Inventarlisten/Inventarliste_Daten_all_Sites_22.12.2016.xls'
sites          = c('Mannheim','Bonn','Berlin')
#path_spss_mnh  = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Mannheim/NGFNplus_Datenmaske_Berlin_05.03.13_KCH.sav'
#path_spss_bnn  = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Bonn/NGFN13_Bonn_alle_Daten.sav'
path_spss_bln  = 'S:/AG/AG-Emotional-Neuroscience-Backup/NGFN/SPSS-Datenmasken/Datenmasken_FINALs/Berlin/NGFNplus_Datenmaske_Berlin_05.03.13_KCH.sav'

# inventories of three places
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

# make one ID variable
inv$ID = paste0('ID__',inv$ID_initials,'__IDN__',inv$ID_initials_NGFN,'__CoID__',inv$`Co-ID`) 

# add scanner variable
inv$scanner = ifelse(inv$scanner == 'Berlin_scanner','Berlin_scanner_PTBnew',inv$scanner)
spss_bln = read.spss(path_spss_bln)

for (ii in 1:length(inv$`Co-ID`)) {
  cur_id      = inv$`Co-ID`[ii]
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

# check if scanner data is there and correct

# REL ABS
# see S:\AG\AG-Emotional-Neuroscience\Restricted\NGFN\eMed\TP6_TU Analysen

tnl$VPPG      = as.character(tnl$VPPG)
tnl           = tnl[!is.na(tnl$VPPG),]
tnl$VPPG      = trimws(tnl$VPPG)
tnl$PhysioVP  = paste0("PhysioVP",tnl$PhysioVP)
tnl$PhysioVP  = trimws(tnl$PhysioVP)
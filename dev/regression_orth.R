# testing if serially orthogonalized design matrix leads to different betas
# than not
library(matlib)

# non orthogonal
cur_dat = data_pdt[data_pdt$subject == unique(data_pdt$subject)[2],]
cur_dat = cur_dat[!is.na(cur_dat$accept_reject),]
mod_n   = glm(accept_reject ~ gain+loss+ed_abs+cat,data = cur_dat,family = needed_family)

# non orthogonal with model matrix
mm       = model.matrix(accept_reject ~ gain+loss+ed_abs+cat,data = cur_dat)
mm       = mm[,-1]
mm       = data.frame(mm, accept_reject=cur_dat$accept_reject)
mod_mn   = glm(accept_reject ~ gain+loss+ed_abs+catgambling+catnegative+catpositive,
               data = mm,family = needed_family)

# orthogonal (have to kill intercept)
# https://cran.r-project.org/web/packages/matlib/vignettes/gramreg.html
mm       = model.matrix(accept_reject ~ gain+loss+ed_abs+cat,data = cur_dat)
mm       = mm[,-1]
mm_o     = gsorth(mm,recenter = T,rescale = T,adjnames = F)
mm_o     = data.frame(mm_o, accept_reject=cur_dat$accept_reject)
mod_o    = glm(accept_reject ~ gain+loss+ed_abs+catgambling+catnegative+catpositive,
               data = mm_o,family = needed_family)

# orthogonal (images come first) (does not make a big difference)
mm       = model.matrix(accept_reject ~ gain+loss+ed_abs+cat,data = cur_dat)
mm       = mm[,-1]
mm       = mm[,c(4,5,6,1,2,3)]
mm_o     = gsorth(mm,recenter = T,rescale = T,adjnames = F)
round(cor(mm_o),3)
mm_o     = data.frame(mm_o, accept_reject=cur_dat$accept_reject)
mod_oc   = glm(accept_reject ~ catgambling+catnegative+catpositive+gain+loss+ed_abs,
               data = mm_o,family = needed_family)

# function lmlist that orthogonalizes the predictors serially
lmList_orth = function(formula, data,family,na.action,pool) {
  # get the subject or ranef variable over which to perform lmlist
  cur_formula = deparse(substitute(formula))
  cur_raneff    = strsplit(cur_formula,'|',fixed = T)
  cur_ranef     = trimws(cur_raneff[[1]][2])
  cur_form_core = cur_raneff[[1]][1]
  y             = strsplit(cur_form_core,'~',fixed=T)
  y             = trimws(y[[1]][1])
  
  # make new data.frame which to use later
  # orthogonalize within each subject
  new_data = list()
  all_subs = unique(data[[cur_ranef]])
  
  for (ii in 1:length(all_subs)) {
    # get the data, kill NA's
    cur_dat = data[data[[cur_ranef]] == all_subs[ii],]
    cur_dat = cur_dat[!is.na(cur_dat$accept_reject),]
    
    # orthogonalize serially
    mm   = model.matrix(as.formula(cur_form_core),data = cur_dat)
    mm   = mm[,-1]
    mm_o = gsorth(mm,recenter = T,rescale = T,adjnames = F)
    mm_o = data.frame(mm_o, y=as.numeric(as.character(cur_dat[[y]])))
    
    # fix name of response variable
    names(mm_o)[length(names(mm_o))] = y
    
    # add the subject variable
    mm_o[[length(mm_o)+1]] = all_subs[ii]
    names(mm_o)[length(names(mm_o))] = cur_ranef
    
    
    new_data[[ii]] = mm_o
  }
  
  # unpack
  new_dat = new_data[[1]]
  for (ii in 2:length(new_data)) {
    new_dat = rbind(new_dat,new_data[[ii]])
  }
  
  # make a new formula
  all_preds = names(new_dat)
  all_preds = all_preds[-grep(y,all_preds)]
  all_preds = all_preds[-grep(cur_ranef,all_preds)]
  new_form  = paste0(y,' ~ ',paste(all_preds,collapse='+'))
  new_form  = as.formula(paste0(new_form,' | ', cur_ranef))
  
  # run the lmList and return
  lml = lmList(formula = new_form, data = new_dat,family = family,na.action = na.action,pool = pool)
  return(lml)
  
}



foo <- function(x) deparse(substitute(x))




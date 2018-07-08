## PREAMBLE ===================================================================
# library of Alexander Genauck
# contains all self-made functions which are useful for future use
# also uses all libraries which are handy to have loaded for my purposes

##    Libraries  ==============================================================
# get user
user              = paste0(as.character(Sys.info()["login"]),"/")

# library call function 
# tries to load library; if not possible, will download and install it
# input is package as string variable
agk.load.ifnot.install <- function(package_name){
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

# libraries
# documentation
# http://lme4.r-forge.r-project.org/
# http://lme4.r-forge.r-project.org/lMMwR/lrgprt.pdf
# http://www.biomedcentral.com/1471-2288/11/77
agk.load.ifnot.install("psych")
agk.load.ifnot.install("pracma")
agk.load.ifnot.install("pls")
agk.load.ifnot.install("Hmisc")
agk.load.ifnot.install("lme4")
agk.load.ifnot.install("reshape2")
agk.load.ifnot.install("R.matlab")
agk.load.ifnot.install("gtools")
agk.load.ifnot.install("plyr")
agk.load.ifnot.install("ggplot2")
agk.load.ifnot.install("plot3D")
agk.load.ifnot.install("rgl")
agk.load.ifnot.install("gridExtra")
agk.load.ifnot.install("kernlab")
agk.load.ifnot.install("boot")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install("fastICA")
agk.load.ifnot.install("corrplot")
agk.load.ifnot.install("glmnet")
# https://cran.r-project.org/src/contrib/Archive/glmnetUtils/
#tmp_path = paste0('C:/Users/',user,'Downloads/glmnetUtils-master.zip')
#install.packages(tmp_path, repos = NULL, type="source")
#glmnetUtils now on CRAN
agk.load.ifnot.install("glmnetUtils")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install("parallel")
agk.load.ifnot.install("foreach")
agk.load.ifnot.install("doSNOW")
agk.load.ifnot.install("simpleboot")
agk.load.ifnot.install("GPArotation")
agk.load.ifnot.install("nnet")
agk.load.ifnot.install("msm")
agk.load.ifnot.install("foreign")
agk.load.ifnot.install("caret")
agk.load.ifnot.install('readxl')
agk.load.ifnot.install("rJava")
agk.load.ifnot.install("xlsx")
agk.load.ifnot.install('e1071')
agk.load.ifnot.install('ptw')
agk.load.ifnot.install('lmPerm')
agk.load.ifnot.install('pROC')
agk.load.ifnot.install('cvTools')
agk.load.ifnot.install('matlib')
agk.load.ifnot.install('robust')
agk.load.ifnot.install('e1071')
#agk.load.ifnot.install('sampling')

library(compiler)

## General Purpose ============================================================
detach.package.rf <- function(pkg, character.only = FALSE)
  # recurrent force detaching of package
  # use package name as string
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE,force=TRUE)
  }
}

agk.clean.intercept.name = function(df) {
  # Function to clean the Intercept name
  names(df)[grep('(Intercept)',names(df))] = 'Intercept'
  return(df)
}

agk.merge.df.by.row.names = function(df1,df2,x1=NULL,x2=NULL) {
  # function to merge two dfs on row.names
  # you can pass additional merge variables by name
  # as char vectors in x1 and x2
  # row.names goes first
  
  # assertions
  stopifnot(length(x1) == length(x2))
  
  df1$'__IDRX__' = row.names(df1) 
  df2$'__IDRX__' = row.names(df2)
  
  # order
  df1 = df1[order(df1$`__IDRX__`),]
  df2 = df2[order(df2$`__IDRX__`),]
  
  mdf = merge(df1,df2,by.x = c('__IDRX__',x1),by.y = c('__IDRX__',x2))
  
  # post assertion
  stopifnot(row.names(df1) == mdf$`__IDRX__`)
  
  # row.names: put them back
  row.names(mdf) = mdf$`__IDRX__`
  mdf$`__IDRX__` = NULL
  
  # return
  return(mdf)
}

agk.mult.fun = function(x) {
  # function that multiplies
  # successively
  # like sum but with product
  res = x[1]
  for (ll in 2:length(x)) {
    res = res*x[ll]
  }
  return(res)
}

agk.in.dupl = function(x,y) {
  # function that works like the %in% operator but returns
  # indices and allows for duplicates
  # x is all subjects
  # y is the sampling of x
  # x is names (e.g. subject names)
  # y is names (e.g. subject names, that name the desired in x)
  # y can name subjects double and triple: the indeces will reflect that
  # x is the whole set
  # y is the subsample desired
  
  # get indeces
  x_ind = 1:length(x)
  
  # map
  return(as.numeric(agk.recode.c(y,x,x_ind)))
}

agk.lme.summary = function(model,type) {
  # will produce no approximation (default),
  # normal distr. approximation
  # to get p-values for lme4 models
  # model: lme4 model
  # type: 'none','norm'
  # 'satter' out cause anoying summary function, takes long
  # and does not allow REML
  # 'ken' also slow
  # https://en.wikipedia.org/wiki/Restricted_maximum_likelihood
  # value: returns the normal summary, and if type
  # is either 'norm', coefficient table with p-value
  
  if (type == 'none') {
    return(summary(model))
  } else if (type == 'norm') {
    cur_sum = summary(model)
    # extract coefficients
    coefs <- data.frame(coef(cur_sum))
    # use normal distribution to approximate p-value
    coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
    cur_sum$coefficients = coefs
    return(cur_sum)
  }
}

which.max.last = function(x) {
  # returns the last maximum
  cur_ind = which(max(x) == x)
  cur_ind = cur_ind[length(cur_ind)]
  return(cur_ind)
}

which.min.last = function(x) {
  # returns the last maximum
  cur_ind = which(min(x) == x)
  cur_ind = cur_ind[length(cur_ind)]
  return(cur_ind)
}

agk.mult.order = function(dms_d) {
  # function to order a data frame according to all
  # variables in the data frame starting at first
  # most left-hand
  all_str = c()
  for (ss in 1:length(dms_d)) {
    all_str[ss] = paste0(deparse(substitute(dms_d)),'[[',ss,']]')
  }
  cmd = paste('order(',paste(all_str,collapse = ','),')')
  cmd = paste0(deparse(substitute(dms_d)),'[',cmd,',]')
  dms_d_o = eval(parse(text=cmd))
  return(dms_d_o)
}

agk.best = function (x, metric, maximize) 
  # for caret train; pick the best metric after training
{
  bestIter <- if (maximize) 
    which.max.last(x[, metric])
  else which.min.last(x[, metric])
  bestIter
}


plot.jpeg = function(path, add=FALSE,r_start = 1, c_start = 1){
  # plot jpeg
  # source:
  # http://stackoverflow.com/questions/9543343/plot-a-jpg-image-using-base-graphics-in-r
  require('jpeg')
  # read the file
  jpg = readJPEG(path, native=T) 
  # get the resolution
  res = dim(jpg)[1:2]
  # initialize an empty plot area if add==FALSE
  if (!add) {plot(1,1,xlim=c(1,res[2]),ylim=c(1,res[1]),asp=1,type='n',xaxs='i',
                  yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')}
  rasterImage(jpg,c_start,r_start,c_start+res[2]-1,r_start+res[1]-1)
}

agk.sample.patch <- function(m,dimension, start_r,start_c) {
  # sample patch of size [x,y] of a 2-D matrix
  # you may indicate a starting point in row and column or put NA there to sample it randomly
  cur_dim <- dim(m)
  if (!is.na(start_r)) {
    start_r <- sample(1:(cur_dim[1]-(dimension[1]-1)),size = 1)
  }
  if (!is.na(start_c)) {
    start_c <- sample(1:(cur_dim[2]-(dimension[2]-1)),size = 1)
  }
  patch   <- m[start_r:(start_r+dimension[1]-1),start_c:(start_c+dimension[2]-1)]
  return(list(patch,(cur_dim[1]-(start_r+dimension[1]-1)+1),start_c))
}

agk.mean.quantile <- function(x,lower,upper) {
  # function to return the x% (given by lower and upper) quantile and the mean
  x_out <- quantile(x,probs = c(lower,upper))
  x_out <- c(mean(x),x_out)
  names(x_out) <- c("mean", "lower", "upper")
  return(x_out)
}

agk.mean.quantile.c = cmpfun(agk.mean.quantile)

agk.density_p = function(x,value,type = 'smaller') {
  # function to return the probability of value
  # given the density estimated based on x
  # default (one-sided 'smaller')/'greater':
  # one-sided test, that we see value or smaller/greater
  # so useful if you wanna test difference measure "x greater y"
  # two.sided:
  # is the one-sided test, that we see abs(value) or smaller
  # so useful if you wanna test difference measure "x different y"
  if(any(is.na(x))) {return(NA)}
  if (type == 'smaller' | type == 'greater') {
    dens_diffs   = density(x,n = 2^15)
    dens_diffs$y = dens_diffs$y/sum(dens_diffs$y)
    dens_diffs$y = cumsum(dens_diffs$y)
    f <- approxfun(dens_diffs, rule=2)
    if (type == 'smaller') {
      return(f(value))
    } else {
      return(1-f(value))
    }
  } else if (type == 'two.sided') {
    if (value == 0) {
      stop('When testing against value = 0 then you probably do not want to do a two-sided test, cause it will be n.s. anyways!')
    }
    x     = abs(x)
    value = abs(value)
    return(agk.density_p(x,value,type = 'greater'))
  } else {
    stop('Wrong input for "type".')
  }
}

agk.density_p.c = cmpfun(agk.density_p)

agk.whiten = function(dat) {
  # function to whiten data
  # provide a data frame with only
  # numeric variables
  # will return list
  # with whitened data
  # eigenvectors
  # eigenvalues
  dat     = scale(dat,center=T,scale=F)
  cur_cen = attr(dat,"scaled:center")
  cur_eig = eigen(cov(dat))
  cur_vec = cur_eig$vectors
  cur_val = diag(cur_eig$values^(-1/2))
  
  # whiten
  Z = dat%*%cur_vec%*%(cur_val)
  return(list(wdata = Z,vectors = cur_vec, values = cur_eig$values,centers = cur_cen,
              eigvalmatsqrt = cur_val))
}

agk.whiten.newdat = function(dat,whitener) {
  # performs a rotation of a vector or matrix (dat)
  # according to the information from a whitening pca
  # provide the output object from agk.whiten (whitener)
  dat               = scale(dat,center = whitener$centers,scale=F)
  cur_rot           = dat %*% whitener$vectors %*% whitener$eigvalmatsqrt
  return(cur_rot)
}

# delete trailing and leading white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# delete trailing white space
trim.trailing <- function (x) sub("\\s+$", "", x)

# delete leading white space
trim.leading <- function (x)  sub("^\\s+", "", x)

# function to plot the cov matrix of variables in df
# cur_range can be used to focus on a specific range of values of values to be color coded
# set NA if not needed
agk.plot.cov <- function(cur_df,cur_range) {
  p1 <- qplot(x=Var1, y=Var2, data=melt(cor(cur_df)), 
              fill=value, geom="tile")
  p1 <- p1 + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
  p1 <- p1 + ggtitle(paste("cov matrix",
                           "\n"," ",sep="")) +xlab('') + ylab('')
  if (!is.na(cur_range)) {
    p1 <- p1 + scale_fill_gradient(limits=cur_range)
  }
  return(p1)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

agk.unique.matrix = function(mat) {
  # function that drops redundant lines
  # line is redundant that has the same
  # set of elements as another one
  # starts from top
  # takes line compares via setequal
  # and if necessary drops the line
  # if a drop occurs is starts from top
  while (T) {
    cur_cut = c()
    for (ii in 1:length(mat[,1])) {
      cur_cmp = mat[ii,]
      cur_chk = apply(mat,FUN=setequal,y=cur_cmp,MARGIN = 1)
      cur_tru = which(cur_chk)
      cur_fix = which(cur_tru == ii)
      cur_cut = cur_tru[-cur_fix]
      
      if (length(cur_cut) != 0) {
        # cutting redundant lines
        mat = mat[-cur_cut,]
        # breaking
        break
      }
    }
    if (ii == length(mat[,1]) & (length(cur_cut) == 0)) {
      # break if ran through all lines and never a cut occurred
      break
    }
  }
  return(mat)
}

agk.scale_range <- function(x,C,D) {
  # function to scale data into a new range
  A <- min(x,na.rm = T)
  B <- max(x,na.rm = T)
  new_x <- C*(1-(x-A)/(B-A))+D*((x-A)/(B-A))
  return(new_x)
}

# function to play sound and scale before to -1 and 1
# assumes rate=8192
# play is commented out in script! 
# to let it run smoothly for knit
playsc <- function(x) {
  new_x <- scale_range(x,-1,1)
  play(new_x,rate=8192)
}

# function to plot scatterplot and marginal histogram (2 variables)
# source: http://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/
scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

# function to plot the cov matrix (?!?!) of variables in df
agk.heatmap <- function(mat) {
  p1 <- qplot(x=Var1, y=Var2, data=melt(mat), 
              fill=value, geom="tile")
  p1 <- p1 + scale_x_continuous(expand=c(0,0),breaks=NULL)
  p1 <- p1 + scale_y_continuous(expand=c(0,0),breaks=NULL)
  p1 <- p1 + xlab('') + ylab('')
  p1 <- p1 + theme(legend.position="none")
  return(p1)
}

# function to scale all variables in a df if numeric
agk.scale.ifpossible <- function(df,docenter = TRUE,doscale=TRUE) {
  for(ii in 1:length(df[1,])) {
    if (is.numeric(df[[ii]])) {
      cur_scale = scale(df[[ii]],center = docenter,scale = doscale)
      cur_attr  = attributes(cur_scale)
      df[ii] = cur_scale
      attributes(df[[ii]]) = cur_attr
    }
  }
  return(df)
}



# function to compute euclidean distance from point Q to line
# line is given by support point (some point on line) and vector indicating direction of line
# in the LA gamble matrix, the point c(26,12) lies on the diagonal, where c(gain, abs(loss))
# vec is c(2,-1)

vec = c(2,1,0)
sp  = c(26,13,0)

agk_get_ed <- function(Q,sp,vec) {
  # compute vec CROSS (Q-sp)
  vcd  = pracma::cross(vec,(Q-sp))
  vcd  = norm(as.matrix(vcd),type="F")
  
  # compute norm of vec
  nvec = norm(as.matrix(vec),type="F")
  
  # compute distance
  d    = vcd/nvec
  return(d)
}

agk_get_ed.c = cmpfun(agk_get_ed)

# function to exclude subjects from a data frame
# df ... data frame, idd_var... what is the sub variabale, subs_exc... char vector indicating the subjects
agk.exclude.subjects = function(df,id_var,subs_excl) {
  for (ii in 1:length(subs_excl)) {
    cur_cmd = paste("df = subset(df,", id_var,"!=subs_excl[ii])" )
    eval(parse(text=cur_cmd))
  }
  return(df)
}

lmList_orth = function(formula, data,family,na.action,pool) {
  # function lmList that orthogonalizes the predictors serially
  # within the subsetting variable (e.g. "subject")
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
    
    # orthogonalize serially, kill intercept if there
    mm   = model.matrix(as.formula(cur_form_core),data = cur_dat)
    if (colnames(mm)[1] == '(Intercept)') {
      no_int = 0
      mm   = mm[,-1]
    } else {
      no_int = 1
    }

    if (length(as.data.frame(mm)) > 1) {
      mm_o = gsorth(mm,recenter = T,rescale = T,adjnames = F)
    } else {
      mm_o = mm
    }
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
  if (no_int) {
    all_preds = c('0', all_preds)
  }
  new_form  = paste0(y,' ~ ',paste(all_preds,collapse='+'))
  new_form  = as.formula(paste0(new_form,' | ', cur_ranef))
  
  # run the lmList and return
  lml = lmList(formula = new_form, data = new_dat,family = family,na.action = na.action,pool = pool)
  return(lml)
  
}

agk.boot.cfint.mermod = function(mermod,num_cpus,num,fun_extract, cur_control,type = 'parametric') {
  # doing a parametric bootstrap of a lmer mod
  # interesting cause bootmer does not seem to accept control arguments
  # function simply bootstraps parametrically the fixef extracts according to
  # fun_extract
  # it returns the naked bootstraps which then can be used to make the cfints
  
  # get the function call
  cur_form                                = mermod@call
  cur_form_n                              = names(cur_form)
  cur_form                                = as.character(mermod@call)
  cur_form[which(cur_form_n %in% 'data')] = 'cur_d'
  new_form                                = paste(cur_form[1],'(')
  for (ff in 2:length(cur_form)) {
    if (ff == length(cur_form)) {final_str = ')'} else {final_str = ','}
    new_form = paste(new_form,paste(cur_form_n[ff],cur_form[ff],sep=' = '), final_str)
  }
  
  mod_call = new_form
  
    # sub function
    agk.boot.cfint.mermod.subfun = function(mermod,mod_call,fun_extract,cur_control,type) {
      cur_d = mermod@frame
      if (type == 'parametric') {
        # simulate under alternative hyothesis
        cur_d$accept_reject <- (simulate(mermod))[[1]]
      } else if (type == 'non-parametric') {
        # resample the data with relacement
        cur_d = cur_d[sample(1:length(cur_d[,1]),replace = T),]
      } else {
        stop('unknown bootstrap type!')
      }
      # fit under alternative hypothesis
      cur_fit <- eval(parse(text=mod_call))
      return(fun_extract(cur_fit))
    }

  ## running the bootstrap (no cov)
  cl<-parallel::makeCluster(num_cpus)
  registerDoSNOW(cl)
  effects_under_ha <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data_pdt","cur_control")) %dopar% {
    agk.boot.cfint.mermod.subfun(mermod,mod_call,fun_extract,cur_control,type)
  }
  return(effects_under_ha)
}

agk.boot.p.mermod = function(mermod,mermod0=NULL,num_cpus,num,fun_extract,cur_control,type,permvars=NULL) {
  # doing a parametric bootstrap of p-values against a 0-hypothesis mermod0
  # interesting cause bootmer does not seem to accept control arguments
  # function simply bootstraps parametrically the fixef extracts according to
  # fun_extract
  # it returns the naked bootstraps which then can be used to make the p-values
  # OR:
  # bootstrapping a distribution under 0 by using permutation
  # it returns the naked bootstraps which then can be used to make the p-values
  # needs a single model for which the permutation should be computed
  # permvars indicates all variables that should be scrambled
  # e.g. c('group','category'); clean nested model comparisons usually call
  # for permvars of length 1
  
  # get the function call
  cur_form                                = mermod@call
  cur_form_n                              = names(cur_form)
  cur_form                                = as.character(mermod@call)
  cur_form[which(cur_form_n %in% 'data')] = 'cur_d'
  new_form                                = paste(cur_form[1],'(')
  for (ff in 2:length(cur_form)) {
    if (ff == length(cur_form)) {final_str = ')'} else {final_str = ','}
    new_form = paste(new_form,paste(cur_form_n[ff],cur_form[ff],sep=' = '), final_str)
  }
  mod_call = new_form
  
  # sub function
  agk.boot.p.mermod.subfun = function(mermod,mermod0,mod_call,fun_extract,cur_control,type,permvars) {
    cur_d = mermod@frame
    if (type == 'parametric') {
      # simulate under 0 hyothesis
      cur_d$accept_reject = (simulate(mermod0))[[1]]
    } else if (type == 'perm') {
      # scramble all variables in permvars
      for (pp in 1:length(permvars)) {
        cur_d[[permvars[pp]]] = sample(cur_d[[permvars[pp]]])
      }
    } else {
      stop('unknown bootstrap type!')
    }
    
    # fit under alternative hypothesis
    cur_fit             = eval(parse(text=mod_call))
    return(fun_extract(cur_fit))
  }
  
  # running the bootstrap
  cl<-parallel::makeCluster(num_cpus)
  registerDoSNOW(cl)
  effects_under_0 <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("cur_control")) %dopar% {
    agk.boot.p.mermod.subfun(mermod,mermod0,mod_call,fun_extract,cur_control,type,permvars)
  }
  return(effects_under_0)
}

## PDT Physio study ===========================================================
## summary and matching function
# needed functions
# LEAVE THIS TO MEAN FOR NOW (Then it aligns with everything);
# ITS USE IS NOT FLEXIBLY and THOROUGHLY ENOUGH
f.summary            = function(x) mean(x, na.rm=TRUE)

# ROC curves
agk.plot.mean.roc = function(rocl,des_ci = 'quant',add = F) {
  # function that takes a list of ROCs
  # smoothes them
  # and computes a mean + CI ROC
  # plots it
  
  # make new ROC objects
  for (ii in 1:length(rocl)) {
    if (var(rocl[[ii]]$original.predictor) == 0) {
      cur_obj = list()
      cur_obj$sensitivities = seq(0,1,length.out = 502)
      cur_obj$specificities = seq(0,1,length.out = 502)
      rocl[[ii]]            = cur_obj
    } else {
      predictor  = rocl[[ii]]$original.predictor
      response   = rocl[[ii]]$original.response
      predictor  = agk.scale_range(predictor,0,1)
      rocl[[ii]] = roc(response,predictor)
      rocl[[ii]] = smooth(rocl[[ii]],n = 500)
    }
  }
  
  # ROC curve (get all the coordinates)
  specificities  = t(as.matrix(rocl[[1]]$specificities))
  sensitivities  = t(as.matrix(rocl[[1]]$sensitivities))
  
  # get all the ROC curves
  for (ii in 2:length(rocl)) {
    if (length(rocl[[ii]]$specificities) != length(specificities[ii-1,])) {
      stop('length not same!')
    }
    
    specificities = rbind(specificities,t(as.matrix(rocl[[ii]]$specificities)))
    sensitivities = rbind(sensitivities,t(as.matrix(rocl[[ii]]$sensitivities)))
  }
  
  # get CI
  specificities  = as.data.frame(specificities)
  sensitivities  = as.data.frame(sensitivities)
  
  if (des_ci == 'ci') {
    spec_mci  = lapply(specificities,FUN = agk.boot.ci.c,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
    sens_mci  = lapply(sensitivities,FUN = agk.boot.ci,R=1000,cur_fun = mean,lower = 0.025,upper=0.975)
  } else {
    # not ci of mean but mean and percentiles over CV rounds
    spec_mci  = lapply(specificities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
    sens_mci  = lapply(sensitivities,FUN = agk.mean.quantile.c,lower = 0.025,upper=0.975)
  }
  
  spec_mci     = as.data.frame(matrix(unlist(spec_mci),ncol = 3,byrow = T))
  sens_mci     = as.data.frame(matrix(unlist(sens_mci),ncol = 3,byrow = T))
  
  # plot
  if (add == F) {
    plot(spec_mci$V1[order(spec_mci$V1)],sens_mci$V1[order(spec_mci$V1)],xlim = c(1,0),type='l',lty=1,
         xlab = 'specificity', ylab = 'sensitivity',col='blue',lwd=4)
    lines(spec_mci$V2[order(spec_mci$V1)],sens_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
    lines(spec_mci$V3[order(spec_mci$V1)],sens_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='blue',lty=1)
  } else if (add == T) {
    # add to existing figure
    lines(spec_mci$V1[order(spec_mci$V1)],sens_mci$V1[order(spec_mci$V1)],xlim = c(1,0),type='l',lty=2,lwd=4,col='red')
    lines(spec_mci$V2[order(spec_mci$V1)],sens_mci$V2[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
    lines(spec_mci$V3[order(spec_mci$V1)],sens_mci$V3[order(spec_mci$V1)],xlim = c(1,0),col='red',lty=2)
  }
  
  abline(a=1,b=-1,lty=4,lwd=2)
}

ttest.group          = function(x,cur_group,cur_alt = 'two.sided') {
  if (!is.numeric(x)){return(NA)}
  level_1 = levels(cur_group)[1]
  level_2 = levels(cur_group)[2]
  x_1     = x[cur_group == level_1]
  x_2     = x[cur_group == level_2]
  return(t.test(x_1,x_2,alternative = cur_alt))
}

ttest.group.boot = function(x,cur_group) {
  # function to perform a two-sample t-test
  # bootstrap the p-value
  # returns the density object which includes the p-value
  
  # get the data
  level_1 = levels(cur_group)[1]
  level_2 = levels(cur_group)[2]
  x_1     = x[cur_group == level_1]
  x_2     = x[cur_group == level_2]
  
  # make two data frames d_1 and d_2
  # missings replaced by group mean
  cur_var = x_1
  cur_var[is.na(cur_var)] = mean(cur_var,na.rm = T)
  ind     = 1:length(cur_var)
  d_1     = data.frame(ind,cur_var)
  cur_var = x_2
  cur_var[is.na(cur_var)] = mean(cur_var,na.rm = T)
  ind     = 1:length(cur_var)
  d_2     = data.frame(ind,cur_var)
  
  cur_fun  = function(d,indices) {
    cur_d  = d[indices,]
    mean_1 = mean.rmna(cur_d$cur_var)
    if (is.na(mean_1)) {
      stop("mean is NA despite using mean.rmna")
    }
    return(mean_1)
  }
  tmp_1  = boot(d_1,cur_fun,R=1500)
  tmp_2  = boot(d_2,cur_fun,R=1500)
  tmp_t0 = tmp_1$t0 - tmp_2$t0
  tmp_t  = tmp_1$t[,1] - tmp_2$t[,1]
  #dens   = density(tmp_t,bw = c("ucv"),n = 1000)
  dens   = density(tmp_t,n = 2000)
  dens$y = dens$y/sum(dens$y)
  dens$y = cumsum(dens$y)
  f <- approxfun(dens, rule=2)
  dens$f = f
  if (tmp_t0 >= 0) {
    dens$p.value = f(0)
  } else {
    dens$p.value = 1-f(0)
  }
  
  # collecting the estimated stats
  dens$med_1 = tmp_1$t0
  dens$med_2 = tmp_2$t0
  
  return(dens)
}

f.difftest    = function(df,cur_var,cur_group) {
  # function to check if sig. group differences there 
  # for checking on group matching
  if (cur_group != "") {
    # different functions to compute it
    if(match_boot == 0) {
      tmp_t       = ttest.group(df[cur_var][[1]]+randn(1,length(df[cur_var][[1]]))*0.000001,df[cur_group][[1]])
      med_1       = round(tmp_t$estimate[[1]],digits = 2)
      med_2       = round(tmp_t$estimate[[2]],digits = 2)
      param_names = c("cur_var", levels(df[cur_group][[1]])[1],"HC_se",levels(df[cur_group][[1]])[2],"PG_se","cur_se","t_p-value","t_matched","k_p-value","k_matched", "matched")
    } else if (match_boot == 1) {
      # boot case
      tmp_t       = ttest.group.boot(df[cur_var][[1]]+randn(1,length(df[cur_var][[1]]))*0.000001,df[cur_group][[1]])
      med_1       = round(tmp_t$med_1,digits = 2)
      med_2       = round(tmp_t$med_2,digits = 2)
      param_names = c("cur_var", levels(df[cur_group][[1]])[1],"HC_se",levels(df[cur_group][[1]])[2],"PG_se","cur_se","t_p-value_boot","t_matched_boot","k_p-value","k_matched", "matched")
    } else if (match_boot == 2) {
      # permutation case
      # next line only for quickly getting the means
      tmp_t         = ttest.group(df[cur_var][[1]]+randn(1,length(df[cur_var][[1]]))*0.000001,df[cur_group][[1]])
      med_1         = round(tmp_t$estimate[[1]],digits = 2)
      med_2         = round(tmp_t$estimate[[2]],digits = 2)
      # permutation test
      cur_form      = as.formula(paste(cur_var,'~',cur_group))
      tmp_t         = summary(lmp(cur_form,data=df,maxIter=1000000,nCycle=1,settings=F,Ca = 0.000001))
      tmp_t$p.value = tmp_t$coefficients[2,3]
      param_names = c("cur_var", levels(df[cur_group][[1]])[1],"HC_se",levels(df[cur_group][[1]])[2],"PG_se","cur_se","t_p-value_perm","t_matched_perm","k_p-value","k_matched", "matched")
    }
    
    # non-metric test
    tmp_k       = kruskal.test(as.numeric(df[cur_var][[1]]+randn(1,length(df[cur_var][[1]]))*0.000001),df[cur_group][[1]])
    
    # boot the se's
    first_bt    = one.boot(df[cur_var][[1]][df[cur_group][[1]] == levels(df[cur_group][[1]])[1]],FUN=f.summary,R=3000)
    secnd_bt    = one.boot(df[cur_var][[1]][df[cur_group][[1]] == levels(df[cur_group][[1]])[2]],FUN=f.summary,R=3000)
    HC_se       = round(sd(first_bt$t,na.rm = T),digits=2)
    PG_se       = round(sd(secnd_bt$t,na.rm = T),digits=2)
    cur_se      = round(mean(c(sd(first_bt$t,na.rm = T),sd(secnd_bt$t,na.rm = T))),digits=2)
    
    # decide if matched
    t_matched   = ifelse(tmp_t$p.value < m_crit,"NO","YES")
    k_matched   = ifelse(tmp_k$p.value < m_crit,"NO","YES")
    
    # different rules for successful matching
    if(match_boot == 0) {
      # normal t.test case (both t and k must be ok)
      matched     = ifelse((t_matched == "YES" & k_matched == "YES"),"YES","NO")
    } else if (match_boot == 1 | match_boot == 2) {
      # boot or perm case
      matched     = ifelse((t_matched == "YES"),"YES","NO")
    }
    
    # pack the stuff
    cur_var     = as.data.frame(cur_var)
    tmp         = data.frame(cur_var,med_1,HC_se,med_2,PG_se,cur_se,tmp_t$p.value,t_matched,tmp_k$p.value,k_matched,matched)
    names(tmp)  = param_names
    return(tmp)
  } else {
    # single group case
    tmp_t       = t.test(df[cur_var][[1]]+randn(1,length(df[cur_var][[1]]))*0.000001)
    med_1       = round(f.summary(df[cur_var][[1]]),digits = 2)
    first_bt    = one.boot(df[cur_var][[1]],FUN=f.summary,R=3000)
    cur_se      = sd(first_bt$t,na.rm = T)
    param_names = c("cur_var","HC","cur_se","t_p-value")
    cur_var     = as.data.frame(cur_var)
    tmp         = data.frame(cur_var,med_1,cur_se,tmp_t$p.value)
    names(tmp)  = param_names
    return(tmp)
  }
}

# agk.boot.t.test = function(x1,x2) {
#   # function to perform a bootstrapped
#   # two sample t-test
# }

## split the data set into pilot study and sanity study
agk_get_group <- function(x) {
  # function to get from the PhysioVP number the correct grouping variable
  
  # additional grouping information (POSTPILOT)
  # HC group
  hcgroup = c(101)
  pggroup = c(102,103,104,105)
  
  if(length(strfind(as.character(x),"Pretest"))!=0) {
    return(99)
  } else {
    tmp_number = strsplit(as.character(x),split = "VP")
    tmp_number = as.numeric(tmp_number[[1]][2])
  }
  
  # specific
  if (sum(tmp_number == hcgroup)) {
    return (1)
  } else if (sum(tmp_number == pggroup)) {
    return (3)
  }
  
  # general
  if (tmp_number > 60) {
    return (3)
  } else if ((tmp_number > 40) & (tmp_number < 60)) {
    return (2)
  } else if (tmp_number < 35) {
    return (1)
  }
  
}

## MACHINE INTELLIGENCE =======================================================
# Oja's rule PCA
# extracts 1st PC via Hebbian learning
# takes in an observation vector of length n;
# first element of obs vector will be assumed as out
# all others as inp
# inp will be prlonged with first element
# init_w must be of same length as x
agk.pca.oja <- function(x, w_init, cur_eta) {
  # linear combination
  outp <- w_init%*%x
  inp  <- x
  w_new <- w_init + (cur_eta*outp)%*%(inp-outp*w_init)
  return(w_new)
}

agk.batch.oja <- function(X, w_init,cur_eta) {
  w_coll <- c()
  w_new <- agk.pca.oja(as.numeric(X[1,]), w_init, cur_eta)
  w_coll <- rbind(w_coll,as.numeric(w_new))
  for(ii in 2:length(X[,1])) {
    w_new <- agk.pca.oja(as.numeric(X[ii,]), as.numeric(w_new), cur_eta)
    w_coll <- rbind(w_coll,as.numeric(w_new))
  }
  return (w_coll)
}

# make a cvTools thingy sub function
agk.get.perm.vec = function(whichv,cur_fold_id) {
  # finds the subset vector to satisfy:
  # whichv[perm_vec] = cur_fold_id
  perm_vec = rep(NA,length(cur_fold_id))
  for (ii in 1:length(perm_vec)) {
    # desired number
    cur_des = cur_fold_id[ii]
    # look for desired number
    for (jj in 1:length(whichv)) {
      if(is.na(whichv[jj])) {next}
      if (whichv[jj] == cur_des) {
        # success
        perm_vec[ii] = jj
        whichv[jj] = NA
        break
      }
    }
  }
  return(perm_vec)
}

agk.cv.glm.force = function(data,glm.fit,cost,K,rec_ct = 1) {
  # function to run cv.glm it again if error occurs
  # rec_ct is the recursion count
  if (rec_ct > 100) {stop('Too many recursions when trying cv.glm!')}
  out = tryCatch(
    # work around for bug
    # https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R 
    {
      out = cv.glm(data,glm.fit,cost,K)
    },
    error=function(cond) {
      message("Had to retry again!")
      # call the function again
      rec_ct = rec_ct + 1
      agk.cv.glm.force(data,glm.fit,cost,K,rec_ct)
    }
  )
  return(out)
}

# function to plot histogram and the kernel
# estimated marginal density of the variables
# takes in name of data.frame to be used 
# from global env
# and character vector with names of variables
# returns a list with the plots
# CODE ADAPTED FROM #
# http://stackoverflow.com/questions/5033240/plot-probability-with-ggplot2-not-density
agk.plot.density <- function(name_df, vars,subtitle) {
  plots <- list()
  for (ii in 1:length(vars)) {
    base_expr <- paste('plots[[ii]] <- ggplot(',
                       name_df,',aes(x=',
                       vars[ii],'))',sep='')
    eval(parse(text=base_expr))
    hist_expr <- paste('plots[[ii]] <- plots[[ii]] +', 
                       'geom_histogram(aes(y = ..density..), binwidth=density(name_df$',
                       vars[ii],',bw="ucv")$bw, color="white")',sep='')
    eval(parse(text=hist_expr))
    plots[[ii]] <- plots[[ii]] + geom_density(fill="red", alpha = 0.2) +
      theme_bw() +
      xlab(vars[ii]) +
      ylab('density')
    plots[[ii]] <- plots[[ii]] + ggtitle(paste("marginal density of", vars[ii],"\n",
                                               subtitle,                                               
                                               " "))
  }
  return(plots)
}

# simpler agk.plot density " NOT DONE"
agk.plot.density.simple <- function(x,kern,width=0) {
  cur_df <- data.frame(x=x)
  plots <- list()
  cur_plot <- ggplot(cur_df,aes(x=x))
  if (width == 0) {
    cur_plot <- cur_plot + geom_histogram(aes(y = ..density..), binwidth=density(cur_df$x,bw="ucv")$bw, color="white")
  } else {
    cur_plot <- cur_plot + geom_histogram(aes(y = ..density..), binwidth=width, color="white")
  }
  cur_plot <- cur_plot + geom_density(fill="red", alpha = 0.2) + theme_bw() + xlab("x") + ylab('density')
  cur_plot <- cur_plot + ggtitle(paste("density of x"))
  return(cur_plot)
}

## Kernel-PCA functions
# function to center a kernel matrix
agk.center.kmatrix <- function(K) {
  Kc <- matrix(NA,nrow=length(K[,1]),ncol=length(K[1,]))
  p  <- length(K[,1])
  last_term <- (1/(p^2))*sum(K[,])
  for (ii in 1:length(K[,1])) {
    for (jj in 1:length(K[1,])) {
      Kc[ii,jj] <- K[ii,jj] - (1/p)*sum(K[,jj]) - (1/p)*sum(K[ii,]) + last_term 
    }
  }
  return(Kc)
}

# function to project new data point
# onto ev (i.e. pc) from kpca
# returns a scalar
agk.proj.on.kpca.ev <- function(ev,dat,new_pt,sigm) {
  rbf <- rbfdot(sigma = sigm)
  # create kernel matrix w.r.t. new data point
  kern_mat_new <- kernelMatrix(rbf,dat,new_pt)
  # center the matrix
  kern_mat_new <- agk.center.kmatrix(kern_mat_new)
  u <- c()
  for (ii in 1:length(kern_mat_new[,1])) {
    u[ii] <- ev[ii]*kern_mat_new[ii,1]
  }
  u_sum <- sum(u)
  return(u_sum)
}

# function to normalize an "a"-eigenvector matrix
# i.e. "a" is a vector which is representings an eigenvector with coeffieciencts 
# for linearly combining data points
# e = sum_k(a_k*datpt_k); to make norm(e) = 1, apply this function to "a"-vector
# takes in cur_eigen; lsit with [[1]] matrix of eigenvalues;
# [[2]] is matrix of "a" - eigenvectors
# returns a matrix with eigv which are normalized
agk.norm.a <- function(cur_eigen) {
  p      <- length(cur_eigen[[2]][1,])
  norm_a <- matrix(nrow=p,ncol=p)
  for (ii in 1:p) {
    cur_a     <- cur_eigen[[2]][,ii]
    cur_l     <- cur_eigen[[1]][ii]
    cur_norma <- (2/(sqrt(cur_l)*sqrt(cur_a%*%cur_a)))*cur_a
    norm_a[,ii] <- as.numeric(cur_norma)
  }
  norm_a <- as.matrix(norm_a)
  return(norm_a)
}

# function for kpca to visualize first n PCs; assumes rbf kernel;
# assumes 2-column data
agk.plot.kpca <- function(res, range,cur_eigen,dat,sigm,n_pcs,dat_kpca,sub_title) {
  # prep grid
  cur_grid <- matrix(nrow=res,ncol=res)
  x1_steps <- seq(range[1],range[2],length.out = res)
  x2_steps <- seq(range[1],range[2],length.out = res)
  
  pc_grids <- list()
  for (ii in 1:n_pcs) {
    for (jj in 1:length(x1_steps)) {
      for (kk in 1:length(x2_steps)) {
        cur_grid[jj,kk] <- agk.proj.on.kpca.ev(cur_eigen[[2]][,ii],
                                               dat,
                                               t(as.matrix(
                                                 c(x1_steps[jj],x2_steps[kk]))),
                                               sigm = sigm)
      }
    }
    pc_grids[[ii]] <- cur_grid
  }
  
  # plotting
  a <- as.matrix(rep(x1_steps,times = res))
  b <- as.matrix(rep(x2_steps,each = res))
  all_grids <- data.frame()
  for (ii in 1:length(pc_grids)) {
    u  <- as.numeric(pc_grids[[ii]])
    pc <- rep(paste("PC",as.character(ii)),length(u))
    tmp <- data.frame(a,b,u,pc)
    all_grids <- rbind(all_grids,tmp)
  }
  
  dat_kpcas <- data.frame(dat,dat_kpca$group)
  names(dat_kpcas) <- c("x1s","x2s","group")
  
  p3 <- ggplot(all_grids,aes(x=a,y=b,fill=u)) + geom_tile()
  p3 <- p3 + scale_fill_gradient(low = "black", high = "white")
  p3 <- p3 + facet_wrap(~ pc,nrow=2)
  p3 <- p3 +  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  p3 <- p3 + geom_point(data = dat_kpcas,
                        mapping = aes(x=x1s,y=x2s,col=group,fill=NULL),
                        alpha = 0.8, size =1)
  
  p3 <- p3 + ggtitle(paste("Visualization of first",
                           n_pcs, "features of data set based on rbf kpca",
                           sub_title))
  return(p3)
}

## ICA ##
## ICA using natural gradient with B-S-source ampl stand.
# obs_dis must be a data.frame/matrix with observations per columns
agk.ICA.nat.grad <- function(obs,obs_dis,W_init,eta_init,w_0_init) {
  W <- W_init
  
  # initializations (learning rate, bias w_0)
  eta   <- eta_init
  w_0   <- w_0_init 
  count <- 0
  ns    <- length(obs_dis[,1])
  
  for (ii in 1:length(obs_dis[1,])) {
    cur_eta <- eta/ii
    cur_pt  <- as.matrix(obs_dis[,ii])
    u       <- W%*%cur_pt + w_0
    y       <- sigmoid(u,a=1,b=0)
    d_W     <- cur_eta*(W + (ones(ns,1) - 2*y)%*%t(cur_pt)%*%(t(W)%*%W))
    d_w_0   <- cur_eta*((ones(ns,1) - 2*y))
    
    # normalize amplitude
    d_W <- d_W-diag(diag(d_W))
    
    # new W and w_0
    W   <- W + d_W
    w_0 <- w_0 + d_w_0
    
  }
  
  # test the estimated unmixing matrix
  est_sources_n <- W%*%obs
  
  return(est_sources_n)
  
}

## simulated annealing

# function to compute local energy
# given a set of nodes which can have states {-1,1} (S)
# given a symmetric weights matrix with diag = 0 (W_0)
# given an index to identify the node in question (cur_i)
agk.local.energy <- function(S_0,W_0,cur_i) {
  loc_en <- -0.5*sum((S_0[cur_i]*S_0)*W_0[cur_i,])
  return(loc_en)
}

# function to compute mean-fields
agk.mn.fields <- function(S_0,W_0,cur_i) {
  mn_fields <- sum(W_0[cur_i,]*S_0)
  return(mn_fields)
}

# function to compute the energy of the complete state
agk.EofS <- function(W_0,S_0) {
  EofS <- -0.5*sum(W_0*(S_0 %*% t(S_0)))
  return(EofS)
}

# The probability that the network is in a state S with energy E(S) 
agk.PofS <- function(E_S,beta, all_E_S) {
  Z    <- sum(exp(-beta*all_E_S))
  PofS <- (1/Z)*exp(-beta*E_S)
  return(PofS)
} 

# The probability that the network is in a state S with energy E(S) 
agk.PofS <- function(E_S,beta, all_E_S) {
  Z    <- sum(exp(-beta*all_E_S))
  PofS <- (1/Z)*exp(-beta*E_S)
  return(PofS)
} 

## K-MEANS clustering ##
# function to assign every data point to nearest prototype
# dat goes in as df with every row as one obs
# prot goes in as every col is one cluster mean
# uses absolute distance
# m is assignment vector (1 is assigned, 0 if not)
# m must sum to 1 and has length of number of clust
agk.kmeans.assign <- function(dat,prot) {
  # compute euclidean distances
  cur_dist <- distmat(t(prot), as.matrix(dat))
  cur_dist <- as.data.frame(cur_dist)
  
  # make Matrix of vectors which indicate assignment
  # also a vector m which indicates the cluster assingment by number
  M <- data.frame(matrix(0,length(prot[1,]),length(dat[,1])))
  m <- c()
  for (ii in 1:length(cur_dist[1,])) {
    M[ii] <- ifelse(cur_dist[ii] == min(cur_dist[ii]),1,0)
    m[ii] <- which.min(as.numeric(cur_dist[[ii]]))
  }
  return(list(M=M,m=m,cur_dist=cur_dist))
} 


# squashes a vector s.t. vector will sum to 1
# uses beta (inverse temperature) as measure of fuzziness
agk.softmax <- function(x,cur_beta) {
  squashed_v <- c()
  sum_v      <- sum(exp(cur_beta*x))
  for(ii in 1:length(x)) {
    squashed_v[ii] = exp(cur_beta*x[ii])/sum_v
  }
  return(squashed_v)
}

## SOM Maps ##
# make a 2D-grid
agk.make.2Dgrid <- function(k){
  x1 <- rep(seq(1:k),k)
  x2 <- rep(seq(1:k),each=k)
  cur_grid <- data.frame(x1,x2)
  return(cur_grid)
}

# get connections of 2D grid (map space)
# returns a data.frame with coord. for line segments
# and the indices data.frame to know which point to connect with which
agk.connect <- function(cur_grid) {
  cur_dist <- distmat(as.matrix(cur_grid),as.matrix(cur_grid))
  smallest <- t(apply(cur_dist, 1,order)[1:3, ])
  smallest <- smallest[,-1]
  
  # get line segments
  ct        <- 0
  line_segm <- matrix(NA,nrow = 4*length(cur_grid[,1]),ncol = 2)
  for(ii in 1:dim(cur_grid)[1]) {
    ct <- ct+1
    line_segm[ct,] <- as.numeric(cur_grid[ii,])
    ct <- ct +1
    line_segm[ct,] <- as.numeric(cur_grid[smallest[ii,1],])
    ct <- ct+1
    line_segm[ct,] <- as.numeric(cur_grid[ii,])
    ct <- ct +1
    line_segm[ct,] <- as.numeric(cur_grid[smallest[ii,2],])
  }
  return(list(line_segm=line_segm,smallest=smallest))
}

# function to connect points in data space
# given the connection rule from map space
# given through "smallest" matrix
# which gives for every point the indices of
# the two points which it needs to connect to
agk.connect.smallest <- function(cur_grid,smallest) {
  # get line segments
  ct        <- 0
  line_segm <- matrix(NA,nrow = 4*length(cur_grid[,1]),ncol = (dim(cur_grid))[2])
  for(ii in 1:dim(cur_grid)[1]) {
    ct <- ct+1
    line_segm[ct,] <- as.numeric(cur_grid[ii,])
    ct <- ct +1
    line_segm[ct,] <- as.numeric(cur_grid[smallest[ii,1],])
    ct <- ct+1
    line_segm[ct,] <- as.numeric(cur_grid[ii,])
    ct <- ct +1
    line_segm[ct,] <- as.numeric(cur_grid[smallest[ii,2],])
  }
  return(line_segm)
}

# initialize 2D-map space
agk.init.2Dmap <- function(cur_k,cex_v) {
  cur_grid <- agk.make.2Dgrid(cur_k)
  plot(cur_grid)
  points(cur_grid,col="red",pch=21,cex=3*cex_v,bg =rgb(1,0,0))
  text(x=cur_grid[,1],cur_grid[,2],labels=as.character(1:length(cur_grid[,1])))
  cur_con <- agk.connect(cur_grid)
  for(kk in seq(1,length(cur_con$line_segm[,1]),by = 2)) {
    lines(cur_con$line_segm[(kk:(kk+1)),])
  }
  title(paste("The map space with k =", cur_k*cur_k))
  return(list(cur_grid=cur_grid,cur_con=cur_con))
}

# init prototypes
# use mean and add noise along two princ components
# get prcomp
agk.init.prot <- function(dat,sigm) {
  pca <- prcomp(dat)
  pca <- pca$rotation
  rapc <- randn()*sigm*pca[,1]+randn()*sigm*pca[,2]
  for (ll in 1:(cur_k*cur_k-1)) {
    rapc <- rbind(rapc,randn()*sigm*pca[,1]+randn()*sigm*pca[,2])
  }
  prot     <- t(rapc)  + matrix(rep(as.numeric(colMeans(dat)),cur_k*cur_k),ncol = cur_k*cur_k)  
  prot <- t(prot)
  row.names(prot) <- NULL
  prot <- t(prot)
  return(prot)
}

## cross-validation ##
# get K-fold CV rule
# source: https://stackoverflow.com/questions/7402313/generate-sets-for-cross-validation-in-r
f_K_fold <- function(Nobs,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

# get K-fold CV rule, but without the negative correlation between training prediction and test truth
f_K_fold_decorr = function(Nobs,K=10) {
  # function that splits the data into two halfs (A,B)
  # then performs f_K_fold on both but with K/2
  # then takes training from A and matches randomly with test from B
  # and vice versa
  # this should reduce correlation between training pred and test truth
  
  # corr K
  if (K < 10) {
    stop('K too small!')
  } else {
    K_adj = round(K/2)
  }
  
  # create A and B Nobs
  Nobs_A = floor(Nobs/2)
  Nobs_B = ceil(Nobs/2)
  
  # create f_K_fold for both cases
  flds_A = f_K_fold(Nobs_A, K_adj)
  flds_B = f_K_fold(Nobs_B, K_adj)
  
  # correct the indices of B
  for (ii in 1:length(flds_B)) {
    flds_B[[ii]]$train = flds_B[[ii]]$train + Nobs_A
    flds_B[[ii]]$test = flds_B[[ii]]$test + Nobs_A
  }
  
  # match training from A to test from B and vice versa
  # A
  cur_order  = gtools::permute(1:K_adj)
  new_flds_A = list()
  for (ii in 1:length(flds_A)) {
    cur_tt = list()
    cur_tt$train = flds_A[[ii]]$train
    cur_tt$test  = flds_B[[cur_order[ii]]]$test
    new_flds_A[[ii]] = cur_tt
  }
  # B
  cur_order  = gtools::permute(1:K_adj)
  new_flds_B = list()
  for (ii in 1:length(flds_B)) {
    cur_tt = list()
    cur_tt$train = flds_B[[ii]]$train
    cur_tt$test  = flds_A[[cur_order[ii]]]$test
    new_flds_B[[ii]] = cur_tt
  }
  
  # make one flds obj
  flds = c(new_flds_A,new_flds_B)
}

f_K_fold <- function(Nobs,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

f_K_fold_b2g = function(cur_data,group,k) {
  # stratified k-fold for 2 groups
  # double k first to get desired k
  k2 = 2*k
  # make a nrow_variable
  cur_data$nrow = 1:length(cur_data[,1])
  # gets balanced k-fold for 2 groups
  cur_glevs   = levels(cur_data[[group]])
  dat_group_A = cur_data[cur_data[[group]] == cur_glevs[1],]
  dat_group_B = cur_data[cur_data[[group]] == cur_glevs[2],]
  
  # do it for each group separately
  # group A
  cur_n_A = length(dat_group_A[,1])
  cur_k_A = (k2/2)
  folds_A = f_K_fold(cur_n_A,cur_k_A)
  
  # group B
  cur_n_B = length(dat_group_B[,1])
  cur_k_B = (k2/2)
  folds_B = f_K_fold(cur_n_B,cur_k_B)
  
  # check if groups balanced
  if(cur_n_A != cur_n_B) {
    warning('Unbalanced group sizes; f_K_fold_b2g will not yield balanced sets on all folds.')
  }
  
  # get the original row.nums back
  # group A
  for(ii in 1:length(folds_A)) {
    folds_A[[ii]]$train = dat_group_A$nrow[folds_A[[ii]]$train]
    folds_A[[ii]]$test  = dat_group_A$nrow[folds_A[[ii]]$test]
  }
  # group B
  for(ii in 1:length(folds_B)) {
    folds_B[[ii]]$train = dat_group_B$nrow[folds_B[[ii]]$train]
    folds_B[[ii]]$test  = dat_group_B$nrow[folds_B[[ii]]$test]
  }
  
  # combine group A and B training and test
  fin_folds=list()
  for(ii in 1:length(folds_A)) {
    cur_folds = list()
    cur_folds$train = c(folds_A[[ii]]$train,folds_B[[ii]]$train)
    cur_folds$test  = c(folds_A[[ii]]$test,folds_B[[ii]]$test)
    fin_folds[[ii]] = cur_folds
  }
  
  # unit test: test if balanced group membership in training and test
  for (gg in 1:length(fin_folds)) {
    check_half = table(cur_data[[group]][fin_folds[[gg]]$train])
    if (check_half[1] != check_half[2]) {
      stop('Unbalanced training set created!')
    }
    check_half = table(cur_data[[group]][fin_folds[[gg]]$test])
    if (check_half[1] != check_half[2]) {
      stop('Unbalanced test set created!')
    }
  }
  
  # return
  return(fin_folds)
}

f_K_fold_b2g_dbl = function(cur_data,group,k,dat_match) {
  # stratified k-fold for 2 groups
  # will deal with duplicates (for strat. bootstr cur_data)
  # and restratifies training sets according to dat_match
  # tests sets are not necessarily balanced
  # training sets are
  
  # double k first to get desired k
  k2 = 2*k
  
  # kill duplicates
  cur_data_bcp = cur_data
  cur_data     = cur_data[duplicated(cur_data$subject) == FALSE,]
  
  # make a nrow_variable
  cur_data$nrow = 1:length(cur_data[,1])
  # gets balanced k-fold for 2 groups
  cur_glevs   = levels(cur_data[[group]])
  dat_group_A = cur_data[cur_data[[group]] == cur_glevs[1],]
  dat_group_B = cur_data[cur_data[[group]] == cur_glevs[2],]
  
  # do it for each group separately
  # group A
  cur_n_A = length(dat_group_A[,1])
  cur_k_A = (k2/2)
  folds_A = f_K_fold(cur_n_A,cur_k_A)
  
  # group B
  cur_n_B = length(dat_group_B[,1])
  cur_k_B = (k2/2)
  folds_B = f_K_fold(cur_n_B,cur_k_B)
  
  # check if groups balanced
  if(cur_n_A != cur_n_B) {
    warning('Unbalanced group sizes; f_K_fold_b2g will not yield balanced sets on all folds.')
  }
  
  # get the original row.nums back
  # group A
  for(ii in 1:length(folds_A)) {
    folds_A[[ii]]$train = dat_group_A$nrow[folds_A[[ii]]$train]
    folds_A[[ii]]$test  = dat_group_A$nrow[folds_A[[ii]]$test]
  }
  # group B
  for(ii in 1:length(folds_B)) {
    folds_B[[ii]]$train = dat_group_B$nrow[folds_B[[ii]]$train]
    folds_B[[ii]]$test  = dat_group_B$nrow[folds_B[[ii]]$test]
  }
  
  # combine group A and B training and test
  fin_folds=list()
  for(ii in 1:length(folds_A)) {
    cur_folds = list()
    cur_folds$train = c(folds_A[[ii]]$train,folds_B[[ii]]$train)
    cur_folds$test  = c(folds_A[[ii]]$test,folds_B[[ii]]$test)
    fin_folds[[ii]] = cur_folds
  }
  
  # fin_folds is now likely unbalanced and maybe unstrat
  # align dat match
  dat_match_aligned = dat_match[agk.in.dupl(dat_match$VPPG,cur_data$subject),]
  for (ff in 1:length(fin_folds)) {
    
    cur_train_data  = dat_match_aligned[fin_folds[[ff]]$train,]
    dms_final_ns_tr = agk.PDT.strat(cur_train_data)
    dms_final_ns_tr = dms_final_ns_tr$dms_final_ns
    
    # coding
    dat_match_aligned$id = 1:length(dat_match_aligned$VPPG)
    strat_tr_id          = as.numeric(agk.recode.c(dms_final_ns_tr$VPPG,dat_match_aligned$VPPG,dat_match_aligned$id))
    
    # strat this fold
    fin_folds[[ff]]$train = strat_tr_id
  }
  
  # unit test: test if balanced group membership in training and test
  for (gg in 1:length(fin_folds)) {
    check_half = table(cur_data[[group]][fin_folds[[gg]]$train])
    if (check_half[1] != check_half[2]) {
      stop('Unbalanced training set created!')
    }
    check_half = table(cur_data[[group]][fin_folds[[gg]]$test])
    if (check_half[1] != check_half[2]) {
      warning('Unbalanced test set created!')
    }
  }
  
  # return
  return(fin_folds)
}

agk.get.foldid = function(flds) {
  # get a foldid object
  # vector of values between 1 and nfold identifying what (test?) fold 
  # each observation is in
  # e.g. for glmnet
  
  # get the unique ids of train and test
  id_vec = c()
  for (ii in 1:length(flds)) {
    id_vec = unique(c(id_vec,unique(flds[[ii]]$train)))
    id_vec = unique(c(id_vec,unique(flds[[ii]]$test)))
  }
  
  # id vector ordered
  id_vec = id_vec[order(id_vec)]
  
  # foldid
  foldid = rep(NA,length(id_vec))
  for (ii in 1:length(flds)) {
    if (all(is.na(foldid[id_vec %in%  flds[[ii]]$test]))) {
      foldid[id_vec %in%  flds[[ii]]$test] = ii
    } else {
      stop('This/These subject(s) has already test fold assigned!')
    }
  }
  return(foldid)
}

agk.get.foldid.dupl = function(flds) {
  # get a foldid object
  # vector of values between 1 and nfold identifying what (test?) fold 
  # each observation is in
  # e.g. for glmnet
  
  # get the unique ids of train and test
  id_vec = c()
  for (ii in 1:length(flds)) {
    id_vec = unique(c(id_vec,unique(flds[[ii]]$train)))
    id_vec = unique(c(id_vec,unique(flds[[ii]]$test)))
  }
  
  # id vector ordered
  id_vec = id_vec[order(id_vec)]
  
  # foldid
  foldid = rep(NA,length(id_vec))
  for (ii in 1:length(flds)) {
    if (all(is.na(foldid[id_vec %in%  flds[[ii]]$test]))) {
      foldid[id_vec %in%  flds[[ii]]$test] = ii
    } else {
      stop('This/These subject(s) has already test fold assigned!')
    }
  }
  return(foldid)
}


# custom kernel for ksvm
# setClass(Class = 'rbfkernel_mahalanobis',contains = c('kernel'))
# 
# rbfdot_mahalanobis = function (sigma = 1, S) {
#   # custom kernel to be used with correlated
#   # features
#   # creates an object of class 'rbfkernel'
#   # works like rbfdot
#   # but with mahalonobis distance
#   # for this the covariance matrix (S) of data is needed
#   rval <- function(x, y = NULL, S) {
#     if (!is(x, "vector"))
#       stop("x must be a vector")
#     if (!is(y, "vector") && !is.null(y))
#       stop("y must a vector")
#     if (is(x, "vector") && is.null(y)) {
#       return(1)
#     }
#     if (is(x, "vector") && is(y, "vector")) {
#       if (!length(x) == length(y))
#         stop("number of dimension must be the same on both data points")
#       # get the mahalonobis distance
#       cur_m_distance = mahalanobis(x,y,S)
#       cur_rbd        = exp(-(cur_m_distance/2*sigma^2))
#       return(cur_rbd)
#     }
#   }
#   return(new('rbfkernel_mahalanobis', .Data = rval, kpar = list(sigma = sigma, S = S)))
# }

rbfdot_mahalanobis = function (x,y=NULL) {
  # custom kernel to be used with correlated
  # features
  # creates an object of class 'rbfkernel'
  # works like rbfdot
  # but with mahalonobis distance
  # for this the covariance matrix (S) of data is needed

  # checks
  if (!is(x, "vector"))
    stop("x must be a vector")
  if (!is(y, "vector") && !is.null(y))
    stop("y must a vector")
  if (is(x, "vector") && is.null(y)) {
    return(1)
  }

  # doing
  if (is(x, "vector") && is(y, "vector")) {
    if (!length(x) == length(y))
      stop("number of dimension must be the same on both data points")
    # get the mahalonobis distance
    cur_m_distance = mahalanobis(x,y,S)
    cur_rbd        = exp(-sigma*cur_m_distance)
    return(cur_rbd)
  }
}

class(rbfdot_mahalanobis) = "kernel"

# myrbfdot # TODO: does not give same results as rbfdot!
myrbfdot = function (x,y=NULL,sigma = 1) {
  # custom kernel to be used with correlated
  # features
  # creates an object of class 'rbfkernel'
  # works like rbfdot
  # but with mahalonobis distance
  # for this the covariance matrix (S) of data is needed

  # checks
  if (!is(x, "vector"))
    stop("x must be a vector")
  if (!is(y, "vector") && !is.null(y))
    stop("y must a vector")
  if (is(x, "vector") && is.null(y)) {
    return(1)
  }

  # doing
  if (is(x, "vector") && is(y, "vector")) {
    if (!length(x) == length(y))
      stop("number of dimension must be the same on both data points")
    # get the mahalonobis distance
    cur_m_distance = (norm(x-y,type = 'F'))^2
    cur_rbd        = exp(-sigma*cur_m_distance)
    return(cur_rbd)
  }
}
class(myrbfdot) = 'kernel'

## FUNCTIONS PDT ==============================================================
agk.get.lambda = function(df,str_gain,str_loss,str_splitter,with_limits = F) {
  # function to compute lambda
  # lambda = -loss/gain
  # uses str_gain, str_loss and str_splitter to find
  # gain, loss in all categories
  # supports only single str_splitters (i.e. simple interaction)
  # CAREFUL: function assumes that categories align in gain and loss
  # and of course both have to be modulated
  # with_limits: logical to say if you want to limit too big or small LAs 
  # NEEDS A DF as input 
  
  # get the df
  if (is.data.frame(df)) {
    df = df
  } else {
    df = coef(df)
  }
  if (!is.data.frame(df)) {
    stop('df is not a data frame nor can I extract coefficients data frame from df object.')
  }
  cur_names    = names(df)
  all_splts    = strsplit(cur_names,str_splitter)
  str_gainloss = c(str_gain,str_loss) 
  
  # getting all the string splitters
  str_spl_cl  = gsub('\\','',str_splitter,fixed=T)
  str_spl_cl  = gsub('(','',str_spl_cl,fixed=T)
  str_spl_cl  = gsub(')','',str_spl_cl,fixed=T)
  str_spl_cl  = strsplit(str_spl_cl,'|',fixed=T)
  str_spl_ind = c()
  for (cc in 1:length(str_spl_cl[[1]])) (
    if (!isempty(grep(str_spl_cl[[1]][cc],cur_names,fixed=T))) {
      str_spl_ind[cc] = T
    } else {
      str_spl_ind[cc] = F
    }
  )
  
  stopifnot(sum(str_spl_ind) == 1 | sum(str_spl_ind) == 0)
  cur_split_str = str_spl_cl[[1]][which(str_spl_ind)]
  
  # getting the gains and losses from the data frame
  cats_res_list = list()
  vars_res_list = list()
  for (ii in 1:length(str_gainloss)) {
    cur_str = str_gainloss[ii]
    vars_res = c()
    cats_res = c()
    for (jj in 1:length(all_splts)) {
      cur_splt = all_splts[[jj]]
      if (length(grep(cur_str,cur_splt)) == 0) {
        # intercepts
        next
      } else {
        if (length(cur_splt) == 1) {
          vars_res[jj] = cur_str[1]
          cats_res[jj] = ''
        } else {
          vars_res[jj] = paste0(cur_splt[1],cur_split_str,cur_splt[2])
          cats_res[jj] = cur_splt[2]
        }
      }
    }
    vars_res            = vars_res[!is.na(vars_res)]
    cats_res            = cats_res[!is.na(cats_res)]
    vars_res_list[[ii]] = vars_res
    cats_res_list[[ii]] = cats_res
  }
  
  # calc lambdas
  LA_names  = c()
  LA_values = list()
  for (ii in 1:length(vars_res_list[[1]])) {
    if(cats_res_list[[1]][ii] == '') {
      # case neutral
      LA_values[[ii]] = -df[vars_res_list[[2]][ii]]/df[vars_res_list[[1]][ii]]
    } else {
      # case categories
      neu_ind = which(cats_res_list[[1]] == '')
      LA_values[[ii]] = -(df[vars_res_list[[2]][ii]] + df[vars_res_list[[2]][neu_ind]])/(df[vars_res_list[[1]][ii]] + df[vars_res_list[[1]][neu_ind]])
    }
    
    if (with_limits) {
      # correcting immediately for NaN
      LA_values[[ii]][is.na(LA_values[[ii]])]         = 1
      LA_values[[ii]][LA_values[[ii]] == Inf]         = 10
      LA_values[[ii]][is.na(LA_values[[ii]] == -Inf)] = -10
    } else {
      if (any(LA_values[[ii]] == Inf) | any(LA_values[[ii]] == Inf)) {
        stop('Inf or -Inf LA_values[[ii]] == Inf when computing LA')
      }
    }


    # naming of LA variable
    if(cats_res_list[[1]][ii] == '') {
      LA_names[ii]    = paste0('LA',cats_res_list[[1]][ii])
    } else {
      LA_names[ii]    = paste0('LA:',cats_res_list[[1]][ii])
    }
  }
  
  if (with_limits) {
    # ensure borders of LA
    message('I am using a [-10 10] border rule for LA; arbitrary; better use ridge?')
    for (ii in 1:length(LA_values)) {
      cur_lav                = LA_values[[ii]]
      cur_lav[cur_lav < -10] = -10
      cur_lav[cur_lav > 10]  = 10
      LA_values[[ii]]        = cur_lav
    }
  }
  
  # baseline correct LA
  if (length(LA_values) > 1) {
    for (ii in 2:length(LA_values)) {
      LA_values[[ii]] = LA_values[[ii]] - LA_values[[1]]
    }
  }
  
  
  # ammend df
  for (ii in 1:length(LA_values)) {
    df[LA_names[ii]] = LA_values[[ii]]
  }
  
  # return
  return(df)
}

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

prob.log.reg <-function(int,b1,b2,x1,x2){
  return (1/(1+exp(-(int+x1*b1+x2*b2))))
}

agk.normalize <- function(x) {
  # normalize into 0-1 range
  return((x-min(x))/(max(x)-min(x)))
}

agk.corboxmetric <- function (cur_metr,cur_box) {
  # check if there is loc and noa in the PCs
  # function to check correlation of boxcar and found component
  # fits multinomial regression and returns predictive value
  cur_loc <- as.factor(cur_box)
  mod <- multinom(cur_loc ~ cur_metr)
  return(mean(as.numeric(predict(mod) == cur_loc)))
}

agk.polym.df = function(preds,degree) {
  # polym expansion of preds using polym
  # preds should be a data frame
  # default is orthogonal polynomials
  names_pred = names(preds)
  cmd        = paste(paste0('preds$',names_pred),collapse = ',')
  cmd        = paste0('polym(',cmd,',degree=',degree,')')
  res        = eval(parse(text=cmd))
  return(res)
}

agk.regress.out = function(x,preds,rob,inf_crit = 'AIC',ro_info) {
  # function to regress out of x what is useful using
  # pred as predictors up to polynomial expansion 
  # defined in des_poly
  # returns the residuals
  
  # pretest
  if (!is.numeric(x)) {
    res = x
    ro  = NA
  } else {
    preds        = agk.polym.df(preds,degree=2)
    preds        = as.data.frame(preds)
    names(preds) = paste0('v_',names(preds))
    names_pred   = names(preds)
    one_df       = data.frame(x,preds)
    
    # formula
    cur_form      = as.formula(paste('x~',paste(names_pred,collapse = '+')))
    
    # modeling
    if (rob == F) {
      cur_mod   = lm(cur_form,data = one_df)
      # model selection
      if (inf_crit == 'BIC') {
        cur_st = step(cur_mod,k=log(length(one_df[,1])),trace = 0)
      } else {
        cur_st = step(cur_mod,trace=0)
      }
    } else {
      cur_mod = lmRob(cur_form,data=one_df)
      # model selection
      cur_st  = step.lmRob(cur_mod) 
    }
    
    # was something regressed out?
    res = resid(cur_st)
    if (cor(res,x) == 1) {
      ro = 0
    } else {
      ro = 1
    }
  }
  
  # returning with or without ro info
  if (ro_info == 1) {
    return(list(res=res,ro=ro))
  } else {
    return(res)
  }
}

agk.regress.out.c = cmpfun(agk.regress.out)

agk.impute.mean <- function(x) {
  # function (for "apply") to impute vecmean for NA
  # give a vector
  cur_mean <- mean(x,na.rm = T)
  cur_fun <- function(y,cur_mean) {
    if (is.na(y)) {
      y<-cur_mean
    } 
    return(y)
  }
  return(sapply(x,cur_fun,cur_mean=cur_mean))
}

# ML with feature elimination
# lasso with many different start foldings
# consensus-based selection

agk.cut.for.high.entrp = function(x) {
  # decides whether entropy is higher
  # than expected entropy under
  # hypothesis of normal distribution
  
  # numeric?
  if(!is.numeric(x)) {
    return(0)
  }
  
  # density of x
  dens_x = density(x)
  
  # get the normal entropy
  make_n = randn(1,length(dens_x$x))*sd(x) + mean(x)
  n_dens = density(make_n)
  n_entr = agk.entr.st(n_dens$y)
  
  # get the observed entropy
  x_entr = agk.entr.st(dens_x$y)
  
  if(is.nan(x_entr))   {return(1)}
  if(x_entr >= n_entr) {return(1)} else {return(0)}
  
}

# all data: split in training and test
feat_sel = function(cur_dat,cur_lab,cnfolds) {
  
  # more preprocessing: killing correlated variables
  # throw out big correlation
  correlationMatrix = cor(cur_dat)
  highlyCorrelated  = findCorrelation(correlationMatrix, cutoff=0.90)
  if (!isempty(highlyCorrelated)) {
    cur_dat           = cur_dat[-highlyCorrelated]
  }
  
  # cut high entropy variables
  if (kill_high_entr) {
    highEntropy = which(as.logical(as.numeric(unlist(lapply(cur_dat,agk.cut.for.high.entrp)))))
    stopifnot(length(cur_dat) != length(highEntropy))
    if (!isempty(highEntropy)) {
      cur_dat           = cur_dat[-highEntropy]
    }
  }

  # repetition of CV for feature elimination
  reps = 30
  
  # make formula
  cur_form = as.formula(paste('HCPG ~ ',paste(names(cur_dat),collapse = ' + ')))
  
  # attach
  cur_data      = cur_dat
  cur_data$HCPG = cur_lab
  
  # surviviving predictors
  surv_list = list()
  
  # cvme vector
  cvmes = c()
  
  # ct
  ct = 0
  for (ii in 1:reps){
    
    disp(ii)
    
    # get new folds
    cur_flds    = f_K_fold_b2g(cur_data,'HCPG',cnfolds)
    cur_fold_id = agk.get.foldid(cur_flds)
    
    # the alphas
    #alphas = agk.norm.range(log(seq(0.001,1,length.out = length(cur_dat))),0,1)
    # this means I am prefering more complex models!
    alphas = seq(1,0,length.out = 4)
    
    # fit the model
    cur_cvmoda         = glmnetUtils::cva.glmnet(cur_form,family=cur_family,
                                                 data=cur_data,type.measure=type_measure,
                                                 alpha = alphas,foldid = cur_fold_id,use.model.frame=useModelFrame,
                                                 grouped=doGrouped,standardize = des_stand_glmnet)
    
    # best model
    winning_model = agk.glmnetUtils.cvalpha.getwinner(cur_cvmoda)
    
    # get the coef
    cur_surv = colnames(winning_model$coef)[abs(winning_model$coef)>0]
    if (!length(cur_surv) == 1) {
      # do not allow no-picks, because this will def mean a chance prediction
      ct = ct + 1
      surv_list[[ct]] = cur_surv
    }
    
    # get the cvme
    cvmes[ii] = winning_model$cvme
  }
  
  return(list(surv_list = surv_list, cvmes = cvmes))
}
feat_sel.c   = cmpfun(feat_sel)

# consensus rating
agk.consensus = function(surv_list, prop_crit,min_var) {
  # TODO: return threshold used
  # function to return which of the surviving features
  # appear in at least prop_crit times of the 
  # length of surv_list
  # will reduce prop_crit if not at least min_var feat surviving
  
  # check first if nested to deeply
  if (prop_crit < 0.2) {
    return(list(surv_feats=c(), prop_crit=prop_crit))
  }
  
  # core code
  all_feat = unlist(surv_list)
  all_feat = all_feat[-which(duplicated(all_feat))]
  
  cnt_fun = function(x,cur_feat) {return(any(cur_feat == x))}
  
  scores = c()
  for (ff in 1:length(all_feat)) {
    scores[ff] = sum(unlist(lapply(surv_list,cnt_fun,cur_feat = all_feat[ff])))
  }
  
  # survival?
  surv_feats = which(scores/length(surv_list) >= prop_crit)
  
  # selection of survival features
  surv_feats = all_feat[surv_feats]
  
  # kill intercept
  surv_feats = surv_feats[-grep('Intercept',surv_feats)]
  
  # check if enough vars survived
  if (length(surv_feats) >= min_var) {
    return(list(surv_feats=surv_feats, prop_crit=prop_crit))
  } else {
    agk.consensus(surv_list,prop_crit-0.01,min_var)
  }
  
}

agk.get.mean.model = function(all_winning_models) {
  the_coefs = list()
  for (ll in 1:length(all_winning_models)) {
    # get coefs
    cur_coefs = all_winning_models[[ll]]$coef
    # kick out intercept only models
    # TODO this maybe drop!
    #if (all(cur_coefs[2:length(cur_coefs)] == 0)) {next}
    # other than that scale the params to make comparable
    # TODO this maybe drop
    #the_coefs[[ll]] = agk.scale_range(cur_coefs,-1,1)
    the_coefs[[ll]] = cur_coefs
  }
  coef_names = colnames(the_coefs[[1]])
  the_coefs = lapply(the_coefs,as.numeric)
  the_coefs_ma = t(as.matrix(the_coefs[[1]]))
  for (ll in 2:length(the_coefs)) {
    the_coefs_ma = rbind(the_coefs_ma,the_coefs[[ll]])
  }
  colnames(the_coefs_ma) = coef_names
  cur_res = (colMeans(the_coefs_ma))
  return(cur_res)
}

## PDT: stratification ========================================================
# bootstrap a sample with stratification and with replacement

# functions
imp_miss = function(x) {
  if (!is.numeric(x)) {return(x)}
  cur_med = median(x,na.rm = T)
  x[is.na(x)] = cur_med
  return(x)
}

agk.mult.order = function(dms_d) {
  # function to order a data frame according to all
  # variables in the data frame starting at first
  # most left-hand
  all_str = c()
  for (ss in 1:length(dms_d)) {
    all_str[ss] = paste0(deparse(substitute(dms_d)),'[[',ss,']]')
  }
  cmd = paste('order(',paste(all_str,collapse = ','),')')
  cmd = paste0(deparse(substitute(dms_d)),'[',cmd,',]')
  dms_d_o = eval(parse(text=cmd))
  return(dms_d_o)
}

agk.droplevels = function(x) {
  # function that drops levels
  # if no levels then returns x
  if (is.null(levels(x))) {return(x)}
  return(droplevels(x))
}

agk.strata.n = function(x) {
  # function that takes x and always expects
  # first is one group, second is other and
  # repeat, takes them by pair and then
  # takes the higher one and replaces lower
  
  # prepare for push to beyond median (HC more PG and vice versa)
  # inds = seq(1,length(x),by=2)
  # dire = c()
  # for (ii in inds) {
  #   if (x[c(ii)] > x[c(ii+1)]) {
  #     dire[ii] = 1
  #   } else {
  #     dire[ii] = -1
  #   }
  # }
  
  # go to median
  inds = seq(1,length(x),by=2)
  for (ii in inds) {
    if(any(x[c(ii,ii+1)] == 0)) {x[c(ii,ii+1)] = 1}
    x[c(ii,ii+1)] = rep(ceil(median(x[c(ii,ii+1)])),2)
  }
  
  # push to beyond median (HC more PG and vice versa)
  # for (ii in inds) {
  #   if(any(x[c(ii,ii+1)] <= 3)) {next}
  #   if (dire[ii] == 1) {
  #     x[ii]   = x[ii]   - 2
  #     x[ii+1] = x[ii+1] + 2 
  #   } else {
  #     x[ii]   = x[ii]   + 2 
  #     x[ii+1] = x[ii+1] - 2 
  #   }
  # 
  # }
  
  # check that it is an even sum
  #if(mod(sum(x),2)) {
  #  x[1] = x[1]-1
  #}
  return(x)
}


agk.PDT.strat.prep = function(dat_match, km_res = NULL) {
  # core function of PDT stratification
  
  # dat_match edu
  dat_match$edu = dat_match$edu_years + dat_match$edu_years_voca
  
  # vars which need to be strat
  vars_strat      = c('HCPG',cur_names_dom)
  vars_strat_cmpl = c('edu','smoking_ftdt','Age','income_personal','audit','dem_gender','handedness','unemployed')
  vars_strat_k    = c('smoking_ftdt','edu','Age','audit','income_personal')
  vars_strat      = c('HCPG',vars_strat_k)
  
  # strat it
  dms            = dat_match
  row.names(dms) = dms$VPPG
  dms_cmpl       = dms[c(vars_strat_cmpl,'VPPG','HCPG')]
  dms            = dms[c(vars_strat,'VPPG')]
  
  # kmeans
  dms_k          = dms[c(vars_strat_k)]
  cur_vp         = row.names(dms_k)
  dms_k          = as.data.frame(lapply(dms_k,as.numeric))
  dms_k          = as.data.frame(lapply(dms_k,imp_miss))
  dms_k          = as.data.frame(lapply(dms_k,agk.scale_range,C=0,D=1))
  dms_k$edu      = dms_k$edu
  #dms_k$HCPG     = dms$HCPG
  #dms_k_PG       = dms_k[dms_k$HCPG == 'PG',]
  #dms_k_PG$HCPG  = NULL
  #dms_k$HCPG     = NULL
  #km_res_PG      = kmeans(dms_k_PG,3) 
  #km_res         = kmeans(dms_k,km_res_PG$centers)
  if (!is.null(km_res)) {
    km_res = kmeans(dms_k,km_res$centers)
  } else {
    km_res         = kmeans(dms_k,3)
  }
  km             = as.data.frame(km_res$cluster)
  row.names(km)  = cur_vp
  names(km)      = 'km_cluster'
  km$VPPG        = row.names(km)
  dms            = merge(dms,km,by='VPPG')
  vars_strat     = c('HCPG','km_cluster')
  
  # using strata
  dms_d              = dms
  
  # replace missings
  dms_d          = as.data.frame(lapply(dms_d,imp_miss))
  
  # with kmeans only few variables
  dms_d = dms_d[c('HCPG','km_cluster','VPPG')]
  
  # order data
  dms_d              = agk.mult.order(dms_d)
  
  # droplevels
  dms_d = as.data.frame(lapply(dms_d,agk.droplevels))
  
  # check table
  current_table = paste('xtabs(~',paste(vars_strat,collapse='+'),',data=dms_d)')
  current_table = eval(parse(text=current_table))
  
  # get useful strata numbers
  strata_n = as.numeric(current_table)
  strata_n = agk.strata.n(strata_n)
  strata_n = as.numeric(t(matrix(strata_n,nrow = 2)))
  
  resf                 = list()
  resf$dms_d           = dms_d 
  resf$dms             = dms 
  resf$dms_cmpl        = dms_cmpl 
  resf$strata_n        = strata_n
  resf$vars_strat      = vars_strat 
  resf$vars_strat_cmpl = vars_strat_cmpl
  resf$dat_match       = dat_match
  resf$km_res          = km_res
  
  return(resf)
}

agk.PDT.strat.core = function(resf) {
  # core function of PDT stratification
  
  disp('unpacking')
  dms_d           = resf$dms_d 
  strata_n        = resf$strata_n
  vars_strat      = resf$vars_strat 
  vars_strat_cmpl = resf$vars_strat_cmpl
  dms             = resf$dms 
  dms_cmpl        = resf$dms_cmpl
  dat_match       = resf$dat_match
  km_res          = resf$km_res
  
  # get the stratum order
  # TODO: find order of strata
  
  # searching a stratification
  searching_strat = 1
  ct              = 0 
  
  while(searching_strat) {
    ct = ct + 1
    # not too often
    if (ct > 20) {
      return(NULL)
    }
    # get a stratification
    disp('getting strat')
    res = strata(dms_d,size=strata_n,stratanames = vars_strat,method='srswr')
    
    # finally sample
    disp('sampling new data')
    new_dat         = getdata(dms_d,res)
    dms_cmpl$id     = 1:length(dms_cmpl$VPPG)
    new_dat$id      = agk.recode.c(new_dat$VPPG,dms_cmpl$VPPG,dms_cmpl$id)
    dms_final_ns    = dms_cmpl[as.numeric(new_dat$id),]
    dms_final_ns$id = NULL
    disp('done sampling new data')
    
    # test all the vars of int (if stratified)
    disp('starting matching tests')
    disp(paste('we have this many tests:', length(vars_strat_cmpl)))
    for (vv in 1:length(vars_strat_cmpl)) {
      disp(paste('I am at test:', vv))
      
      # check if data frame is not too long
      #if(length(dms_final_ns[,1]) > 70) {break}
      
      # check if HCPG is strat
      HCPG_tab = table(dms_final_ns$HCPG)
      if(HCPG_tab[1] != HCPG_tab[2]) {
        disp('unequal number HC and PG')
        return(NULL)
        }
      if (vars_strat_cmpl[vv] == 'HCPG') {next}
      
      disp(vars_strat_cmpl[vv])
      
      if (is.numeric(dms_final_ns[[vars_strat_cmpl[vv]]])) {
        cur_sum = summary(lm(dms_final_ns[[vars_strat_cmpl[vv]]]~dms_final_ns$HCPG))
        if (cur_sum$coefficients[2,4] < 0.1) {
          disp('a no-match happened')
          break
        }
      } else {
        cur_x_tab = xtabs(~dms_final_ns$HCPG+dms_final_ns[[vars_strat_cmpl[vv]]])
        cur_sum = fisher.test(cur_x_tab)
        if (cur_sum$p.value < 0.1) {
          disp('a no-match happened')
          break
        }
      }
    }
    
    # found a strat
    if (vv == length(vars_strat_cmpl)) {
      disp('Found a stratification!')
      return(list(dms_final_ns = dms_final_ns,km_res = km_res))
    }
  }
}

agk.PDT.strat = function(dat_match,km_res = NULL) {
  # PDT stratification
  try_strat = 1
  while(try_strat) {
    resf         = agk.PDT.strat.prep(dat_match,km_res = km_res)
    dms_final_ns = agk.PDT.strat.core(resf)
    if (!is.null(dms_final_ns)) {break}
  }
  return(dms_final_ns)
}

## ADDITIONAL FUNCTIONS =======================================================

moving.average <- function(x,n=5) {
  # moving average function
  filter(x,rep(1/n,n), sides=2,circular = T)
}

mylogit_predict <- function (cur.data,cur.coeff) {
  # predict function
  cur.pred <- c()
  for (ii in 1:length(cur.data[,1])) {
    tmp <- prob.log.reg(cur.coeff$Gewinn, cur.coeff$Verlust, cur.coeff$intercept, cur.data$Gewinn[ii], cur.data$Verlust[ii])
    cur.pred <- rbind(cur.pred, as.numeric(tmp))    
  }
  return (cur.pred)
}

mean.rmna <- function(x) {
  # mean function that per default removes NA
  mean(x,na.rm=T)
}

sum.rmna <- function(x) {
  # sum function that per default removes NA
  sum(x,na.rm=T)
}


mod<-function(x,m){
  # modulus
  t1<-floor(x/m)
  return(x-t1*m)
}

get.log <- function(v1){
  # this get_log function deals with negative values in a 2step approach:
  # the vector will firstly be transformed such that
  # it will be will always have 1 as its smallest value
  # then the log function is applied
  # i.e. the constant to be added is custom-made to the vector;
  # source: http://blogs.sas.com/content/iml/2011/04/27/log-transformations-how-to-handle-negative-data-values.html
  a = 1-min(v1)
  return (log(v1 + 1 - min(v1)))
}

get.log.base <- function(v1, base){
  # this get_log function deals with negative values in a 2step approach:
  # the vector will firstly be transformed such that
  # it will be will always have 1 as its smallest value
  # then the log function is applied
  # i.e. the constant to be added is custom-made to the vector;
  # source: http://blogs.sas.com/content/iml/2011/04/27/log-transformations-how-to-handle-negative-data-values.html
  # base can be other than e
  a = 1-min(v1)
  return (log((v1 + 1 - min(v1)),base=base))
}

agk.normality.check = function(mod) {
  # normality check of the resid of a model
  # shapiro only if n is appropriate
  cr  = resid(mod)
  crt = rnorm(length(cr), mean = mean(cr), sd=sd(cr))
  print(ks.test(cr,crt))
  hist(cr)
  qqnorm(cr)
  if ((length(cr) >= 3) & (length(cr) <= 5000)) {
    print(shapiro.test(resid(mod)))
  } else {
    disp(paste0("N of observations (",length(cr),
               ") is inappropriate for shapiro test"))
  }
}

agk.check.mvnorm <- function(v1,v2) {
  # calc the mv normality
  require(mvnormtest)
  v1 <- v1[!is.na(v1)]
  v1 <- v1[!is.na(v2)]
  v2 <- v2[!is.na(v1)]
  v2 <- v2[!is.na(v2)]
  return (mshapiro.test(t(as.matrix(data.frame(v1,v2)))))
}

spearman_brown <- function(r,f) {
  # Spearman-Brown-Formula
  # r is current correlation
  # f is factor by which test gets prolonged (e.g. 2 in split-half case)
  rel = f*r/(1+(f-1)*r)
  return(rel)
}

v.ran <- function(v1,perc) {
  # a function that asks for a 0,1 vector and will randomly pick perc% elements
  # and replace them with random 0 or 1
  how_many <- round((perc/100)*length(v1)) # how many elements should be replaced
  where    <- sample(1:length(v1), how_many, replace=F) # get indices of where to replace
  v1[where]<- sample(0:1,how_many, replace=T) # fill in random 0' and 1's in v1
  return (v1)
}

calc.lambda<-function(bg,bl){
  # loss aversion: calculate lambda
  return(bl*(-1)/bg)
}

## gets a subsample of trials (rows) split by "variable" (e.g. "subject")

make_sample_of_trials <- function(variable, data.la, n,check,same,verbosity) {
  got_rand <- 0
  subs <- levels(eval(parse(text=paste("data.la$",variable,sep = ""))))
  new_data <- data.frame()
  for (ii in 1:length(subs)) {
    if (verbosity == 1) {
      print(paste("now running...", subs[ii]))
    }
    cur.data <- data.la[eval(parse(text=paste("data.la$", variable,sep = ""))) == subs[ii],]
    # save this data frame as baseline for later represent. check
    cur.d.sv <- cur.data
    cur.ind  <- seq(1,length(cur.data[,1]))
    repr_result <- 0
    # make a new sample as long as not representative
    count <- 0
    while (repr_result == 0 & count < 3000){
      if (verbosity == 1) {
        print(count)
      }
      count <- count + 1
      if (got_rand == 0) {
        cc.ind   <- sample(cur.ind, n,replace = F)
        # same random sample for all subjects
        if (same == 1) {got_rand <- 1}
      }
      cur.data <- cur.d.sv[cc.ind,]
      if (check == 0) {
        repr_result <- 1
      } else {repr_result <- check_representativeness(cur.data,cur.d.sv,0.05)}
    }
    if (count >= 3000) {print(paste("attention, did not find adequate sample in...", subs[ii]))}
    new_data <- rbind(new_data,cur.data)
  }
  return(new_data)
}

## make sample of trials but based on gambles; sample of certain gambles is picked and then this is been looked for in the data and taken

make.gamble.sample <- function(variable, data.la, n,verbosity) {
  subs <- levels(eval(parse(text=paste("data.la$",variable,sep = ""))))
  new_data <- data.frame()
  count <- 0
  for (ii in 1:length(subs)) {
    count <- count +1
    if (verbosity == 1) {
      print(paste("now running...", subs[ii]))
    }
    cur.data <- data.la[eval(parse(text=paste("data.la$", variable,sep = ""))) == subs[ii],]
    # save this data frame as baseline for later represent. check
    cur.d.sv <- cur.data
    # here now get the random gamble sample
    if (count == 1) {
      gamble_sample <- sample.gambles(n,data.la)
    }
    # here we pull the desired gamble sample from the current subject's data
    cur.data <- get.sample.trials(cur.d.sv,gamble_sample)
    new_data <- rbind(new_data,cur.data)
  }
  return(list(new_data,gamble_sample))
}

## function to decide for sample of gambles

sample.gambles <- function(n_gam,data.la) {
  pos_gain <- as.numeric(levels(as.factor(as.character(data.la$Gewinn))))
  pos_loss <- as.numeric(levels(as.factor(as.character(data.la$Verlust))))
  comvec   <- agk.combvec(pos_gain,pos_loss)
  which_tr <- comvec[sample(length(1:length(comvec[,1])),n_gam),]
  return(which_tr)
}

## function that takes the gamble_sample and takes the certain trials out of the data.la
get.sample.trials <- function(data.la,gamble_sample) {
  new_data <- data.frame()
  for(ii in 1:length(gamble_sample[,1])) {
    cur_new <- data.la[data.la$Gewinn == gamble_sample[ii,1] & data.la$Verlust == gamble_sample[ii,2],]
    new_data <- rbind(new_data,cur_new)
  }
  return(new_data)
}

## comvec function
agk.combvec <- function(x,y){
  out_m <- matrix(0,nrow = length(x)*length(y),ncol=2)
  count <- 0
  for (ii in 1:length(x)) {
    for (jj in 1:length(y)) {
      count <- count+1
      out_m[count,1] <- x[ii]
      out_m[count,2] <- y[jj]
    }
  }
  return(out_m)
}

## function to check whether the drawn sample is representative

check_representativeness <- function(cur.data.la, data.la,alpha) {
  check_result <- 1
  
  desired_variables <- list("ev","ed", "ed.abs","ratio","diff","RiskMar")
  
  for (ii in 1:length(desired_variables)) {
    # testing a against b baseline
    cur.a <- eval(parse(text=paste("cur.data.la$",desired_variables[[ii]],sep = "")))
    cur.b <- eval(parse(text=paste("data.la$",desired_variables[[ii]],sep = "")))
    
    # test
    cur.t <-t.test(cur.data.la$ev,data.la$ev)
    cur.t <- cur.t$p.value
    if (cur.t < alpha) {
      check_result <- 0
      return(check_result)}
    
    cur.t <- cur.t <-bartlett.test(list(cur.data.la$ev,data.la$ev))
    cur.t <- cur.t$p.value
    if (cur.t < alpha) {
      check_result <- 0
      return(check_result)}
  }
  # all tests passed
  return(check_result)
}

## for beginning plotting
prep_ggplot_la <- function() {
  
  library(ggplot2)
  library(Hmisc)
  
  # set the theme
  theme_update(panel.background=element_rect(fill="white", colour="black", size=2.0),
               panel.grid = element_blank(),
               plot.background = element_rect(fill="white", color = "white"),
               axis.text = element_text(colour = "black",size=23,face="bold"),
               axis.ticks = element_line(colour ="black", size=1),
               axis.title = element_text(colour = "black",size=23,face="bold"),
               legend.position = c(0.7, 0.85),
               legend.background = element_rect(fill="white", color = "white"),
               legend.title = element_blank(),
               legend.text = element_text(colour = "black",size=18,face="bold"))
  print("Feedback: updated theme; if you want to change, change the function in 'La_functions.R' and run again.")
}

## function to paste a short variable to a long data.frame given an extending variable
# e.g. age is a variable per subject; does not vary per trials, but in long format i need to have it per trial
# in this version long and short var must be ordered by extend var already!

paste_short_on_long <- function (long, short_var, extend_var, new_name) {
  sub_var_expr <- parse(text=paste("data.la$",extend_var,sep = ""))
  subs <- levels(eval(sub_var_expr))
  new.var <- c()
  for (ii in 1:length(subs)) {
    cur.trials <- long[eval(sub_var_expr) == subs[ii],]
    cur.newvar <- rep(short_var[ii],length(cur.trials[,1]))
    new.var <- c(new.var,as.numeric(cur.newvar))
  }
  new.var <- as.data.frame(new.var)
  names(new.var) <- new_name
  new_long <- data.frame(long,new.var)
}

## payout (BGG style)

agk.calc.payout <- function(cur.data, extend_var, base_value,n) {
  sub_var_expr <- parse(text=paste("cur.data$",extend_var,sep = ""))
  eval(parse(text=paste("cur.data$",extend_var,"<-as.factor(as.character(cur.data$",extend_var,"))",sep = "")))
  subs <- levels(eval(sub_var_expr))
  sub.payout <- c()
  for (ii in 1:length(subs)) {
    cur.trials <- cur.data[eval(sub_var_expr) == subs[ii],]
    cur.acc    <- subset(cur.trials,Button==1)
    cur.ind    <- seq(1,length(cur.acc[,1]))
    cur.ind    <- sample(cur.ind,n)
    cur.gam    <- cur.acc[cur.ind,]
    cur.res    <- c()
    for (kk in 1:length(cur.gam[,1])) {
      cur.res[kk] <- agk.local.payout(cur.gam$Gewinn[kk], cur.gam$Verlust[kk]) 
    }
    sub.payout[ii] <- sum(as.numeric(cur.res)) + base_value
  }
  return(sub.payout)
}

agk.local.payout <- function(gain,loss) {
  if (runif(1)>0.5) {
    return(gain)} else {return(loss)}
}

## function to plot the stability of lambda per trial
agk.lambda.stab <- function(cur.data, extend_var,start_point,end_point) {
  lambda_evolution <- c()
  for (ii in start_point:end_point) {
    cur.data.loop <- agk.sel.trials(cur.data,extend_var,ii)
    cur.lambdas   <- agk.lmlist.data.la(cur.data.loop)
    cur.lambdas   <- cur.lambdas$Lambda_ed
    lambda_evolution <- cbind(lambda_evolution,cur.lambdas)
  }
  for (kk in 1: length(tmp[,1])) {
    plot(as.numeric(tmp[kk,]),type="line") 
  }
  return(lambda_evolution)
}

## function to assess stability of model by prediciting next vali_trials on every step
agk.lambda.stab.pred <- function(cur.data, extend_var,start_point,end_point,vali_trials) {
  pred_evolution   <- c()
  collect_vali     <- list()
  collect_train    <- list()
  collect_lmlist   <- list()
  count            <- 0
  for (ii in start_point:end_point) {
    count <- count +1
    cur.data.loop <- agk.sel.trials(cur.data,extend_var,1,ii)
    collect_train[[count]] <- cur.data.loop
    cur.data.loop <- subset(cur.data.loop,!is.na(accept.reject))
    collect_vali[[count]]  <- agk.sel.trials(cur.data,extend_var,ii+1,ii+vali_trials)
    collect_lmlist[[count]] <- agk.lmlist.data.la.list(cur.data.loop)
  }
  # calc and plot pred evolution
  # go in every trial step
  for (kk in 1: length(collect_lmlist)) {
    # at trial step, go in every subject
    cur.pred <- c()
    for(ll in 1:length(collect_lmlist[[kk]])) {
      # make subject variable generic
      collect_vali[[kk]]$subject <- as.factor(as.numeric(collect_vali[[kk]]$subject))
      tmp.newdata <- collect_vali[[kk]][as.numeric(collect_vali[[kk]]$subject) ==ll,]
      tmp <- table((round(predict(collect_lmlist[[kk]][[ll]],
                                  newdata=tmp.newdata,"response"))-as.numeric(as.character(tmp.newdata$accept.reject)))==0)
      cur.pred[ll] <- as.numeric(tmp["TRUE"]/(tmp["TRUE"]+tmp["FALSE"]))
    }
    pred_evolution <- cbind(pred_evolution,cur.pred)
  }
  # return the calculations
  return(list(pred_evolution))
}

## function to get a data.la and get via lm.list the lambda vector out of it
## watch out; in "cur.data" the extend_variable must be called "subject"
agk.lmlist.data.la <- function(cur.data) {
  lmlist_ed <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = cur.data, na.action = na.omit, family = "binomial")
  crm_lmlist <- c()
  for (ii in 1:length(lmlist_ed)) {
    crm_lmlist <- rbind(crm_lmlist,t(as.matrix(as.numeric(lmlist_ed[[ii]][[1]]))))
  }
  crm_lmlist <- as.data.frame(crm_lmlist)
  names(crm_lmlist) <- c("intercept", "ed.abs", "Gewinn", "Verlust")
  crm_lmlist$id       <- as.factor(as.numeric(as.matrix(names(lmlist_ed))))
  crm <- crm_lmlist
  crm$Lambda_ed <- calc_lambda(crm$Gewinn,crm$Verlust)
  return(crm$Lambda_ed)
}

## function to get a data.la and get via lm.list the whole lmlist
agk.lmlist.data.la.list <- function(cur.data) {
  lmlist_ed <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = cur.data, na.action = na.omit, family = "binomial")
  return(lmlist_ed)
}


## function to just select x trials from data.la (every row is a trial)
agk.sel.trials <- function(cur.data, extend_var,start_n,end_n) {
  sub_var_expr <- parse(text=paste("cur.data$",extend_var,sep = ""))
  eval(parse(text=paste("cur.data$",extend_var,"<-as.factor(as.character(cur.data$",extend_var,"))",sep = "")))
  subs <- levels(eval(sub_var_expr))
  data.out <- data.frame()
  for (ii in 1:length(subs)) {
    cur.trials <- cur.data[eval(sub_var_expr) == subs[ii],]
    cur.trials <- cur.trials[start_n:end_n,]
    data.out   <- rbind(data.out,cur.trials)
  }
  return(data.out)
}

# function to measure the amount of match
agk.perc.match <- function(x,y) {
  if (length(x) != length(y)) {
    print("x and y must have same length!")
    return} else {
      tmp <- sum(as.numeric(x==y))/length(x)
      return(tmp)
    }
}

# function to measure the amount of match
# allow to reset one level in the second argument
# e.g. when there are two neutral categories which are
# just neutral in the first argument
agk.perc.match.reval <- function(x,y,reval_code) {
  y <- agk.recode(y,reval_code[1],reval_code[2])
  y <- as.factor(as.numeric(y))
  if (length(x) != length(y)) {
    print("x and y must have same length!")
    return} else {
      tmp <- sum(as.numeric(x==y))/length(x)
      return(tmp)
    }
}

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

# a recode function
# function to recode x, given a source vector y and a translated vector z
# source vector is giving the centers on a scale; minimal distance will lead to its translation value
# for recoding
# e.g. given a vector with values ranging from 1 to 100; source vector gives centers 10, 40, 60, 80
# then the value 9 will be put into category "1" with center 10;
# distance is euclidean distance
agk.recode.dist <- function(x,y,z) {
  for (ii in 1:length(x)) {
    cur_dist <- c()
    for(jj in 1:length(y)) {
      cur_dist[jj] <- dist(c(x[ii],y[jj]))
    }
    x[ii] <- which(cur_dist == min(cur_dist))
  }
  return(x)
}

# a recode function
# same as agk.recode.dist but for data.frame; does it for all the variables in data frame
agk.recode.dist.df <- function(df,y,z) {
  sv_names <- colnames(df)
  for (kk in 1:length(df[1,])) {
    x <- df[,kk]
    x <- agk.recode.dist(x,y,z)
    df[kk] <- x
  }
  colnames(df) <- sv_names
  return(df)
}

# scale if exists
# ifnot: does not do anything
agk.scale.ifexists <-function(df,var) {
  tmp_var <- eval(parse(text=paste(df,"$",var,sep="")))
  # does the variable exist?
  if(length(eval(parse(text=paste(df,"$",var,sep=""))))!=0) {
    tmp_var <- scale(tmp_var)
    return (tmp_var)
  }
  # if does not exist then return NA-variable
  tmp <- as.matrix(rep(NA,length(eval(parse(text=paste(df,"[,1]",sep=""))))))
  colnames(tmp) <- var
  return(tmp)
}

# varselect from a data frame
# takes variables of a df and makes a new df of this subselection
# if variable does not exist then NAs will be produced
agk.varselect <- function(cur_df, vars) {
  # check if exists
  sel_df <- data.frame(rep(NA,length(cur_df[,1])))
  for (ii in 1:length(vars)) {
    cur_cmd <- paste("length(","cur_df","$",vars[ii],")==0",sep="")
    if (eval(parse(text=cur_cmd))) {
      sel_df <- cbind(sel_df,rep(NA,length(cur_df[,1])))
    } else {
      cur_var <- eval(parse(text=paste("cur_df","$",vars[ii],sep="")))
      sel_df <- cbind(sel_df,cur_var)
    }
  }
  sel_df <- sel_df[,-1]
  names(sel_df) <- vars
  return(sel_df)
}

# batch scaling; puts together all variables as if one; gets var and sd and then scales every cell with that
agk.batch.scale <- function(int_df,fixed_mean) {
  all_nums <- na.omit(as.numeric(as.matrix(int_df)))
  cur_sd <- sd(all_nums)
  if (!is.null(fixed_mean)) {
    cur_mn = fixed_mean
  } else {
    cur_mn <- mean(all_nums)
  }
  int_df <- (int_df-cur_mn)/cur_sd
  names(int_df) <- paste(names(int_df),"s",sep="")
  return(int_df)
}

# computes the glmer model and bootstraps confintervals of desired parameters 
agk.la.by.defcat <- function(cur_cat,cur_fun,ncpus,nsim,run_boot,study_name,ranef_cat) {
  
  ## make the name
  cur_date  <- date()
  cur_date  <- strsplit(cur_date," ")
  cur_date  <- paste(cur_date[[1]][5],cur_date[[1]][2],cur_date[[1]][3],sep="_")
  # get the n per cat
  tmp       <- subset(data_pdt, !is.na(accept_reject))
  eval(parse(text=paste("tmp <- subset(tmp,!is.na(",cur_cat,"))",sep="")))
  n_levels  <- eval(parse(text=paste("length(levels(data_pdt$",cur_cat,"))",sep="")))
  cur_n     <- round(length(tmp[,1])/n_levels)
  
  cur_name  <- paste("boot_pdt_", cur_cat, "_",cur_n,"_",cur_date,".RData",sep="")
  
  ## make the model command and run
  if (ranef_cat == 0) {
    first         <- "pdt_bobyqua_cat <- glmer(accept_reject ~ (gain + loss)*"
    third         <- '+(gain + loss|subject), data=data_pdt, family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))'
    cur_model     <- paste(first,cur_cat,third,sep="")
  } else if (ranef_cat == 1) {
    first         <- "pdt_bobyqua_cat <- glmer(accept_reject ~ (gain + loss)*"
    third         <- '+(((gain + loss)*'
    fifth         <- ')|subject), data=data_pdt, family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10000)))'
    cur_model     <- paste(first,cur_cat,third,sep="")  
  }
  
  eval(parse(text=cur_model))
  
  ## BIC whether category as additional fixed effect is useful
  tmp     <- anova(pdt_bobyqua_noc,pdt_bobyqua_cat)
  cur_bic <- tmp$BIC[2]
  
  ## bootstrap  
  if (run_boot == 1) {
    # bootsrapping CIs
    cur_boot <- bootMer(pdt_bobyqua_cat,
                        FUN = cur_fun,
                        nsim = nsim, seed = NULL, use.u = FALSE,
                        type = "parametric",
                        verbose = TRUE, .progress = "txt",
                        parallel = "snow", ncpus = ncpus)
  }
  
  # save the bootstrap and the model
  setwd(path_res)
  if (run_boot == 1) {
    save(list=c("cur_boot","pdt_bobyqua_cat"), file=cur_name) } else {
      save(list=c("pdt_bobyqua_cat"),file=cur_name)
    }
  
  # plot the boot object
  if (run_boot == 1) {
    # set results path
    setwd(path_res)
    cur_plot <- agk.plot.boot(cur_name = cur_name,cur_n = cur_n, cur_cat = cur_cat,nsim = nsim,study_name)
    return(list(cur_bic=cur_bic,cur_plot=cur_plot))
  } else {
    return(list(cur_bic=cur_bic,model=pdt_bobyqua_cat))
  }
}

# function to plot a boot object and save the plot
agk.plot.boot <- function(cur_name, cur_n, cur_cat, nsim,study_name) {
  
  #   # check if cur_boot is already there and delete if yes
  #   if (exists("cur_boot",envir = .GlobalEnv)) {rm(list = c("cur_boot"),envir = .GlobalEnv)}
  
  ne <- new.env()
  load(file=cur_name, env=ne)
  cur_boot <- eval(parse(text=ls(env=ne,pattern = "*boot*")),envir = ne)
  
  # how many failed and nsim
  nfail <- attr(cur_boot,"bootFail")
  nsim  <- as.character(nsim)
  
  # plot bootstrapped distribution of fixed effects
  boot_pdt_df <- cur_boot
  boot_pdt_df <- data.frame(boot_pdt_df$t)
  boot_pdt_df <- melt(boot_pdt_df)
  boot_pdt_df$variable <- as.factor(boot_pdt_df$variable)
  
  # define whiskers
  f <- function(x) {
    r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  # do it
  p1 <- ggplot(boot_pdt_df, aes(variable, value)) + stat_summary(fun.data = f, geom="boxplot")
  p1 <- p1 + geom_abline(intercept=0,slope=0)
  p1 <- p1 + theme(axis.text.x=element_text(angle=-45, size =13, color = "black"))
  tmp_title <- paste("bootstrapped distribution of",
                     "fixed effects\n",
                     "study name: ", study_name,"\n",
                     "n per cat level =",cur_n,"\n",
                     "n of bootstrap =", nsim, "\n whiskers repr. 95% intervals\n",
                     "n of failed sim =", nfail,"\n",
                     "BY:", cur_cat)
  p1 <- p1 + ggtitle(tmp_title)
  p1
  return(p1)
}

## not my function
# for exporting plots
ExportPlot <- function(gplot, filename, width=2, height=1.5) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
  postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
  print(gplot)
  dev.off()
  png(file = paste(filename, '_.png', sep=""), width = width * 100, height = height * 100)
  print(gplot)
  dev.off()
}

## Kernel PCA functions
# function to center a kernel matrix
agk.center.kmatrix <- function(K) {
  Kc <- matrix(NA,nrow=length(K[,1]),ncol=length(K[1,]))
  p  <- length(K[,1])
  last_term <- (1/(p^2))*sum(K[,])
  for (ii in 1:length(K[,1])) {
    for (jj in 1:length(K[1,])) {
      Kc[ii,jj] <- K[ii,jj] - (1/p)*sum(K[,jj]) - (1/p)*sum(K[ii,]) + last_term 
    }
  }
  return(Kc)
}

# function to project new data point
# onto ev (i.e. pc) from kpca
# returns a scalar
agk.proj.on.kpca.ev <- function(ev,dat,new_pt,sigm) {
  rbf <- rbfdot(sigma = sigm)
  # create kernel matrix w.r.t. new data point
  kern_mat_new <- kernelMatrix(rbf,dat,new_pt)
  # center the matrix
  kern_mat_new <- agk.center.kmatrix(kern_mat_new)
  u <- c()
  for (ii in 1:length(kern_mat_new[,1])) {
    u[ii] <- ev[ii]*kern_mat_new[ii,1]
  }
  u_sum <- sum(u)
  return(u_sum)
}

# function to normalize an "a"-eigenvector matrix
# i.e. "a" is a vector which is representings an eigenvector with coeffieciencts 
# for linearly combining data points
# e = sum_k(a_k*datpt_k); to make norm(e) = 1, apply this function to "a"-vector
# takes in cur_eigen; lsit with [[1]] matrix of eigenvalues;
# [[2]] is matrix of "a" - eigenvectors
# returns a matrix with eigv which are normalized
agk.norm.a <- function(cur_eigen) {
  p      <- length(cur_eigen[[2]][1,])
  norm_a <- matrix(nrow=p,ncol=p)
  for (ii in 1:p) {
    cur_a     <- cur_eigen[[2]][,ii]
    cur_l     <- cur_eigen[[1]][ii]
    cur_norma <- (2/(sqrt(cur_l)*sqrt(cur_a%*%cur_a)))*cur_a
    norm_a[,ii] <- as.numeric(cur_norma)
  }
  norm_a <- as.matrix(norm_a)
  return(norm_a)
}

## Prune function - Milan, 19.05.2015
## This function takes a dataset and a z score related parameter and then decides: 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# (A) which pictures have to go, because they are not seen as members of their category according to the ratings?
#       pos: val ratings (> 0.5)
#       neg: val ratings (< 0.5)
#       neu: val aro dom ratings (should stay within a cube of [+0.5 - 0.5] of val,aro,dom
#   CUTOFF: 90% of the subjects should agree
#   example: picture X has "> 0.5" on positive rating in 9 out of 10 subjects, so it is indeed a pos picture and can stay
# (B) which pictures have to go, because they are not seen as members of their category according to the qualitative ratings?
#       gam,pos,neg,neu: check out fit_cluster variable
#   CUTOFF: 90% of the subjects should agree
#   example: picture X is a "gam" pic by definition but only 7 out of 10 subjects think so too: so it must go;
# (C) which pictures have to go because they don't meet either A OR B; 
#       only relevant for the pos, neg, neu pics
#
#the function returns the pruned data frames "data_pdt_pa", "data_pdt_pb", "data_pdt_pab";
# and lists of images to throw out: "throwOutA", "throwOutB", "throwOutAB";
#
# 22.06.2015 Update: the output shoould yield pictures that have to go by category and percentage of pictures
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



prune.df <- function(data_pdt,param, agree) {
  
  throwOutA=list()
  throwOutA$pos=c()
  throwOutA$neg=c()
  throwOutA$neu=c()
  throwOutB=list()
  throwOutB$pos=c()
  throwOutB$neg=c()
  throwOutB$neu=c()
  throwOutB$gam=c()
  data_pdt$pruneA=NA
  data_pdt$pruneB=NA
  for (i in unique(data_pdt$stim)){
    ## was before:
    image_locations=eval(which(data_pdt$stim==i))
    
    # should it not be this...?
    #image_locations=eval(which(data_pdt$stim==(unique(data_pdt$stim))[i]))
    ######################
    ######  PRUNE (A)  ###
    ######################
    
    valRes=data_pdt$valence[image_locations][!is.na(data_pdt$valence[image_locations])]
    
    if (!is.na(data_pdt$imageGroup[eval(image_locations[1])])) {
      
      if (data_pdt$imageGroup[eval(image_locations[1])]=='pos'){
        # if the number of subjects who rated the image higher then given parameter (correclty) / number of subjects who rated this image
        # is higher than agree
        if (length(valRes[valRes>param])/length(valRes)>agree){
          data_pdt$pruneA[image_locations]=1
        }
        else {
          throwOutA$pos=c(throwOutA$pos, i)
          data_pdt$pruneA[image_locations]=0
        }
      }
      if (data_pdt$imageGroup[eval(image_locations[1])]=='neg'){    
        # if the number of subjects who rated the image lower then given parameter (correclty) / number of subjects who rated this image
        # is higher than agree
        if (length(valRes[valRes<param*-1])/length(valRes)>agree){
          data_pdt$pruneA[image_locations]=1
        }
        else {
          throwOutA$neg=c(throwOutA$neg, i)
          data_pdt$pruneA[image_locations]=0
        }
      }
      if (data_pdt$imageGroup[eval(image_locations[1])]=='neu'){
        arousRes=data_pdt$arousal[image_locations][!is.na(data_pdt$arousal[image_locations])]
        domiRes=data_pdt$dominance[image_locations][!is.na(data_pdt$dominance[image_locations])]
        # if the number of subjects who rated the image within the limints of +-parameter (correclty) / number of subjects who rated this image
        # is higher than agree
        if ((length(valRes[valRes<param & valRes>param*-1])/length(valRes)>agree) &
              (length(arousRes[arousRes<param & arousRes>param*-1])/length(valRes)>agree) &
              (length(domiRes[domiRes<param & domiRes>param*-1])/length(valRes)>agree)) {
          data_pdt$pruneA[image_locations]=1
        }
        else {
          throwOutA$neu=c(throwOutA$neu, i)
          data_pdt$pruneA[image_locations]=0
        }
      }
      ######################
      ######  PRUNE (B)  ###
      ######################
      group=as.character(eval(data_pdt$imageGroup[eval(image_locations[1])]))
      if (group=='nap') {
        group='neu'
      }
      cluster=data_pdt$fit_cluster[image_locations][!is.na(data_pdt$fit_cluster[image_locations])]
      if (length(cluster)!=0) {
        if (length(cluster[cluster==group])/length(cluster) >agree) {
          data_pdt$pruneB[image_locations]=1
        }
        else {
          nam=paste0(group)
          throwOutB[[nam]]= c(eval(parse(text=paste0('throwOutB$',group))), i)
          data_pdt$pruneB[image_locations]=0
        }
      }   
    }
  }
  data_pdt_a <- subset(data_pdt, pruneA !=0 | is.na(pruneA))
  data_pdt_b <- subset(data_pdt, pruneB ==0 | is.na(pruneB))
  data_pdt_ab <- subset(data_pdt, pruneA !=0 | is.na(pruneA) & pruneB ==0 | is.na(pruneB))
  
  ## A function that concatinates two lists with some common elements appending the elements
  appendList <- function (x, val) 
  {
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
      x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
        appendList(x[[v]], val[[v]])
      else c(x[[v]], val[[v]])
    }
    x
  }
  throwOutAB=appendList(throwOutA, throwOutB)
  throwOutAB=lapply( throwOutAB, unique, simplify=FALSE)
  
  #To calculate the portion of images that were thrown out by each procedure I am dividing the number of images that we are throwing out
  #by the number of images that were processed, i.e. fullfill the criteria that:
  #for A their imageGroup variable is not an NA and that they are either pos neg or neu
  #for B that they got a fit_cluster value for at least one participant
  throwOutAportion=list()
  throwOutBportion=list()
  throwOutABportion=list()
  
  throwOutAportion$overall=Reduce('+',lapply(throwOutA, length))/length(unique(subset(data_pdt, pruneA !=0 | pruneA !=1)$stim))
  throwOutBportion$overall=Reduce('+',lapply(throwOutA, length))/length(unique(subset(data_pdt, pruneB !=0 | pruneB !=1)$stim))
  throwOutABportion$overall=Reduce('+',lapply(throwOutA, length))/length(unique(subset(data_pdt, pruneA !=0 | pruneA !=1 | pruneB !=0 | pruneB !=1)$stim))
  
  Subset=subset(data_pdt, pruneA !=0 | pruneA !=1)
  for (i in unique(Subset$imageGroup)){
    throwOutAportion[[i]]= length(throwOutA[[i]])/length(unique(Subset$stim[Subset$imageGroup==i]))
  }
  Subset=subset(data_pdt, pruneB !=0 | pruneB !=1)
  for (i in Subset$imageGroup){
    throwOutBportion[[i]]= length(throwOutB[[i]])/length(unique(Subset$stim[Subset$imageGroup==i]))
  }
  Subset=subset(data_pdt, pruneA !=0 | pruneA !=1 | pruneB !=0 | pruneB !=1)
  for (i in Subset$imageGroup){
    throwOutABportion[[i]]= length(throwOutAB[[i]])/length(unique(Subset$stim[Subset$imageGroup==i]))
  }
  
  out <- list(data_pdt_a =data_pdt_a, data_pdt_b = data_pdt_b, data_pdt_ab = data_pdt_ab, 
              throwOutA=throwOutA, throwOutB=throwOutB, throwOutAB=throwOutAB, 
              throwOutAportion=throwOutAportion, throwOutBportion=throwOutBportion, throwOutABportion=throwOutABportion)
  
  return (out)
}

#################################################################################

##################
## functions LA ##
##################

## all the functions and libraries I need
## functions
# get K-fold CV rule
f_K_fold <- function(Nobs,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

agk.plot.density <- function(name_df, vars,subtitle) {
  plots <- list()
  for (ii in 1:length(vars)) {
    base_expr <- paste('plots[[ii]] <- ggplot(',
                       name_df,',aes(x=',
                       vars[ii],'))',sep='')
    eval(parse(text=base_expr))
    hist_expr <- paste('plots[[ii]] <- plots[[ii]] +', 
                       'geom_histogram(aes(y = ..density..), binwidth=density(', name_df,'$',
                       vars[ii],',bw="ucv")$bw, color="white")',sep='')
    eval(parse(text=hist_expr))
    plots[[ii]] <- plots[[ii]] + geom_density(fill="red", alpha = 0.2) +
      theme_bw() +
      xlab(vars[ii]) +
      ylab('density')
    plots[[ii]] <- plots[[ii]] + ggtitle(paste("marginal density of", vars[ii],"\n",
                                               subtitle,                                               
                                               " "))
  }
  return(plots)
}

# return first element
first <- function(x) {return(x[1])}

agk.get.compl.coef <- function(x,int=NULL) {
  ## function to extract all fixef+ranef from glmer fit
  # int is a variable which was used for interaction in fixed effect
  ## THIS FUNCTION TO BE EXTENDED FOR MULTIPLE VARIABLES IN INTERACTION IN FIXED EFFECT AND RANEF
  # will only extract subject specific random effects; even if there are, e.g. stimulus specific ones
  # also 
  no_intercept = 0
  if(is.null(int)) {
    disp("Warning: no interaction variable given, returning coef")
    return(coef(x))
  }
  cur_coef   <- coef(x)
  cur_coef_c = cur_coef$subject
  if (isempty(grep('(Intercept)',names(cur_coef_c)))) {
    no_intercept = 1
    # in case no intercept then add and set to 0
    cur_coef_c$'(Intercept)' = 0
    cur_coef_c = cur_coef_c[c(length(cur_coef_c),1:(length(cur_coef_c)-1))]
    cur_coef$subject = cur_coef_c
  }
  cur_raneff <- names(cur_coef)
  cur_dat    <- x@frame
  cur_raneff_var <- eval(parse(text=paste("cur_dat$",cur_raneff,sep="")))
  cur_int  <- eval(parse(text=paste("cur_dat$",int)))
  lev_int  <- levels(cur_int)
  cur_con  <- contrasts(cur_int)
  f <- function(x) {return(!any(x == 1))}
  for (ii in 1:length(cur_con[,1])) {
    if (f(cur_con[ii,])==T) {
      cur_bl <- ii
      break
    }
  }
  cur_bl <- (row.names(cur_con))[cur_bl]
  cur_subs <- row.names(cur_coef$subject)
  
  grp_by_subs <- aggregate(cur_int,by = list(cur_raneff_var),FUN = "first")
  
  # which are the fixed effects variables which are not involved with interactions
  cur_fix  <- fixef(x)
  if (no_intercept == 1) {
    cur_fix = c(0,cur_fix)
    names(cur_fix)[1] = '(Intercept)'
  }
  ind_int  <- grep(int,names(cur_fix))
  nint     <- (names(cur_fix))[-ind_int]
  
  # which are the fixef which fit to which none-int fixef?
  int_nint_match <- list()
  for (kk in 1:length(nint)) {
    int_nint_match[[kk]] <- grep(nint[kk],names(cur_fix))
  }
  # fill the intercept (1) indices
  ul_int_nint_match <- unlist(int_nint_match)
  
  intercept_int <- c()
  for (ll in 1:length(cur_fix)) {
    if (any(ll == ul_int_nint_match)) {
      intercept_int[ll] <- 0
    } else {
      intercept_int[ll] <- 1
    }
  }
  add_ind <- which(intercept_int == 1)
  int_nint_match[[1]] <- c(int_nint_match[[1]],add_ind) 
  
  # which are the fixef which fit to which group?
  fix_int_match <- list()
  for (kk in 1:length(lev_int)) {
    fix_int_match[[kk]] <- grep(lev_int[kk],names(cur_fix))
  }
  
  out_coef <- matrix(0,nrow = length(grp_by_subs$x),ncol=length(nint))
  colnames(out_coef) <- nint
  for (ii in 1:length(cur_coef$subject[,1])) {
    cur_grp <- grp_by_subs$x[ii]
    # is it a control group subject?
    if (cur_grp == cur_bl) {
      out_coef[ii,] <- as.numeric(cur_coef$subject[ii,-ind_int])
    } else {
      # not a control subject
      acute_coef    <- cur_coef$subject[ii,-ind_int]
      all_cur_fixef <- cur_coef$subject[ii,]
      for (jj in 1:length(acute_coef)) {
        cur_nint     <- (names(acute_coef))[jj]
        cur_ind_nint <- which(cur_nint == nint)
        
        # which are the right additives according to nint?
        nint_add <- int_nint_match[[cur_ind_nint]]
        # which are the right additives according to int?
        int_add  <- fix_int_match[[which(cur_grp == lev_int)]]
        
        # adding
        tmp_ind <- Reduce(intersect, list(nint_add,int_add))
        # check if it is a non-interactive covariate
        if (length(tmp_ind)==0) {
          next
        }
        acute_coef[jj] <- acute_coef[jj] + all_cur_fixef[tmp_ind]
      }
      out_coef[ii,] <- as.numeric(acute_coef)
    }
  }
  
  out_coef       <- as.data.frame(out_coef)
  eval(parse(text=paste("out_coef$",int,"<-grp_by_subs$x",sep="")))
  out_coef$subject <- grp_by_subs$Group.1
  
  return(out_coef)
}

# normalize into 0-1 range
agk.normalize <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

# normalize generally to any range
agk.norm.range <- function (array, x, y) {
  # Normalize to [0, 1]:
  m = min(array);
  range = max(array) - m;
  array = (array - m) / range;
  
  # Then scale to [x,y]:
  range2 = y - x;
  normalized = (array*range2) + x;
  return(normalized)
}

# predict function for lmList object (lme4)
agk.predperc.lmList.logit <- function(x) {
  
  subs <- names(x)
  pred_perc <- c()
  
  for(ii in 1:length(subs)) {
    cur_call      <- cur_call <- x@.Data[[ii]]
    cur_resp      <- round(predict.glm(cur_call,type="response"))
    cur_dat       <- cur_call$data
    real_resp     <- na.omit(cur_dat$accept.reject)
    pred_perc[ii] <- mean(as.numeric(cur_resp == real_resp))
  } 
  
  cur_df <- data.frame(subs,pred_perc)
  
  return(cur_df)
  
}

# moving average
mav <- function(x,n=5){filter(x,rep(1/n,n), sides=2,circular = T)}

# function to assign every data point to nearest prototype
# dat goes in as df with every row as one obs
# prot goes in as every col is one cluster mean
# uses absolute distance
# m is assignment vector (1 is assigned, 0 if not)
# m must sum to 1 and has length of number of clust
agk.kmeans.assign <- function(dat,prot) {
  # compute euclidean distances
  cur_dist <- distmat(t(prot), as.matrix(dat))
  cur_dist <- as.data.frame(cur_dist)
  
  # make Matrix of vectors which indicate assignment
  # also a vector m which indicates the cluster assingment by number
  M <- data.frame(matrix(0,length(prot[1,]),length(dat[,1])))
  m <- c()
  for (ii in 1:length(cur_dist[1,])) {
    M[ii] <- ifelse(cur_dist[ii] == min(cur_dist[ii]),1,0)
    m[ii] <- which.min(as.numeric(cur_dist[[ii]]))
  }
  return(list(M=M,m=m,cur_dist=cur_dist))
} 



# check if there is loc and noa in the PCs
# function to check correlation of boxcar and found component
# fits multinomial regression and returns predictive value
agk.corboxmetric <- function (cur_metr,cur_box) {
  cur_loc <- as.factor(cur_box)
  mod <- multinom(cur_loc ~ cur_metr)
  return(mean(as.numeric(predict(mod) == cur_loc)))
}

agk.entr <- function(x){
  # entropy of vector
  # x is a vector of probabilities
  x <- agk.normalize(x)
  x <- x+0.01
  x <- x/sum(x)
  entr <- -sum(x*log(x))
  entr_st <- -sum(x*log(x))
  return(entr)
}

agk.entr.st <- function(x){
  # normalized entropy (standardized by maximal entropy)
  # x is a vector of probabilities
  x <- agk.normalize(x)
  x <- x+0.01
  x <- x/sum(x)
  entr <- -sum(x*log(x,base=length(x)))
  return(entr)
}

agk.entr.x = function(emp_x) {
  # function to return the entropy of x
  # x is just a sample of a distribution
  # density estimation will be don on fly
  
  # numeric?
  if(!is.numeric(emp_x)) {
    return(NA)
  }
  
  # density of x
  dens_x_emp = density(emp_x)
  
  # get the observed entropy
  x_entr = agk.entr(dens_x_emp$y)
  
  return(x_entr)
}

agk.entr.x.st = function(emp_x) {
  # function to return the stand. entropy of x
  # x is just a sample of a distribution
  # density estimation will be don on fly
  
  # numeric?
  if(!is.numeric(emp_x)) {
    return(NA)
  }
  
  # density of x
  dens_x_emp = density(emp_x)
  
  # get the observed entropy
  x_entr = agk.entr.st(dens_x_emp$y)
  
  return(x_entr)
}



# function (for "apply") to impute vecmean for NA
# give a vector
agk.impute.mean <- function(x) {
  cur_mean <- mean(x,na.rm = T)
  cur_fun <- function(y,cur_mean) {
    if (is.na(y)) {
      y<-cur_mean
    } 
    return(y)
  }
  return(sapply(x,cur_fun,cur_mean=cur_mean))
}

# function to make the loss of control and negelct of other areas component vector
agk.make.locnoa <- function(var_names,loc_vecs,noa_vecs) {
  
  zeroone_loc <- list()
  zeroone_noa <- list()
  
  for (kk in 1:length(loc_vecs)) {
    ind      <- loc_vecs[[kk]]
    cur_locn <- rep(0,length(var_names[[kk]]))
    for(ii in 1:length(cur_locn)) {
      if(any(ii == abs(ind))) {
        cur_var <- which(ii==abs(ind))
        cur_locn[ii] <- (1*sign(ind[cur_var]))}
    }
    zeroone_loc[[kk]] <- cur_locn 
  }
  
  for (kk in 1:length(noa_vecs)) {
    ind      <- noa_vecs[[kk]]
    cur_noan <- rep(0,length(var_names[[kk]]))
    for(ii in 1:length(cur_noan)) {
      if(any(ii == abs(ind))) {
        cur_var <- which(ii==abs(ind))
        cur_noan[ii] <- (1*sign(ind[cur_var]))}
    }
    zeroone_noa[[kk]] <- cur_noan 
  }
  
  loc <- unlist(zeroone_loc)
  noa <- unlist(zeroone_noa)
  
  return(list(loc=loc,noa=noa))
  
}

# log regr function for two pred vars
prob.log.reg <-function(int,b1,b2,x1,x2){
  
  return (1/(1+exp(-(int+x1*b1+x2*b2))))
}

# predict function
mylogit_predict <- function (cur.data,cur.coeff) {
  cur.pred <- c()
  for (ii in 1:length(cur.data[,1])) {
    tmp <- prob.log.reg(cur.coeff$Gewinn, cur.coeff$Verlust, cur.coeff$intercept, cur.data$Gewinn[ii], cur.data$Verlust[ii])
    cur.pred <- rbind(cur.pred, as.numeric(tmp))    
  }
  return (cur.pred)
}

# mean function that per default removes NA
mean.rmna <- function(x) {
  mean(x,na.rm=T)
}

# modulus
mod<-function(x,m){
  t1<-floor(x/m)
  return(x-t1*m)
}

# calc the mv normality

check_mvnorm <- function(v1,v2) {
  require(mvnormtest)
  v1 <- v1[!is.na(v1)]
  v1 <- v1[!is.na(v2)]
  v2 <- v2[!is.na(v1)]
  v2 <- v2[!is.na(v2)]
  return (mshapiro.test(t(as.matrix(data.frame(v1,v2)))))
}

# Spearman-Brown-Formula
# r is current correlation
# f is factor by which test gets prolonged (e.g. 2 in split-half case)
spearman_brown <- function(r,f) {
  rel = f*r/(1+(f-1)*r)
  return(rel)
}

# functions for check effect of randomness

## a function that asks for a 0,1 vector and will randomly pick perc% elements
## and replace them with random 0 or 1
v_ran <- function(v1,perc) {
  
  how_many <- round((perc/100)*length(v1)) # how many elements should be replaced
  where    <- sample(1:length(v1), how_many, replace=F) # get indices of where to replace
  v1[where]<- sample(0:1,how_many, replace=T) # fill in random 0' and 1's in v1
  
  return (v1)
  
}

## loss aversion: calculate lambda
calc_lambda<-function(bg,bl){
  return(bl*(-1)/bg)
}

# cross validation function for a glmer model
# given a single random effect unit (subject)
# k-fold (folding over single random effect unit)
# 10% subjects will be left out and the data fit
# but during prediction we only use the fixed effects;
agk.cv.glmer <- function(model,byvar,k) {
  
  # get data
  dat <- model@frame
  
  # get the levels of the byvar
  lev <- eval(parse(text=paste("levels(dat$",byvar,")",sep="")))
  
  # create folds
  flds <- f_K_fold(length(lev),k)
  f = function(.,x) {any(.==x)}
  
  # get the column that has the byvar that the folds will be done over ("subject")
  bycol <- eval(parse(text=(paste("as.matrix(dat$",byvar,")",sep=""))))
  
  pred_score <- c()
  for (ii in 1:k){
    print(paste("fold...",ii))
    
    # create train data
    ind_train  <- flds[[ii]]$train
    ind_train  <- apply(bycol,1,f,x=lev[ind_train])
    data_train <- dat[ind_train,]
    
    # create test data
    ind_test   <- flds[[ii]]$test
    ind_test   <- apply(bycol,1,f,x=lev[ind_test])
    data_test  <- dat[ind_test,]
    
    # training
    curr_model <- update(model,data = data_train)
    
    # testing
    curr_targ  <- data_test$accept.reject
    preds      <- predict(curr_model,newdata = data_test,type="response",re.form=~0)
    preds      <- round(preds)
    cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
    cor.preds  <- ifelse(cor.preds==0,1,0)
    
    pred_score[ii] <- mean(cor.preds,na.rm = T)
  }
  EG    <- mean(pred_score)
  EG_sd <- sd(pred_score)
  EG_list <- list(EG=EG,EG_sd=EG_sd,pred_score=pred_score)
  return(EG_list)
}

agk.getflds <- function(k,byvar,model) {
  
  # get data
  dat <- model@frame
  
  # get the levels of the byvar
  if (!is.null(byvar)) {
    lev <- eval(parse(text=paste("levels(dat$",byvar,")",sep="")))
    
    # create folds
    flds <- f_K_fold(length(lev),k)
    
  } else {
    flds <- f_K_fold(length(dat[,1]),k)
  }
  return(flds)
}

# cross validation function for a glmer model
# given a single random effect unit (e.g. subject)
# this unit is the natural unit of the experiment
# k-fold CV can now be done using only fixed effects (
# predicting a completely new subject)
# or by using k-fold observations CV; all subs 
# are always used, just observations are split into
# training and test (set byvar to "observations")
# you have to provide the flds
# which is a list having on every entry the indices
# of byvar for training and test
agk.setcv.glmer <- function(model,byvar,flds) {
  
  # get data
  dat <- model@frame
  
  # get k
  k <- length(flds)
  
  if (byvar == "observations") {
    dat$observations <- as.factor(1:length(dat[,1]))
  }
  
  # get the levels of the byvar
  lev <- eval(parse(text=paste("unique(dat$",byvar,")",sep="")))
  
  # prep function for train/test set creation
  f = function(.,x) {any(.==x)}
  
  # get the column that has the byvar that the folds will be done over ("subject")
  bycol <- eval(parse(text=(paste("as.matrix(dat$",byvar,")",sep=""))))
  
  pred_score <- c()
  for (ii in 1:k){
    print(paste("fold...",ii))
    
    # create train data
    ind_train  <- flds[[ii]]$train
    ind_train  <- apply(bycol,1,f,x=lev[ind_train])
    data_train <- dat[ind_train,]
    
    # create test data
    ind_test   <- flds[[ii]]$test
    ind_test   <- apply(bycol,1,f,x=lev[ind_test])
    data_test  <- dat[ind_test,]
    
    # training
    curr_model <- update(model,data = data_train)
    
    # target
    cur_call    = curr_model@call
    cur_formula = cur_call$formula
    cc          = as.character(cur_formula)
    ct          = cc[2]
    
    # testing
    curr_targ  <- data_test[[ct]]
    if (byvar == "observations") {
      preds      <- predict(curr_model,newdata = data_test,type="response",re.form=NULL)
    } else {
      preds      <- predict(curr_model,newdata = data_test,type="response",re.form=~0)
    }
    preds      <- round(preds)
    cor.preds  <- curr_targ == preds
    
    # pred score (the higher the better)
    pred_score[ii] <- mean(cor.preds,na.rm = T)
  }
  EG    <- mean(pred_score)
  EG_sd <- sd(pred_score)
  EG_list <- list(EG=EG,EG_sd=EG_sd,pred_score=pred_score)
  return(EG_list)
}



## gets a subsample of trials (rows) split by "variable" (e.g. "subject")

make_sample_of_trials <- function(variable, data.la, n,check,same,verbosity) {
  got_rand <- 0
  subs <- levels(eval(parse(text=paste("data.la$",variable,sep = ""))))
  new_data <- data.frame()
  for (ii in 1:length(subs)) {
    if (verbosity == 1) {
      print(paste("now running...", subs[ii]))
    }
    cur.data <- data.la[eval(parse(text=paste("data.la$", variable,sep = ""))) == subs[ii],]
    # save this data frame as baseline for later represent. check
    cur.d.sv <- cur.data
    cur.ind  <- seq(1,length(cur.data[,1]))
    repr_result <- 0
    # make a new sample as long as not representative
    count <- 0
    while (repr_result == 0 & count < 3000){
      if (verbosity == 1) {
        print(count)
      }
      count <- count + 1
      if (got_rand == 0) {
        cc.ind   <- sample(cur.ind, n,replace = F)
        # same random sample for all subjects
        if (same == 1) {got_rand <- 1}
      }
      cur.data <- cur.d.sv[cc.ind,]
      if (check == 0) {
        repr_result <- 1
      } else {repr_result <- check_representativeness(cur.data,cur.d.sv,0.05)}
    }
    if (count >= 3000) {print(paste("attention, did not find adequate sample in...", subs[ii]))}
    new_data <- rbind(new_data,cur.data)
  }
  return(new_data)
}

## make sample of trials but based on gambles; sample of certain gambles is picked and then this is been looked for in the data and taken

make.gamble.sample <- function(variable, data.la, n,verbosity) {
  subs <- levels(eval(parse(text=paste("data.la$",variable,sep = ""))))
  new_data <- data.frame()
  count <- 0
  for (ii in 1:length(subs)) {
    count <- count +1
    if (verbosity == 1) {
      print(paste("now running...", subs[ii]))
    }
    cur.data <- data.la[eval(parse(text=paste("data.la$", variable,sep = ""))) == subs[ii],]
    # save this data frame as baseline for later represent. check
    cur.d.sv <- cur.data
    # here now get the random gamble sample
    if (count == 1) {
      gamble_sample <- sample.gambles(n,data.la)
    }
    # here we pull the desired gamble sample from the current subject's data
    cur.data <- get.sample.trials(cur.d.sv,gamble_sample)
    new_data <- rbind(new_data,cur.data)
  }
  return(list(new_data,gamble_sample))
}

## function to decide for sample of gambles

sample.gambles <- function(n_gam,data.la) {
  pos_gain <- as.numeric(levels(as.factor(as.character(data.la$Gewinn))))
  pos_loss <- as.numeric(levels(as.factor(as.character(data.la$Verlust))))
  comvec   <- agk.combvec(pos_gain,pos_loss)
  which_tr <- comvec[sample(length(1:length(comvec[,1])),n_gam),]
  return(which_tr)
}

## function that takes the gamble_sample and takes the certain trials out of the data.la
get.sample.trials <- function(data.la,gamble_sample) {
  new_data <- data.frame()
  for(ii in 1:length(gamble_sample[,1])) {
    cur_new <- data.la[data.la$Gewinn == gamble_sample[ii,1] & data.la$Verlust == gamble_sample[ii,2],]
    new_data <- rbind(new_data,cur_new)
  }
  return(new_data)
}

## comvec function
agk.combvec <- function(x,y){
  out_m <- matrix(0,nrow = length(x)*length(y),ncol=2)
  count <- 0
  for (ii in 1:length(x)) {
    for (jj in 1:length(y)) {
      count <- count+1
      out_m[count,1] <- x[ii]
      out_m[count,2] <- y[jj]
    }
  }
  return(out_m)
}

## function to check whether the drawn sample is representative

check_representativeness <- function(cur.data.la, data.la,alpha) {
  check_result <- 1
  
  desired_variables <- list("ev","ed", "ed.abs","ratio","diff","RiskMar")
  
  for (ii in 1:length(desired_variables)) {
    # testing a against b baseline
    cur.a <- eval(parse(text=paste("cur.data.la$",desired_variables[[ii]],sep = "")))
    cur.b <- eval(parse(text=paste("data.la$",desired_variables[[ii]],sep = "")))
    
    # test
    cur.t <-t.test(cur.data.la$ev,data.la$ev)
    cur.t <- cur.t$p.value
    if (cur.t < alpha) {
      check_result <- 0
      return(check_result)}
    
    cur.t <- cur.t <-bartlett.test(list(cur.data.la$ev,data.la$ev))
    cur.t <- cur.t$p.value
    if (cur.t < alpha) {
      check_result <- 0
      return(check_result)}
  }
  # all tests passed
  return(check_result)
}

## for beginning plotting

## setting the colour friendly pallette
# The palette with black:
#cbbPalette   = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette   = grey.colors(3, start = 0.3, end = 0.8, gamma = 2.2, alpha = NULL)
cbbPalette   = c(cbbPalette,"#000000")
Palette_heat = grey.colors(3, start = 0.0, end = 1.0, gamma = 2.2, alpha = NULL)


# To use for line and point colors, add
scale_colour_manual(values=cbbPalette)

# make the theme (based on theme_bw)
theme_la <- function(base_size = 8, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background=element_rect(fill="white", colour="black", size=2.0),
      panel.grid = element_blank(),
      plot.background = element_rect(fill="white", color = "white"),
      plot.title = element_text(colour = "black",size=25,face="bold",vjust=2),
      axis.text = element_text(colour = "black",size=25,face="bold"),
      axis.ticks = element_line(colour ="black", size=1.5),
      axis.title.x = element_text(colour = "black",size=25,face="bold",vjust=0.2),
      axis.title.y = element_text(colour = "black",size=25,face="bold",vjust=2.0,angle = 90),
      plot.margin = unit(c(1,0.5,1,1),"cm"),
      #legend.position = c(0.92, 0.83),
      legend.background = element_rect(fill="white", color = "white"),
      legend.text = element_text(colour = "black",size=20,face="bold"),
      legend.title = element_text(colour = "black",size=20,face="bold"),
      #legend.margin=unit(.85,"cm"),
      legend.text.align=0,
      strip.text.x = element_text(size=25,face="bold",vjust=100),
      strip.text.y = element_text(size=25,face="bold",vjust=100),
      strip.background = element_rect(colour="white", fill="white")
    )
}

## function to paste a short variable to a long data.frame given an extending variable
# e.g. age is a variable per subject; does not vary per trials, but in long format i need to have it per trial
# in this version long and short var must be ordered by extend var already!
paste_short_on_long <- function (long, short_var, extend_var, new_name) {
  sub_var_expr <- parse(text=paste("data.la$",extend_var,sep = ""))
  subs <- levels(eval(sub_var_expr))
  new.var <- c()
  for (ii in 1:length(subs)) {
    cur.trials <- long[eval(sub_var_expr) == subs[ii],]
    cur.newvar <- rep(short_var[ii],length(cur.trials[,1]))
    new.var <- c(new.var,as.numeric(cur.newvar))
  }
  new.var <- as.data.frame(new.var)
  names(new.var) <- new_name
  new_long <- data.frame(long,new.var)
}

# # function that takes an lmLIST object and sets all its coefs zero which are not significant
# agk.prune.lmList = function(lml) {
#   cur_coef = coef(lml)
#   cur_p    = summary(lml)
#   cur_p    = cur_p$coefficients
#   for (ii = 1:length(cur_p)) {
#     
#   }
# }

# subfunction that does stepwise variable selection using a lm or glm model
agk.step.lm = function(lml,fam_arg) {
  cur_mod=lml
  # if model is empty
  if (is.null(lml)) {
    return(NA)
  }
  cur_data = cur_mod$data
  n = dim(cur_data)[1];
  cur_form = cur_mod$formula
  cur_form = as.character(cur_form)
  cur_form = paste(cur_form[2],cur_form[1],cur_form[3])
  if (fam_arg == "binomial") {
    cur_cmd  = paste("step(glm(",cur_form,",data=cur_data,family='binomial'),direction='both',k=log(n))")
  } else {
    cur_cmd  = paste("step(lm(",cur_form,",data=cur_data),direction='both')")
  }
  mod_pr   = eval(parse(text=cur_cmd))
  return(mod_pr$coefficients)
}

agk.glmnet.lm = function(lml,fam_arg,type_measure,alpha,lambda, verbosity = 0,
                         nfolds = 10) {
  # subfunction that does variable selection by using glmnet
  # on the lm or glm (via CV for lambda; alpha is also high-level cv'd)
  # nfolds is an optional argument for cv.glmnet
  # to determine how many folds should be used for lambda CV
  # set to "LOOCV" if leave one out CV is desired
  
  # the current model
  cur_mod  = lml
  # if model is empty
  if (is.null(lml)) {
    return(NULL)
  }
  
  # verbosity
  if (verbosity) {
    # display which subject this is
    disp(as.character(lml$data$subject[1]))
  }
  
  # extract data, formula, etc.
  cur_data = cur_mod$data
  cur_form = cur_mod$formula
  cur_form = as.character(cur_form)
  cur_form = paste(cur_form[2],cur_form[1],cur_form[3])
  
  # get the desired nfolds
  if (missing(nfolds)) {
    cnfolds = 10
  } else if (nfolds == "LOOCV") {
    cnfolds = length(cur_data[,1])
  } else {
    cnfolds = nfolds
  }
  
  # special case LOOCV
  if (nfolds == "LOOCV") {
    # perform fixed LOOCV, meaning foldid is always set as seq(1:cnfolds)
    # useful to reproduce results but slow
    folds_string = paste0("',nfolds=cnfolds",",foldid=seq(1:cnfolds)")
  } else {
    folds_string = "',nfolds=cnfolds"
  }
  
  # different model families
  # if (fam_arg == "gaussian") {
  #   cur_data = cur_mod$model
  #   cur_form = eval(parse(text=as.character(cur_mod$call$formula)))
  #   cur_form = as.character(cur_form)
  #   cur_form = paste(cur_form[2],cur_form[1],cur_form[3])
  # }
  
  # do cross validation on the model
  cur_cmd  = paste("glmnetUtils::cv.glmnet(formula=",cur_form,",data=cur_data,grouped=FALSE,",
                   "lambda=lambda,use.model.frame=T,family=fam_arg,alpha=",alpha,
                   ",type.measure='",type_measure,folds_string,")",sep="")
  tmp <-eval(parse(text=cur_cmd))
  
  # return the resulting models and the winning coefficients
  cur_coef = t(as.matrix(coef(tmp,s="lambda.min")))
  row.names(cur_coef) = NULL
  return(list(cur_coef = cur_coef, cur_cvres = tmp))
}

## OLD
# agk.glmnet.lm.cvalpha = function(lml,fam_arg,type_measure,lambda,
#                                  alphas = c(0,0.05,0.3,1),verbosity = 0,
#                                  nfolds = 10) {
#   # grid search alpha in agk.glmnet.lm
#   # wrapper function to do 
#   # grid search for best alpha value
#   # alpha == 0 is ridge
#   # alpha == 1 is lasso
#   # nfolds is an optional argument for
#   # agk.glmnet.lm to say how many
#   # CV folds should be used to 
#   # decide which lambda to use
#   # set to "LOOCV" if leave one out CV
#   # is desired
#   
#   # get the desired nfolds
#   if (missing(nfolds)) {
#     cnfolds = 10
#     } else {
#     cnfolds = nfolds
#     }
#   
#   # if model is empty
#   if (is.null(lml)) {
#     return(NULL)
#   }
#   
#   # run cv with different alphas
#   scores_m = c()
#   res      = list()
#   for (jj in 1:length(alphas)) {
#     if (verbosity) {
#       disp(paste("I am testing alpha = ",alphas[jj]))
#     }
#     cur_alpha    = alphas[jj]
#     res[[jj]]    = agk.glmnet.lm(lml,fam_arg,type_measure,cur_alpha,
#                                  lambda,nfolds = cnfolds)
#     tmp          = res[[jj]]$cur_cvres
#     scores_m[jj] = tmp$cvm[which(tmp$lambda.min == tmp$lambda)]
#   }
#   
#   # getting the smallest error across all alphas
#   # if multiple, then the last (preferring sparse models)
#   # WRONG: cause late lambdas are low and then the model is less sparse
#   cur_pick = which(min(scores_m)==scores_m)
#   return(list(coef = res[[cur_pick[length(cur_pick)]]]$cur_coef, cvme = scores_m,
#               winning_model = res[[cur_pick[1]]]$cur_cvres,winning_alpha=alphas[cur_pick]))
# }
## OLD

agk.glmnet.lm.cvalpha = function(lml,fam_arg,type_measure,lambda,
                                 alphas = c(0,0.05,0.3,1),verbosity = 0,
                                 nfolds = 10) {
  # grid search alpha in agk.glmnet.lm allowing to input a fit model
  # using cva.glmnet
  # wrapper function to do 
  # alpha == 0 is ridge
  # alpha == 1 is lasso
  # nfolds is an optional argument for
  # agk.glmnet.lm to say how many
  # CV folds should be used to 
  # decide which lambda to use
  # set to "LOOCV" if leave one out CV
  # is desired
  
  # the current model
  cur_mod  = lml
  # if model is empty
  if (is.null(lml)) {
    return(NULL)
  }
  
  # verbosity
  if (verbosity) {
    # display which subject this is
    disp(as.character(lml$data$subject[1]))
  }
  
  # extract data, formula, etc. x,y
  cur_data = cur_mod$data
  cur_y    = cur_mod$y
  cur_form = cur_mod$formula
  cur_form = as.character(cur_form)
  if (!isempty(grep('(0 +|0+|-1|- 1)',cur_form[3]))) {
    # glmnet needs a different way to suppress intercept
    cur_form[3] = gsub('(0 +|0+|-1|- 1)','', cur_form[3])
    cur_intercept = F
  } else {
    cur_intercept = T
  }
  cur_form = as.formula(paste(cur_form[2],cur_form[1],cur_form[3]))
  cur_x    = model.matrix.lm(cur_form,cur_data)

  # get the desired nfolds
  if (missing(nfolds)) {
    cnfolds = 10
  } else if (nfolds == "LOOCV") {
    cnfolds = length(cur_data[,1])
  } else {
    cnfolds = nfolds
  }
  
  # run cv with different alphas
  # use model frame F otherwise two-predictor without intercept wont work
  scores_m  = c()
  res       = list()
  cur_alpha = alphas
  all_res   = list()
  for (rr in 1:ridge_bv_reps) {
    res           = glmnetUtils::cva.glmnet(cur_form,family=fam_arg,
                                        type.measure=type_measure,data= cur_data,
                                        alpha = alphas,nfolds = cnfolds,
                                        use.model.frame=T,grouped=doGrouped,lambda = exp(seq(log(0.001), log(10), length.out=100)),
                                        standardize = des_stand_glmnet,intercept = cur_intercept)
    res           = agk.glmnetUtils.cvalpha.getwinner(res)
    all_res[[rr]] = res
  }
  
  cur_coefs = all_res[[1]]$coef
  all_cvmes = all_res[[1]]$cvme
  
  if  (ridge_bv_reps > 1) {
    for (rr in 2:ridge_bv_reps) {
      cur_coefs = rbind(cur_coefs,all_res[[rr]]$coef)
      all_cvmes = c(all_cvmes,all_res[[rr]]$cvme)
    }
    cur_coefs = colMedians(cur_coefs)
    all_cvmes = median(all_cvmes)
  }
  
  # drop useless intercept column
  if (cur_intercept == F) {
    cur_coefs = cur_coefs[-which(colnames(cur_coefs) == '(Intercept)')]
  }
  
  return(list(coef = cur_coefs, cvme = all_cvmes,
              winning_model = NULL,winning_alpha=NULL))
}

agk.glmnetUtils.cvalpha.getwinner = function(cur_cvmoda) {
  # extract winner from cvalpha (glmnetUtils) mod
  cur_res = list()
  # get winning alpha
  if (cur_cvmoda$modlist[[1]]$name == 'AUC') {
    agk.glmnetUtils.cvalpha.getwinner.getalpha = function(cur_cvmoda) {
      return(max(cur_cvmoda$cvm))
    }
  } else {
    agk.glmnetUtils.cvalpha.getwinner.getalpha = function(cur_cvmoda) {
      return(min(cur_cvmoda$cvm))
    }
  }

  cur_mins               = unlist(lapply(cur_cvmoda$modlist,FUN=agk.glmnetUtils.cvalpha.getwinner.getalpha))
  # preferring higher alphas (more sparse models)
  if (cur_cvmoda$modlist[[1]]$name == 'AUC') {
    ind_alpha              = which(max(cur_mins) == cur_mins)[length(which(max(cur_mins) == cur_mins))]
  } else {
    ind_alpha              = which(min(cur_mins) == cur_mins)[length(which(min(cur_mins) == cur_mins))]
  }
  cur_res$winning_alpha  = cur_cvmoda$alpha[ind_alpha]
  
  
  # now get the best lambda
  if (cur_cvmoda$modlist[[1]]$name == 'AUC') {
    cur_res$cvme           = min(1-cur_cvmoda$modlist[[ind_alpha]]$cvm)
  } else {
    cur_res$cvme           = min(cur_cvmoda$modlist[[ind_alpha]]$cvm)
  }
  cur_res$coef           = t(as.matrix(coef(cur_cvmoda$modlist[[ind_alpha]],s='lambda.min')))
  
  # get the winning model (winning alpha cvmod)
  cur_res$winning_model  = cur_cvmoda$modlist[[ind_alpha]]
  # return
  return(cur_res)
}


# mean function that returns 0 for numeric(0) 
# and replaces in a vector all NA's with 0 before mean
agk.mean.zero <- function (x) {
  if (length(x) == 0) {return(0)} 
  x[is.na(x)] = 0; x[is.nan(x)] = 0;
  return(mean(x))
}


## payout (BGG style)

agk.calc.payout <- function(cur.data, extend_var, base_value,n) {
  sub_var_expr <- parse(text=paste("cur.data$",extend_var,sep = ""))
  eval(parse(text=paste("cur.data$",extend_var,"<-as.factor(as.character(cur.data$",extend_var,"))",sep = "")))
  subs <- levels(eval(sub_var_expr))
  sub.payout <- c()
  for (ii in 1:length(subs)) {
    cur.trials <- cur.data[eval(sub_var_expr) == subs[ii],]
    cur.acc    <- subset(cur.trials,Button==1)
    cur.ind    <- seq(1,length(cur.acc[,1]))
    cur.ind    <- sample(cur.ind,n)
    cur.gam    <- cur.acc[cur.ind,]
    cur.res    <- c()
    for (kk in 1:length(cur.gam[,1])) {
      cur.res[kk] <- agk.local.payout(cur.gam$Gewinn[kk], cur.gam$Verlust[kk]) 
    }
    sub.payout[ii] <- sum(as.numeric(cur.res)) + base_value
  }
  return(sub.payout)
}

agk.local.payout <- function(gain,loss) {
  if (runif(1)>0.5) {
    return(gain)} else {return(loss)}
}

## function to plot the stability of lambda per trial
agk.lambda.stab <- function(cur.data, extend_var,start_point,end_point) {
  lambda_evolution <- c()
  for (ii in start_point:end_point) {
    cur.data.loop <- agk.sel.trials(cur.data,extend_var,ii)
    cur.lambdas   <- agk.lmlist.data.la(cur.data.loop)
    cur.lambdas   <- cur.lambdas$Lambda_ed
    lambda_evolution <- cbind(lambda_evolution,cur.lambdas)
  }
  for (kk in 1: length(tmp[,1])) {
    plot(as.numeric(tmp[kk,]),type="line") 
  }
  return(lambda_evolution)
}

## function to assess stability of model by prediciting next vali_trials on every step
agk.lambda.stab.pred <- function(cur.data, extend_var,start_point,end_point,vali_trials) {
  pred_evolution   <- c()
  collect_vali     <- list()
  collect_train    <- list()
  collect_lmlist   <- list()
  count            <- 0
  for (ii in start_point:end_point) {
    count <- count +1
    cur.data.loop <- agk.sel.trials(cur.data,extend_var,1,ii)
    collect_train[[count]] <- cur.data.loop
    cur.data.loop <- subset(cur.data.loop,!is.na(accept.reject))
    collect_vali[[count]]  <- agk.sel.trials(cur.data,extend_var,ii+1,ii+vali_trials)
    collect_lmlist[[count]] <- agk.lmlist.data.la.list(cur.data.loop)
  }
  # calc and plot pred evolution
  # go in every trial step
  for (kk in 1: length(collect_lmlist)) {
    # at trial step, go in every subject
    cur.pred <- c()
    for(ll in 1:length(collect_lmlist[[kk]])) {
      # make subject variable generic
      collect_vali[[kk]]$subject <- as.factor(as.numeric(collect_vali[[kk]]$subject))
      tmp.newdata <- collect_vali[[kk]][as.numeric(collect_vali[[kk]]$subject) ==ll,]
      tmp <- table((round(predict(collect_lmlist[[kk]][[ll]],
                                  newdata=tmp.newdata,"response"))-as.numeric(as.character(tmp.newdata$accept.reject)))==0)
      cur.pred[ll] <- as.numeric(tmp["TRUE"]/(tmp["TRUE"]+tmp["FALSE"]))
    }
    pred_evolution <- cbind(pred_evolution,cur.pred)
  }
  # return the calculations
  return(list(pred_evolution))
}

## function to get a data.la and get via lm.list the lambda vector out of it
## watch out; in "cur.data" the extend_variable must be called "subject"
agk.lmlist.data.la <- function(cur.data) {
  lmlist_ed <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = cur.data, na.action = na.omit, family = "binomial")
  crm_lmlist <- c()
  for (ii in 1:length(lmlist_ed)) {
    crm_lmlist <- rbind(crm_lmlist,t(as.matrix(as.numeric(lmlist_ed[[ii]][[1]]))))
  }
  crm_lmlist <- as.data.frame(crm_lmlist)
  names(crm_lmlist) <- c("intercept", "ed.abs", "Gewinn", "Verlust")
  crm_lmlist$id       <- as.factor(as.numeric(as.matrix(names(lmlist_ed))))
  crm <- crm_lmlist
  crm$Lambda_ed <- calc_lambda(crm$Gewinn,crm$Verlust)
  return(crm$Lambda_ed)
}

## function to get a data.la and get via lm.list the whole lmlist
agk.lmlist.data.la.list <- function(cur.data) {
  lmlist_ed <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = cur.data, na.action = na.omit, family = "binomial")
  return(lmlist_ed)
}


## function to just select x trials from data.la (every row is a trial)
agk.sel.trials <- function(cur.data, extend_var,start_n,end_n) {
  sub_var_expr <- parse(text=paste("cur.data$",extend_var,sep = ""))
  eval(parse(text=paste("cur.data$",extend_var,"<-as.factor(as.character(cur.data$",extend_var,"))",sep = "")))
  subs <- levels(eval(sub_var_expr))
  data.out <- data.frame()
  for (ii in 1:length(subs)) {
    cur.trials <- cur.data[eval(sub_var_expr) == subs[ii],]
    cur.trials <- cur.trials[start_n:end_n,]
    data.out   <- rbind(data.out,cur.trials)
  }
  return(data.out)
}

# a recode function
# function to recode x, given a source vector y and a translated vector z
agk.recode <- function(x,y,z) {
  for (ii in 1:length(x)) {
    done <- 0
    for (jj in 1:length(y)) {
      # NA in x will be NA
      if(is.na(x[ii])) {
        x[ii] <- NA
        next
      }
      if (x[ii] == y[jj]) {
        x[ii] <- z[jj]
        done <- 1
      }
      if (done == 1) {next}
    }
  }
  return(x)
}

# function to downsample x, keeping every n-th position
# careful: the signal will be shortened by this!
agk.downsample <- function (x,n) {
  core <- c(1, rep(0,n-1))
  which_samples <- rep(core,length.out = length(x))
  all_samples   <- 1:length(x)
  x_out         <- x[as.logical(which_samples)]
  return(x_out)
}


# function to aggregate data.la (gain and loss)
# levels: 1,2,3 
agk.aggregate.data.la.gain.loss <- function(v,aggr) {
  steps  = sort(unique(v))
  new_st = agk.downsample(steps, aggr)
  n      = length(steps)/length(new_st)
  new_st = rep(new_st,each=n)
  
  # now recode
  v_ds  = as.numeric(agk.recode.c(v,steps,new_st))
  return(v_ds)
}

agk.aggregate.data.la.gain.loss.c = cmpfun(agk.aggregate.data.la.gain.loss)

# function for extracting the estimates per subject
agk.get.subef <- function(la05,data.la) {
  
  # make group vector
  tmp <- suppressWarnings(aggregate(data.la, by=list(data.la$subject, data.la$group), FUN=mean))
  HC<-rep("HC",as.numeric(table(tmp$Group.2)[2]))
  AD<-rep("AD",as.numeric(table(tmp$Group.2)[1]))
  PG<-rep("PG",as.numeric(table(tmp$Group.2)[3]))
  group <- rbind(as.matrix(HC),as.matrix(PG),as.matrix(AD))
  
  # get intercept, gain, loss of HC fixef
  intercept <- matrix(rep(fixef(la05)["(Intercept)"], length.out=nrow(tmp)))
  gain      <- matrix(rep(fixef(la05)["Gewinn"], length.out=nrow(tmp)))
  loss      <- matrix(rep(fixef(la05)["Verlust"], length.out=nrow(tmp)))
  
  # get group specific int
  int.AD    <- matrix(rep(fixef(la05)["groupAD>HC"], as.numeric(table(group)["AD"])))
  int.PG    <- matrix(rep(fixef(la05)["groupPG>HC"], as.numeric(table(group)["PG"])))
  int.HC    <- matrix(rep(0,as.numeric(table(group)["HC"])))
  int.group <- rbind(int.HC, int.PG, int.AD)
  intercept <- intercept + int.group
  
  # get group specific gain
  gain.HC   <- matrix(rep(0,as.numeric(table(group)["HC"])))
  gain.PG   <- matrix(rep(fixef(la05)["Gewinn:groupPG>HC"], as.numeric(table(group)["PG"])))
  gain.AD   <- matrix(rep(fixef(la05)["Gewinn:groupAD>HC"], as.numeric(table(group)["AD"])))
  gain.group<- rbind(gain.HC,gain.PG,gain.AD)
  gain      <- gain + gain.group
  
  # get group specific loss
  loss.HC   <- matrix(rep(0,as.numeric(table(group)["HC"])))
  loss.PG   <- matrix(rep(fixef(la05)["Verlust:groupPG>HC"], as.numeric(table(group)["PG"])))
  loss.AD   <- matrix(rep(fixef(la05)["Verlust:groupAD>HC"], as.numeric(table(group)["AD"])))
  loss.group<- rbind(loss.HC,loss.PG,loss.AD)
  loss      <- loss + loss.group
  
  # fixef as a whole
  fixef.la  <- cbind(intercept, gain, loss)
  # add ranef
  tmp <- ranef(la05)$subject
  test <- data.frame(tmp[,1],tmp$Gewinn,tmp$Verlust)
  test      <- fixef.la + test
  
  crm       <- data.frame(id=rownames(ranef(la05)$subject),
                          group =group,
                          intercept=test["tmp...1."],
                          gain=test["tmp.Gewinn"],
                          loss=test["tmp.Verlust"])
  crm$Lambda<- crm$tmp.Verlust*(-1)/crm$tmp.Gewinn
  
  names(crm) <- c("id", "group", "intercept", "Gewinn", "Verlust", "Lambda")
  
  if(do_lmlist == 1){
    crm_lmlist <- c()
    for (ii in 1:length(lmlist_ed)) {
      crm_lmlist <- rbind(crm_lmlist,t(as.matrix(as.numeric(lmlist_ed[[ii]][[1]]))))
    }
    crm_lmlist <- as.data.frame(crm_lmlist)
    names(crm_lmlist) <- c("intercept", "ed.abs", "Gewinn", "Verlust")
    crm_lmlist$id       <- as.factor(as.numeric(as.matrix(names(lmlist_ed))))
    crm <- crm_lmlist
  }
  
  return(crm)
}

# function for extracting fixefs of interest
get_la_fixef <- function(.) {
  
  # intercepts
  an_AD <- as.numeric((fixef(.)["(Intercept)"] + fixef(.)["groupAD>HC"]))
  an_PG <- as.numeric((fixef(.)["(Intercept)"] + fixef(.)["groupPG>HC"]))
  an_HC <- as.numeric(fixef(.)["(Intercept)"])
  
  la_AD <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupAD>HC"])*(-1)/(fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupAD>HC"]))
  la_PG <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupPG>HC"])*(-1)/(fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupPG>HC"]))
  la_HC <- as.numeric((fixef(.)["Verlust"])*(-1)/(fixef(.)["Gewinn"]))
  
  bg_AD <- as.numeric((fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupAD>HC"]))
  bg_PG <- as.numeric((fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupPG>HC"]))
  bg_HC <- as.numeric(fixef(.)["Gewinn"])
  
  bl_AD <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupAD>HC"]))
  bl_PG <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupPG>HC"])) 
  bl_HC <- as.numeric(fixef(.)["Verlust"])
  
  # contrasts
  x_la_HCgrAD <- la_HC - la_AD
  x_la_HCgrPG <- la_HC - la_PG
  x_la_PGgrAD <- la_PG - la_AD
  
  # return the coefficients
  tmp   <- ls()
  out_v <- c() 
  for (ii in 1:length(tmp)) {
    out_v[ii] <- eval(parse(text=tmp[ii]))
  }
  names(out_v) <- tmp
  return(out_v)
}

# function for extracting fixefs of interest
get_la_fixef_pdt <- function(.) {
  
  # intercepts
  an_PG <- as.numeric((fixef(.)["(Intercept)"] + fixef(.)["HCPGPG"]))
  an_HC <- as.numeric(fixef(.)["(Intercept)"])
  
  la_PG <- as.numeric((fixef(.)["loss"] + fixef(.)["loss:HCPGPG"])*(-1)/(fixef(.)["gain"] + fixef(.)["gain:HCPGPG"]))
  la_HC <- as.numeric((fixef(.)["loss"])*(-1)/(fixef(.)["gain"]))
  
  bg_PG <- as.numeric((fixef(.)["gain"] + fixef(.)["gain:HCPGPG"]))
  bg_HC <- as.numeric(fixef(.)["gain"])
  
  bl_PG <- as.numeric((fixef(.)["loss"] + fixef(.)["loss:HCPGPG"])) 
  bl_HC <- as.numeric(fixef(.)["loss"])
  
  # contrasts
  x_la_HCgrPG <- la_HC - la_PG
  
  # return the coefficients
  tmp   <- ls()
  out_v <- c() 
  for (ii in 1:length(tmp)) {
    out_v[ii] <- eval(parse(text=tmp[ii]))
  }
  names(out_v) <- tmp
  return(out_v)
}



# function to plot a boot object and save the plot (with CI)
agk.barplot.boot <- function(boot_pdt_df, cur_n, cur_cat, nsim,study_name) {
  
  temp        <- 1:length(boot_pdt_df[,1])
  boot_pdt_df <- data.frame(boot_pdt_df,temp)
  boot_pdt_df <- melt(boot_pdt_df,id.vars = c("temp"))
  boot_pdt_df$variable <- factor(boot_pdt_df$variable,ordered = T)
  
  # now make a grouping variable
  # first discard all the contrasts; won't plot them here
  # also intercepts
  to_discard  <- grep("gr",as.character(boot_pdt_df$variable))
  boot_pdt_df <- boot_pdt_df[-to_discard,]
  to_discard  <- grep("in_",as.character(boot_pdt_df$variable))
  boot_pdt_df <- boot_pdt_df[-to_discard,]
  # make the grouping var
  ADgroup <- grep("AD",as.character(boot_pdt_df$variable))
  PGgroup <- grep("PG",as.character(boot_pdt_df$variable))
  HCgroup <- grep("HC",as.character(boot_pdt_df$variable))
  grouping_var <- c()
  for (kk in 1:length(boot_pdt_df$variable)) {
    if (any(ADgroup==kk)) {
      grouping_var[kk] <- "AD" 
    } else if (any(PGgroup==kk)) {
      grouping_var[kk] <- "PG" 
    } else if (any(HCgroup == kk)) {
      grouping_var[kk] = "HC"
    }
  }
  
  boot_pdt_df <- data.frame(boot_pdt_df,as.factor(grouping_var))
  names(boot_pdt_df)[4] <- "group"
  boot_pdt_df$group = factor(boot_pdt_df$group, levels=c("HC","PG","AD"), labels = c("HC","PG","AD"))
  
  # rename the variables
  boot_pdt_df$variable <- as.character(boot_pdt_df$variable)
  for (kk in 1:length(boot_pdt_df[,1])) {
    boot_pdt_df$variable[kk] <- strsplit(as.character(boot_pdt_df$variable[kk]),"_")[[1]][1]
  }
  
  boot_pdt_df <- as.data.frame(acast(boot_pdt_df, group + temp  ~ variable, value.var = "value"))
  boot_pdt_df$group <- row.names(boot_pdt_df)
  
  boot_pdt_df$group <- as.character(boot_pdt_df$group)
  for (kk in 1:length(boot_pdt_df$group)) {
    boot_pdt_df$group[kk] <- strsplit(as.character(boot_pdt_df$group[kk]),"_")[[1]][1]
    
  }
  
  # get the mean and errorbar values
  # reorder
  #boot_pdt_df <- boot_pdt_df[,c(3,1,2,4,5)]
  plot.dat <- data.frame()
  var.names <- c()
  for (ii in 1:(length(boot_pdt_df[1,])-1)){
    # bootstrapped 95% CI
    bt.df     <- as.data.frame(as.list(aggregate(boot_pdt_df[ii],by = list(group = boot_pdt_df$group), agk.mean.quantile,lower=0.025,upper=0.975)))
    names(bt.df) <- c("group", "mean","lower", "upper")
    var.names <- c(var.names,rep(names(boot_pdt_df[ii]),length(levels(bt.df$group))))
    plot.dat <- rbind(plot.dat,bt.df)
  }
  plot.dat$variable <- var.names
  plot.dat$variable <- as.factor(plot.dat$variable)
  #plot.dat$variable <- factor(plot.dat$variable,levels(plot.dat$variable)[c(3,1,2,4)])
  #plot.dat$group    <- factor(plot.dat$group,levels(plot.dat$group)[c(2,3,1)])
  plot.dat$group = factor(plot.dat$group, levels=c("HC","PG","AD"), labels = c("HC","PG","AD"))
  
  # plot
  p <- ggplot(data = plot.dat, aes(group,mean,fill=variable))
  p <- p+geom_bar(position="dodge",stat="identity")
  p
  p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                         position=position_dodge(width=0.9), color=cbbPalette[4],
                         width=0.1) + ylab("mean (95% boots. CI)\n")
  
  p <- p + ggtitle("Fixed effects \n")
  
  return(p)
}

# bootstrap confint of a parameter, given a vector x
agk.boot.ci = function(x,cur_fun,lower,upper,R) {
  if (!is.numeric(x)) {return(c(NA,NA,NA))}
  if (all(is.na(x)))  {return(c(NA,NA,NA))}
  tmp = one.boot(x,FUN=cur_fun,R)
  return(as.numeric(agk.mean.quantile(tmp$t,lower,upper)))
}

agk.boot.ci.c = compiler::cmpfun(agk.boot.ci)

# function that computes mean with disregarding NA's
agk.mean.narm <- function(x) {
  return(mean(x,na.rm =T ))
}

# function to do the correlations in crm given the desired variables in data.la
# does it for one group in data.la
# given in cur_group
# thresh is the p-threshold used for the fdr correction
# returns a list with matrices of always two variables which show sig. corr.
agk.la.cor <- function(desired_vars, data.la, crm,fmri_dat, thresh, cur_group) {
  data_la_covs <- subset(data.la[desired_vars], group == cur_group)
  sub_grp      <- data_la_covs[,c(1,2)]
  data_la_covs <- data_la_covs[,-c(1,2)]
  covs <- aggregate(data_la_covs,list(sub_grp$subject),agk.mean.narm)
  covs <- covs[,-1]
  
  # kick out vars which have missings
  f= function(.) {!any(is.na(.))}
  covs <- covs[,as.logical(sapply(covs,FUN = f))]
  covs <- scale(covs)
  
  # use first PCA of severity
  cur_pca <- prcomp(covs)
  covs    <- as.data.frame(cur_pca$x[,1])
  names(covs) <- "sev_PC1"
  
  crmc <- subset(crm,group==cur_group)
  crmc <- crmc[,-c(1,2)]
  
  # get the fMRI data attached
  fmri_dat_c <- subset(fmri_dat,crm$group == cur_group)
  crmc <- cbind(crmc,fmri_dat_c)
  covs <- cbind(covs)
  
  cur_cor <- corr.test(crmc, covs, adjust="fdr", method = "spearman")
  cur_cor_r <- cur_cor[1]$r
  cur_cor_p <- cur_cor[4]$p
  
  # extract the sig correlations and prep for plotting
  plot_cors <- list()
  count <- 0
  for (kk in 1:length(cur_cor_r[,1])) {
    for (ll in 1:length(cur_cor_r[1,])) {
      if (cur_cor_p[kk,ll] < thresh) {
        count <- count + 1
        v1 <- crmc[row.names(cur_cor_p)[kk]]
        v2 <- covs[names(covs)[ll]]
        plt_out <- data.frame(v1,v2)
        names(plt_out) <- c(row.names(cur_cor_p)[kk],names(covs)[ll])
        plot_cors[[count]] <- plt_out
      }
    }
  }
  return(plot_cors)
}

# extracts of a cor object the sig corelltions and returns th respective bivariate dataframe for plotting
agk.extract.sig.corr <- function(cur_cor,crmc,covs,thresh=0.05) {
  cur_cor_r <- cur_cor[1]$r
  cur_cor_p <- cur_cor[4]$p
  
  # extract the sig correlations and prep for plotting
  plot_cors <- list()
  count <- 0
  for (kk in 1:length(cur_cor_r[,1])) {
    for (ll in 1:length(cur_cor_r[1,])) {
      if (cur_cor_p[kk,ll] < thresh) {
        count <- count + 1
        v1 <- crmc[row.names(cur_cor_p)[kk]]
        v2 <- covs[names(covs)[ll]]
        plt_out <- data.frame(v1,v2)
        names(plt_out) <- c(row.names(cur_cor_p)[kk],names(covs)[ll])
        plot_cors[[count]] <- plt_out
      }
    }
  }
  return(plot_cors)
}

# agk.plot.cor <- function(cur_group,dat_list) {
#   compl_df <- data.frame()
#   for (ii in 1:length(dat:list)) {
#     cur_df <- dat_list[[ii]]
#     compl_df <- cbind(compl_df, cur_df)
#   }
#   compl_df <- scale(compl_df)
#   # add sub id
#   sub <- 1:length(compl_df[,1])
#   compl_df <- data.frame(sub,compl_df)
#   compl_df <- melt(compl_df,id.vars=list("sub"))
#   
#   p1 <- ggplot(compl_df, aes(x=))
# }

agk.plot.cor <- function(cur_group,dat_list) {
  plot_list <- list()
  for (ii in 1:length(dat_list)) {
    cur_df <- dat_list[[ii]]
    plot_list[[ii]] <- ggplot(cur_df,aes_string(x=names(cur_df)[1],y=names(cur_df)[2])) + geom_point()
  }
  return(plot_list)
}

# function to control fmri for gmd and age according to eigv name
agk.control.fmri <- function(fmri_dat, gmd_dat,nms_anly,age) {
  sv_names <- names(fmri_dat)
  v_init <- rep(0,length(fmri_dat[,1]))
  dat_out <- data.frame(v_init)
  for (ii in 1:length(fmri_dat[1,])) {
    cur_fmri  <- as.numeric(fmri_dat[[ii]])
    cur_gmd   <- gmd_dat[ii]
    cur_name  <- as.character(nms_anly$V1[1])
    
    if (as.logical(length((grep(c("bpm"),cur_name)) !=0))) {
      if (as.logical(length((grep(c("Age"),cur_name)) !=0))) {
        cur_v <- resid(lm(cur_fmri ~ cur_gmd + age))
      } else {cur_v <- resid(lm(cur_fmri ~ cur_gmd))}
    } else {
      if (as.logical(length((grep(c("Age"),cur_name)) !=0))) {
        cur_v <- resid(lm(cur_fmri ~ age))
      } else {cur_v <- cur_fmri}
    }
    dat_out <- data.frame(dat_out,cur_v)  
  }
  dat_out <- dat_out[,-1]
  dat_out <- as.data.frame(dat_out)
  names(dat_out) <- sv_names
  return(dat_out)
}

# function to plot bar plot of fmri_vector eigenvariate
agk.fmri.mean.plot <- function(fmri_con,group,nms_anly,ind) {
  
  fmri_v <- fmri_con[ind]
  cur_df <- data.frame(fmri_v,group)
  p <- ggplot(cur_df, aes_string(x=names(cur_df)[2],y=names(cur_df)[1]))
  
  p <- p + stat_summary(fun.y = mean, geom = "bar", width=0.2, fill="#333333")
  p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
  p <- p + stat_summary(fun.data = mean_cl_normal, geom = "linerange",
                        color="#FF3333", size =4.5)
  p <- p + ylab(names(cur_df)[1])
  return(p)
}

##   COMPILING ================================================================
join = cmpfun(join)

# # checking for unusual behavior towards gain and loss
# 
# check_gain_loss <- function (data.la, break_var) {
#   prep_ggplot_la()
#   weird_subs <- c()
#   subs <- levels(eval(parse(text=paste("data.la$",break_var,sep = ""))))
#   for (ii in 1:length(subs)) {
#     print(paste("now running ", subs[ii], sep=""))
#     cur.data <- subset(data.la,subset = eval(parse(text=paste("data.la$",break_var,sep = ""))) == subs[ii])
#     cur.data <- subset(cur.data, !is.na(accept.reject))
#     
#     # first test
#     cur.res  <- describeBy(cur.data$Gewinn, cur.data$accept.reject)
#     if (cur.res$'0'$mean > cur.res$'1'$mean) {
#       weird_subs[ii] <- subs[ii]
#       p <- ggplot(cur.data, aes(x=accept.reject, y=Gewinn)) + stat_summary(fun.y = mean, geom = "bar")
#       p
#       next
#     }
#     cur.res  <- describeBy(cur.data$Verlust, cur.data$accept.reject)
#     if (cur.res$'0'$mean < cur.res$'1'$mean) {
#       weird_subs[ii] <- subs[ii]
#       p <- ggplot(cur.data, aes(x=accept.reject, y=Verlust)) + stat_summary(fun.y = mean, geom = "bar")
#       p
#       next
#     }
#     
#     # other test
#     cur.res <- 
#   }
# }



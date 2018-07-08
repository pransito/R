## PREAMBLE ===================================================================
# replicating the LA paper behavioral analyses using the lae model and
# lae: loss aversion model with euclidean distance
# bootstrapping p-values
# bootstrap the p-values?!
# function to compute a p-value for a fixed effect of interest
# in a glmer model
# mod_bl ist das baseline model gefittet
# mod ist ein string der das fitting des complete models angbit
# fun_extract ist die funktion, die fiexedeffect und contrasts of interests extrahiert
# mod_compl_fit ist das komplette model fit on the original data

## PARAMETERS =================================================================
num           = 1000 # bootstrapping reps
mod_preps     = 0    # estimate the models?
with_ed       = 0    # lae or la?
boots_p_nocov = 1    # bootstrap p-values no covariate?
boots_p_cov   = 1    # bootstrap p-values with covariate?
boots_ci      = 1    # bootstrap ci's ?
eval_res      = 1    # evaluate and plot results?
with_cov      = 0    # only for evaluation: with or without covariate?

## FUNCTIONS ==================================================================
# fixef extraction function
fun_extract = function(.) {
  
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

agk.boot.p.lmer.subfun = function(mod_bl,mod,fun_extract,mod_compl_fit) {
  data_pdt = mod_compl_fit@frame
  # simulate under 0 hyothesis
  data_pdt$accept_reject <- (simulate(mod_bl))[[1]]
  # fit under alternative hypothesis
  cur_fit <- eval(parse(text=mod))
  return(fun_extract(cur_fit))
}

# function to plot a boot object and save the plot (with CI)
agk.barplot.boot <- function(boot_pdt_la, cur_n, cur_cat, nsim,study_name) {
  
  temp        <- 1:length(boot_pdt_la$t[,1])
  boot_pdt_df <- data.frame(boot_pdt_la$t,temp)
  boot_pdt_df <- melt(boot_pdt_df,id.vars = c("temp"))
  boot_pdt_df$variable <- factor(boot_pdt_df$variable,ordered = T)
  
  # now make a grouping variable
  # first discard all the contrasts; won't plot them here
  # also intercepts
  to_discard  <- c(grep("gr",as.character(boot_pdt_df$variable)),grep("in_",as.character(boot_pdt_df$variable)))
  if (!isempty(to_discard)) {
    boot_pdt_df <- boot_pdt_df[-to_discard,]
  }
  #to_discard  <- grep("AD",as.character(boot_pdt_df$variable))
  #boot_pdt_df <- boot_pdt_df[-to_discard,]
  # make the grouping var
  PGgroup <- grep("PG",as.character(boot_pdt_df$variable))
  HCgroup <- grep("HC",as.character(boot_pdt_df$variable))
  grouping_var <- c()
  for (kk in 1:length(boot_pdt_df$variable)) {
    if (any(PGgroup==kk)) {
      grouping_var[kk] <- "PG" 
    } else if (any(HCgroup == kk)) {
      grouping_var[kk] = "HC"
    }
  }
  
  boot_pdt_df <- data.frame(boot_pdt_df,as.factor(grouping_var))
  names(boot_pdt_df)[4] <- "group"
  boot_pdt_df$group = factor(boot_pdt_df$group, levels=c("HC","PG"), labels = c("HC","PG"))
  
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
    names(bt.df) <- c("group", "mean", "lower","upper")
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
                         position=position_dodge(width=0.9), color=cbbPalette[5],
                         width=0.1) + ylab("mean (95% boots. CI)\n")
  
  p <- p + ggtitle("Fixed effects \n")
  
  return(p)
}

cbbPalette   = grey.colors(4, start = 0.3, end = 0.8, gamma = 2.2, alpha = NULL)
cbbPalette   = c(cbbPalette,"#000000")

## PREPARATIONS ===============================================================
# path
#path_repl_LA = 'C:/Users/genaucka/Google Drive/Library/R/PDT/analyses/in_progress/replicate_LA'
path_repl_LA = 'E:/Google Drive/Library/R/PDT/analyses/in_progress/replicate_LA'

# taking out smoking vars
data_pdt = data_pdt[-grep('smoking',names(data_pdt))]

if (mod_preps) {
  # replicating reduced LA?
  cur_control = glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=400))
  if (with_ed) {
    mod_lae     = glmer(accept_reject ~ gain + loss + ed_abs + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_laeg    = glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_laeg_s   = 'glmer(accept_reject ~ (gain + loss + ed_abs)*HCPG + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)'
  } else {
    mod_la      = glmer(accept_reject ~ gain + loss + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_lag     = glmer(accept_reject ~ (gain + loss)*HCPG + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_lag_s   = 'glmer(accept_reject ~ (gain + loss)*HCPG + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)'
  }
}

if (with_ed) {
  mod_compl_fit = mod_laeg
  mod_bl        = mod_lae
  mod           = mod_laeg_s
} else {
  mod_compl_fit = mod_lag
  mod_bl        = mod_la
  mod           = mod_lag_s
}

anova(mod_bl,mod_compl_fit)

# la_coefs = agk.get.compl.coef(mod_compl_fit,'HCPG')
# la_coefs$lambda = -la_coefs$loss/la_coefs$gain
# plot(la_coefs$lambda ~ la_coefs$HCPG)
# lagrp = lmp(lambda ~ HCPG, data = la_coefs,Ca=0.000001, maxIter = 1000000,nCycle=1)
# summary(lagrp)


## RUNNING THE BOOTSTRAP (NO COV) =============================================
if (boots_p_nocov) {
  cl<-makeCluster(6)  
  registerDoSNOW(cl)
  effects_under_0 <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data_pdt",'cur_control')) %dopar% {
    agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
  } 
  stopCluster(cl)
  
  setwd(path_repl_LA)
  if (with_ed) {
    save(file="effects_under_0_PDT_LA.RData",list=c("effects_under_0"))
  } else {
    save(file="effects_under_0_PDT_LA_noed.RData",list=c("effects_under_0"))
  }
}

## PREPARATIONS (WITH COV) ====================================================
smoking_df = dat_match[c('VPPG','smoking_ftdt')]
data_pdt = merge(data_pdt,smoking_df,by.x = 'subject',by.y='VPPG')
if (mod_preps) {
  if (with_ed) {
    mod_lae_cv     = glmer(accept_reject ~ (gain + loss + ed_abs)*smoking_ftdt + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_laeg_cv    = glmer(accept_reject ~ (gain + loss + ed_abs)*smoking_ftdt*HCPG + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_laeg_cv_s  = 'glmer(accept_reject ~ (gain + loss + ed_abs)*smoking_ftdt*HCPG + (gain + loss + ed_abs|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)'
  } else {
    mod_la_cv      = glmer(accept_reject ~ (gain + loss)*smoking_ftdt + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_lag_cv     = glmer(accept_reject ~ (gain + loss)*smoking_ftdt*HCPG + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)
    mod_lag_cv_s   = 'glmer(accept_reject ~ (gain + loss)*smoking_ftdt*HCPG + (gain + loss|subject), data = data_pdt, family = binomial,nAGQ = 0, control = cur_control)'
  }
}

anova(mod_la,mod_lag)

if (with_ed) {
  mod_compl_fit = mod_laeg_cv
  mod_bl        = mod_lae_cv
  mod           = mod_laeg_cv_s
} else {
  mod_compl_fit = mod_lag_cv
  mod_bl        = mod_la_cv
  mod           = mod_lag_cv_s
}


## RUNNING THE BOOTSTRAP (COV) ================================================
if (boots_p_cov) {
  cl<-makeCluster(6)  
  registerDoSNOW(cl)
  effects_under_0 <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data_pdt",'cur_control')) %dopar% {
    agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
  } 
  stopCluster(cl)
  
  setwd(path_repl_LA)
  if (with_ed) {
    save(file="effects_under_0_smoking_PDT_LA.RData",list=c("effects_under_0"))
  } else {
    save(file="effects_under_0_smoking_PDT_LA_noed.RData",list=c("effects_under_0"))
  }
}

## BOOTSRAP CI's ==============================================================
if (boots_ci) {
  if (with_ed) {
    mod_compl_fit = mod_laeg
  } else {
    mod_compl_fit = mod_lag
  }
  # always without cov
  boot_pdt_la = bootMer(mod_compl_fit,use.u = F,
                        FUN = fun_extract,
                        nsim = num,
                        type = "parametric",
                        verbose = TRUE,
                        parallel = "snow",
                        ncpus = 7)
  
  setwd(path_repl_LA)
  if (with_ed) {
    save(file="boot_pdt_la_compl_no_cov.RData", list=c("boot_pdt_la"))
  } else {
    save(file="boot_pdt_la_compl_noed_no_cov.RData", list=c("boot_pdt_la"))
  }
}

## EVALUATION =================================================================
if (eval_res) {
  setwd(path_repl_LA)
  if (with_ed & with_cov == 0) {
    disp('evaluating with ed and no cov')
    load("effects_under_0_PDT_LA.RData")
    mod_compl_fit = mod_laeg
  } else if (with_ed & with_cov == 1) {
    disp('evaluating with ed and with cov')
    load("effects_under_0_smoking_PDT_LA.RData")
    mod_compl_fit = mod_laeg_cv
  } else if (with_ed == 0 & with_cov == 1) {
    disp('evaluating without ed and with cov')
    load("effects_under_0_smoking_PDT_LA_noed.RData")
    mod_compl_fit = mod_lag_cv
  } else if (with_ed == 0 & with_cov == 0) {
    disp('evaluating without ed and with cov')
    load("effects_under_0_PDT_LA_noed.RData")
    mod_compl_fit = mod_lag
  }
  
  cur_vars    <- c("an_HC","an_PG","bg_HC","bg_PG","bl_HC", "bl_PG", "x_la_HCgrPG")
  for (ii in 1:length(cur_vars)) {
    cur_fun     <- function(l) {return(l[cur_vars[ii]])}
    eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0,FUN = cur_fun))",sep="")))
  }
  
  con_bl_PG_HC <- bl_PG-bl_HC 
  con_bg_PG_HC <- bg_PG-bg_HC
  con_an_PG_HC <- an_PG-an_HC
  
  ## p-values for HC gr PG
  obs_fixef <- fun_extract(mod_compl_fit)
  disp('LA HC gr PG')
  disp(1-agk.density_p(x_la_HCgrPG,obs_fixef[9]))
  
  # p-value for PG gr HC: beta gain (one-sided)
  disp('bg PG gr HC')
  disp(agk.density_p(con_bg_PG_HC,obs_fixef[4]-obs_fixef[3]))
  
  # p-value for PG gr HC: beta loss (one-sided)
  disp('bl PG gr HC')
  disp(1-agk.density_p(con_bl_PG_HC,obs_fixef[6]-obs_fixef[5]))
  
  # p-value for PG gr HC: beta loss (one-sided)
  disp('intercept PG gr HC')
  disp(1-agk.density_p(con_an_PG_HC,obs_fixef[2]-obs_fixef[1]))
}

## PLOT RESULTS ===============================================================
if (eval_res) {
  if (with_ed) {
    load('boot_pdt_la_compl_no_cov.RData')
  } else {
    load('boot_pdt_la_compl_noed_no_cov.RData')
  }
  p1 = agk.barplot.boot(boot_pdt_la,cur_n = 1000,nsim = 1,study_name = 'LA')
  p1 <- p1 + scale_fill_manual(labels=c(expression(beta[int],beta[gain]), 
                                        expression(beta[loss]), 
                                        expression(lambda)),
                               values=cbbPalette,
                               name="fixed effect")
  p1 <- p1 + ggtitle("Differences in loss aversion by group\n")
  print(p1 + theme_la())
}

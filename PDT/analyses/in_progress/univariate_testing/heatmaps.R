## heatmap
cur_med = median(data_pdt$v_imageRating1s)
data_pdt$craving = as.factor(ifelse(data_pdt$v_imageRating1s <=cur_med,"low","high"))

# get a group specific heatmap by agregating over subs
hmdf <- aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$gain_bcp,data_pdt$loss_bcp),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","gain","loss","PoA")

tmp <- aggregate(hmdf$PoA,by = list(hmdf$HCPG,hmdf$gain,hmdf$loss),FUN=mean)
names(tmp) <- c("HCPG","gain", "loss","PoA")
tmp$loss <- as.factor(tmp$loss)
tmp$gain <- as.factor(tmp$gain)
tmp$HCPG = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 

# heatmap
p_hm <- qplot(x=gain, y=loss, data=tmp, fill=PoA, geom="tile")
#p_hm <- p_hm + scale_fill_gradientn(colours=rainbow(4))
p_hm <- p_hm + scale_fill_gradientn(colours=c(Palette_heat[1], Palette_heat[3]))
p_hm <- p_hm + facet_grid(.~HCPG)
p_hm <- p_hm + scale_x_discrete(expand=c(0,0),breaks=c(14,36))+ scale_y_discrete(expand=c(0,0),breaks=c(-18,-7))
p_hm <- p_hm + labs(x = "possible gain (euros)", y = "possible loss (euros)")
p_hm  = p_hm + theme_la()
p_neu = (p_hm + theme_update(strip.text.x = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.text.y = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.background = element_rect(fill="white", color = "white"),
                             strip.switch.pad.grid = unit(c(3),"cm"),
                             strip.switch.pad.wrap = unit(c(3),"cm"),
                             panel.margin.x =unit(0.5,"cm"),
                             panel.margin.y =unit(1,"cm"),
                             plot.title = element_text(colour = "black",size=25),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.ticks = element_line(colour ="black", size=1.5),
                             axis.title.x = element_text(colour = "black",size=25,face="bold",vjust=0.2),
                             axis.title.y = element_text(colour = "black",size=25,face="bold",vjust=2.0,angle = 90),
                             plot.margin = unit(c(1,0.5,1,1),"cm"),
                             legend.background = element_rect(fill="white", color = "white"),
                             legend.text = element_text(colour = "black",size=20,face="bold"),
                             legend.title = element_text(colour = "black",size=20,face="bold"),
                             #legend.margin=unit(.85,"cm"),
                             legend.text.align=0,
                             legend.margin=unit(0.2,"cm")))
p_neu + guides(fill = guide_legend(title.theme = element_text(size=25, face="bold", colour = "black", angle = 0)))


# get a group and craving specific heatmap by agregating over subs
hmdf <- aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$craving,data_pdt$gain_bcp,data_pdt$loss_bcp),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","craving","gain","loss","PoA")

tmp <- aggregate(hmdf$PoA,by = list(hmdf$HCPG,hmdf$craving,hmdf$gain,hmdf$loss),FUN=mean)
names(tmp) <- c("HCPG","craving","gain", "loss","PoA")
tmp$loss <- as.factor(tmp$loss)
tmp$gain <- as.factor(tmp$gain)
tmp$HCPG = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 

# heatmap
p_hm <- qplot(x=gain, y=loss, data=tmp, fill=PoA, geom="tile")
#p_hm <- p_hm + scale_fill_gradientn(colours=rainbow(4))
p_hm <- p_hm + scale_fill_gradientn(colours=c(Palette_heat[1], Palette_heat[3]))
p_hm <- p_hm + facet_grid(HCPG~craving)
p_hm <- p_hm + scale_x_discrete(expand=c(0,0),breaks=c(14,36))+ scale_y_discrete(expand=c(0,0),breaks=c(-18,-7))
p_hm <- p_hm + labs(x = "possible gain (euros)", y = "possible loss (euros)")
p_hm  = p_hm + theme_la()
p_neu = (p_hm + theme_update(strip.text.x = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.text.y = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.background = element_rect(fill="white", color = "white"),
                             strip.switch.pad.grid = unit(c(3),"cm"),
                             strip.switch.pad.wrap = unit(c(3),"cm"),
                             panel.margin.x =unit(0.5,"cm"),
                             panel.margin.y =unit(1,"cm"),
                             plot.title = element_text(colour = "black",size=25),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.ticks = element_line(colour ="black", size=1.5),
                             axis.title.x = element_text(colour = "black",size=25,face="bold",vjust=0.2),
                             axis.title.y = element_text(colour = "black",size=25,face="bold",vjust=2.0,angle = 90),
                             plot.margin = unit(c(1,0.5,1,1),"cm"),
                             legend.background = element_rect(fill="white", color = "white"),
                             legend.text = element_text(colour = "black",size=20,face="bold"),
                             legend.title = element_text(colour = "black",size=20,face="bold"),
                             #legend.margin=unit(.85,"cm"),
                             legend.text.align=0,
                             legend.margin=unit(0.2,"cm")))
p_neu + guides(fill = guide_legend(title.theme = element_text(size=25, face="bold", colour = "black", angle = 0)))

# barplot acceptance rate by craving
hmdf <- aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$craving,data_pdt$gain_bcp,data_pdt$loss_bcp),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","craving","gain","loss","PoA")

tmp <- aggregate(hmdf$PoA,by = list(hmdf$subject, hmdf$HCPG,hmdf$craving),FUN=mean)
names(tmp) <- c("subject","HCPG","craving","PoA")
tmp$HCPG = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 
wide = reshape(tmp,direction = "wide",idvar=c("subject","HCPG"),timevar = c("craving"))
wide$'PoA.high-low' = wide$PoA.high-wide$PoA.low
wide = wide[-1]
wide = wide[-c(2,3)]
boot_pdt_df = wide

# get the mean and errorbar values
one.boot.mean = function(x) {
  tmp = one.boot(x,FUN=mean,R=2000)
  return(agk.mean.quantile(tmp$t))
  }
plot.dat <- data.frame()
var.names <- c()
for (ii in 2:length(boot_pdt_df[1,])){
  # bootstrapped 95% CI
  bt.df     <- as.data.frame(as.list(aggregate(boot_pdt_df[ii],by = list(group = boot_pdt_df$HCPG), FUN=one.boot.mean)))
  names(bt.df) <- c("group", "lower", "upper","mean")
  var.names <- c(var.names,rep(names(boot_pdt_df[ii]),length(levels(bt.df$group))))
  plot.dat <- rbind(plot.dat,bt.df)
}
plot.dat$variable <- var.names
plot.dat$variable <- as.factor(plot.dat$variable)
plot.dat$group = factor(plot.dat$group, levels=c("HC","PG"), labels = c("HC","PG"))

# plot
p <- ggplot(data = plot.dat, aes(group,mean,fill=variable))
p <- p+geom_bar(position="dodge",stat="identity",fill="#333333")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                       position=position_dodge(width=0.9),
                       color="#FF3333", size =2.5) + ylab("mean shift in acc. rate \n(95% boots. CI)\n")

p <- p + ggtitle("Shift in acceptance rate \n depending on craving value of pictures\n")
p

# get a group and category specific heatmap by agregating over subs
hmdf <- aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$cat,data_pdt$gain_bcp,data_pdt$loss_bcp),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","craving","gain","loss","PoA")

tmp <- aggregate(hmdf$PoA,by = list(hmdf$HCPG,hmdf$craving,hmdf$gain,hmdf$loss),FUN=mean)
names(tmp) <- c("HCPG","craving","gain", "loss","PoA")
tmp$loss <- as.factor(tmp$loss)
tmp$gain <- as.factor(tmp$gain)
tmp$HCPG = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 

# heatmap
p_hm <- qplot(x=gain, y=loss, data=tmp, fill=PoA, geom="tile")
#p_hm <- p_hm + scale_fill_gradientn(colours=rainbow(4))
p_hm <- p_hm + scale_fill_gradientn(colours=c(Palette_heat[1], Palette_heat[3]))
p_hm <- p_hm + facet_grid(HCPG~craving)
p_hm <- p_hm + scale_x_discrete(expand=c(0,0),breaks=c(14,36))+ scale_y_discrete(expand=c(0,0),breaks=c(-18,-7))
p_hm <- p_hm + labs(x = "possible gain (euros)", y = "possible loss (euros)")
p_hm  = p_hm + theme_la()
p_neu = (p_hm + theme_update(strip.text.x = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.text.y = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                             strip.background = element_rect(fill="white", color = "white"),
                             strip.switch.pad.grid = unit(c(3),"cm"),
                             strip.switch.pad.wrap = unit(c(3),"cm"),
                             panel.margin.x =unit(0.5,"cm"),
                             panel.margin.y =unit(1,"cm"),
                             plot.title = element_text(colour = "black",size=25),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.text = element_text(colour = "black",size=25,face="bold"),
                             axis.ticks = element_line(colour ="black", size=1.5),
                             axis.title.x = element_text(colour = "black",size=25,face="bold",vjust=0.2),
                             axis.title.y = element_text(colour = "black",size=25,face="bold",vjust=2.0,angle = 90),
                             plot.margin = unit(c(1,0.5,1,1),"cm"),
                             legend.background = element_rect(fill="white", color = "white"),
                             legend.text = element_text(colour = "black",size=20,face="bold"),
                             legend.title = element_text(colour = "black",size=20,face="bold"),
                             #legend.margin=unit(.85,"cm"),
                             legend.text.align=0,
                             legend.margin=unit(0.2,"cm")))
p_neu + guides(fill = guide_legend(title.theme = element_text(size=25, face="bold", colour = "black", angle = 0)))

# barplot acceptance rate by cat
hmdf <- aggregate(as.numeric(as.character(data_pdt$accept_reject)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$cat,data_pdt$gain_bcp,data_pdt$loss_bcp),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","category","gain","loss","PoA")

tmp <- aggregate(hmdf$PoA,by = list(hmdf$subject, hmdf$HCPG,hmdf$category),FUN=mean)
names(tmp) <- c("subject","HCPG","category","PoA")
tmp$HCPG = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 
wide = reshape(tmp,direction = "wide",idvar=c("subject","HCPG"),timevar = c("category"))
wide$'PoA.high-low' = wide$PoA.gambling-wide$PoA.neutral
wide = wide[-1]
wide = wide[-c(2,3,4,5)]
boot_pdt_df = wide

boot_pdt_df_behav$study = "Behav"
boot_pdt_df_MRI$study = "MRI"

boot_pdt_df = rbind(boot_pdt_df_behav,boot_pdt_df_MRI)
boot_pdt_df = boot_pdt_df[c(1,3,2)]

# get the mean and errorbar values
one.boot.mean = function(x) {
  tmp = one.boot(x,FUN=mean,R=2000)
  return(agk.mean.quantile(tmp$t))
}
plot.dat <- data.frame()
var.names <- c()
for (ii in 3:length(boot_pdt_df[1,])){
  # bootstrapped 95% CI
  bt.df     <- as.data.frame(as.list(aggregate(boot_pdt_df[ii],by = list(group = boot_pdt_df$HCPG, study = boot_pdt_df$study), FUN=one.boot.mean)))
  names(bt.df) <- c("group","study","lower", "upper","mean")
  var.names <- c(var.names,rep(names(boot_pdt_df[ii]),length(levels(bt.df$group))))
  plot.dat <- rbind(plot.dat,bt.df)
}
plot.dat$variable <- var.names
plot.dat$variable <- as.factor(plot.dat$variable)
plot.dat$group = factor(plot.dat$group, levels=c("HC","PG"), labels = c("HC","PG"))
plot.dat$study = as.factor(plot.dat$study)


# plot
p <- ggplot(data = plot.dat, aes(study,mean,fill=group))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                       position=position_dodge(width=0.9),
                       color="#FF3333", size =2.5) + ylab("mean shift in acc. rate \n(95% boots. CI)\n")

p <- p + ggtitle("Shift in acceptance rate \n moving from neutral to gambling pictures\n") + coord_cartesian(ylim = c(-0.3, 0.5)) 
p + theme_la() + scale_fill_brewer(palette = "Paired")



## plotting the means and CIs for beta_gain, beta_loss and LA (using delta_method)
# behav HC
fm_behav = fm
cur_coef = fixef(fm_behav)
cur_cov  = vcov(fm_behav)
se_LA_behav_HC = deltamethod(~(-x3/x2),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3])/(cur_coef[2]))
LA_behav_HC    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)

# behav PG
cur_coef       = fixef(fm_behav)
cur_cov        = vcov(fm_behav)
se_LA_behav_PG = deltamethod(~(-(x3+x7)/(x2+x6)),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3]+cur_coef[7])/(cur_coef[2]+cur_coef[6]))
LA_behav_PG    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)
se_bg_behav_PG = deltamethod(~(x2+x6),cur_coef,cur_cov,ses=T)
cur_bg         = cur_coef[2]+cur_coef[6]
bg_behav_PG    = as.numeric(c(cur_bg,cur_bg-1.96*se_bg_behav_PG,cur_bg+1.96*se_bg_behav_PG))
se_bl_behav_PG = deltamethod(~(x3+x7),cur_coef,cur_cov,ses=T)
cur_bl         = cur_coef[3]+cur_coef[7]
bl_behav_PG    = as.numeric(c(cur_bl,cur_bl-1.96*se_bl_behav_PG,cur_bl+1.96*se_bl_behav_PG))

# CIs
nranpars <- length(getME(fm_behav,"theta"))
nfixpars <- length(fixef(fm_behav))
c1 <- confint(fm_behav,method="Wald", parm=(nranpars+1):(nranpars+nfixpars))
bg_behav_HC    = c(cur_coef[2],c1[2,])
bl_behav_HC    = c(cur_coef[3],c1[3,])

ci_behav = rbind(bg_behav_HC,bl_behav_HC,LA_behav_HC,bg_behav_PG,bl_behav_PG,LA_behav_PG)


# MRI HC
fm_behav = fm
cur_coef = fixef(fm_behav)
cur_cov  = vcov(fm_behav)
se_LA_behav_HC = deltamethod(~(-x3/x2),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3])/(cur_coef[2]))
LA_behav_HC    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)

# MRI PG
cur_coef       = fixef(fm_behav)
cur_cov        = vcov(fm_behav)
se_LA_behav_PG = deltamethod(~(-(x3+x7)/(x2+x6)),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3]+cur_coef[7])/(cur_coef[2]+cur_coef[6]))
LA_behav_PG    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)
se_bg_behav_PG = deltamethod(~(x2+x6),cur_coef,cur_cov,ses=T)
cur_bg         = cur_coef[2]+cur_coef[6]
bg_behav_PG    = as.numeric(c(cur_bg,cur_bg-1.96*se_bg_behav_PG,cur_bg+1.96*se_bg_behav_PG))
se_bl_behav_PG = deltamethod(~(x3+x7),cur_coef,cur_cov,ses=T)
cur_bl         = cur_coef[3]+cur_coef[7]
bl_behav_PG    = as.numeric(c(cur_bl,cur_bl-1.96*se_bl_behav_PG,cur_bl+1.96*se_bl_behav_PG))

# CIs
nranpars <- length(getME(fm_behav,"theta"))
nfixpars <- length(fixef(fm_behav))
c1 <- confint(fm_behav,method="Wald", parm=(nranpars+1):(nranpars+nfixpars))
bg_behav_HC    = c(cur_coef[2],c1[2,])
bl_behav_HC    = c(cur_coef[3],c1[3,])

ci_MRI = rbind(bg_behav_HC,bl_behav_HC,LA_behav_HC,bg_behav_PG,bl_behav_PG,LA_behav_PG)
ci_MRI$study   = "MRI"
ci_behav$study = "Behav"

ci_all = rbind(ci_MRI,ci_behav)
colnames(ci_all) = c("mean","lower","upper")
ci_all = (as.matrix(ci_all))
ci_df  = as.data.frame(matrix(as.numeric(ci_all),ncol = 3,byrow = F))
names(ci_df) = colnames(ci_all)
ci_df$param = as.factor(c("bg","bl","LA","bg","bl","LA","bg","bl","LA","bg","bl","LA"))
ci_df$group = as.factor(c("HC","HC","HC","PG","PG","PG","HC","HC","HC","PG","PG","PG"))
ci_df$study = as.factor(c("Behav","Behav","Behav","Behav","Behav","Behav","MRI","MRI","MRI","MRI","MRI","MRI"))

# plotting
p <- ggplot(data = ci_df, aes(group,mean,fill=param))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                       position=position_dodge(width=0.9),
                       color="#FF3333", size =2.5) + ylab("param est. \n(95% Wald or delta method CI)\n")

p <- p + ggtitle("Gain and loss sensitivity \n and loss aversion in study 2\n")
p <- p + theme_la() + scale_fill_brewer(palette = "Paired")
p + facet_grid(~study) + scale_fill_manual(labels=c(expression(beta[gain]), expression(beta[loss]), expression(lambda)),
                                           values= RColorBrewer::brewer.pal(3,"Paired"),
                                           name="fixed effect") +coord_cartesian(ylim = c(-0.8, 2.7))


## plotting the means and CIs for beta_gain, beta_loss and LA (using delta_method) [USING THE SIMPLE LA MODEL WITHOUT CATEGORY]
# behav HC
fm_behav = la_modg_ranef
cur_coef = fixef(fm_behav)
cur_cov  = vcov(fm_behav)
se_LA_behav_HC = deltamethod(~(-x3/x2),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3])/(cur_coef[2]))
LA_behav_HC    = c(cur_LA,cur_LA-1.96*se_LA_behav_HC,cur_LA+1.96*se_LA_behav_HC)

# behav PG
cur_coef       = fixef(fm_behav)
cur_cov        = vcov(fm_behav)
se_LA_behav_PG = deltamethod(~(-(x3+x5)/(x2+x4)),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3]+cur_coef[5])/(cur_coef[2]+cur_coef[4]))
LA_behav_PG    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)
se_bg_behav_PG = deltamethod(~(x2+x4),cur_coef,cur_cov,ses=T)
cur_bg         = cur_coef[2]+cur_coef[4]
bg_behav_PG    = as.numeric(c(cur_bg,cur_bg-1.96*se_bg_behav_PG,cur_bg+1.96*se_bg_behav_PG))
se_bl_behav_PG = deltamethod(~(x3+x5),cur_coef,cur_cov,ses=T)
cur_bl         = cur_coef[3]+cur_coef[5]
bl_behav_PG    = as.numeric(c(cur_bl,cur_bl-1.96*se_bl_behav_PG,cur_bl+1.96*se_bl_behav_PG))

# CIs
nranpars <- length(getME(fm_behav,"theta"))
nfixpars <- length(fixef(fm_behav))
c1 <- confint(fm_behav,method="Wald", parm=(nranpars+1):(nranpars+nfixpars))
bg_behav_HC    = c(cur_coef[2],c1[2,])
bl_behav_HC    = c(cur_coef[3],c1[3,])

ci_behav = rbind(bg_behav_HC,bl_behav_HC,LA_behav_HC,bg_behav_PG,bl_behav_PG,LA_behav_PG)


# MRI HC
fm_behav = fm
cur_coef = fixef(fm_behav)
cur_cov  = vcov(fm_behav)
se_LA_behav_HC = deltamethod(~(-x3/x2),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3])/(cur_coef[2]))
LA_behav_HC    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)

# MRI PG
cur_coef       = fixef(fm_behav)
cur_cov        = vcov(fm_behav)
se_LA_behav_PG = deltamethod(~(-(x3+x7)/(x2+x6)),cur_coef,cur_cov,ses=T)
cur_LA         = as.numeric(-(cur_coef[3]+cur_coef[7])/(cur_coef[2]+cur_coef[6]))
LA_behav_PG    = c(cur_LA,cur_LA-1.96*se_LA_behav_PG,cur_LA+1.96*se_LA_behav_PG)
se_bg_behav_PG = deltamethod(~(x2+x6),cur_coef,cur_cov,ses=T)
cur_bg         = cur_coef[2]+cur_coef[6]
bg_behav_PG    = as.numeric(c(cur_bg,cur_bg-1.96*se_bg_behav_PG,cur_bg+1.96*se_bg_behav_PG))
se_bl_behav_PG = deltamethod(~(x3+x7),cur_coef,cur_cov,ses=T)
cur_bl         = cur_coef[3]+cur_coef[7]
bl_behav_PG    = as.numeric(c(cur_bl,cur_bl-1.96*se_bl_behav_PG,cur_bl+1.96*se_bl_behav_PG))

# CIs
nranpars <- length(getME(fm_behav,"theta"))
nfixpars <- length(fixef(fm_behav))
c1 <- confint(fm_behav,method="Wald", parm=(nranpars+1):(nranpars+nfixpars))
bg_behav_HC    = c(cur_coef[2],c1[2,])
bl_behav_HC    = c(cur_coef[3],c1[3,])

ci_MRI = rbind(bg_behav_HC,bl_behav_HC,LA_behav_HC,bg_behav_PG,bl_behav_PG,LA_behav_PG)
ci_MRI$study   = "MRI"
ci_behav$study = "Behav"

ci_all = rbind(ci_MRI,ci_behav)
colnames(ci_all) = c("mean","lower","upper")
ci_all = (as.matrix(ci_all))
ci_df  = as.data.frame(matrix(as.numeric(ci_all),ncol = 3,byrow = F))
names(ci_df) = colnames(ci_all)
ci_df$param = as.factor(c("bg","bl","LA","bg","bl","LA","bg","bl","LA","bg","bl","LA"))
ci_df$group = as.factor(c("HC","HC","HC","PG","PG","PG","HC","HC","HC","PG","PG","PG"))
ci_df$study = as.factor(c("Behav","Behav","Behav","Behav","Behav","Behav","MRI","MRI","MRI","MRI","MRI","MRI"))

# plotting
p <- ggplot(data = ci_df, aes(group,mean,fill=param))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                       position=position_dodge(width=0.9),
                       color="#FF3333", size =2.5) + ylab("param est. \n(95% Wald or delta method CI)\n")

p <- p + ggtitle("Gain and loss sensitivity \n and loss aversion in study 2\n")
p <- p + theme_la() + scale_fill_brewer(palette = "Paired")
p + facet_grid(~study) + scale_fill_manual(labels=c(expression(beta[gain]), expression(beta[loss]), expression(lambda)),
                                           values= RColorBrewer::brewer.pal(3,"Paired"),
                                           name="fixed effect") +coord_cartesian(ylim = c(-0.8, 2.7))

## What's behind the gambling images?
hmdf <- aggregate(as.numeric(as.character(data_pdt$eda)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$cat),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","category","eda")

tmp        = aggregate(hmdf$eda,by = list(hmdf$subject, hmdf$HCPG,hmdf$category),FUN=mean)
names(tmp) = c("subject","HCPG","category","value")
tmp$HCPG   = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 
wide = reshape(tmp,direction = "wide",idvar=c("subject","HCPG"),timevar = c("category"))
wide = wide[-1]
wide = wide[-c(2,4,5)]
boot_pdt_df     = wide
boot_pdt_df_eda = boot_pdt_df
boot_pdt_df_eda$measure = "eda"

hmdf <- aggregate(as.numeric(as.character(data_pdt$imageRating1s)),by=list(data_pdt$subject, data_pdt$HCPG,data_pdt$cat),FUN="mean.rmna")
names(hmdf) <- c("subject", "HCPG","category","craving")
tmp        = aggregate(hmdf$craving,by = list(hmdf$subject, hmdf$HCPG,hmdf$category),FUN=mean)

names(tmp) = c("subject","HCPG","category","value")
tmp$HCPG   = factor(tmp$HCPG,levels=c("HC","PG"),labels=c("HC","PG")) 
wide = reshape(tmp,direction = "wide",idvar=c("subject","HCPG"),timevar = c("category"))
wide = wide[-1]
wide = wide[-c(2,4,5)]
boot_pdt_df     = wide
boot_pdt_df_cra = boot_pdt_df
boot_pdt_df_cra$measure = "craving"

boot_pdt_df = rbind(boot_pdt_df_eda,boot_pdt_df_cra)
boot_pdt_df = boot_pdt_df[c(1,3,2)]
boot_pdt_df = boot_pdt_df_cra

# get the mean and errorbar values
one.boot.mean = function(x) {
  tmp = one.boot(x,FUN=mean.rmna,R=2000)
  return(agk.mean.quantile(tmp$t))
}
plot.dat <- data.frame()
var.names <- c()
for (ii in 3:length(boot_pdt_df[1,])){
  # bootstrapped 95% CI
  bt.df     <- as.data.frame(as.list(aggregate(boot_pdt_df[ii],by = list(group = boot_pdt_df$HCPG, measure = boot_pdt_df$measure), FUN=one.boot.mean)))
  names(bt.df) <- c("group","measure","lower", "upper","mean")
  var.names <- c(var.names,rep(names(boot_pdt_df[ii]),length(levels(bt.df$group))))
  plot.dat <- rbind(plot.dat,bt.df)
}
plot.dat$variable <- var.names
plot.dat$variable <- as.factor(plot.dat$variable)
plot.dat$group = factor(plot.dat$group, levels=c("HC","PG"), labels = c("HC","PG"))
plot.dat$measure = as.factor(plot.dat$measure)


# plot
p <- ggplot(data = plot.dat, aes(measure,mean,fill=group))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                       position=position_dodge(width=0.9),
                       color="#FF3333", size =2.5) + ylab("mean shift in subjective rating \n(95% boots. CI)\n")

p <- p + ggtitle("Effect of gambling pictures \n on craving\n") + coord_cartesian(ylim = c(-0.3, 1.2)) 
p + theme_la() + scale_fill_brewer(palette = "Paired")



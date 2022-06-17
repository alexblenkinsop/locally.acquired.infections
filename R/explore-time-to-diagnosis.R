

do <- subset(dsubgraphtaxa,REP=='000' & SELECT!='Ams')

#do <- subset(do, HIV1_POS_D<2015 & HIV1_POS_D>2010)

# number infected in 2015
# HSX 34
# MSM 108

length(unique(do$ID[do$INF_D>=2015 & do$INF_D<2016]))
# 81

length(unique(do$ID[do$INF_D>=2015 & do$INF_D<2016 & do$HIV1_POS_D<=2020 & do$SELECT=='AmsHSX']))
length(unique(do$ID[do$INF_D>=2015 & do$INF_D<2016 & do$HIV1_POS_D<=2020 & do$SELECT=='AmsMSM']))

# get time to diagnosis distribution
do[, time:=HIV1_POS_D - INF_D]

dt <- subset(do,SELECT=='AmsHSX' )
dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=2010 & do$INF_D<2011 & do$HIV1_POS_D<=2020)
dt <- subset(do,SELECT=='AmsMSM' & do$INF_D>=2015 & do$INF_D<2016 & do$HIV1_POS_D<=2020)
n <- 34
n <- 108
require(fitdistrplus)
do.exp <- fitdist(dt$time, "exp")

par(mfrow = c(1,1))
denscomp(list(do.exp),
				  legendtext = c("exponential"), fitlty = 1)

gofstat(list(do.exp),
				fitnames = c("exponential"))

rexp(n,rate=do.exp$estimate)

# proportion undiagnosed by 2019 (i.e. within 3 years)
table(rexp(n,rate=do.exp$estimate)<=3)


rnorm(n,mean=mean(dt$time),sd=sd(dt$time))
table(rexp(n,rate=1/mean(dt$time))<=3)


## HSX
# overll
di <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(37,35,34,33,34,34,35,37,39))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=2010)
	di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
mean(di$und,na.rm=T)

# Non-b = 20-23%
di <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(37,35,34,33,34,34,35,37,39))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsHSX' & do$ST!='B' & do$INF_D>=2010 & do$INF_D<2015)
	di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
mean(di$und[di$year>=2015],na.rm=T)

#B = 10-12%
di <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(37,35,34,33,34,34,35,37,39))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsHSX' & do$ST=='B' & do$INF_D>=2010)
	di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
mean(di$und[di$year>=2015],na.rm=T)

#### MSM
## overall
dm <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(160,152,143,132,121,108,94,78,61))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsMSM' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsMSM' & do$INF_D>=2010)
	dm[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
dm[, weight:=N_inf/sum(N_inf)]
dm[, un_w:=und*weight]
mean(dm$und[di$year>=2015],na.rm=T)
sum(dm$un_w[di$year>=2015])

#B
dm <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(160,152,143,132,121,108,94,78,61))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsMSM' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsMSM' & do$ST=='B' & do$INF_D>=2010)
	dm[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
dm <- subset(dm,year>=2015)
dm[, weight:=N_inf/sum(N_inf)]
dm[, un_w:=und*weight]
mean(dm$und[dm$year>=2015],na.rm=T)
sum(dm$un_w[dm$year>=2015])


#non-B
dm <- data.table(year=c(2010,2011,2012,2013,2014,2015,2016,2017,2018),N_inf=c(160,152,143,132,121,108,94,78,61))
for(i in 2010:2018){
	#dt <- subset(do,SELECT=='AmsMSM' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsMSM' & do$ST!='B' & do$INF_D>=2010)
	dm[year==i,und:=(N_inf - sum(rexp(N_inf,rate=1/mean(dt$time))<=2019-i))/N_inf]
}
dm <- subset(dm,year>=2015)
dm[, weight:=N_inf/sum(N_inf)]
dm[, un_w:=und*weight]
mean(dm$und[dm$year>=2015],na.rm=T)
sum(dm$un_w[dm$year>=2015])

#########

## HSX
# overall
di <- data.table(year=c(2015,2016,2017,2018),N_inf=c(34,35,37,39))
for(i in 2015:2018){
	#dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=i & do$INF_D<(i+1) & do$HIV1_POS_D<=2020)
	dt <- subset(do,SELECT=='AmsHSX' & do$INF_D>=2010 & do$INF_D<2013)
	do.exp <- fitdist(dt$time, "weibull")
	
	par(mfrow = c(1,1))
	denscomp(list(do.exp),
					 legendtext = c("exponential"), fitlty = 1)
	di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=do.exp$estimate)<=2019-i))/N_inf]
}
mean(di$und,na.rm=T)

###########################
# BY risk group and subtype
###########################
require(data.table)
require(fitdistrplus)
require(truncdist)
dsubgraphtaxa <- readRDS('/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/branching_process_model/subgraphs_withmetadata.RDS')
do <- subset(dsubgraphtaxa,REP=='000' & SELECT!='Ams')
do[, time:=HIV1_POS_D - INF_D]

# HSX B
di <- data.table(year=c(2015,2016,2017,2018),N_inf=c(37,37,37,37))
dt <- subset(do,SELECT=='AmsHSX' & do$ST=='B' & do$INF_D>=2010 & do$INF_D<2013)

do.hsx.b.exp <- fitdist(dt$time[dt$time>0], distr = "exp")
do.hsx.b.weibull <- fitdist(dt$time[dt$time>0], distr = "weibull", method = "mle", lower = c(0, 0))
do.hsx.b.gamma <- fitdist(dt$time[dt$time>0], distr = "gamma")
do.hsx.b.lnorm <- fitdist(dt$time[dt$time>0], distr = "lnorm")
do.hsx.b.llogis <- fitdist(dt$time[dt$time>0], distr = "llogis")

par(mfrow = c(1,1))
denscomp(list(do.hsx.b.exp,do.hsx.b.weibull,do.hsx.b.gamma,do.hsx.b.lnorm),
				 legendtext = c("exponential","weibull","gamma","lnorm"), fitlty = 1)

gofstat(list(do.hsx.b.exp,do.hsx.b.weibull,do.hsx.b.gamma,do.hsx.b.lnorm),
				fitnames = c("exponential","weibull","gamma","lnorm"))

di <-  data.table(year=seq(2010,2018,1))
for(i in 2010:2018){
	#di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=do.hsx.b$estimate)<=2019-i))/N_inf]
	#di[year==i,und:=1 - pgamma(2019-i,shape=do.hsx.b$estimate[1],rate=do.hsx.b$estimate[2])]
	#di[year==i,und:=1 - pweibull(2019-i,shape=do.hsx.b.weibull$estimate[1],scale=do.hsx.b.weibull$estimate[2])]
	di[year==i,und:=1 - ptrunc(2019-i,"weibull",shape=do.hsx.b.weibull$estimate[1],scale=do.hsx.b.weibull$estimate[2], a=0, b=8)]
	#di[year==i,und:=1 - pexp(2019-i,rate=do.hsx.b.exp$estimate)]
}
g_hsx_b <- ggplot() + geom_line(data=di,aes(x=year,y=und)) +
	theme_bw() +
	scale_y_continuous(label = scales::percent) +
	labs(x='year of infection',y='% undiagnosed by 2019')

## average undiagnosed 2015-2019
# no weights
av_p1 <- mean(di$und[di$year %in% c(2010:2014)])
av_p2 <- mean(di$und[di$year %in% c(2015:2019)])
# weight by infected
dt <- subset(do,SELECT=='AmsHSX' & do$ST=='B' & do$INF_D>=2015 & do$INF_D<2020)
dt[, INF_Y:=floor(INF_D)]
inf <- dt[, list(inf=length(ID)),by=c('INF_Y')]

# hsx nonb
#di <- data.table(year=c(2015,2016,2017,2018),N_inf=c(34,35,37,39))
di <- data.table(year=c(2015,2016,2017,2018),N_inf=c(37,37,37,37))
dt <- subset(do,SELECT=='AmsHSX' & do$ST!='B' & do$INF_D>=2010 & do$INF_D<2013)

do.hsx.nb.exp <- fitdist(dt$time[dt$time>0], "exp")
do.hsx.nb.gamma <- fitdist(dt$time[dt$time>0], distr = "gamma")
do.hsx.nb.weibull <- fitdist(dt$time[dt$time>0], distr = "weibull", method = "mle", lower = c(0, 0))
#do.hsx.nb.invweibull <- fitdist(dt$time[dt$time>0], distr = "invweibull")
#do.hsx.nb.invgamma <- fitdist(dt$time[dt$time>0], distr = "invgamma")
#do.hsx.nb.lnorm <- fitdist(dt$time[dt$time>0], distr = "lnorm")
par(mfrow = c(1,1))
denscomp(list(do.hsx.nb.exp,do.hsx.nb.weibull,do.hsx.nb.gamma),
				 legendtext = c("exponential","weibull","gamma"), fitlty = 1)

gofstat(list(do.hsx.nb.exp,do.hsx.nb.weibull,do.hsx.nb.gamma),
				fitnames = c("exponential","weibull","gamma"))

di <-  data.table(year=seq(2010,2018,1))
for(i in 2010:2018){
	#di[year==i,und:=(N_inf - sum(rexp(N_inf,rate=do.hsx.nb$estimate)<=2019-i))/N_inf]
	#di[year==i,und:=1 - pgamma(2019-i,shape=do.hsx.nb$estimate[1],rate=do.hsx.nb$estimate[2])]
	#di[year==i,und:=1 - pexp(2019-i,rate=do.hsx.nb.exp$estimate)]
	di[year==i,und:=1 - ptrunc(2019-i,"exp",rate=do.hsx.nb.exp$estimate, a=0, b=8)]
}
g_hsx_nb <- ggplot() + geom_line(data=di,aes(x=year,y=und)) +
	theme_bw() +
	scale_y_continuous(label = scales::percent) +
	labs(x='year of infection',y='% undiagnosed by 2019')
## average undiagnosed 2015-2019
# no weights
av_p1 <- mean(di$und[di$year %in% c(2010:2014)])
av_p2 <- mean(di$und[di$year %in% c(2015:2019)])

# MSM
#B
#dm <- data.table(year=c(2015,2016,2017,2018),N_inf=c(108,94,78,61))
dm <- data.table(year=c(2015,2016,2017,2018),N_inf=c(80,80,80,80))
dt <- subset(do,SELECT=='AmsMSM' & do$ST=='B' & do$INF_D>=2010 & do$INF_D<2013)
do.msm.b.exp <- fitdist(dt$time[dt$time>0], "exp")
do.msm.b.gamma <- fitdist(dt$time[dt$time>0], distr = "gamma")
do.msm.b.weibull <- fitdist(dt$time[dt$time>0], distr = "weibull", method = "mle", lower = c(0, 0))
gofstat(list(do.msm.b.exp,do.msm.b.weibull,do.msm.b.gamma),
				fitnames = c("exponential","weibull","gamma"))
dm <-  data.table(year=seq(2010,2018,1))
for(i in 2010:2018){
	#dm[year==i,und:=(N_inf - sum(rexp(N_inf,rate=do.msm.b$estimate)<=2019-i))/N_inf]
	#dm[year==i,und:=1 - pgamma(2019-i,shape=do.msm.b$estimate[1],rate=do.msm.b$estimate[2])]
	#dm[year==i,und:=1 - pweibull(2019-i,shape=do.msm.b.weibull$estimate[1],scale=do.msm.b.weibull$estimate[2])]
	dm[year==i,und:=1 - ptrunc(2019-i,"weibull",shape=do.msm.b.weibull$estimate[1],scale=do.msm.b.weibull$estimate[2], a=0, b=8)]
}
av_p1 <- mean(dm$und[dm$year %in% c(2010:2014)])
av_p2 <- mean(dm$und[dm$year %in% c(2015:2019)])

g_msm_b <- ggplot() + geom_line(data=dm,aes(x=year,y=und)) +
	theme_bw() +
	scale_y_continuous(label = scales::percent) +
	labs(x='year of infection',y='% undiagnosed by 2019')
## average undiagnosed 2015-2019
# no weights
av <- mean(dm$und)


#non-B
dt <- subset(do,SELECT=='AmsMSM' & do$ST!='B' & do$INF_D>=2010 & do$INF_D<2013)
do.msm.nb.exp <- fitdist(dt$time[dt$time>0], "exp")
do.msm.nb.gamma <- fitdist(dt$time[dt$time>0], distr = "gamma")
do.msm.nb.weibull <- fitdist(dt$time[dt$time>0], distr = "weibull", method = "mle", lower = c(0, 0))
gofstat(list(do.msm.nb.exp,do.msm.nb.weibull,do.msm.nb.gamma),
				fitnames = c("exponential","weibull","gamma"))
#dm <- data.table(year=c(2015,2016,2017,2018),N_inf=c(108,94,78,61))
dm <- data.table(year=c(2015,2016,2017,2018),N_inf=c(80,80,80,80))
dp <- data.table(x=seq(1,8,1))
dp[,y:=pexp(x,rate=do.msm.nb.exp$estimate)]
dm <-  data.table(year=seq(2010,2018,1))
for(i in 2010:2018){
	#dm[year==i,und:=(N_inf - sum(rexp(N_inf,rate=do.msm.nb$estimate)<=2019-i))/N_inf]
	#dm[year==i,und:=1 - pgamma(2019-i,shape=do.msm.nb$estimate[1],rate=do.msm.nb$estimate[2])]
	#dm[year==i,und:=1 - pexp(2019-i,rate=do.msm.nb.exp$estimate)]
	dm[year==i,und:=1 - ptrunc(2019-i,"exp",rate=do.msm.nb.exp$estimate, a=0, b=8)]
}
g_msm_nb <- ggplot() + geom_line(data=dm,aes(x=year,y=und)) +
	theme_bw() +
	scale_y_continuous(label = scales::percent) +
	labs(x='year of infection',y='% undiagnosed by 2019')
## average undiagnosed 2015-2019
# no weights
av_p1 <- mean(dm$und[dm$year %in% c(2010:2014)])
av_p2 <- mean(dm$und[dm$year %in% c(2015:2019)])


#plotdist(dt$time, histo = TRUE, demp = TRUE)
png(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/inf-diag-time-distribution-gamma.png')

par(mfrow = c(2,2))
#denscomp(list(do.hsx.b,do.hsx.nb,do.msm.b,do.msm.nb),
#											legendtext = c("HSX B",'HSX non-B','MSM B','MSM non-B'), fitlty = 1)
denscomp(list(do.hsx.b),
				 legendtext = c("HSX B"), fitlty = 1)
denscomp(list(do.hsx.nb),
				 legendtext = c('HSX non-B'), fitlty = 1)
denscomp(list(do.msm.b),
				 legendtext = c('MSM B'), fitlty = 1)
denscomp(list(do.msm.nb),
				 legendtext = c('MSM non-B'), fitlty = 1)
dev.off()

png(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/inf-diag-time-distribution-comparedist-hsx2007.png')

par(mfrow = c(2,2))
#denscomp(list(do.hsx.b,do.hsx.nb,do.msm.b,do.msm.nb),
#											legendtext = c("HSX B",'HSX non-B','MSM B','MSM non-B'), fitlty = 1)
denscomp(list(do.hsx.b.exp,do.hsx.b.weibull,do.hsx.b.gamma),
				 legendtext = c("exponential","weibull","gamma"), main= "HSX B", fitlty = 1)
denscomp(list(do.hsx.nb.exp,do.hsx.nb.weibull,do.hsx.nb.gamma),
				 legendtext = c("exponential","weibull","gamma"), main= "HSX non-B", fitlty = 1)
denscomp(list(do.msm.b.exp,do.msm.b.weibull,do.msm.b.gamma),
				 legendtext = c("exponential","weibull","gamma"), main= "MSM B", fitlty = 1)
denscomp(list(do.msm.nb.exp,do.msm.nb.weibull,do.msm.nb.gamma),
				 legendtext = c("exponential","weibull","gamma"), main= "MSM non-B", fitlty = 1)
dev.off()

g_hsx_b_l <- annotate_figure(g_hsx_b,	top = text_grob("HSX B",size=20))
g_hsx_nb_l <- annotate_figure(g_hsx_nb,	top = text_grob("HSX non-B",size=20))
g_msm_b_l <- annotate_figure(g_msm_b,	top = text_grob("MSM B",size=20))
g_msm_nb_l <- annotate_figure(g_msm_nb,	top = text_grob("MSM non-B",size=20))

g <- ggarrange(g_hsx_b_l,g_hsx_nb_l,g_msm_b_l,g_msm_nb_l,nrow=2,ncol=2,align="hv")
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/pr_undiagnosed_hsx2007.png',g,h=10,w=12)


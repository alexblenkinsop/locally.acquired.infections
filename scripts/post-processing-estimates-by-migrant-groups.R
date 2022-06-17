# post-processing-make-mwmb-figure.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running post-processing-make-mwmb-figures.R\n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(bayesplot, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggpubr, quietly = TRUE))
suppressMessages(library(viridis, quietly = TRUE))
suppressMessages(library(forcats, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

args_dir <- list()
args_dir[['stanModelFileMSM']] <- 'branching_process_210810b_cmdstan'
args_dir[['stanModelFileHSX']] <- 'branching_process_210810m_cmdstan'
args_dir[['period']] <- '2014-2018'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['job_name']] <- 'undiagnosed_untilmay'
args_dir[['infdate']] <- 1
args_dir[['start_d']] <- 2014
args_dir[['end_d']] <- 2019


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFileMSM')	
	stopifnot(args_line[[3]]=='-stanModelFileHSX')	
	stopifnot(args_line[[5]]=='-in_dir')
	stopifnot(args_line[[7]]=='-job_name')
	stopifnot(args_line[[9]]=='-period')
	stopifnot(args_line[[11]]=='-start_d')
	stopifnot(args_line[[13]]=='-end_d')
	stopifnot(args_line[[15]]=='-source_dir')
	
	args_dir <- list()
	args_dir[['stanModelFileMSM']] <- args_line[[2]]
	args_dir[['stanModelFileHSX']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['job_name']] <- args_line[[8]]
	args_dir[['period']] <- args_line[[10]]
	args_dir[['start_d']] <- args_line[[12]]
	args_dir[['end_d']] <- args_line[[14]]
	args_dir[['source_dir']] <- args_line[[16]]
	args_dir[['infdate']] <- 1
	
} 

## load functions
source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'post-processing-summary-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileMSM']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_MSM')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_MSM')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileMSM , "-", args_dir$job_tag)
cat(outfile.base)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.MSM.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_mwmb.RDS')
cat("\n read RDS:", file)
msm_st_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_mwmb.RDS')
msm_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
local.MSM.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
st.MSM.p1 <- readRDS(file)
st.MSM.p1 <- merge(st.MSM.p1, plot.pars.basic.MSM.p1$ds,by.x='subtype',by.y='subtypes')
file <- paste0(outfile.base,'-stanout-bplace-gqs.RDS')
cat("\n read RDS:", file)
fit.bplace.MSM <- readRDS(file)

msm_st_p1[, trsm:='MSM']
msm_p1[, trsm:='MSM']
msm_st_p1[, time:='2014-2019']
msm_p1[, time:='2014-2019']

args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFileHSX']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_HSX')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_HSX')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileHSX , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.HSX.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_mwmb.RDS')
cat("\n read RDS:", file)
hsx_st_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_mwmb.RDS')
hsx_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
local.HSX.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
st.HSX.p1 <- readRDS(file)
st.HSX.p1 <- merge(st.HSX.p1, plot.pars.basic.HSX.p1$ds,by.x='subtype',by.y='subtypes')
file <- paste0(outfile.base,'-stanout-bplace-gqs.RDS')
cat("\n read RDS:", file)
fit.bplace.HSX <- readRDS(file)

hsx_st_p1[, trsm:='HSX']
hsx_p1[, trsm:='HSX']
hsx_st_p1[, time:='2014-2019']
hsx_p1[, time:='2014-2019']

cat(" \n -------------------------------- define plotting functions -------------------------------- \n")

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

quantiles_95 <- function(x) {
	r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}

col_pal = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

cat(" \n -------------------------------- make overall local infections -------------------------------- \n")

st.MSM.p1[, time:='2014-2019']
st.HSX.p1[, time:='2014-2019']
st.MSM.p1[, trsm:='MSM']
st.HSX.p1[, trsm:='HSX']
st <- rbind(st.MSM.p1,st.HSX.p1)
setnames(st,'subtypes_name','ST')

st$ST <- factor(st$ST,levels=c("B","nonB"),labels=c("B","Non-B"))
st[, mwmb:='All']
#st[, period:=NULL]

#st$TRANSM <- factor(st$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
#setnames(st,'TRANSM','trsm')
setnames(st,'prop','p_st')
st_all <- st

local.MSM.p1[, time:='2014-2019']
local.HSX.p1[, time:='2014-2019']
dat <- rbind(local.MSM.p1,local.HSX.p1)
#dat$TRANSM <- factor(dat$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
#setnames(dat,'TRANSM','trsm')
setnames(dat,'TRM','trsm')
dat[, mwmb:='All']
dat_all <- dat

##########
st <- rbind(msm_st_p1,hsx_st_p1)
st <- merge(st,st_all,by=c('mwmb','ST','p_st','trsm','time'),all=T)
st[ST=='nonB', ST:='Non-B']

st$ST <- factor(st$ST,levels=c("Non-B", "B"))
#st$ST <- factor(st$ST,levels=c("06cpx","G","02AG","A1","D","C","01AE", "B"))

st$mwmb_lab <- factor(st$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','S. America &\n Caribbean','Other'))



dat <- rbind(msm_p1,hsx_p1)
dat <- merge(dat,dat_all,by=c('mwmb','trsm','time','iteration','inf_Ams'),all=T)

dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','S. America &\n Caribbean','Other'))

dat_s <- copy(dat)

cat(" \n -------------------------------- make labels -------------------------------- \n")

infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')

dsubgraphtaxa <- readRDS(infile.subgraphs)

if(args_dir$infdate==1){
	dsubgraphtaxa[INF_D>=args_dir$start_d & INF_D<args_dir$end_d, keep:=1]
}else{
	dsubgraphtaxa[HIV1_POS_D>=args_dir$start_d & HIV1_POS_D<args_dir$end_d, keep:=1]
}
dsubgraphtaxa <- subset(dsubgraphtaxa, REP=='000' & SELECT %in% c('AmsHSX','AmsMSM') & keep==1)


geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']

load(file='/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200917/misc/200917_sequence_labels.rda')
dseq <- merge(dseq,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
setnames(dseq,'ORIGIN','BIRTH_COUNTRY_ISO')
dsubgraphtaxa <- merge(dsubgraphtaxa,unique(subset(dseq,select=c('PATIENT','BIRTH_COUNTRY_ISO','WRLD'))),by.x='ID',by.y='PATIENT',all.x=T)

## msm
dsubgraphtaxa[TRANSM=='MSM', mwmb:="Other"]
dsubgraphtaxa[TRANSM=='MSM' & BIRTH_COUNTRY %in% c("Netherlands"), mwmb:="NL"]
# western countires (non-NL)
dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("WEurope","NorthAm","Oceania") & BIRTH_COUNTRY!='Netherlands', mwmb:="G1"]
# eastern and central europe
dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("EEurope", "CEurope"), mwmb:="G2"]
# caribbean and south america
dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]

## hsx
dsubgraphtaxa[TRANSM=='HSX', mwmb:="Other"]
dsubgraphtaxa[TRANSM=='HSX' & BIRTH_COUNTRY %in% c("Netherlands"), mwmb:="NL"]
# sub-saharan africa
dsubgraphtaxa[TRANSM=='HSX' & LOC_BIRTH %in% c("Africa"), mwmb:="G4"]
# caribbean and south america
dsubgraphtaxa[TRANSM=='HSX' & LOC_BIRTH %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]

cases <- dsubgraphtaxa[, list(N=length(unique(ID))),by=c('TRANSM','mwmb')]
cases_all <- dsubgraphtaxa[, list(N=length(unique(ID))),by=c('TRANSM')]
cases_all[,mwmb:='All']
cases <- merge(cases,cases_all,by=c('TRANSM','mwmb','N'),all=T)
setnames(cases,'TRANSM','trsm')

`%notin%` <- Negate(`%in%`)

## load birth places from GQs
bp.MSM <- as.data.table( reshape2::melt( fit.bplace.MSM ) )
setnames(bp.MSM, 1:6, c('iteration','subtype','bplace','index_cases','N','chaintype'))
bp.MSM[, trsm:='MSM']
locs <- unique(dat$mwmb[dat$trsm=='MSM'])
labs <- data.table(mwmb=sort(locs[locs!='All']),bplace=seq(1:length(locs[locs!='All'])))
bp.MSM <- merge(bp.MSM, labs,by='bplace',all=T)

bp.HSX <- as.data.table( reshape2::melt( fit.bplace.HSX ) )
setnames(bp.HSX, 1:6, c('iteration','subtype','bplace','index_cases','N','chaintype'))
bp.HSX[, trsm:='HSX']
locs <- unique(dat$mwmb[dat$trsm=='HSX'])
labs <- data.table(mwmb=sort(locs[locs!='All']),bplace=seq(1:length(locs[locs!='All'])))
bp.HSX <- merge(bp.HSX, labs,by='bplace',all=T)

bp <- rbind(bp.MSM,bp.HSX)
bp <- subset(bp,is.finite(N))
tmp <- bp[, list(N=sum(N)),by=c('trsm','iteration','bplace','mwmb')]
bp_s <- copy(tmp)
tmp <- tmp[, list(N = quantile(N, prob=ps,na.rm=T),
									 q_label=p_labs),by=c('trsm','bplace','mwmb')]		

tmp2 <- bp[, list(N=sum(N)),by=c('trsm','iteration')]
# save samples by migrant group and totals
bp_all_s <- copy(tmp2)
bp_all_s[, mwmb:='All']
bp_s <- merge(bp_s,bp_all_s,by=c('trsm','mwmb','iteration','N'),all=T)
tmp2 <- tmp2[, list(mwmb='All',N = quantile(N, prob=ps,na.rm=T),
									q_label=p_labs),by=c('trsm')]		
bp <- merge(tmp,tmp2,by=c('trsm','mwmb','N','q_label'),all=T)
setnames(bp,'N','N_star')
bp <- dcast.data.table(bp, trsm+mwmb~q_label, value.var='N_star')
bp[, L:= paste0( round(M, d=0), ' [',  round(CL, d=0),'-', round(CU, d=0),']')]
set(bp, NULL, c('CL','CU','M'),NULL)
# add subtype labs
# bp <- merge(bp, pars.basic$ds,by.x='subtype',by.y='subtypes')
# add birthplace labs (groups are ordered alphabetically)
# labs <- data.table(mwmb=sort(unique(st_mwmb$mwmb)),bplace=seq(1:length(unique(st_mwmb$mwmb))))
# bp <- merge(bp, labs,by='bplace',all=T)

dat <- merge(dat,cases,by=c('trsm','mwmb'),all=T)
dat <- merge(dat,bp,by=c('trsm','mwmb'),all=T)
#dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','S. America &\n Caribbean','Other'))
dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,Oceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan Africa','S. America &\n Caribbean','Other'))
#dat[, xlab:=paste0(mwmb_lab,' \n(N=',N,',\nN*=',L,')')]
dat[, xlab:=paste0(mwmb_lab,' \n(N=',N,', N*=',L,')')]

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

order <- unique(subset(dat,select=c('trsm','mwmb_lab','xlab')))
order <- order[order(trsm,mwmb_lab),]
#dat$xlab <- factor(dat$xlab,levels=c("All \n(N=304)","All \n(N=44)","NL \n(N=147)","NL \n(N=13)",
#																		 "S. America &\n Caribbean \n(N=54)","S. America &\n Caribbean \n(N=14)",
#																		 "W.Europe,\nN.America,\nOceania \n(N=42)","E. & C. Europe \n(N=18)" ,
#																		 "Sub-Saharan\nAfrica \n(N=12)","Other \n(N=43)","Other \n(N=5)"))
dat$xlab <- factor(dat$xlab,levels=unique(order$xlab))

cat(" \n -------------------------------- plot -------------------------------- \n")

subt_msm <- ggplot(subset(st,trsm=='MSM'), aes(x=mwmb_lab)) +
	geom_bar(aes(y=p_st,fill=ST), stat='identity',position = "fill",width=0.9) +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#facet_grid(time ~ .) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Amsterdam MSM \n Proportion", x='',fill='Subtype') +
	ggsci::scale_fill_npg(alpha=1)


subt_hsx <- ggplot(subset(st,trsm=='HSX'), aes(x=mwmb_lab)) +
	geom_bar(aes(y=p_st,fill=ST), stat='identity',position = "fill",width=0.9) +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#facet_grid(time ~ .) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Amsterdam heterosexuals \n Proportion", x='Birthplace',fill='Subtype') +
	ggsci::scale_fill_npg(alpha=1)

leg <- get_legend(subt_hsx +
										guides(fill = guide_legend(nrow = 1,
																							 #title.position="top",
																							 #title.theme = element_text(
																							 #size = 30,
																							 #face = "bold"
																							 #)
										)) +
										theme(legend.position = "bottom"))


g_mwmb_msm <- ggplot(subset(dat,trsm=='MSM'), aes(x=mwmb_lab, y=inf_Ams,fill=mwmb_lab)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#facet_grid(time ~ .) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)
	#ggsci::scale_color_lancet() +
	#ggsci::scale_fill_lancet()

g_mwmb_hsx <- ggplot(subset(dat,trsm=='HSX'), aes(x=mwmb_lab, y=inf_Ams,fill=mwmb_lab)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#facet_grid(time ~ .) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside",
				legend.text = element_blank()
				, legend.title = element_blank()) +
	#scale_fill_manual(guide=FALSE) +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='Birthplace') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)
	#ggsci::scale_color_lancet() +
	#ggsci::scale_fill_lancet()

g_st <- ggarrange(subt_msm,subt_hsx,ncol=1,
									 align=c("hv"))
g_inf <- ggarrange(g_mwmb_msm,g_mwmb_hsx,ncol=1,
									align=c("hv"))
g_all <- ggarrange(g_st,g_inf,ncol=2,align=c('hv'))
#g_all <- ggarrange(g_msm,g_hsx,ncol=1,align=c('hv'))
g_l <- ggarrange(g_all,leg,ncol=1,heights=c(1,0.05))
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX','.png'), g_l, w=12, h=10)


g_mwmb <- ggplot(subset(dat, time=='2014-2019'), aes(x=xlab, y=inf_Ams,fill=xlab)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	facet_grid(. ~ trsm,scales="free_x") +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text(angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX_2014-2019_cases','.png'), g_mwmb, w=20, h=10)

g_st <- ggarrange(subt_msm+theme(legend.position="bottom"),subt_hsx+theme(legend.position="bottom"),ncol=2,align=c("hv"),common.legend=T)
g_all <- ggarrange(g_st,g_mwmb ,ncol=1,align=c('hv'))
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX','.png'), g_all, w=14, h=15)

## make plot with bars instead

tab <- dat[, list(qs= quantile(inf_Ams, prob=ps), qlab=p_labs), by=c('trsm','mwmb_lab','xlab','time')]
ds <- dcast(tab,trsm+mwmb_lab+time+xlab~qlab,value.var='qs')
g_mwmb <- ggplot(subset(ds), aes(x=xlab, y=M,fill=xlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=xlab,ymin=CL, ymax=CU,fill=xlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_wrap(. ~ trsm,scales="free_y") +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	coord_flip() +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text( size=26,vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal) +
  scale_colour_manual(values=col_pal)
	#ggsci::scale_color_lancet() +
	#ggsci::scale_fill_lancet()
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2019_cases_barplot_big','.png'), g_mwmb, w=22, h=10)
annotate_figure(g_mwmb,fig.lab = "A", fig.lab.face = "bold")
g_mwmb <- ggarrange(g_mwmb,ncol=1,labels=c('A'),font.label=list(size=40))
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2019_cases_barplot_lab','.png'), g_mwmb, w=22, h=10)

## save numbers

tab <- dat[, list(qs= quantile(inf_Ams, prob=ps), qlab=p_labs), by=c('trsm','mwmb_lab','time')]
par <- dcast.data.table(tab, trsm+mwmb_lab+time~qlab, value.var='qs')
par[, L:= paste0( round(M*100, d=0), '% [',  round(CL*100, d=0),'-', round(CU*100, d=0),'%]')]
write.csv(par,file=paste0(outfile.base,'-','local_infections_migrantgroup.csv'))

cat(" \n -------------------------------- calculate number of locally acquired infections -------------------------------- \n")

dat <- merge(dat_s,bp_s,by=c('trsm','mwmb','iteration'),all=T)
dat[, N_local:=floor(inf_Ams*N)]

dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. \nEurope','S. America \n& Caribbean','Sub-Saharan \nAfrica','S. America \n& Caribbean','Other'))
dat[, fb:='Foreign-born']
dat[mwmb=='NL', fb:='NL']
dat_c <- subset(dat,mwmb_lab!='All')

tab <- dat[, list(pct_local= quantile(inf_Ams, prob=ps), 
									N_local= quantile(N_local, prob=ps),
									qlab=p_labs), by=c('trsm','mwmb_lab','time')]
dat <- dcast(tab,trsm+mwmb_lab+time~qlab,value.var=c('pct_local','N_local'))
dat[, labs:=paste0(round(pct_local_M*100,1),"%")]
dat[, pct_local_L:=paste0(round(pct_local_M*100,1),'% [',round(pct_local_CL*100,1),'-',round(pct_local_CU*100,1),'%]')]
dat[, N_local_L:=paste0(N_local_M,' [',N_local_CL,'-',round(N_local_CU,0),']')]

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
dat[, mwmb_lab_trsm:=paste0(trsm,' ',mwmb_lab)]
dat <- dat[order(trsm,mwmb_lab),]
write.csv(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup.csv'))
saveRDS(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup.RDS'))

order <- unique(subset(dat,select=c('trsm','mwmb_lab','mwmb_lab_trsm')))
order <- order[order(trsm,mwmb_lab),]
dat$mwmb_lab_trsm <- factor(dat$mwmb_lab_trsm,levels=unique(order$mwmb_lab_trsm))
dat_m <- copy(dat)

g_mwmb_N <- ggplot(subset(dat), aes(x=mwmb_lab, y=N_local_M,fill=mwmb_lab_trsm)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mwmb_lab,ymin=N_local_CL, ymax=N_local_CU,fill=mwmb_lab_trsm),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(aes(label=labs), position=position_dodge(width=0.9), vjust=-2, hjust=0,size=8) +
	facet_wrap(. ~ trsm,scales="free_x") +
	scale_y_continuous(expand=expansion(mult = c(0, .05))) +
										 #, breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 #limits = c(0, 1),labels=scales::percent) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text( size=24, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Number of infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2018_N_infections','.png'), g_mwmb_N, w=22, h=10)
annotate_figure(g_mwmb,fig.lab = "B", fig.lab.face = "bold")
g_mwmb_N <- ggarrange(g_mwmb_N,ncol=1,labels=c('B'),font.label=list(size=40))
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2018_N_infections','.png'), g_mwmb_N, w=22, h=10)

dat_c <- dat_c[, list(N_local=sum(N_local)),
										by=c('trsm','fb','time','iteration')]
tab <- dat_c[, list(N_local= quantile(N_local, prob=ps),
									qlab=p_labs), by=c('trsm','fb','time')]
dat <- dcast(tab,trsm+fb+time~qlab,value.var=c('N_local'))
dat[, N_local_L:=paste0(M,' [',CL,'-',round(CU,0),']')]

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
dat$fb <- factor(dat$fb,levels=c('NL','Foreign-born'),labels=c('Dutch-born','Foreign-born'))
dat <- dat[order(trsm,fb),]

write.csv(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup_agg-foreign-born.csv'))

saveRDS(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup_agg-foreign-born.RDS'))

## get % of all local infections in FB by trm group
dl <- dcast(dat_c,trsm+iteration+time~fb,value.var=c('N_local'))
dl[, pct_fb:=`Foreign-born`/(`Foreign-born`+NL)*100]
tab <- dl[, list(pct_fb= quantile(pct_fb, prob=ps),
										qlab=p_labs), by=c('trsm','time')]
dat <- dcast(tab,trsm+time~qlab,value.var=c('pct_fb'))
dat[, pct_fb_L:=paste0(round(M,0),' [',round(CL,0),'-',round(CU,0),']')]
write.csv(dat,file=paste0(outfile.base,'-','pct_local_infections_in_foreign-born.csv'))

## get % of all local infections in FB
dl <- dcast(dat_c,iteration+time~trsm+fb,value.var=c('N_local'))
dl[, pct_fb:=(`MSM_Foreign-born`)/(`HSX_Foreign-born` + `MSM_Foreign-born` + HSX_NL + MSM_NL)*100]
tab <- dl[, list(pct_fb= quantile(pct_fb, prob=ps),
								 qlab=p_labs), by=c('time')]
dat <- dcast(tab,time~qlab,value.var=c('pct_fb'))
dat[, pct_fb_L:=paste0(round(M,0),' [',round(CL,0),'-',round(CU,0),']')]
write.csv(dat,file=paste0(outfile.base,'-','pct_local_infections_in_foreign-born-MSM.csv'))


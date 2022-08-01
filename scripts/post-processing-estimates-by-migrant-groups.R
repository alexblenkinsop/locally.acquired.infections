# post-processing-make-mwmb-figure.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running post-processing-estimates-by-migrant-groups.R\n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggpubr, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

args_dir <- list()
args_dir[['stanModelFileMSM']] <- 'branching_process_210810b_cmdstan'
args_dir[['stanModelFileHSX']] <- 'branching_process_210810m_cmdstan'
args_dir[['period']] <- '2014-2018'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['job_name']] <- 'undiagnosed_weighted_inf_rate'
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
cat("\n read RDS:", file)
msm_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
cat("\n read RDS:", file)
local.MSM.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
cat("\n read RDS:", file)
st.MSM.p1 <- readRDS(file)
st.MSM.p1 <- merge(st.MSM.p1, plot.pars.basic.MSM.p1$ds,by.x='subtype',by.y='subtypes')
file <- paste0(outfile.base,'-stanout-bplace-gqs.RDS')
cat("\n read RDS:", file)
fit.bplace.MSM <- readRDS(file)

msm_st_p1[, trsm:='MSM']
msm_p1[, trsm:='MSM']
msm_st_p1[, time:=args_dir[['period']]]
msm_p1[, time:=args_dir[['period']]]

args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileHSX']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_HSX')
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
cat("\n read RDS:", file)
hsx_p1 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
cat("\n read RDS:", file)
local.HSX.p1 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
cat("\n read RDS:", file)
st.HSX.p1 <- readRDS(file)
st.HSX.p1 <- merge(st.HSX.p1, plot.pars.basic.HSX.p1$ds,by.x='subtype',by.y='subtypes')
file <- paste0(outfile.base,'-stanout-bplace-gqs.RDS')
cat("\n read RDS:", file)
fit.bplace.HSX <- readRDS(file)

hsx_st_p1[, trsm:='HSX']
hsx_p1[, trsm:='HSX']
hsx_st_p1[, time:=args_dir[['period']]]
hsx_p1[, time:=args_dir[['period']]]

infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')
infile.chaintype.infd <- file.path(args_dir$source_dir,'data','subgraph_metadata',paste0('subgraph_classification_infdate_',args_dir$start_d,'.csv'))
infile.chaintype.diagd <- file.path(args_dir$source_dir,'data','subgraph_metadata',paste0('subgraph_classification_diagdate_',args_dir$start_d,'.csv'))
infile.bplaces.trm.mwmb <- file.path(args_dir$source_dir,'data','patient_data',paste0('birthplaces_chaintype_',args_dir$start_d,'.csv'))

cat(" \n -------------------------------- define plotting functions -------------------------------- \n")

`%notin%` <- Negate(`%in%`)

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

quantiles_95 <- function(x) {
	r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}

col_pal = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

cat(" \n -------------------------------- make overall local infections -------------------------------- \n")

st.MSM.p1[, time:=args_dir[['period']]]
st.HSX.p1[, time:=args_dir[['period']]]
st.MSM.p1[, trsm:='MSM']
st.HSX.p1[, trsm:='HSX']
st <- rbind(st.MSM.p1,st.HSX.p1)
setnames(st,'subtypes_name','ST')

st$ST <- factor(st$ST,levels=c("B","nonB"),labels=c("B","Non-B"))
st[, mwmb:='All']
setnames(st,'prop','p_st')
st_all <- st

local.MSM.p1[, time:=args_dir[['period']]]
local.HSX.p1[, time:=args_dir[['period']]]
dat <- rbind(local.MSM.p1,local.HSX.p1)
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

dsubgraphtaxa <- readRDS(infile.subgraphs)

cat(paste0('\n Load subgraph classifications by ', args_dir$start_d, '\n'))

if(args_dir$infdate==1){
	dc <- data.table(read.csv(infile.chaintype.infd, header=T))
}else{
	dc <- data.table(read.csv(infile.chaintype.diagd, header=T))
}
dc[, ST2:= 'nonB']
dc[ST=='B', ST2:= 'B']
setnames(dc,c('ST','ST2'),c('ST2','ST'))

dsubgraphtaxa <- merge(dsubgraphtaxa,dc,by=c('SELECT','NAME','FULL_NAME','ST','ST_CLADE'),all.x=T)

# count chains which started before startd
dsubgraphtaxa[period=='pre-startd', chaintype := "pre-existing"]
dsubgraphtaxa[period=='post-startd', chaintype := "emergent"]

# select patients in study window
dsubgraphtaxa[inf_after_startd==1 & inf_after_endd==0, keep:=1]

dsubgraphtaxa <- subset(dsubgraphtaxa, keep==1)
dsubgraphtaxa[, TRANSM:= gsub('Ams','',SELECT)]

# read in number of individuals by risk group and chaintype
cases <- data.table(read.csv(infile.bplaces.trm.mwmb,header=T))
cases <- cases[, list(N=sum(N)),by=c('TRANSM','mwmb')]

cases_all <- dsubgraphtaxa[, list(N=length(unique(ID))),by=c('TRANSM')]
cases_all[,mwmb:='All']
cases <- merge(cases,cases_all,by=c('TRANSM','mwmb','N'),all=T)
setnames(cases,'TRANSM','trsm')

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

# merge with observed data
dat <- merge(dat,cases,by=c('trsm','mwmb'),all=T)
dat <- merge(dat,bp,by=c('trsm','mwmb'),all=T)

# make labels
dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,Oceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan Africa','S. America &\n Caribbean','Other'))
dat[, xlab:=paste0(mwmb_lab,' \n(N=',N,', N*=',L,')')]

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

order <- unique(subset(dat,select=c('trsm','mwmb_lab','xlab')))
order <- order[order(trsm,mwmb_lab),]
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
										guides(fill = guide_legend(nrow = 1)) +
										theme(legend.position = "bottom"))


g_mwmb_msm <- ggplot(subset(dat,trsm=='MSM'), aes(x=mwmb_lab, y=inf_Ams,fill=mwmb_lab)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)

g_mwmb_hsx <- ggplot(subset(dat,trsm=='HSX'), aes(x=mwmb_lab, y=inf_Ams,fill=mwmb_lab)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside",
				legend.text = element_blank()
				, legend.title = element_blank()) +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='Birthplace') +
	scale_fill_manual(values=col_pal) +
	scale_colour_manual(values=col_pal)

g_st <- ggarrange(subt_msm,subt_hsx,ncol=1,
									 align=c("hv"))
g_inf <- ggarrange(g_mwmb_msm,g_mwmb_hsx,ncol=1,
									align=c("hv"))
g_all <- ggarrange(g_st,g_inf,ncol=2,align=c('hv'))
g_l <- ggarrange(g_all,leg,ncol=1,heights=c(1,0.05))
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX','.png'), g_l, w=12, h=10)


g_mwmb <- ggplot(dat, aes(x=xlab, y=inf_Ams,fill=xlab)) +
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
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX_cases','.png'), g_mwmb, w=20, h=10)

g_st <- ggarrange(subt_msm+theme(legend.position="bottom"),subt_hsx+theme(legend.position="bottom"),ncol=2,align=c("hv"),common.legend=T)
g_all <- ggarrange(g_st,g_mwmb ,ncol=1,align=c('hv'))
ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_birthplace_MSM-HSX','.png'), g_all, w=14, h=15)

## save numbers

tab <- dat[, list(qs= quantile(inf_Ams, prob=ps), qlab=p_labs), by=c('trsm','mwmb_lab','xlab','time')]
ds <- dcast(tab,trsm+mwmb_lab+time+xlab~qlab,value.var='qs')
saveRDS(ds,file=paste0(outfile.base,'-','locally_acquired_bplace_MSM-HSX_cases.RDS'))

tab <- dat[, list(qs= quantile(inf_Ams, prob=ps), qlab=p_labs), by=c('trsm','mwmb_lab','time')]
par <- dcast.data.table(tab, trsm+mwmb_lab+time~qlab, value.var='qs')
par[, L:= paste0( round(M*100, d=0), '% [',  round(CL*100, d=0),'-', round(CU*100, d=0),'%]')]
write.csv(par,file=paste0(outfile.base,'-','local_infections_migrantgroup.csv'))

saveRDS(dat_s,file=paste0(outfile.base,'-','p_local_infections_migrantgroup_samples.RDS'))
saveRDS(bp_s,file=paste0(outfile.base,'-','N_infections_migrantgroup_samples.RDS'))

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
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_N_infections','.png'), g_mwmb_N, w=22, h=10)
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_N_infections_nolab','.pdf'), g_mwmb_N, w=22, h=10)

# Number of locally acquired infections with % in text
dat[trsm=='Amsterdam MSM', all_lab:= 'Amsterdam \nMSM\n']
dat[trsm=='Amsterdam heterosexual', all_lab:= 'Amsterdam \nheterosexual\n']
dat$all_lab <- factor(dat$all_lab,levels=c('Amsterdam \nMSM\n','Amsterdam \nheterosexual\n'))

col_pal = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])
dat[mwmb_lab=='All', flab:='Overall']

col_pal_all <- c(col_pal[1],col_pal[1])
g_all <- ggplot(subset(dat,mwmb_lab=='All'), aes(x=all_lab, y=N_local_M,fill=mwmb_lab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=all_lab,ymin=N_local_CL, ymax=N_local_CU,fill=mwmb_lab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(aes(label=labs), position=position_dodge(width=0.9), vjust=-2, hjust=0,size=8) +
	facet_wrap(. ~ flab,scales="free_x") +
	scale_y_continuous(expand=expansion(mult = c(0, .05)),breaks=seq(0,450,100),labels=seq(0,400,100)) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text( size=24, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Number and proportion of infections in Amsterdam \nin 2014-2018 that were acquired in Amsterdam", x='') +
	scale_fill_manual(values=col_pal_all) #+

col_pal_mwmmb  <- c(pal_npg("nrc")(9)[c(1:5)],pal_npg("nrc")(9)[c(8,6,7,9)])

dat[, trsm2:= paste0(trsm,' stratified by place of birth')]
dat[mwmb_lab=='NL',mwmb_lab:='Dutch-born']
dat[,mwmb_lab:=factor(mwmb_lab,levels=c("All","Dutch-born", "W.Europe,\nN.America,\nOceania","E. & C. \nEurope" ,             
																				"S. America \n& Caribbean","Sub-Saharan \nAfrica",
																				"Other"))]
g_mwmb <- ggplot(subset(dat,mwmb_lab!='All'), aes(x=mwmb_lab, y=N_local_M,fill=mwmb_lab_trsm)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mwmb_lab,ymin=N_local_CL, ymax=N_local_CU,fill=mwmb_lab_trsm),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(aes(label=labs), position=position_dodge(width=0.9), vjust=-2, hjust=0,size=8) +
	facet_wrap(. ~ trsm,scales="free_x") +
	scale_y_continuous(expand=expansion(mult = c(0, .05)),breaks=seq(0,470,100),labels=seq(0,400,100)) +
	coord_cartesian(ylim=c(0,470)) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text( size=24, vjust = 0.5, hjust=0.5),
				axis.text.y=element_blank(),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="", x='') +
	scale_fill_manual(values=col_pal_mwmmb) +
	scale_colour_manual(values=col_pal_mwmmb)

g_mwmb_N <- ggarrange(g_all,g_mwmb,ncol=2,widths=c(0.25,0.8),font.label=list(size=40))
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2018_N_infections_all_final','.png'), g_mwmb_N, w=24, h=10)
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_MSM-HSX_2014-2018_N_infections_all_final','.pdf'), g_mwmb_N, w=24, h=10)


# calculate % of all local infections in each risk group
dat_c <- dat_c[, list(N_local=sum(N_local)),
										by=c('trsm','fb','time','iteration')]
N_local <- copy(dat_c)
tab <- dat_c[, list(N_local= quantile(N_local, prob=ps),
									qlab=p_labs), by=c('trsm','fb','time')]
dat <- dcast(tab,trsm+fb+time~qlab,value.var=c('N_local'))
dat[, N_local_L:=paste0(M,' [',CL,'-',round(CU,0),']')]

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
dat$fb <- factor(dat$fb,levels=c('NL','Foreign-born'),labels=c('Dutch-born','Foreign-born'))
dat <- dat[order(trsm,fb),]

write.csv(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup_agg-foreign-born.csv'))
saveRDS(dat,file=paste0(outfile.base,'-','local_infections_pct_N_migrantgroup_agg-foreign-born.RDS'))

## % of all local infections in FB by trm group
dl <- dcast(N_local,trsm+iteration+time~fb,value.var=c('N_local'))
dl[, pct_fb:=`Foreign-born`/(`Foreign-born`+NL)*100]
tab <- dl[, list(pct_fb= quantile(pct_fb, prob=ps),
										qlab=p_labs), by=c('trsm','time')]
dat <- dcast(tab,trsm+time~qlab,value.var=c('pct_fb'))
dat[, pct_fb_L:=paste0(round(M,0),' [',round(CL,0),'-',round(CU,0),']')]
write.csv(dat,file=paste0(outfile.base,'-','pct_local_infections_in_foreign-born.csv'))

# total local inf
tot <- N_local[, list(N_local=sum(N_local)),
									by=c('time','iteration')]
tot <- tot[, list(N_local= quantile(N_local, prob=ps),
										qlab=p_labs), by=c('time')]
tot <- dcast(tot,time~qlab,value.var=c('N_local'))
tot[, N_local_L:=paste0(M,' [',CL,'-',round(CU,0),']')]
write.csv(dat,file=paste0(outfile.base,'-','total_local_infections.csv'))

## get % of all local infections by risk grp & place of birth
dat_c <- N_local[, list(time=time,trsm=trsm,fb=fb,
											N_local_pct=N_local/sum(N_local)),
							 by=c('iteration')]
tab <- dat_c[, list(N_local_pct= quantile(N_local_pct, prob=ps),
										qlab=p_labs), by=c('trsm','fb','time')]
dat <- dcast(tab,trsm+fb+time~qlab,value.var=c('N_local_pct'))
dat[, N_local_pct_L:=paste0(round(M*100,2),'% [',round(CL*100,2),'-',round(CU*100,2),'%]')]
dat[, N_local_pct_L_round:=paste0(round(M*100,0),'% [',round(CL*100,0),'-',round(CU*100,0),'%]')]
write.csv(dat,file=paste0(outfile.base,'-','pct_total_local_infections_trsm_fb.csv'))

## get % of all local infections in FB MSM
dl <- dcast(N_local,iteration+time~trsm+fb,value.var=c('N_local'))
dl[, pct_fb:=(`MSM_Foreign-born`)/(`HSX_Foreign-born` + `MSM_Foreign-born` + HSX_NL + MSM_NL)*100]
tab <- dl[, list(pct_fb= quantile(pct_fb, prob=ps),
								 qlab=p_labs), by=c('time')]
dat <- dcast(tab,time~qlab,value.var=c('pct_fb'))
dat[, pct_fb_L:=paste0(round(M,0),' [',round(CL,0),'-',round(CU,0),']')]
write.csv(dat,file=paste0(outfile.base,'-','pct_local_infections_in_foreign-born-MSM.csv'))


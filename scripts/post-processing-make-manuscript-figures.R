# post-processing-make-report-figures.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running post-processing-make-report-figures.R\n \n -------------------------------- \n")

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
args_dir[['infdate']] <- 1
args_dir[['period']] <- '2014-2018'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
#args_dir[['job_name']] <- 'rho_infd_pre_infd_em_diagd_infd_1000sg'
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
	
} 

## load functions
source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'post-processing-summary-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

#### MSM

#args_dir[['job_name']] <- 'infd_undiagnosed_invscale'
#args_dir[['period']] <- '2015-2019'
args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFileMSM']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_MSM')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_MSM')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileMSM , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.MSM.p2 <- readRDS(file)
file <- paste0(outfile.base,'-stanout-transmission_pars.RDS')
cat("\n read RDS:", file)
plot.pars.trmspars.MSM.p2 <- readRDS(file)
file <- paste0(outfile.base,'-external_importations_sbt_samples.RDS')
ex.im.MSM.p2 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
local.MSM.p2 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
st.MSM.p2 <- readRDS(file)
file <- paste0(outfile.base,'-','unobs_emergent_chains.RDS')
N_em.MSM.p2 <- readRDS(file)

#args_dir[['job_name']] <- 'rho_infd_pre_infd_em_infd_1000sg'
#args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFileHSX']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_HSX')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_HSX')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileHSX , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.HSX.p2 <- readRDS(file)
file <- paste0(outfile.base,'-stanout-transmission_pars.RDS')
cat("\n read RDS:", file)
plot.pars.trmspars.HSX.p2 <- readRDS(file)
file <- paste0(outfile.base,'-external_importations_sbt_samples.RDS')
ex.im.HSX.p2 <- readRDS(file)
file <- paste0(outfile.base,'-local_infections_samples_TRANSM.RDS')
local.HSX.p2 <- readRDS(file)
file <- paste0(outfile.base,'-subtypes_prop_TRANSM.RDS')
st.HSX.p2 <- readRDS(file)
file <- paste0(outfile.base,'-','unobs_emergent_chains.RDS')
N_em.HSX.p2 <- readRDS(file)

cat(" \n -------------------------------- Order subtypes by number of individuals -------------------------------- \n")

outfile <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')
dsubgraphtaxa <- readRDS(file=outfile)

cat(paste('\n Subset data to individuals diagnosed up until ',2019,' \n'))
# Remove individuals with unknown HIV positive date
#if(args_dir$infdate==1){
if(1){
		dsubgraphtaxa[INF_D>=args_dir$start_d & INF_D<args_dir$end_d, keep:=1]
}else{
	dsubgraphtaxa[HIV1_POS_D>=args_dir$start_d & HIV1_POS_D<args_dir$end_d, keep:=1]
}

#if(args_dir$infdate==1){
if(1){
	dsubgraphsize <- dsubgraphtaxa[, list(icases=length(ID[INF_D<2010 & RECART_D>=2010]),jcases=length(ID[INF_D>=2010 & INF_D<=2019])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
}else{
	dsubgraphsize <- dsubgraphtaxa[, list(icases=length(ID[HIV1_POS_D<2010 & RECART_D>=2010]),jcases=length(ID[HIV1_POS_D>=2010 & HIV1_POS_D<=2019])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
}
cases <- subset(dsubgraphsize,REP=='000',select=c('NAME','SELECT','ST','ST_CLADE','icases','jcases')) # size calculated from adj_icases and j
cases$SIZE <- cases$icases + cases$jcases
freqs <- cases[, list(N=length(NAME)), by=c('SELECT','ST','SIZE')]
freqs <- subset(freqs,SELECT=='Ams')
freqs[, f:= SIZE*N]
freqs <- freqs[, list(N=sum(f)), by=c('ST')]
freqs <- freqs[order(-N),]
st_order <- unique(freqs$ST)


#cat('\n Correct misrecorded dates of diagnosis/ART start \n')
#infile.indinfo <- file.path(args_dir$in_dir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
#dsubgraphtaxa <- correct_misrecorded_dates(dsubgraphtaxa,infile.indinfo)


# summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
#dsubgraphsize <- dsubgraphtaxa[, list(icases=length(ID[HIV1_POS_D<2010 & RECART_D>=2010]),jcases=length(ID[HIV1_POS_D>=2010 & HIV1_POS_D<=2019])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]

dsubgraphsize <- subset(dsubgraphsize,REP=='000' & SELECT!='Ams')
sizes <- dsubgraphsize[, list(cases=sum(jcases)),by=c('ST')]
sizes <- sizes[order(-cases),]
st_order <- unique(sizes$ST)

sizes <- dsubgraphsize[, list(cases=sum(jcases)),by=c('ST','SELECT')]
sizes[SELECT=='AmsHSX',group:='Amsterdam heterosexual']
sizes[SELECT=='AmsMSM',group:='Amsterdam MSM']

cat(" \n -------------------------------- Plot all Rts -------------------------------- \n")

colnames(plot.pars.trmspars.MSM.p2$log_r0_sbts) <- plot.pars.basic.MSM.p2$ds$subtypes_name
colnames(plot.pars.trmspars.HSX.p2$log_r0_sbts) <- plot.pars.basic.HSX.p2$ds$subtypes_name

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

MSM.2 <- as.data.table( reshape2::melt( plot.pars.trmspars.MSM.p2$log_r0_sbts ) )
setnames(MSM.2, 1:3, c('iteration','subtypes','value'))
MSM.2[,trsm:='MSM']
MSM.2[,period:='2015-2019']

HSX.2 <- as.data.table( reshape2::melt( plot.pars.trmspars.HSX.p2$log_r0_sbts ) )
setnames(HSX.2, 1:3, c('iteration','subtypes','value'))
HSX.2[,trsm:='HSX']
HSX.2[,period:='2015-2019']

Rt <- rbind(MSM.2,HSX.2)

Rt[,value:=exp(value)]

Rt[,group:="Amsterdam heterosexual"]
Rt[trsm=='MSM',group:="Amsterdam MSM"]

Rt$subtypes <- factor(Rt$subtypes,levels=c(st_order,'All'))
Rt$group <- factor(Rt$group,levels=c('Amsterdam MSM','Amsterdam heterosexual'))
Rt[,par:='Effective reproduction\n number']


cat(" \n -------------------------------- summarise vmr -------------------------------- \n")

colnames(plot.pars.trmspars.MSM.p2$vmr_minus_one) <- plot.pars.basic.MSM.p2$ds$subtypes_name
colnames(plot.pars.trmspars.HSX.p2$vmr_minus_one) <- plot.pars.basic.HSX.p2$ds$subtypes_name

MSM.2 <- as.data.table( reshape2::melt( plot.pars.trmspars.MSM.p2$vmr_minus_one ) )
#setnames(MSM.2, 1:2, c('iteration','value'))
setnames(MSM.2, 1:3, c('iteration','subtypes','value'))
MSM.2[,trsm:='MSM']
MSM.2[,period:='2014-2019']

HSX.2 <- as.data.table( reshape2::melt( plot.pars.trmspars.HSX.p2$vmr_minus_one ) )
#setnames(HSX.2, 1:2, c('iteration','value'))
setnames(HSX.2, 1:3, c('iteration','subtypes','value'))
HSX.2[,trsm:='HSX']
HSX.2[,period:='2014-2019']

vmr <- rbind(MSM.2,HSX.2)

vmr[,group:="Amsterdam heterosexual"]
vmr[trsm=='MSM',group:="Amsterdam MSM"]
#vmr[, subtypes:='All']
#vmr$subtypes <- factor(vmr$subtypes,levels=c(st_order,'All'))
vmr$subtypes <- factor(vmr$subtypes,levels=c(st_order))
vmr$group <- factor(vmr$group,levels=c('Amsterdam MSM','Amsterdam heterosexual'))
vmr[,par:='Variance-to-mean \nratio']
# add one since samples are for vmr-1
vmr[, value:=value+1]

cat(" \n -------------------------------- local infections -------------------------------- \n")

ex.im.MSM.p2[,period:='2015-2019']
ex.im.MSM.p2 <- merge(ex.im.MSM.p2,plot.pars.basic.MSM.p2$ds,by.x='subtype',by.y='subtypes')
ex.im.HSX.p2[,period:='2015-2019']
ex.im.HSX.p2 <- merge(ex.im.HSX.p2,plot.pars.basic.HSX.p2$ds,by.x='subtype',by.y='subtypes')

do <- rbind(ex.im.MSM.p2,ex.im.HSX.p2)

setnames(do,'subtypes_name.x','subtypes_name')

do[,group:="Amsterdam heterosexual"]
do[TRM=='MSM',group:="Amsterdam MSM"]
do$subtypes_name <- factor(do$subtypes_name,levels=c(st_order,'All'))
do$group <- factor(do$group,levels=c('Amsterdam MSM','Amsterdam heterosexual'))
do[,value:=inf_Ams]
do[,par:='local_inf']
do$subtypes <- factor(do$subtypes_name,levels=c(st_order,'All'))

cat(" \n -------------------------------- combine estimates -------------------------------- \n")

par <- merge(Rt,vmr,by=c('iteration','subtypes','value','trsm','period','group','par'),all=T)
par <- merge(par,subset(do,select=c('iteration','subtypes','group','period','par','value')),by=c('iteration','subtypes','group','period','par','value'),all=T)
par[par=='Effective reproduction\n number',par:='Rt']
par[par=='Variance-to-mean \nratio',par:='vmr']

par <- merge(par,sizes,by.x=c('subtypes','group'),by.y=c('ST','group'),all.x=T)
par[, st_lab:=paste0(subtypes,' (n=',cases,')')]
par$subtypes <- factor(par$subtypes,levels=c(st_order,'All'))
st <- unique(subset(par,select=c('subtypes','st_lab')))
st <- st[order(subtypes),]
par$st_lab <- factor(par$st_lab,levels=st$st_lab)
par$group <- factor(par$group,levels=c('Amsterdam MSM','Amsterdam heterosexual'))

# don't plot subtypes with no observed subgraph data (or less than 10 cases 2010-2019)
par <- par[!(par=='Rt' & cases==0),]
par <- par[!(par=='local_inf' & cases==0),]
par <- par[!(par!='vmr' & cases<10),]
cat(" \n -------------------------------- make figure -------------------------------- \n")

quantiles_95 <- function(x) {
	r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}

g1 <- ggplot(subset(par,par=='Rt'), aes(x=st_lab, y=value,col=period,fill=period)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),colour="black") +
	scale_y_continuous(expand=c(0,0)) +
	coord_cartesian(ylim=c(0,1.1)) +
	geom_hline(yintercept=1) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_blank(),
				#axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_text(size=22),
				strip.background=element_blank(),
				strip.placement = "outside",
				plot.margin=unit(c(0,0,0,0), "cm")) +
	labs(y='Effective reproduction\n number', x=NULL,fill='Time period') +
	facet_wrap(.~group,scales="free") +
	ggsci::scale_fill_npg(alpha=1) +
	ggsci::scale_colour_npg()

g2 <- ggplot(subset(par,par=='vmr'), aes(x=st_lab, y=value,col=period,fill=period)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),colour="black") +
	scale_y_continuous(expand=c(0,0)) +
	#coord_cartesian(ylim=c(0,30)) +
	geom_hline(yintercept=1) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_blank(),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside",
				plot.margin=unit(c(0,0,0,0), "cm")) +
	labs(y='Variance-to-mean \nratio', x=NULL,fill='Time period') +
	#facet_wrap(.~group,scales="free_y") +
	facet_wrap(.~group,scales="free") +
	ggsci::scale_fill_npg(alpha=1) +
	ggsci::scale_colour_npg()

g3 <- ggplot(subset(par,par=='local_inf'), aes(x=st_lab, y=value,col=period,fill=period)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),colour="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="bottom",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside",
				plot.margin=unit(c(0,0,0,0), "cm")) +
	labs(y='Locally acquired \ninfections', x='',fill='Time period') +
	facet_wrap(.~group,scales="free") +
	ggsci::scale_fill_npg(alpha=1) +
	ggsci::scale_colour_npg()

g <- ggarrange(g1+ theme(
								axis.title.x = element_blank(),
								axis.text.x = element_blank(),
								axis.ticks.x = element_blank()),
							g2+ theme(
								axis.title.x = element_blank(),
								axis.text.x = element_blank(),
								axis.ticks.x = element_blank()),
							g3,ncol=1,common.legend=T,legend="bottom",
							 #heights=c(1.5,1.2,1.5),
							 heights=c(1.5,1.5,1.5),
							 labels="auto",font.label = list(size = 20, face = "bold"),
							 align=c("hv"))
ggsave(paste0(outfile.base,'-parameters_sbt_trsm_boxplots.png'), g, w=18, h=17)

cat(" \n -------------------------------- make overall local infections -------------------------------- \n")

st.MSM.p2[, period:='2015-2020']
st.HSX.p2[, period:='2015-2020']
st <- rbind(st.MSM.p2,st.HSX.p2)

st$ST <- factor(st$ST,levels=c("nonB","06cpx","G","02AG","A1","D","C","01AE", "B"))

st$TRANSM <- factor(st$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

subt <- ggplot(subset(st, TRANSM %in% c('Amsterdam MSM','Amsterdam heterosexual')), aes(x=period)) +
	geom_bar(aes(y=p_st,fill=ST), stat='identity',position = "fill",width=0.9) +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="right",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	facet_grid(TRANSM~.) +
	labs(y="Proportion", x='',fill='Subtype') +
	ggsci::scale_fill_npg(alpha=1)

local.MSM.p2[, period:='2015-2020']
local.HSX.p2[, period:='2015-2020']
dat <- rbind(local.MSM.p2,local.HSX.p2)
dat$TRANSM <- factor(dat$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

g_trsm <- ggplot(subset(dat, TRANSM %in% c('Amsterdam MSM','Amsterdam heterosexual')), aes(x=period, y=inf_Ams,fill=period)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 15) +
	theme(axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5),
				legend.position="none",
				strip.text.x=element_blank(),
				strip.background=element_blank(),
				strip.placement = "outside") +
	facet_grid(TRANSM~.) +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	ggsci::scale_color_lancet() +
	ggsci::scale_fill_lancet()

#g <- ggarrange(subt,g_trsm,ncol=2,
#							 align=c("h"),widths=c(1.3,1))
#ggsave(paste0(outfile.base,'-','dist_subtypes_locally_acquired_trsm','.png'), g, w=12, h=8)

g_trsm <- ggplot(subset(dat, TRANSM %in% c('Amsterdam MSM','Amsterdam heterosexual')), aes(x=TRANSM, y=inf_Ams,fill=TRANSM)) +
	stat_summary(fun.data = quantiles_95, geom="boxplot",position = position_dodge2(preserve = "single"),color="black") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text(size=20,angle=15, vjust = 0.5, hjust=0.5),
				legend.position="none",
			  strip.text.x=element_text(size=20),
				strip.background=element_blank(),
				strip.placement = "outside") +
	facet_grid(.~period) +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	ggsci::scale_color_lancet() +
	ggsci::scale_fill_lancet()

#ggsave(paste0(outfile.base,'-','overall_locally_acquired_trsm','.png'), g_trsm, w=12, h=8)

cat(" \n -------------------------------- summarise inputs to % local infections -------------------------------- \n")

do[, ratio:=NT/NI]
#par <- subset(do, subtypes %in% c('B','02AG','A1','01AE','All'))

# infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')
# 
# dsubgraphtaxa <- readRDS(infile.subgraphs)
# 
# if(infdate==1){
# 	dsubgraphtaxa[INF_D>=start_d & INF_D<end_d, keep:=1]
# }else{
# 	dsubgraphtaxa[HIV1_POS_D>=start_d & HIV1_POS_D<end_d, keep:=1]
# }
# dsubgraphtaxa <- subset(dsubgraphtaxa, REP=='000' & SELECT==paste0('Ams',args$trsm) & keep==1)
# 
# dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
# dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
# do <- subset(dsubgraphtaxa,TRANSM==args$trsm)
# do <- do[,list(N=length(FULL_NAME)),by=c('ORIGIN','TRANSM','ST','ST_CLADE')]
# setnames(do,c('TRANSM'),c('TRM_GROUP'))
# 
# do <- subset(do, ORIGIN!='Unknown')
# 
# da <- do[, list(N=sum(N)), by=c('ORIGIN','TRM_GROUP','ST')]
# da <- merge(da, da[, list(TOTAL=sum(N)), by=c('TRM_GROUP','ST')], by=c('TRM_GROUP','ST'))
# tmp <- dcast.data.table(da,ORIGIN~ST,value.var='N')
# 
# loc <- data.table(location=seq(1:length(unique(tmp$ORIGIN))),
# 									loc_label=unique(tmp$ORIGIN))
N_em.HSX.p2[, period:='2015-2019']
N_em.MSM.p2[, period:='2015-2019']
N_em.HSX.p2[, TRM:='HSX']
N_em.MSM.p2[, TRM:='MSM']
N_em.HSX.p2 <- merge(N_em.HSX.p2,plot.pars.basic.HSX.p2$ds,by.x='subtype',by.y='subtypes')
N_em.MSM.p2 <- merge(N_em.MSM.p2,plot.pars.basic.MSM.p2$ds,by.x='subtype',by.y='subtypes')

em_c <- rbind(N_em.HSX.p2,N_em.MSM.p2)
em_c[, NT:=NULL]

do <- merge(do,em_c,by=c('TRM','period','subtypes_name','iteration'))

table1 <- do[, list(Ams_origin=quantile(1-p,prob=ps,na.rm=T),
										ext_origin=quantile(p,prob=ps,na.rm=T),
										N_obs_em_chains=quantile(N_o,prob=ps,na.rm=T),
										N_unobs_em_chains=quantile(N_u,prob=ps,na.rm=T),
										N_em_chains=quantile(NT,prob=ps,na.rm=T),
										N_cases=quantile(NI,prob=ps,na.rm=T),
										ratio=quantile(NT/NI,prob=ps,na.rm=T),
										ext_imp=quantile(p*(NT/NI),prob=ps,na.rm=T),
										inf_Ams=quantile(inf_Ams,prob=ps,na.rm=T),
										q_label=p_labs),
						 by=c('TRM','period','subtypes_name')]

ans <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='Ams_origin')
ans[, Ams_origin:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
set(ans,NULL,c('CL','CU','M'),NULL)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='ext_origin')
ans2[, ext_origin:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','ext_origin')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='N_obs_em_chains')
ans2[, N_obs_em_chains:= paste0( round(M, d=1), ' [',  round(CL, d=1),'-', round(CU, d=1),']')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','N_obs_em_chains')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='N_unobs_em_chains')
ans2[, N_unobs_em_chains:= paste0( round(M, d=1), ' [',  round(CL, d=1),'-', round(CU, d=1),']')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','N_unobs_em_chains')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='N_em_chains')
ans2[, N_em_chains:= paste0( round(M, d=1), ' [',  round(CL, d=1),'-', round(CU, d=1),']')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','N_em_chains')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='N_cases')
ans2[, N_cases:= paste0( round(M, d=1), ' [',  round(CL, d=1),'-', round(CU, d=1),']')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','N_cases')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='ratio')
ans2[, ratio:= paste0( round(M, d=2), ' [',  round(CL, d=2),'-', round(CU, d=2),']')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','ratio')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='ext_imp')
ans2[, ext_imp:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','ext_imp')),by=c('TRM','period','subtypes_name'),all=T)

ans2 <- dcast.data.table(table1, TRM+period+subtypes_name~q_label, value.var='inf_Ams')
ans2[, inf_Ams:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,subset(ans2,select=c('TRM','period','subtypes_name','inf_Ams')),by=c('TRM','period','subtypes_name'),all=T)

ans[subtypes_name=='nonB',subtypes_name:='Non-B']
ans$TRM <- factor(ans$TRM,levels=c('MSM','HSX'),labels=c('Amsterdam \nMSM','Amsterdam \nhetersexual'))

saveRDS(ans,file=paste0(outfile.base,'-','local_inf_inputs.RDS'))
write.csv(ans,paste0(outfile.base,'-','local_inf_inputs.csv'))

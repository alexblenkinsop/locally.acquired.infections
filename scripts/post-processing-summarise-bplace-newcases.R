
cat(" \n -------------------------------- \n \n Running post-processing-summarise-bplace-newcases.R\n \n -------------------------------- \n")

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

#### MSM
args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileMSM']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_MSM')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_MSM')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileMSM , "-", args_dir$job_tag)

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

#### HSX
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


cat(" \n -------------------------------- define functions -------------------------------- \n")

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

quantiles_95 <- function(x) {
	r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
	names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
	r
}


cat(" \n -------------------------------- summarise observed birth places -------------------------------- \n")

infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')

dsubgraphtaxa <- readRDS(infile.subgraphs)

dsubgraphtaxa[,SG:= FULL_NAME]
dsubgraphtaxa$SG <- as.factor(dsubgraphtaxa$SG)
dsubgraphtaxa[!is.na(ST_CLADE),SG:= paste0(SG,'_',ST_CLADE)]
# get min/max date of diagnosis date in subgraphs
mind <- subset(dsubgraphtaxa,REP=='000')[,list(MIND=min(INF_D,na.rm=T),MAXD=max(INF_D,na.rm=T)),by=c('TRANSM','SG')]
dsubgraphtaxa <- merge(dsubgraphtaxa,mind,by=c('TRANSM','SG'))
start_y <- as.integer(substr(args_dir$period,1,4))

# count chains which started before 2015
dsubgraphtaxa$chaintype <- "pre-existing"
dsubgraphtaxa$chaintype[dsubgraphtaxa$MIND>=start_y] <- "emergent"

if(args_dir$infdate==1){
	dsubgraphtaxa[INF_D>=args_dir$start_d & INF_D<args_dir$end_d, keep:=1]
}else{
	dsubgraphtaxa[HIV1_POS_D>=args_dir$start_d & HIV1_POS_D<args_dir$end_d, keep:=1]
}
dsubgraphtaxa <- subset(dsubgraphtaxa, REP=='000' & SELECT %in% c('AmsHSX','AmsMSM') & keep==1)


geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']

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

cases <- dsubgraphtaxa[, list(N=length(unique(ID))),by=c('TRANSM','chaintype','mwmb')]
cases_all <- dsubgraphtaxa[, list(N=length(unique(ID))),by=c('TRANSM','chaintype')]
cases_all[,mwmb:='All']
cases <- merge(cases,cases_all,by=c('TRANSM','chaintype','mwmb','N'),all=T)
setnames(cases,'TRANSM','trsm')
cases$chaintype <- factor(cases$chaintype,levels=c('pre-existing','emergent'))
cases <- dcast.data.table(cases, trsm+mwmb~chaintype, value.var='N')
cases[, Total:=`pre-existing` + emergent]
cases[, p_pre:=paste0(round((`pre-existing`/Total)*100,1),"%")]
cases[, p_em:=paste0(round((emergent/Total)*100,1),"%")]
setcolorder(cases, c('trsm' , 'mwmb' ,'Total','pre-existing','p_pre', 'emergent','p_em' ))
setnames(cases,c('Total','pre-existing','p_pre', 'emergent','p_em'),
				 paste0(c('Total','pre-existing','p_pre', 'emergent','p_em'),"_o"))

cat(" \n -------------------------------- summarise predicted birth places -------------------------------- \n")
## load birth places from GQs
bp.MSM <- as.data.table( reshape2::melt( fit.bplace.MSM ) )
setnames(bp.MSM, 1:6, c('iteration','subtype','bplace','index_cases','N','chaintype'))
bp.MSM[, trsm:='MSM']
locs <- unique(dsubgraphtaxa$mwmb[dsubgraphtaxa$TRANSM=='MSM'])
labs <- data.table(mwmb=sort(locs[locs!='All']),bplace=seq(1:length(locs[locs!='All'])))
bp.MSM <- merge(bp.MSM, labs,by='bplace',all=T)

bp.HSX <- as.data.table( reshape2::melt( fit.bplace.HSX ) )
setnames(bp.HSX, 1:6, c('iteration','subtype','bplace','index_cases','N','chaintype'))
bp.HSX[, trsm:='HSX']
locs <- unique(dsubgraphtaxa$mwmb[dsubgraphtaxa$TRANSM=='HSX'])
labs <- data.table(mwmb=sort(locs[locs!='All']),bplace=seq(1:length(locs[locs!='All'])))
bp.HSX <- merge(bp.HSX, labs,by='bplace',all=T)

bp <- rbind(bp.MSM,bp.HSX)
bp <- subset(bp,is.finite(N))
bp <- dcast(bp,trsm+iteration+subtype+index_cases+bplace+mwmb~chaintype,value.var='N')
bp[, Total:= origins_ind_x + origins_ind_e]
setnames(bp,c('origins_ind_x','origins_ind_e'),c('pre-existing','emergent'))

# sum over index cases and subtypes, report by migrant group
tmp <- bp[, list(`pre-existing`=sum(`pre-existing`),
									emergent=sum(emergent),
								  Total=sum(Total)),by=c('trsm','iteration','bplace','mwmb')]

tmp[, p_pre:=`pre-existing`/Total]
tmp[, p_em:=emergent/Total]

# summarise all MSM/HSX
tmp2 <- bp[, list(`pre-existing`=sum(`pre-existing`),
									emergent=sum(emergent),
									Total=sum(Total)),by=c('trsm','iteration')]

tmp2[, p_pre:=`pre-existing`/Total]
tmp2[, p_em:=emergent/Total]
tmp3 <- tmp2
tmp3[,bplace:= 0]
tmp3[,mwmb:= 'All']
tmp3 <- merge(tmp,tmp3,by=c('trsm','iteration','bplace','mwmb','pre-existing','p_pre','emergent','p_em','Total'),all=T)
saveRDS(tmp3,file=paste0(outfile.base,'-','birthplaces_cases_samples.RDS'))

tmp <- tmp[, list(`pre-existing` = quantile(`pre-existing`, prob=ps,na.rm=T),
									p_pre = quantile(p_pre, prob=ps,na.rm=T),
									emergent = quantile(emergent, prob=ps,na.rm=T),
									p_em = quantile(p_em, prob=ps,na.rm=T),
									Total = quantile(Total, prob=ps,na.rm=T),
									q_label=p_labs),by=c('trsm','bplace','mwmb')]		

tmp2 <- tmp2[, list(mwmb='All',
									`pre-existing` = quantile(`pre-existing`, prob=ps,na.rm=T),
										p_pre = quantile(p_pre, prob=ps,na.rm=T),
										emergent = quantile(emergent, prob=ps,na.rm=T),
										p_em = quantile(p_em, prob=ps,na.rm=T),
										Total = quantile(Total, prob=ps,na.rm=T),
									  q_label=p_labs),by=c('trsm')]		

bp <- merge(tmp,tmp2,by=c('trsm','mwmb','pre-existing','p_pre','emergent','p_em','Total','q_label'),all=T)
#setnames(bp,'N','N_star')
tmp <- copy(bp) # copy for figure

bp <- dcast.data.table(bp, trsm+mwmb~q_label, value.var=c('pre-existing','p_pre','emergent','p_em','Total'))

bp[, Total:= paste0( round(Total_M, d=0), ' [',  round(Total_CL, d=0),'-', round(Total_CU, d=0),']')]
bp[, `pre-existing`:= paste0( round(`pre-existing_M`, d=0), ' [',  round(`pre-existing_CL`, d=0),'-', round(`pre-existing_CU`, d=0),']')]
bp[, p_pre:= paste0( round(p_pre_M*100, d=1), '% [',  round(p_pre_CL*100, d=1),'-', round(p_pre_CU*100, d=1),'%]')]
bp[, emergent:= paste0( round(emergent_M, d=0), ' [',  round(emergent_CL, d=0),'-', round(emergent_CU, d=0),']')]
bp[, p_em:= paste0( round(p_em_M*100, d=1), '% [',  round(p_em_CL*100, d=1),'-', round(p_em_CU*100, d=1),'%]')]

bp <- subset(bp,select=c('trsm','mwmb','Total','pre-existing','p_pre','emergent','p_em'))
bp <- merge(cases,bp,by=c('trsm','mwmb'),all=T)

write.csv(bp,file=paste0(outfile.base,'-','birthplaces_cases.csv'))


### plot
tmp <- melt(subset(tmp,select=c('trsm','mwmb','pre-existing','emergent','q_label')), value.var=c('pre-existing','emergent'))
tmp <- dcast.data.table(tmp, trsm+mwmb+variable~q_label, value.var=c('value'))
tmp[, analysis:='predicted']

cases <- subset(cases,select=c('trsm','mwmb','pre-existing_o','emergent_o'))
#cases <- melt(subset(cases,select=c('trsm','mwmb','pre-existing_o','emergent_o')), value.var=c('pre-existing_o','emergent_o'))
#cases <- dcast(cases,trsm+mwmb~variable,value.var='value')
cases <- melt(subset(cases,select=c('trsm','mwmb','pre-existing_o','emergent_o')), value.var=c('pre-existing_o','emergent_o'))

cases[variable=='pre-existing_o',variable:='pre-existing']
cases[variable=='emergent_o',variable:='emergent']
cases[, analysis:='observed']
setnames(cases,'value','M')

bp <- merge(cases,tmp,by=c('trsm','mwmb','analysis','variable','M'),all=T)

bp[, analysis:=factor(analysis,levels=c('observed','predicted'))]

bp[mwmb=='G1' & trsm=='MSM',mgid:=1]
bp[mwmb=='G2' & trsm=='MSM',mgid:=2]
bp[mwmb=='G3' & trsm=='MSM',mgid:=3]
bp[mwmb=='NL' & trsm=='MSM',mgid:=4]
bp[mwmb=='Other' & trsm=='MSM',mgid:=5]
bp[mwmb=='G4' & trsm=='HSX',mgid:=1]
bp[mwmb=='G5'& trsm=='HSX',mgid:=2]
bp[mwmb=='NL'& trsm=='HSX',mgid:=3]
bp[mwmb=='Other'& trsm=='HSX',mgid:=4]

## relabel migrant groups
bp[mwmb=='All' & trsm=='MSM', mlab:='All']
bp[mwmb=='G1' & trsm=='MSM', mlab:='W.Europe,\nN.America,\nOceania']
bp[mwmb=='G2' & trsm=='MSM', mlab:='E. & C. Europe']
bp[mwmb=='G3' & trsm=='MSM', mlab:='S. America &\n Caribbean']
bp[mwmb=='NL' & trsm=='MSM', mlab:='NL']
bp[mwmb=='Other' & trsm=='MSM', mlab:='Other']

bp[mwmb=='All' & trsm=='HSX', mlab:='All']
bp[mwmb=='G4' & trsm=='HSX', mlab:='Sub-Saharan\nAfrica']
bp[mwmb=='G5' & trsm=='HSX', mlab:='S. America &\n Caribbean']
bp[mwmb=='NL' & trsm=='HSX', mlab:='NL']
bp[mwmb=='Other' & trsm=='HSX', mlab:='Other']

bp$mlab <- factor(bp$mlab,levels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe',
																	 'S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'),
									labels=c('All','NL',' W.Europe,\nN.America,Oceania','E. & C. Europe',
													 'S. America &\n Caribbean','Sub-Saharan Africa','Other'))
bp[, mwmb_lab_trsm:=paste0(trsm,' ',mlab)]
order <- unique(subset(bp,select=c('trsm','mlab','mwmb_lab_trsm')))
order <- order[order(trsm,mlab),]
bp$mwmb_lab_trsm <- factor(bp$mwmb_lab_trsm,levels=unique(order$mwmb_lab_trsm))

#do$mlab <- factor(do$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))
bp$trsm <- factor(bp$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

# make plot of where cases were from in pre-existing/emergent subgraphs
col_pal_c = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

plot <- ggplot(data=bp) +
	geom_bar(aes(x=mlab,y=M,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=CL, ymax=CU,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_wrap(trsm~variable,scales="free") +
	labs(x='\nPlace of birth',y="Number of cases between 2014-2018",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				axis.text.x=element_text(angle=40, vjust = 0.5, hjust=0.5)) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-observed_predicted_cases_origins-2014_2018.png'),plot,w=15, h=15)


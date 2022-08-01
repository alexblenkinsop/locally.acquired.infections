require(scales)
require(ggsci)
require(tidyr)
require(reshape)
require(patchwork)
require(ggpubr)

## set up ----
if(0){
	# setup
	home <- '/Users/alexb/Box Sync/Roadmap'
	home <- file.path(args$in_dir,args$analysis,'Data')
	job_tag <- "undiag_untilmay2019_weights"
	#job_tag <- "2010-2012_inf_NL_sens_midpoint_SC"
	sens <- F # whether to run sensitivity analysis
}

args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	dir.ecdc='/rds/general/project/ratmann_roadmap_data_analysis/live/Data/Undiagnosed/ECDC_model',
	#source_dir= '~/Documents/GitHub/bpm',
	#indir='~/Box Sync/Roadmap',
	#outdir= '~/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/hierarchical_model',
	dir.ecdc='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/ECDC_model',
	stanModelFile= 'undiagnosed_211102',
	analysis= 'analysis_211101',
	#analysis= 'analysis_200917',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	#job_tag= 'undiag_untilmay2019',
	job_tag= 'undiag_untilmay2019_weights',
	#job_tag= '2010-2012_inf_NL_sens_midpoint_SC',
	#job_tag= 'sens_q40',
	#sens='mp',
	sens=F,
	weights='ECDC'
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-analysis')
	stopifnot(args_line[[7]]=='-seed', !is.na(as.integer(args_line[[8]])))	
	stopifnot(args_line[[9]]=='-indir')
	stopifnot(args_line[[11]]=='-outdir')
	stopifnot(args_line[[13]]=='-jobtag')
	stopifnot(args_line[[15]]=='-sens')
	stopifnot(args_line[[17]]=='-dir.ecdc')
	stopifnot(args_line[[19]]=='-weights')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['seed']] <- as.integer(args_line[[8]])
	args[['indir']] <- args_line[[10]]
	args[['outdir']] <- args_line[[12]]
	args[['job_tag']] <- args_line[[14]]
	args[['sens']] <- args_line[[16]]
	args[['dir.ecdc']] <- args_line[[18]]
	args[['weights']] <- args_line[[18]]
} 

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag)) 
#home <- file.path(args$indir,args$analysis,'Data')
#outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outpath <- file.path(args$outdir)
#outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model',job_tag)
outdir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag))

args$trsm <- 'MSM'
job_tag <- args$job_tag

## make stan data ----
cat(" \n -------------------------------- \n Load data \n -------------------------------- \n")

# load data
file.seqlabels <- file.path(args$indir,args$analysis,'misc/200917_sequence_labels.rda')
infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','roadmap_cd4_vl_est.csv')
if(args$sens!=F){
	infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','roadmap_cd4_vl_est-quantiles.csv')
}
geo.file <- file.path(args$indir,'misc/NEWGEO.csv')

infile.ecdc.msm <- file.path(args$dir.ecdc,'AMS_MSM_LOCAL_MR_Result_main.csv')
infile.ecdc.nonmsm <- file.path(args$dir.ecdc,'AMS_NO_MSM_LOCAL_MR_Result_main.csv')

load(file.seqlabels)
dinf <- read.csv(infile.inftime,header=T)
geo <- data.table(read.csv(geo.file))

# load ECDC estimates
dw <- data.table(read.csv(infile.ecdc.msm))
dw[, trsm:='MSM']
tmp <- data.table(read.csv(infile.ecdc.nonmsm))
tmp[, trsm:='HSX']
dw <- rbind(dw,tmp)

cat(" \n -------------------------------- \n Define migrant groups \n -------------------------------- \n")

geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

dind <- data.table(dind)
if(args$sens==F){
	dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
}else if(args$sens=='mp'){
	dinf <- subset(dinf,select=c('id','u','estsctodiagMedian','hiv_pos_d','hiv1_neg_d'))
}else{
	dinf <- subset(dinf,select=c('id','estsctodiagMedian','estsctodiag30','estsctodiag40','hiv_pos_d'))
}
dinf <- unique(dinf)
if(args$sens=='mp'){
	dinf <- subset(dinf,!is.na(hiv1_neg_d)) # just keep patients with a last negative test for sensitivity analysis
}
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D','INF_CNTRY','MIG_D')),by.x='id',by.y='PATIENT',all.x=T)
do <- data.table(dinf)
if(args$sens==30){
	do[, time:=estsctodiag30]
}else if(args$sens==40){
	do[, time:=estsctodiag40]
}else if(args$sens=='mp'){
	do[, estsctodiagMedian:=u/2] # use midpoint between time at risk (last neg test and first pos test)	
	do[, time:=estsctodiagMedian]
}else{
	do[, time:=estsctodiagMedian]
}

do[, INF_D:=HIV1_POS_D - time]
do <- merge(do,subset(dind,select=c('PATIENT','ORIGIN')),by.x='id', by.y='PATIENT',all.x=T)
do <- merge(do,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)

# reclassify by main migrant groups
## msm
do[TRANSM=='MSM', mwmb:="Other"]
do[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
# western countires (non-NL)
do[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania") & ORIGIN!='NL', mwmb:="G1"]
# eastern and central europe
do[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="G2"]
# caribbean and south america
do[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]

## hsx
do[TRANSM=='HSX', mwmb:="Other"]
do[TRANSM=='HSX' & ORIGIN %in% c("NL"), mwmb:="NL"]
# sub-saharan africa
do[TRANSM=='HSX' & WRLD_born %in% c("Africa"), mwmb:="G4"]
# caribbean and south america
do[TRANSM=='HSX' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]

# exclude individuals without an infection time estimate
do <- subset(do,!is.na(time))

# summarise number diagnosed to adjust for undiagnosed
do[, infdate:=INF_D]
do[is.na(INF_D) & HIV1_POS_D>=2014, infdate:=HIV1_POS_D] ## in model we use HIV pos date where there is no infection date - check this
n_diag <- do[, list(N_diag=length(unique(id[infdate>=2014]))),by=c('TRANSM','mwmb')]
n_diag <- n_diag[order(mwmb),]

da <- subset(do,do$INF_D>=2010 & do$INF_D<2013 & TRANSM %in% c('HSX','MSM'))

cat(" \n -------------------------------- \n Make summary table of characteristics \n -------------------------------- \n")

tab <- data.table(var1='TRANSM',var2='TOTAL',table(da$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(da$mwmb,da$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_all <- copy(tab)

# make synthetic subset of data
de <- subset(do,do$INF_D>=2010 & do$INF_D<2013 & TRANSM %in% c('HSX','MSM'))
dexcl <- da

tab <- data.table(var1='TRANSM',var2='TOTAL',table(dexcl$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(dexcl$mwmb,dexcl$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]
tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_excl <- copy(tab)

tab <- data.table(var1='TRANSM',var2='TOTAL',table(de$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(de$mwmb,de$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

time <- de[, list(var='INFTIME',
									V1='Estimated time to diagnosis (years)',
									q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
													 " [",
													 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
													 "-",
													 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
													 "]")),
					 by=c('TRANSM','mwmb')]
time_all <- de[, list(mwmb='TOTAL',
											var='INFTIME',
											V1='Estimated time to diagnosis (years)',
											q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
															 " [",
															 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
															 "-",
															 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
															 "]")),
							 by=c('TRANSM')]
time <- rbind(time,time_all)
time <- dcast(time,TRANSM+mwmb~V1,value.var='q')
tab[var2=='TOTAL',V1:=var2]
tab <- merge(tab,time,by.y=c('TRANSM','mwmb'),by.x=c('V2','V1'),all=T)
tab <- subset(tab,!(N==0 & is.na(`Estimated time to diagnosis (years)`)))
tab_incl <- copy(tab)

setnames(tab_all,'N','N_all')
setnames(tab_excl,'N','N_excluded')
setnames(tab_incl,'N','N_included')
tab <- merge(tab_all,tab_excl,by=c('var1','var2','V1','V2'),all=T)
tab <- merge(tab,tab_incl,by=c('var1','var2','V1','V2'),all=T)

tab[V1=='G1', V1:='W.Europe, N.America, Oceania']
tab[V1=='G2', V1:='E. & C. Europe']
tab[V1=='G3', V1:='S. America & Caribbean']
tab[V1=='G4', V1:='Sub-Saharan Africa']
tab[V1=='G5', V1:='S. America & Caribbean']
tab[V1=='NL', V1:='Netherlands']
tab[V1=='Other', V1:='Other']
tab[V1=='TOTAL', V1:='All']

tab$V2 <- factor(tab$V2,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
tab <- tab[order(V2),]
tab <- subset(tab,is.na(N_excluded)) # drop the pts infected outside NL
saveRDS(tab,file=file.path(outdir, paste0("characteristics_patients_undiagnosed.RDS")))


## summarise empirical quantiles of time to diagnosis for table 1 ----

ps <- c(0.5,0.025,0.975)
p_labs <- c('M','CL','CU')

time_m_ci <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
											 q_label=p_labs),by=c('TRANSM','mwmb')]
time_m_ci <- subset(time_m_ci,TRANSM %in% c('MSM','HSX'))
time_m_ci <- dcast(time_m_ci,TRANSM+mwmb~q_label,value.var='q')
time_m_ci <- subset(time_m_ci,select=c('TRANSM','mwmb','M','CL','CU'))

time_m_ci_all <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
													 q_label=p_labs),by=c('TRANSM')]
time_m_ci_all <- subset(time_m_ci_all,TRANSM %in% c('MSM','HSX'))
time_m_ci_all <- dcast(time_m_ci_all,TRANSM~q_label,value.var='q')
time_m_ci_all <- subset(time_m_ci_all,select=c('TRANSM','M','CL','CU'))

time_m_ci_total <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
														 q_label=p_labs)]
time_m_ci_total <- dcast(time_m_ci_total,.~q_label,value.var="q")

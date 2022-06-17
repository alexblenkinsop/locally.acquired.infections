
cat(" \n -------------------------------- \n \n estimate-total-infections.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['trsm']] <- 'HSX'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'undiagnosed_weighted_inf_rate'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['infdate']] <- 1
	args_dir[['start_y']] <- 2014
	args_dir[['p_undiag_gps']]='G1,G2,G3,NL,Other'
	args_dir[['p_undiag_probs']]='0.1139,0.1729,0.1812,0.13880,0.2268'
	#args_dir[['undiag_job']]= 'undiag_untilmay2019'
	#args_dir[['undiag_job']]= 'sens_q30'
	args_dir[['undiag_job']]= 'weights_diagrates'
}

# save args for report before loading those from running session 
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFile')	
	stopifnot(args_line[[3]]=='-analysis')
	stopifnot(args_line[[5]]=='-in_dir')
	stopifnot(args_line[[7]]=='-out_dir')
	stopifnot(args_line[[9]]=='-job_tag')
	stopifnot(args_line[[11]]=='-trsm')
	stopifnot(args_line[[13]]=='-source_dir')
	stopifnot(args_line[[15]]=='-start_y')
	stopifnot(args_line[[17]]=='-undiag_job')
	
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['analysis']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['out_dir']] <- args_line[[8]]
	args_dir[['job_tag']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
	args_dir[['source_dir']] <- args_line[[14]]
	args_dir[['start_y']] <- args_line[[16]]
	args_dir[['undiag_job']] <- args_line[[18]]
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)

#	load all input variables for this analysis run
do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes
start_d = args$start_d
end_d = args$end_d

p_undiagnosed <- data.table(mwmb=args$p_undiag_gps,du=args$p_undiag_probs)

dsubgraphtaxa <- readRDS(file.path(args_dir[['out_dir']],'subgraphs_withmetadata.RDS'))

infile.meta <- file.path(args$indir, args$analysis, 'misc', '200917_sequence_labels.rda')
infile.seq <-	file.path(args$indir, 'Data', 'data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
infile.cd4 <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_CD4.csv')
infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')

cat("\nRead patient metadata \n")
load(infile.seq)
load(infile.meta)
setnames(dind, gsub("CENS_D", "RECART_D", names(dind)))

dind_diag <- copy(dind)
dind_diag$SEQ <- dind_diag$PATIENT %in% ds$PATIENT
dind_diag <- merge(dind_diag,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)
dind_diag <- data.table(dind_diag)

cat('\n Add birth place region \n')
dind <- data.table(dind)
dind <- add_georeg_bplace_ind_data(dind,infile.geo)
dind <- map_mwmb_regions(dind)

cat("\nLoad time to diagnosis estimates \n")
dt_m <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_MSM.RDS'))
dt_h <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_HSX.RDS'))
dt_m <- readRDS(dt_m)
dt_h <- readRDS(dt_h)
dt_inf <- rbind(dt_m,dt_h)
setnames(dt_inf,'trsm','TRANSM')

cat("\nUpdate col names if using infection date \n")
if(args$infdate==1){
	infile.inftime <- file.path(args$indir,'Data','infection_time_estimates','roadmap_cd4_vl_est.csv')
	#dsubgraphtaxa <- add_infection_time(dsubgraphtaxa,infile.inftime)
	dind <- data.table(dind)
	dind[, ID:=PATIENT]
	dind <- add_infection_time(dind,infile.inftime,dt_inf)
	dsubgraphtaxa[, HIV1_POS_D:=INF_D]
	dind[, HIV1_POS_D:=INF_D]
}

# flag which diagnosed individuals have a sequence
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
dind <- merge(dind,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)
#ams <- dind[dind$CITY=='Amsterdam' & dind$HIV1_POS_D>=2010,]
nl <- dind[dind$HIV1_POS_D>=args_dir$start_y,]
ams <- dind[dind$CITY=='Amsterdam' & dind$HIV1_POS_D>=args_dir$start_y & dind$HIV1_POS_D<args$end_d,]

## add CD4 counts at diagnosis
ams <- add_cd4_counts(ams,infile.cd4,args_dir$start_y,NA)
dind_diag[, ID:=PATIENT]
dind_diag <- add_cd4_counts(dind_diag,infile.cd4,args_dir$start_y,NA)
dind[, ID:=PATIENT]
dind <- add_cd4_counts(dind,infile.cd4,args_dir$start_y,NA)

## add georegions/migrant groups

geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

#ams <- merge(ams,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
dind_diag <- merge(dind_diag,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
ams[BIRTH_CNTRY=='Tajikistan',WRLD_born:='Asia']
ams[WRLD_born=='PacificOceania',WRLD_born:='Asia']
ams[WRLD_born=='Oceania',WRLD_born:='Australia & New Zealand']
ams[BIRTH_CNTRY=='Netherlands',WRLD_born:='Netherlands']

#ams <- map_mwmb_regions(ams)
dind_diag <- map_mwmb_regions(dind_diag)
cat(" \n -------------------------------- load prob(undiagnosed) -------------------------------- \n")

msm_und <- readRDS(file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args_dir$undiag_job),paste0('p_undiagnosed_samples_',args_dir$undiag_job,'_MSM.RDS')))
hsx_und <- readRDS(file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args_dir$undiag_job),paste0('p_undiagnosed_samples_',args_dir$undiag_job,'_HSX.RDS')))
dd <- rbind(msm_und,hsx_und)

summ <- data.table(
	N_diagnoses_MSM_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM'])),
	N_sequences_MSM_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_MSM_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_MSM_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_NL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN!='NL']))),
	N_diagnoses_G1_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G1']))),
	N_diagnoses_G2_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G2']))),
	N_diagnoses_G3_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G3']))),
	N_diagnoses_Oth_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='Other']))),
	N_diagnoses_NL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL'])),
	N_diagnoses_G1_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G1'])),
	N_diagnoses_G2_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G2'])),
	N_diagnoses_G3_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G3'])),
	N_diagnoses_Oth_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='Other'])),
	N_diagnoses_HSX_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX'])),
	N_sequences_HSX_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_HSX_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_HSX_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN!='NL'])),
	N_diagnoses_G4_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G4'])),
	N_diagnoses_G5_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G5'])),
	N_diagnoses_Oth_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='Other'])),
	N_diagnoses_NL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN!='NL'])))
)

diag <- reshape2::melt(summ)
setnames(diag,'value','diag')
diag <- subset(diag,variable %in% c('N_diagnoses_NL_MSM','N_diagnoses_G1_MSM','N_diagnoses_G2_MSM','N_diagnoses_G3_MSM',
																		'N_diagnoses_Oth_MSM','N_diagnoses_NL_HSX','N_diagnoses_G4_HSX','N_diagnoses_G5_HSX',
																		'N_diagnoses_Oth_HSX'))
diag[grepl('_MSM',variable),trsm:='MSM']
diag[grepl('_HSX',variable),trsm:='HSX']
diag[grepl('_NL_',variable) & trsm=='MSM',mg:=1]
diag$mg <- factor(diag$variable,levels=c('N_diagnoses_NL_MSM','N_diagnoses_G1_MSM','N_diagnoses_G2_MSM','N_diagnoses_G3_MSM','N_diagnoses_Oth_MSM',
																				 'N_diagnoses_NL_HSX','N_diagnoses_G4_HSX','N_diagnoses_G5_HSX','N_diagnoses_Oth_HSX'),
									labels=c(4,1,2,3,5,3,1,2,4))
diag$mg <- as.integer(as.character(diag$mg))
dd <- merge(dd,diag,by=c('trsm','mg'),all=T)
dd[, N_inf:=round(diag/(1-av_undiagnosed),digits=0)]
dd[, mg_fb:=1]
dd[trsm=='MSM' & mg==4, mg_fb:=0]
dd[trsm=='HSX' & mg==3, mg_fb:=0]
dd[, tot:='All']
dd_fb <- dd[, list(diag=sum(diag),
									 N_inf=sum(N_inf)),
						by=c('trsm','mg','iter')]
dd_all <- dd[, list(diag=sum(diag),
										N_inf=sum(N_inf)),
						 by=c('trsm','tot','iter')]
setnames(dd_all,'tot','mg')
dat <- rbind(dd_fb,dd_all)
dat[, undiag:=1-(diag/N_inf)]
dat <- dat[, list(undiag=quantile(undiag,prob=c(0.025,0.5,0.975)),
									N_diag=quantile(diag,prob=c(0.025,0.5,0.975)),
									N_inf=quantile(N_inf,prob=c(0.025,0.5,0.975)),
									#undiag=quantile(undiag,prob=c(0.025,0.5,0.975)),
									qlabel=c('L','M','U')),
					 by=c('trsm','mg')] # summarise quantiles for Jan of each year

dat <- dcast(dat,trsm+mg~qlabel,value.var=c('undiag','N_diag','N_inf'))

dat[, undiag:= paste0(round(undiag_M*100,0),"% [",round(undiag_L*100,0),"-",round(undiag_U*100,0),"%]")]
dat[, N_inf:= paste0(N_inf_M," [",N_inf_L,"-",N_inf_U,"]")]
write.csv(dat,file=paste0(outfile.base,'-','number_infections.csv'))

dat <- rbind(dd_fb,dd_all)
dat <- dat[, list(diag=sum(diag),
									 N_inf=sum(N_inf)),
						by=c('mg','iter')]
dat[, undiag:=1-(diag/N_inf)]
dat <- dat[, list(undiag=quantile(undiag,prob=c(0.025,0.5,0.975)),
									N_diag=quantile(diag,prob=c(0.025,0.5,0.975)),
									N_inf=quantile(N_inf,prob=c(0.025,0.5,0.975)),
									#undiag=quantile(undiag,prob=c(0.025,0.5,0.975)),
									qlabel=c('L','M','U')),
					 by=c('mg')] # summarise quantiles for Jan of each year

dat <- dcast(dat,mg~qlabel,value.var=c('undiag','N_diag','N_inf'))

dat[, undiag:= paste0(round(undiag_M*100,0),"% [",round(undiag_L*100,0),"-",round(undiag_U*100,0),"%]")]
dat[, N_inf:= paste0(N_inf_M," [",N_inf_L,"-",N_inf_U,"]")]
write.csv(dat,file=paste0(outfile.base,'-','number_infections_MSM_HSX.csv'))

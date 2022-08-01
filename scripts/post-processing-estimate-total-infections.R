
cat(" \n -------------------------------- \n \n estimate-total-infections.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_211101'
	#	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['in_dir']] <- '/Users/alexb/Documents/Roadmap/refactor_code'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'test_refactor_gqs'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	#args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['out_dir']] <- paste0('/Users/alexb/Documents/Roadmap/refactor_code/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['with_subtypes']] <- 1
	#args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['source_dir']] <- '~/Documents/GitHub/locally.acquired.infections-private'
	args_dir[['infdate']] <- 1
	args_dir[['start_y']] <- 2014
	args_dir[['undiag_job']]= 'undiag_untilmay2019_weights'
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

source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

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

infile.bplaces <- file.path(args$source_dir,'data','patient_data',paste0('birthplaces_subtype_',args$start_d,'.csv'))
infile.diagnoses.i <- file.path(args$source_dir,'data','patient_data',paste0('N_diagnosed_by_mg_infdate_since_',args$start_d,'.csv'))

cat("\nLoad time to diagnosis estimates \n")
dt_m <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_MSM.RDS'))
dt_h <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_HSX.RDS'))
dt_m <- readRDS(dt_m)
dt_h <- readRDS(dt_h)
dt_inf <- rbind(dt_m,dt_h)
setnames(dt_inf,'trsm','TRANSM')


cat(" \n -------------------------------- load prob(undiagnosed) -------------------------------- \n")

msm_und <- readRDS(file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args_dir$undiag_job),paste0('p_undiagnosed_samples_',args_dir$undiag_job,'_MSM.RDS')))
hsx_und <- readRDS(file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args_dir$undiag_job),paste0('p_undiagnosed_samples_',args_dir$undiag_job,'_HSX.RDS')))
dd <- rbind(msm_und,hsx_und)


cat("\nRead patient metadata \n")

dseq <- data.table(read.csv(infile.bplaces, header=T))
dt   <- data.table(read.csv(infile.diagnoses.i, header=T))

dseq <- dseq[, list(N_seq=sum(N)),by =c('TRM_GROUP','mwmb')]
setnames(dseq,'TRM_GROUP','TRANSM')

diag <- merge(dseq,dt,by=c('TRANSM','mwmb'),all=T)
diag[, mg_lab:= paste0(TRANSM,'_',mwmb)]
diag[, mg:= as.integer(factor(mg_lab,levels=c('MSM_NL','MSM_G1','MSM_G2','MSM_G3','MSM_Other',
																	 'HSX_NL','HSX_G4','HSX_G5','HSX_Other'),
									labels=c(4,1,2,3,5,3,1,2,4)))]
setnames(diag,'TRANSM','trsm')

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
									qlabel=c('L','M','U')),
					 by=c('trsm','mg')] 

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
									qlabel=c('L','M','U')),
					 by=c('mg')]

dat <- dcast(dat,mg~qlabel,value.var=c('undiag','N_diag','N_inf'))

dat[, undiag:= paste0(round(undiag_M*100,0),"% [",round(undiag_L*100,0),"-",round(undiag_U*100,0),"%]")]
dat[, N_inf:= paste0(N_inf_M," [",N_inf_L,"-",N_inf_U,"]")]
write.csv(dat,file=paste0(outfile.base,'-','number_infections_MSM_HSX.csv'))

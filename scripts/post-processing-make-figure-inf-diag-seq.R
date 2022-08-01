
cat(" \n -------------------------------- \n \n make-figure-inf-diag-seq.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
	args_dir[['analysis']] <- 'analysis_211101'
		args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	#args_dir[['in_dir']] <- '/Users/alexb/Documents/Roadmap/refactor_code'
	args_dir[['trsm']] <- 'HSX'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'test_refactor_gqs'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	#args_dir[['out_dir']] <- paste0('/Users/alexb/Documents/Roadmap/refactor_code/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/locally.acquired.infections-private'
	#args_dir[['source_dir']] <- '~/Documents/GitHub/locally.acquired.infections-private'
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

cat("\nRead patient metadata \n")

dseq <- data.table(read.csv(infile.bplaces, header=T))
dt   <- data.table(read.csv(infile.diagnoses.i, header=T))

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

dseq <- dseq[, list(N_seq=sum(N)),by =c('TRM_GROUP','mwmb')]
setnames(dseq,'TRM_GROUP','TRANSM')
diag <- merge(dseq,dt,by=c('TRANSM','mwmb'),all=T)
diag[, mg_lab:= paste0(TRANSM,'_',mwmb)]
diag[, mg:= as.integer(as.character(factor(mg_lab,levels=c('MSM_NL','MSM_G1','MSM_G2','MSM_G3','MSM_Other',
																							'HSX_NL','HSX_G4','HSX_G5','HSX_Other'),
															labels=c(4,1,2,3,5,3,1,2,4))))]
setnames(diag,'TRANSM','trsm')

dd <- merge(dd,diag,by=c('trsm','mg'),all=T)
dd[, N_inf:=round(diag/(1-av_undiagnosed),digits=0)]
dd[, mg_fb:=1]
dd[trsm=='MSM' & mg==4, mg_fb:=0]
dd[trsm=='HSX' & mg==3, mg_fb:=0]
dd[, tot:='All']
dd_fb <- dd[, list(seq=sum(N_seq),
									 diag=sum(diag),
									 N_inf=sum(N_inf)),
						by=c('trsm','mg_fb','iter')]
dd_all <- dd[, list(seq=sum(N_seq),
										diag=sum(diag),
										N_inf=sum(N_inf)),
						 by=c('trsm','tot','iter')]
setnames(dd_all,'tot','mg_fb')

dat <- rbind(dd_fb,dd_all)
dat <- dat[, list(seq=sum(seq),
									diag=sum(diag),
									N_inf=sum(N_inf)),
					 by=c('trsm','mg_fb','iter')]
dat[, undiag:=1-(diag/N_inf)]
dat <- dat[, list(N_seq=quantile(seq,prob=c(0.025,0.5,0.975)),
									N_diag=quantile(diag,prob=c(0.025,0.5,0.975)),
									N_inf=quantile(N_inf,prob=c(0.025,0.5,0.975)),
									qlabel=c('L','M','U')),
					 by=c('trsm','mg_fb')]

cat(" \n -------------------------------- plot -------------------------------- \n")

summ <- copy(dat)

summ$trsm <- factor(summ$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

summ[,bplace:='Dutch-born']
summ[mg_fb=='1',bplace:='Foreign-born']
summ[mg_fb=="All",bplace:='All']
summ$bplace <- factor(summ$bplace,levels=c('All','Dutch-born','Foreign-born'))

setnames(summ,c('N_seq','N_diag','N_inf'),c('Sequenced','Diagnosed','Infected'))

pct <- melt(subset(summ,qlabel=='M'),id.vars=c('trsm','mg_fb','bplace','qlabel'))
seqs <- subset(pct,variable=='Sequenced')
setnames(seqs,'value','seqs')
pct <- merge(pct,subset(seqs,select=c('trsm','bplace','seqs')),by=c('trsm','bplace'),all.x=T)
pct[, value:= round(value,digits=0)]
pct[, pct:=round(seqs/value*100,digits=0)]
pct[, lab:=paste0(value,'\n(',pct,'%)')]

dat <- melt(summ,id.vars=c('trsm','mg_fb','bplace','qlabel'))
dat <- dcast(dat,trsm+mg_fb+bplace+variable~qlabel,value.var='value')
dat[variable=='Sequenced' | variable=='Diagnosed', L:=NA]
dat[variable=='Sequenced' | variable=='Diagnosed', U:=NA]
dat$variable <- factor(dat$variable,levels=c('Infected','Diagnosed','Sequenced'))
pct$variable <- factor(pct$variable,levels=c('Infected','Diagnosed','Sequenced'))

plot <- ggplot(data=subset(dat)) +
	geom_bar(aes(x=bplace,y=M,fill=variable),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=bplace,ymin=L, ymax=U,fill=variable),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(data=pct,aes(label=lab,x=bplace,y=value,fill=variable), position=position_dodge(width=0.9), vjust=-1.35) +
	facet_grid(trsm~.,scales="free_y") +
	scale_y_continuous(expand = expansion(mult = c(0, .3)))  +
	labs(x='',y="Number of new infections",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-seq_diag_inf_2014-2018_CIs.png'),plot,w=8, h=9)

saveRDS(pct,file=paste0(outfile.base,'-seq_diag_inf_2014-2018_pct.RDS'))

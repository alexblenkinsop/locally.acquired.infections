# plot-subgraphs-diagdate-infdate.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running plot-subgraphs-diagdate-infdate.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['job_name']] <- 'undiagnosed_untilmay'
	args_dir[['period']] <- '2014-2018'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
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
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['analysis']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['out_dir']] <- args_line[[8]]
	args_dir[['job_tag']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
	args_dir[['source_dir']] <- args_line[[14]]
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))

#	load all input variables for this analysis run
do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes
start_d = args$start_d
end_d = args$end_d

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic.msm.i <- readRDS(file)
file <- paste0(outfile.base,'-obs_actual_cs_distribution_long.rds')
cat("\n read RDS:", file)
msm_i <- readRDS(file)
msm_i[, lab:='MSM 2014-2019']
msm_i[, date:='infdate']

args_dir[['job_name_1']] <- args_dir[['job_name']]
args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
args_dir[['trsm']] <- 'HSX'
args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic.hsx.i <- readRDS(file)
file <- paste0(outfile.base,'-obs_actual_cs_distribution_long.rds')
cat("\n read RDS:", file)
hsx_i <- readRDS(file)
hsx_i[, lab:='HSX 2014-2019']
hsx_i[, date:='infdate']

args_dir[['job_name']] <- 'undiagnosed_untilmay_diagd'
args_dir[['trsm']] <- 'MSM'
args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)
# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic.msm.d <- readRDS(file)
file <- paste0(outfile.base,'-obs_actual_cs_distribution_long.rds')
cat("\n read RDS:", file)
msm_d <- readRDS(file)
msm_d[, lab:='MSM 2014-2019']
msm_d[, date:='diagdate']

args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
args_dir[['trsm']] <- 'HSX'
args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic.hsx.d <- readRDS(file)
file <- paste0(outfile.base,'-obs_actual_cs_distribution_long.rds')
cat("\n read RDS:", file)
hsx_d <- readRDS(file)
hsx_d[, lab:='HSX 2014-2019']
hsx_d[, date:='diagdate']

args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name_1']],'_',args_dir[['period']],'_',args_dir[['trsm']])
args_dir[['job_tag']] <- paste0(args_dir[['job_name_1']],'_',args_dir[['period']],'_',args_dir[['trsm']])

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

cat(" \n -------------------------------- define plotting functions -------------------------------- \n")

all <- rbind(msm_i,hsx_i)

all[,N:=paste0( QN, ' (',round(p*100, d=1), '%)')]
QN_L <- subset(all,P=='p0.025',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_L,c('QN','p'),c('QN_L','p_L'))
QN_U <- subset(all,P=='p0.975',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_U,c('QN','p'),c('QN_U','p_U'))

all <- merge(subset(all,P=='p0.5'),QN_L,by=c('lab','analysis','chains','new_cases'),all=T)
all <- merge(all,QN_U,by=c('lab','analysis','chains','new_cases'),all=T)
all[analysis=='predicted',N:=paste0( QN, ' [', QN_L,'-',QN_U,']',' (',round(p*100, d=1), '%' , ' [', round(p_L*100, d=1),'-',round(p_U*100, d=1),'%',']' ,')')]

obs <- subset(all,analysis=='observed')

inftime <- obs



all <- rbind(msm_d,hsx_d)

all[,N:=paste0( QN, ' (',round(p*100, d=1), '%)')]
QN_L <- subset(all,P=='p0.025',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_L,c('QN','p'),c('QN_L','p_L'))
QN_U <- subset(all,P=='p0.975',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_U,c('QN','p'),c('QN_U','p_U'))

all <- merge(subset(all,P=='p0.5'),QN_L,by=c('lab','analysis','chains','new_cases'),all=T)
all <- merge(all,QN_U,by=c('lab','analysis','chains','new_cases'),all=T)
all[analysis=='predicted',N:=paste0( QN, ' [', QN_L,'-',QN_U,']',' (',round(p*100, d=1), '%' , ' [', round(p_L*100, d=1),'-',round(p_U*100, d=1),'%',']' ,')')]

obs <- subset(all,analysis=='observed')
diagtime <- obs

inftime[,date:='infection date']
diagtime[,date:='diagnosis date']
dat <- rbind(inftime,diagtime)


dat <- subset(dat,P=='p0.5' & new_cases %in% c('0','1','2','3','4','5','6','7+'))
dat[,trsm:='Amsterdam MSM']
dat[lab %in% c('HSX 2014-2019'),trsm:='Amsterdam heterosexual']

dat[chains=='pre-existing', chains:='pre-existing*']
dat$chains <- factor(dat$chains,levels=c('pre-existing*','emergent'))

dat$trsm <- factor(dat$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexual'))

dat[, time:= '2014-2019']

plot <- ggplot(data=subset(dat)) +
	geom_bar(aes(x=new_cases,y=QN,fill=date),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=QN_L, ymax=QN_U,fill=date),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#scale_y_continuous(trans="sqrt",labels = scales::comma) +
	scale_y_sqrt(breaks=seq(0,600,50)) +
	#scale_y_continuous(trans='log2') +
	facet_grid(trsm~chains,scales="free_y") +
	labs(x='\nNumber of cases',y="Number of phylogenetically observed subgraphs",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_nejm()
ggsave(file=paste0(outfile.base,'-observed_chains_diagd_infd.png'),plot,w=15, h=12)

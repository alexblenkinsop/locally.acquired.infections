# post-processing-posterior-predictive-check-trsm-group.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n post-processing-posterior-predictive-check-trsm-group.R\n \n -------------------------------- \n")

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
args_dir[['stanModelFileHSX']] <- 'branching_process_210810i_cmdstan'
args_dir[['period']] <- '2014-2019'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['job_name']] <- 'weibull_est_undiag'
args_dir[['infdate']] <- 1

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFileMSM')	
	stopifnot(args_line[[3]]=='-stanModelFileHSX')	
	stopifnot(args_line[[5]]=='-in_dir')
	stopifnot(args_line[[7]]=='-job_name')
	stopifnot(args_line[[9]]=='-period')
	stopifnot(args_line[[11]]=='-source_dir')
	
	args_dir <- list()
	args_dir[['stanModelFileMSM']] <- args_line[[2]]
	args_dir[['stanModelFileHSX']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['job_name']] <- args_line[[8]]
	args_dir[['period']] <- args_line[[10]]
	args_dir[['source_dir']] <- args_line[[12]]
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

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.MSM.p1 <- readRDS(file)
file <- paste0(outfile.base,'-posteriorpredictivecheck.RDS')
#file <- paste0(outfile.base,'-posteriorpredictivecheck_allsbts.RDS')
cat("\n read RDS:", file)
ppc.MSM <- readRDS(file)

ppc.MSM[, trsm:='MSM']
ppc.MSM[, time:='2014-2019']

args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileHSX']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_HSX')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_HSX')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileHSX , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic.HSX.p1 <- readRDS(file)
#file <- paste0(outfile.base,'-posteriorpredictivecheck_allsbts.RDS')
file <- paste0(outfile.base,'-posteriorpredictivecheck.RDS')
cat("\n read RDS:", file)
ppc.HSX <- readRDS(file)

ppc.HSX[, trsm:='HSX']
ppc.HSX[, time:='2014-2019']

cat(" \n -------------------------------- define plotting functions -------------------------------- \n")

dat <- rbind(ppc.MSM,ppc.HSX)
dat$subtypes_name <- factor(dat$subtypes_name, levels=c('B','nonB'),labels=c('B','Non-B'))
dat$trsm <- factor(dat$trsm, levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

plot <- ggplot(data=dat) +
	geom_bar(aes(x=new_cases,y=p0.5,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=p0.025, ymax=p0.975,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#geom_text(aes(12, N,label=paste("N= ", N, "p=", round(p,2)))) +
	xlim(-1,15) +
	scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
	facet_grid(trsm~subtypes_name+chains,scales="free_y") +
	labs(x='\nNumber of cases since 2015',y="Number of chains",fill="") +
	theme_bw(base_size=20) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_newcases_chaintype_trmgroup_allsbts.png'),plot,w=15, h=10)


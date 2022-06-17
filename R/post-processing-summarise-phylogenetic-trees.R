
cat(" \n -------------------------------- \n \n Summarise-origins-chains.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	#args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['in_dir']] <- '/Users/alexb/Documents/Roadmap'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'undiagnosed_weighted_inf_rate'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	#args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['out_dir']] <- paste0('/Users/alexb/Documents/Roadmap/elife_results_resubmission/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['infdate']] <- 1
	args_dir[['start_d']] <- 2014
	args_dir[['end_d']] <- 2019
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
	stopifnot(args_line[[13]]=='-start_d')
	stopifnot(args_line[[15]]=='-end_d')
	stopifnot(args_line[[17]]=='-source_dir')
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['analysis']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['out_dir']] <- args_line[[8]]
	args_dir[['job_tag']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
	args_dir[['start_d']] <- args_line[[14]]
	args_dir[['end_d']] <- args_line[[16]]
	args_dir[['source_dir']] <- args_line[[18]]
	args_dir[['infdate']] <- 1
} 


#	load all input variables for this analysis run
do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic <- readRDS(file)

file <- paste0(outfile.base,'-stanout-origins-gqs.RDS')
cat("\n read RDS:", file)
fit.origins <- readRDS(file)

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

cat(" \n -------------------------------- Reading data -------------------------------- \n")

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)


ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

`%notin%` <- Negate(`%in%`)

cat(" \n -------------------------------- summarise observed origins -------------------------------- \n")

cat('\nSummarise predicted origins...')

infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')

dsubgraphtaxa <- readRDS(infile.subgraphs)

dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
do <- dsubgraphtaxa
# count unique subgraphs with origin in each region (by clade too)
do <- do[,list(N=length(unique(FULL_NAME))),by=c('REP','ORIGIN','TRANSM','ST','ST_CLADE')]
# aggregate over clades for B
do <- do[,list(N=sum(N)),by=c('REP','ORIGIN','TRANSM','ST')]
do <- subset(do, ORIGIN!='Unknown')
do <- do[, list(ORIGIN=ORIGIN,N=N,pct=N/sum(N)),by=c('REP','TRANSM','ST')]
da <- do[, list(N= round(quantile(N, prob=ps,na.rm=T),0),
								pct= round(quantile(pct, prob=ps,na.rm=T)*100,1),
								q_label=p_labs),by=c('TRANSM','ST','ORIGIN')]		
dn <- dcast(da,TRANSM+ST+ORIGIN~q_label,value.var="N")
dp <- dcast(da,TRANSM+ST+ORIGIN~q_label,value.var="pct")
dn[, N_CI:= paste0(M, " [",CL,"-",CU,"]")]
dp[, pct_CI:= paste0(M, "% [",CL,"-",CU,"%]")]

dt <- merge(subset(dn,select=c('TRANSM','ST','ORIGIN','N_CI')),
						subset(dp,select=c('TRANSM','ST','ORIGIN','pct_CI')),
						by=c('TRANSM','ST','ORIGIN'),all=T)

write.csv(da,file=paste0(outfile.base,'-','observed_phylo_origins.csv'))
cat(" \n -------------------------------- \n \n Running post-processing-combine-MSM-HSX-results.R\n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))

args <- list()
args[['source_dir']] <- '~/git/locally.acquired.infections'
args[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args[['stanModelFileMSM']] <- 'branching_process_210810b_cmdstan'
args[['stanModelFileHSX']] <- 'branching_process_210810m_cmdstan'
args[['job_name']] <- 'elife_paper'
args[['period']] <- '2014-2018'
args[['start_d']] <- 2014
args[['end_d']] <- 2019

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-in_dir')
	stopifnot(args_line[[5]]=='-stanModelFileMSM')	
	stopifnot(args_line[[7]]=='-stanModelFileHSX')	
	stopifnot(args_line[[9]]=='-job_name')

	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['in_dir']] <- args_line[[4]]
	args[['stanModelFileMSM']] <- args_line[[6]]
	args[['stanModelFileHSX']] <- args_line[[8]]
	args[['job_name']] <- args_line[[10]]
} 

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=23, hpc.select=1, hpc.nproc=8, hpc.mem= "80gb", hpc.load= "module load anaconda3/personal\nsource activate bpm", hpc.q=NA, hpc.array=1 )
{	
	pbshead <- "#!/bin/sh"
	tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
	pbshead <- paste(pbshead, tmp, sep = "\n")
	tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
	pbshead <- paste(pbshead, tmp, sep = "\n")
	pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
	if(hpc.array>1)
	{
		pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
	}				
	if(!is.na(hpc.q))
	{
		pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
	}		
	pbshead	<- paste(pbshead, hpc.load, sep = "\n")
	pbshead
}

cmd <- paste0('echo "----------- Combining MSM/HSX results: ------------"\n')
cmd2 <- make.PBS.header( hpc.walltime=23, 
												 hpc.select=1, 
												 hpc.nproc=48, 
												 hpc.mem= "30gb", 
												 hpc.load= "module load anaconda3/personal\nsource activate bpm", 
												 hpc.q=NA,
												 hpc.array= 1)
cmd2 <- paste0(cmd2,'\n')
# set up env variables
cmd2 <- paste0( cmd2, 'SCRIPT_DIR=',args$source_dir,'\n')

tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-summarise-bplace-newcases.R'),
							' -stanModelFileMSM "', args$stanModelFileMSM,'" -stanModelFileHSX "', args$stanModelFileHSX, '" -in_dir "', args$in_dir,	'" -job_name "', args$job_name,				
							'" -period "', args$period ,'" -start_d ', args$start_d ,' -end_d ', args$end_d ,' -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')

tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-estimates-by-migrant-groups.R'),
							' -stanModelFileMSM "', args$stanModelFileMSM,'" -stanModelFileHSX "', args$stanModelFileHSX, '" -in_dir "', args$in_dir, '" -job_name "', args$job_name,	
							'" -period "', args$period ,'" -start_d ', args$start_d ,' -end_d ', args$end_d ,' -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')

tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-make-observed-actual-chain-sizes-figure-table.R'),
							' -stanModelFileMSM "', args$stanModelFileMSM,'" -stanModelFileHSX "', args$stanModelFileHSX, '" -in_dir "', args$in_dir, '" -job_name "', args$job_name,				
							'" -period "', args$period ,'" -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')

tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-posterior_predictive_check_trsm_group.R'),
							' -stanModelFileMSM "', args$stanModelFileMSM,'" -stanModelFileHSX "', args$stanModelFileHSX, '" -in_dir "', args$in_dir, '" -job_name "', args$job_name,				
							'" -period "', args$period ,'" -source_dir $SCRIPT_DIR')
cmd2 <- paste0(cmd2,tmp,'\n')

# write submission file	
tmpdir <- paste0(args[['in_dir']],'/branching_process_model/',args$stanModelFileHSX,'-',args$job_name,'_',args$period,'_HSX')
post.processing.file <- file.path(tmpdir, 'post_processing_combine_analyses.sh')
cat(cmd2, file=post.processing.file)
# set permissions
Sys.chmod(post.processing.file, mode='644')	

cmd 		<- paste("qsub", post.processing.file)
cat(cmd)
cat(system(cmd, intern= TRUE))	


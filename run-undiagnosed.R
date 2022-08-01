require(data.table)

hmc_chains_n <- 3

## important note:
# the combination of stanModelFile and job_tag should be unique for each analysis
# all outputs with stanModelFile-job_tag are assumed to be several HMC chains run in parallel

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

# input args alex
if(1)
{
	hpc.nproc.cmdstan <- 5	
	args <- data.table(
		source_dir= '/rds/general/user/ablenkin/home/git/locally.acquired.infections',
		cmdstan_dir = '/apps/cmdstan/2.33.0',
		in_dir='/rds/general/project/ratmann_roadmap_data_analysis/live',
		out_dir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
    script_file='scripts/undiagnosed_bytrmgroup.R',
		post_script_file='scripts/post-processing-make-figures-undiagnosed.R',
		stanModelFile= 'undiagnosed_211102',
		analysis= 'analysis_211101',
		hmc_stepsize= 0.25,
		hmc_num_samples= 2000,
		hmc_num_warmup= 500,			
		seed= 42,
		chain= 1,
		job_tag= 'test_bash',
		cmdstan = 1L,
		sens=F,
		weights='ECDC'
	)	
}


if(1)
{
	tmp <- data.table(chain=1:hmc_chains_n)		
	#tmp[, seed:= round(runif(seq_len(nrow(tmp)))*1e6)]		
	set(args, NULL, colnames(tmp), NULL)
	tmp[, dummy:= 1L]
	args[, dummy:= 1L]
	args <- merge(args, tmp, by='dummy')
	set(args, NULL, 'dummy', NULL)
}

i <- 1
# make commands
	cmd       <- ''			
	#	general housekeeping
	cmd       <- paste0(cmd,"CWD=$(pwd)\n")
	cmd       <- paste0(cmd,"echo $CWD\n")	
	tmpdir.prefix	<- paste0('undgsd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
	tmpdir    <- paste0("$CWD/",tmpdir.prefix)
	cmd       <- paste0(cmd,"mkdir -p ",tmpdir,'\n')	
	#	generate data set and run if not using cmdstan
	cmd       <- paste0( cmd, 'echo "----------- Generating input data: ------------"\n')
	tmp       <- paste0('Rscript ', file.path(args$source_dir[i],args$script_file[i]), 
									 ' -source_dir "', args$source_dir[i],'"',
									 ' -stanModelFile "', args$stanModelFile[i],'"',
									 ' -analysis "', args$analysis[i],'"',
									 ' -seed ', args$seed[i],
									 #' -chain ', args$chain[i],							
									 ' -indir ', args$in_dir[i],'',
									 ' -outdir ', tmpdir,'',
									 ' -jobtag "', args$job_tag[i],'"',
									 ' -sens "', args$sens[i],'"',
									 ' -weights "', args$weights[i],'"'
	)
	cmd				<- paste0(cmd, tmp, '\n')
	#	generate plots
	cmd       <- paste0( cmd, 'echo "----------- Post-processing: ------------"\n')
	tmp       <- paste0('Rscript ', file.path(args$source_dir[i],args$post_script_file[i]), 
											' -source_dir "', args$source_dir[i],'"',
											' -stanModelFile "', args$stanModelFile[i],'"',
											' -analysis "', args$analysis[i],'"',
											' -indir ', args$in_dir[i],'',
											' -outdir ', tmpdir,'',
											' -jobtag "', args$job_tag[i],'"',
											' -sens "', args$sens[i],'"'
	)
	cmd				<- paste0(cmd, tmp, '\n')
	
	
	#	general housekeeping
	cmd 	<- paste0( cmd, 'echo "----------- Copy files to out directory: ------------"\n')
	tmpdir2	<- file.path(args$out_dir[i], paste0(args$stanModelFile[i],'-',args$job_tag[i]))
	dir.create(tmpdir2)		  		
	cmd		<- paste0(cmd,"mkdir -p ",tmpdir2,'\n')
	cmd		<- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', tmpdir2,'\n')
	cmd		<- paste0(cmd, 'chmod -R g+rw ', tmpdir2,'\n')

	if(args$cmdstan[1]==0)
	{
		pbshead <- make.PBS.header(	hpc.walltime=23, 
																hpc.select=1, 
																hpc.nproc=1, 
																hpc.mem= "30gb", 
																hpc.load= "module load anaconda3/personal\nsource activate bpm", 
																hpc.q=NA,
																hpc.array= length(cmd) )
	}
	if(args$cmdstan[1]==1)
	{ 
		pbshead <- make.PBS.header(	hpc.walltime=23, 
																hpc.select=1, 
																hpc.nproc=hpc.nproc.cmdstan, 
																hpc.mem= paste0(hpc.nproc.cmdstan*9,'gb'), 
																hpc.load= paste0("module load cmdstan/2.33.0 anaconda3/personal\nsource activate bpm\nexport STAN_NUM_THREADS=",hpc.nproc.cmdstan,"\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"), 
																hpc.q=NA,
																hpc.array= length(cmd) )
	}
	
cmd		<- paste(pbshead,cmd,sep='\n')

#	submit job
outfile		<- gsub(':','',paste("bpm",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out_dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))	

require(data.table)

hmc_chains_n <- 3

## important note:
# the combination of stanModelFile and job_tag should be unique for each analysis. job_tag should be the same for MSM and HSX model
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

# input args
if(1)
{
  hpc.nproc.cmdstan <- 5	
  args <- data.table(
    source_dir= '/rds/general/user/ablenkin/home/git/locally.acquired.infections',
    cmdstan_dir = '/apps/cmdstan/2.33.0',
    in_dir='/rds/general/project/ratmann_roadmap_data_analysis/live',
    out_dir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
    report_dir = '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/reports',
    script_file= 'scripts/stan-make-data.r',
    script_converting_file = "scripts/stan-convert-csv-to-rda.r",
    script_generate_quantities_file = "scripts/generate-quantities.R",
    script_rmd_file = "scripts/post-processing-make-report.Rmd",
    stanModelFile= 'branching_process_210810m_cmdstan',
    analysis= 'analysis_211101',
    hmc_stepsize= 0.25,
    hmc_num_samples= 2000,
    hmc_num_warmup= 500,			
    seed= 42,
    chain= 1,
    job_tag= 'elife_paper',
    trsm= 'HSX',
    cmdstan = 1L,
    max_index_cases = 22L,
    start_d = 2014,
    end_d = 2019,
    index_flag= 1,
    infdate=1,
    nonB=1,
    sensitivity_infdate=0,
    undiag_job='undiagnosed_weighted_ECDC'
  )	
}

if(1)
{
  tmp <- data.table(chain=1:hmc_chains_n)		
  tmp[, seed:= c(82525,735379,877954)]
  set(args, NULL, colnames(tmp), NULL)
  tmp[, dummy:= 1L]
  args[, dummy:= 1L]
  args <- merge(args, tmp, by='dummy')
  set(args, NULL, 'dummy', NULL)
  args$job_tag <- paste0(args$job_tag,'_',args$start_d,'-',args$end_d - 1,'_',args$trsm)
}

# make commands
cmds <- vector('list', nrow(args))
for(i in seq_len(nrow(args)))
{
  cmd				<- ''			
  #	general housekeeping
  cmd				<- paste0(cmd,"CWD=$(pwd)\n")
  cmd				<- paste0(cmd,"echo $CWD\n")	
  tmpdir.prefix	<- paste0('bpm_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir			<- paste0("$CWD/",tmpdir.prefix)
  cmd				<- paste0(cmd,"mkdir -p ",tmpdir,'\n')	
  #	generate data set and run if not using cmdstan
  cmd 			<- paste0( cmd, 'echo "----------- Generating input data: ------------"\n')
  tmp 			<- paste0('Rscript ', file.path(args$source_dir[i],args$script_file[i]), 
                   ' -source_dir "', args$source_dir[i],'"',
                   ' -stanModelFile "', args$stanModelFile[i],'"',
                   ' -seed ', args$seed[i],
                   ' -chain ', args$chain[i],							
                   ' -indir ', args$in_dir[i],'',
                   ' -outdir ', tmpdir,'',
                   ' -jobtag "', args$job_tag[i],'"',
                   ' -trsm "', args$trsm[i],'"',
                   ' -cmdstan ', args$cmdstan[i],
                   ' -max_index_cases ', args$max_index_cases[i],
                   ' -start_d ', args$start_d[i],
                   ' -end_d ', args$end_d[i],
                   ' -index_flag ', args$index_flag[i],
                   ' -infdate ', args$infdate[i],
                   ' -nonB ', args$nonB[i],
                   ' -sensitivity_infdate ', args$sensitivity_infdate[i],
                   ' -undiag_job ', args$undiag_job[i]
  )
  cmd				<- paste0(cmd, tmp, '\n')
  #	if using cmdstan  
  if(args$cmdstan[i]==1) 
  {
    cmd <- paste0(cmd, 'echo "----------- Building Stan model file: ------------"\n')
    #	clean up any existing model code
    cmd <- paste0(cmd, 'rm ', file.path('$CWD',paste0(args$stanModelFile[i],'.*')), ' \n')
    #	copy stan model file
    cmd	<- paste0(cmd, 'cp -R ',file.path(args$source_dir[i], 'stan-models',paste0(args$stanModelFile[i],'.stan')),' .\n')
    #	build model		
    cmd <- paste0(cmd, 'cd ', args$cmdstan_dir[i], '\n')
    cmd <- paste0(cmd, 'make STAN_THREADS=TRUE ', file.path('$CWD',args$stanModelFile[i]), ' \n')
    cmd <- paste0(cmd, 'cd $CWD\n')
    #	set up env variables
    cmd <- paste0( cmd, 'JOB_DIR=$(ls -d "',tmpdir,'"/*/)\n')
    cmd <- paste0( cmd, 'JOB_DIR=${JOB_DIR%?}\n')
    cmd <- paste0( cmd, 'JOB_DIR_NAME=${JOB_DIR##*/}\n')
    cmd <- paste0( cmd, 'SCRIPT_DIR=',args$source_dir[i],'\n')
    cmd <- paste0( cmd, 'IN_DIR=',args$in_dir[i],'\n')
    cmd <- paste0( cmd, 'ANALYSIS=',args$analysis[i],'\n')
    cmd <- paste0( cmd, 'TRSM=',args$trsm[i],'\n')
    cmd <- paste0( cmd, 'START_Y=',args$start_d[i],'\n')
    cmd <- paste0( cmd, 'UNDIAG_JOB=',args$undiag_job[i],'\n')
    cmd <- paste0( cmd, 'STAN_MODEL_FILE=',args$stanModelFile[i],'\n')
    cmd <- paste0( cmd, 'STAN_DATA_FILE=$(find ', tmpdir, ' -name "*cmdstanin.R")\n')
    cmd <- paste0( cmd, 'STAN_INIT_FILE=$(find ', tmpdir, ' -name "*cmdstaninit.R")\n')
    cmd <- paste0( cmd, 'STAN_OUT_FILE=', file.path('$JOB_DIR','${JOB_DIR##*/}_stanout.csv'),' \n')
    #	run model
    cmd <- paste0( cmd, 'echo "----------- env variables are: ------------"\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR_NAME\n')
    cmd <- paste0( cmd, 'echo $STAN_DATA_FILE\n')
    cmd <- paste0( cmd, 'echo $STAN_OUT_FILE\n')
    cmd <- paste0( cmd, 'echo $IN_DIR\n')
    cmd <- paste0( cmd, 'echo $ANALYSIS\n')
    cmd <- paste0( cmd, 'echo "----------- Starting Stan sampling: ------------"\n')
    #	
    tmp <- paste0( './',args$stanModelFile[i],' ',
                   'sample num_samples=',args$hmc_num_samples[i],' num_warmup=',args$hmc_num_warmup[i],' save_warmup=0 thin=1 ',
                   'adapt delta=0.95 ',
                   'algorithm=hmc engine=nuts max_depth=15 stepsize=',args$hmc_stepsize[i],' ',
                   'data file=$STAN_DATA_FILE ',
                   'init=$STAN_INIT_FILE ',
                   'random seed=',args$seed[i],' ',
                   'output file=$STAN_OUT_FILE' )
    cmd <- paste0(cmd, tmp, '\n')
    # convert csv to rdata
    cmd		<- paste0( cmd, 'echo "----------- Converting Stan output to RDA file: ------------"\n')
    tmp		<- paste0('Rscript ', file.path(args$source_dir[i],args$script_converting_file[i]), 
                   ' -csv_file "', "$STAN_OUT_FILE",'"',
                   ' -rda_file "', file.path('$JOB_DIR','${JOB_DIR##*/}_stanout.RData'),'"'
    )
    cmd		<- paste0(cmd, tmp, '\n')		
  }			
  
  #	general housekeeping
  cmd 	<- paste0( cmd, 'echo "----------- Copy files to out directory: ------------"\n')
  tmpdir2	<- file.path(args$out_dir[i], paste0(args$stanModelFile[i],'-',args$job_tag[i]))
  if(i==1)
  {
    dir.create(tmpdir2)		  		
  }
  cmd		<- paste0(cmd,"mkdir -p ",tmpdir2,'\n')
  cmd		<- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ', tmpdir2,'\n')
  cmd		<- paste0(cmd, 'chmod -R g+rw ', tmpdir2,'\n')
  
  #	generate quantities 
  cmd <- paste0( cmd, 'echo "----------- Generating quantities: ------------"\n')
  cmd <- paste0( cmd, 'JOB_DIR2="',tmpdir2,'"/"$JOB_DIR_NAME" \n')
  cmd <- paste0( cmd, 'OUT_DIR="',tmpdir2,'" \n')
  cmd <- paste0( cmd, 'echo $JOB_DIR2\n')
  cmd <- paste0( cmd, 'echo $OUT_DIR\n')
  cmd <- paste0( cmd, paste0("echo {1..",args$max_index_cases[i],"} | tr ' ' '\\n' | ") )
  tmp <- ifelse(args$cmdstan[i]==1, hpc.nproc.cmdstan, 1)
  stopifnot(is.numeric(tmp))
  cmd <- paste0( cmd, paste0('xargs -P ',tmp,' -n 1 -I {} ') )
  tmp <- paste0('Rscript ', file.path(args$source_dir[i],args$script_generate_quantities_file[i]),
                ' -source_dir "$SCRIPT_DIR"',
                ' -indir "$IN_DIR"',
                ' -outdir "$OUT_DIR"',
                ' -analysis "$ANALYSIS"',
                ' -indir.results "$JOB_DIR2"',
                ' -trsm "$TRSM"',
                ' -index.case.index {}')		
  cmd <- paste0(cmd, tmp,'\n')
  
  # create post-processing shell script for central analyses
  if(i==1)
  {
    cmd2 <- make.PBS.header( hpc.walltime=23, 
                             hpc.select=1, 
                             hpc.nproc=48, 
                             hpc.mem= "124gb", 
                             hpc.load= "module load anaconda3/personal\nsource activate bpm", 
                             hpc.q=NA,
                             hpc.array= 1)
    cmd2 <- paste0(cmd2,'\n')
    # set up env variables
    cmd2 <- paste0(cmd2,'SCRIPT_DIR=',args$source_dir[i],'\n',			
                   'IN_DIR=',args$in_dir[i],'\n',
                   'OUT_DIR=',tmpdir2,'\n',
                   'JOB_TAG=',args$job_tag[i],'\n',
                   'STAN_MODEL_FILE=',args$stanModelFile[i],'\n',
                   'ANALYSIS=',args$analysis[i],'\n',
                   'NUMB_CHAINS=', max(args$chain),'\n',
                   'OVERWRITE=0\n',
                   'TRSM=', args$trsm[i],'\n',
                   'START_Y=', args$start_d[i],'\n',
                   'END_Y=', args$end_d[i],'\n',
                   'UNDIAG_JOB=', args$undiag_job[i],'\n'
    )
    # save posterior samples
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-save-posterior-samples.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG -numb_chains $NUMB_CHAINS -source_dir $SCRIPT_DIR -trsm $TRSM')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing assess HMC mixing
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-assess-mixing.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing make plot  of transmission pars
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-transmission-pars.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing importations
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-locally-acquired-infections.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing posterior predictive check
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-posterior_predictive_check.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing make tables
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-make-summary-tables.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing summarise chain sizes
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-summarise-chain-sizes.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing summarise chain sizes
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-summarise-subgraphs-transmission-chains.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR -start_y $START_Y -undiag_job $UNDIAG_JOB')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing plot seq diag inf
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-make-figure-inf-diag-seq.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR -start_y $START_Y -undiag_job $UNDIAG_JOB')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing summarise origins
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-summarise-origins-chains.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -start_d $START_Y -end_d $END_Y -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing make phylogenetic subgraphs plot
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-figure-phylogenetic-subgraphs.R'),
    							' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -start_d $START_Y -end_d $END_Y -source_dir $SCRIPT_DIR')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing make avidity assay plot
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-make-plot-recent-infections.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM')
    cmd2 <- paste0(cmd2,tmp,'\n')
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-estimate-total-infections.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -analysis $ANALYSIS -in_dir $IN_DIR -out_dir $OUT_DIR -job_tag $JOB_TAG -trsm $TRSM -source_dir $SCRIPT_DIR -start_y $START_Y -undiag_job $UNDIAG_JOB')
    cmd2 <- paste0(cmd2,tmp,'\n')
    # postprocessing knit report
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','post-processing-knit-report.R'),
                  ' -rmd_file "', args$script_rmd_file[i],'" -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -job_tag $JOB_TAG',				
                  ' -report_dir "', args$report_dir[i],'" -source_dir $SCRIPT_DIR'					
    )
    cmd2 <- paste0(cmd2,tmp,'\n')
    # write submission file	
    post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
    cat(cmd2, file=post.processing.file)
    # set permissions
    Sys.chmod(post.processing.file, mode='644')	
  }
  #	schedule post-processing	
  cmd		<- paste0( cmd, 'echo "----------- Post-processing: ------------"\n')
  tmp		<- paste("if [ $(find ",tmpdir2," -name '*_stanout.RData' | wc -l) -ge ",max( args$chain )," ]; then\n",sep='')
  cmd		<- paste(cmd,tmp,sep='')	
  post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
  cmd 	<- paste0(cmd, '\tcd ', dirname(post.processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(post.processing.file),'\n')
  cmd		<- paste0(cmd,"fi\n")
  cmd		<- paste(cmd, "rm -rf $CWD/", basename(args$source_dir[i]),'\n',sep='')
  cat(cmd)
  cmds[[i]]	<- cmd	
}	
if(args$cmdstan[1]==0)
{
  pbshead <- make.PBS.header(	hpc.walltime=23, 
                              hpc.select=1, 
                              hpc.nproc=1, 
                              hpc.mem= "30gb", 
                              hpc.load= "module load anaconda3/personal\nsource activate bpm", 
                              hpc.q=NA,
                              hpc.array= length(cmds) )
}
if(args$cmdstan[1]==1)
{ 
  pbshead <- make.PBS.header(	hpc.walltime=23, 
                              hpc.select=1, 
                              hpc.nproc=hpc.nproc.cmdstan, 
                              hpc.mem= paste0(hpc.nproc.cmdstan*9,'gb'), 
                              hpc.load= paste0("module load cmdstan/2.33.0 anaconda3/personal\nsource activate bpm\nexport STAN_NUM_THREADS=",hpc.nproc.cmdstan,"\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"), 
                              hpc.q=NA,
                              hpc.array= length(cmds) )
}

#	make array job
for(i in seq_len(nrow(args)))
{
  cmds[[i]] <- paste0(i,')\n',cmds[[i]],';;\n')
}
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')			
cmd		<- paste(pbshead,cmd,sep='\n')

#	submit job
outfile		<- gsub(':','',paste("bpm",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out_dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))	

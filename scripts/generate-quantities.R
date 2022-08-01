require(data.table)
require(rstan)
require(EnvStats)

## command line parsing if any
# local
if(0){
	args_dir = list()
	args_dir[['index.case.index']] = 1
	args_dir[['indir.results']] = '~/Box Sync/Roadmap/RQ1 Estimating introductions/branching_process_model/branching_process_210810a_cmdstan-undiagnosed_B_nonB_vmrre_2015-2019_HSX/branching_process_210810b_cmdstan-undiagnosed_B_nonB_vmrre_2015-2019_HSX-4268449[1].pbs'
	args_dir[['outdir']] = '~/Box Sync/Roadmap/RQ1 Estimating introductions/branching_process_model/branching_process_210810a_cmdstan-undiagnosed_B_nonB_vmrre_2015-2019_HSX'
	args_dir[['source_dir']] <- '~/Documents/GitHub/bpm'
	args_dir[['indir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['trsm']] <- 'HSX'
}

# hpc
if(0){
  args_dir = list()
  args_dir[['index.case.index']] = 1
  args_dir[['indir.results']] = '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810b_cmdstan-test_refactor_2014-2018_MSM/branching_process_210810b_cmdstan-test_refactor_2014-2018_MSM-5978946[1].pbs'
  args_dir[['outdir']] = '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810b_cmdstan-test_refactor_2014-2018_MSM'
  args_dir[['source_dir']] <- '~/git/locally.acquired.infections-private'
  args_dir[['indir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  args_dir[['analysis']] <- 'analysis_211101'
  args_dir[['trsm']] <- 'MSM'
  args_dir[['infdate']] <- 1
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{  
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-indir')
	stopifnot(args_line[[5]]=='-outdir')
	stopifnot(args_line[[7]]=='-analysis')
	stopifnot(args_line[[9]]=='-indir.results')
	stopifnot(args_line[[11]]=='-trsm')
	stopifnot(args_line[[13]]=='-index.case.index')

	args_dir <- list()    
	args_dir[['source_dir']] <- args_line[[2]]  
	args_dir[['indir']] <- args_line[[4]]  
	args_dir[['outdir']] <- args_line[[6]]  
	args_dir[['analysis']] <- args_line[[8]]   
	args_dir[['indir.results']] <- args_line[[10]]  
	args_dir[['trsm']] <- args_line[[12]]  
	args_dir[['index.case.index']] <- args_line[[14]]  
}

#	print args
str(args_dir)

#	read args
indir <- args_dir$indir
analysis <- args_dir$analysis
indir.results <- args_dir$indir.results
index.case.index <- as.integer(args_dir$index.case.index)

## load functions
source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

#	read Stan input data and add index.case.index to parallelise computations
cat('\nReading Stan input data...')
infile.stanin <- list.files(indir.results, pattern='stanin.RData$', recursive=TRUE)

stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(indir.results, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)

# add origins & birthplaces to stan.data
cat("\nAdd origins of subgraphs and birthplaces of cases \n")
infile.bplaces <- file.path(args_dir$source_dir,'data','patient_data',paste0('birthplaces_subtype_',args$start_d,'.csv'))
stan.data = stan_data_add_origins_sbt(args_dir$trsm,stan.data,args_dir$outdir,args_dir$indir.results,args$start_d,args$end_d,args$infdate)
stan.data = stan_data_add_bplace_sbt(args_dir$trsm,stan.data,args_dir$outdir,args_dir$indir.results,infile.bplaces,args$start_d,args$end_d,args$infdate)

# add number of subgraphs to predict
cat("\nAdd number of subgraphs to predict\n")

stan.data$N_sgs_m <- sum(stan.data$cs_obs)
stan.data$z <- seq(1:max(stan.data$N_cs_actual))

rstan::stan_rdump( names(stan.data), file=file.path(args_dir$indir.results, paste0(basename(args_dir$indir.results), '_gqs_cmdstanin.R')), envir=list2env(stan.data))  	

# generate sequence of chain sizes for prediction

stan.data$INDEX_CASE_IDX <- index.case.index

#	reset args
args[['work_dir']] <- getwd()

#	read stanfits
cat('\nReading Stanfit ...')
infile.stanfits <- list.files(indir.results, pattern='*_stanout.RData$', recursive=TRUE)
stopifnot(length(infile.stanfits)>0)
stopifnot(length(infile.stanfits)<=1)
tmp <- load(file.path(indir.results, infile.stanfits[1]))

#	find stan model and stan model for generating quantities

file_stanModel <- gsub('.*locally.acquired.infections-private/(.*)','\\1',args$file_stanModel)
file_stanModel <- file.path(args_dir$source_dir, file_stanModel)
file_stanModel_gqs <- gsub('cmdstan','gqs',file_stanModel)
stopifnot( file.exists(file_stanModel_gqs) )

cat('\nCompiling gqs model file ...')
m2 <- rstan::stan_model(file_stanModel_gqs)

cat('\nGenerating quantities ...')
draws <- as.matrix(fit)
fit2 <- rstan::gqs(m2, data=stan.data, draws=draws)
fit.gqs <- rstan::extract(fit2)

file <- file.path(indir.results, paste0(basename(args$job_dir), '_indexcase',index.case.index,'_stangqs.RDS'))
cat('\nSave ', file, '\n')
if(!is.null(fit.gqs)) saveRDS(fit.gqs, file=file)

cat('\nFinished base-ages-generate-quantities.r ...')

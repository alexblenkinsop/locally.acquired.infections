library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
library(ggsci)


args <- list( 
	#source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	#indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	#outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	source_dir= '~/Documents/GitHub/locally.acquired.infections-private',
	indir='~/Documents/Roadmap/refactor_code/branching_process_model',
	outdir= '~/Documents/Roadmap/refactor_code/branching_process_model',
	stanModelFile= 'branching_process_210810b_cmdstan',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	job_tag= 'test',
	trsm= 'MSM',
	cmdstan = 1L,
	max_index_cases = 22,
	start_d = 2014,
	end_d = 2019,
	index_flag = 1,
	upper.bound.multiplier = 10,
	infdate=1,
	nonB=1,
	sensitivity_infdate=0,
	undiag_job='undiag_untilmay2019_weights'
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-seed', !is.na(as.integer(args_line[[6]])))	
	stopifnot(args_line[[7]]=='-chain', !is.na(as.integer(args_line[[8]])))
	stopifnot(args_line[[9]]=='-indir')
	stopifnot(args_line[[11]]=='-outdir')
	stopifnot(args_line[[13]]=='-jobtag')
	stopifnot(args_line[[15]]=='-trsm')
	stopifnot(args_line[[17]]=='-cmdstan')
	stopifnot(args_line[[19]]=='-max_index_cases')
	stopifnot(args_line[[21]]=='-start_d')
	stopifnot(args_line[[23]]=='-end_d')
	stopifnot(args_line[[25]]=='-index_flag')
	stopifnot(args_line[[27]]=='-infdate')
	stopifnot(args_line[[29]]=='-nonB')
	stopifnot(args_line[[31]]=='-sensitivity_infdate')
	stopifnot(args_line[[33]]=='-undiag_job')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['seed']] <- as.integer(args_line[[6]])
	args[['chain']] <- as.integer(args_line[8])
	args[['indir']] <- args_line[[10]]
	args[['outdir']] <- args_line[[12]]
	args[['job_tag']] <- args_line[[14]]
	args[['trsm']] <- args_line[[16]]
	args[['cmdstan']] <- as.integer(args_line[[18]])  
	args[['max_index_cases']] <- as.integer(args_line[[20]])  
	args[['start_d']] <- as.integer(args_line[[22]])  
	args[['end_d']] <- as.integer(args_line[[24]])  
	args[['index_flag']] <- as.integer(args_line[[26]])  
	args[['infdate']] <- as.integer(args_line[[28]])  
	args[['nonB']] <- args_line[[30]]  
	args[['sensitivity_infdate']] <- args_line[[32]]  
	args[['undiag_job']] <- args_line[[34]]  
} 

## load functions
source(file.path(args$source_dir, 'R', 'stan-functions.r' ))

## set other args
args$file_stanModel <- file.path(args$source_dir, 'stan-models',paste0(args$stanModelFile,'.stan'))
args$upper.bound.multiplier <- 10
tmp <- Sys.getenv("PBS_JOBID")
args$job_id <- ifelse(tmp!='', tmp, as.character(abs(round(rnorm(1) * 1e6))) )
args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag,'-',args$job_id)) 
args$savedata <- TRUE
if(grepl("\\[1\\]", args$job_id)) args$savedata = TRUE

## determine other args
args$index_in_pre = 0

str(args)
set.seed(args$seed)

## make job dir
dir.create( args$job_dir )

## save input args
saveRDS( args, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_args.RDS')))

cat('\n Set input file paths \n')

indir.phsc <- file.path(args$source_dir, 'data', 'subgraphs') 
infile.suppressed <- file.path(args$source_dir,'data','patient_data',paste0('suppressed_',args$start_d,'.csv'))
infile.inftime <- file.path(args$indir,'Data','infection_time_estimates','roadmap_cd4_vl_est.csv')

# input files using infection date
infile.infdate <- file.path(args$source_dir,'data','patient_data',paste0('infdate_after_',args$start_d,'.csv'))
infile.chaintype.infd <- file.path(args$source_dir,'data','subgraph_metadata',paste0('subgraph_classification_infdate_',args$start_d,'.csv'))
infile.diagnoses.i <- file.path(args$source_dir,'data','patient_data',paste0('N_diagnosed_by_mg_infdate_since_',args$start_d,'.csv'))

# input files using diagnosis date
infile.diagdate <- file.path(args$source_dir,'data','patient_data',paste0('diagdate_after_',args$start_d,'.csv'))
infile.chaintype.diagd <- file.path(args$source_dir,'data','subgraph_metadata',paste0('subgraph_classification_diagdate_',args$start_d,'.csv'))
infile.diagnoses.d <- file.path(args$source_dir,'data','patient_data',paste0('N_diagnosed_by_mg_diagdate_since_',args$start_d,'.csv'))

cat('\n Extracting taxa from phylogenetic subgraphs \n')
dsubgraphtaxa <- extract_subgraphs(indir.phsc)

cat('\n Adding suppression status from viral load data \n')

dsupp <- read.csv(infile.suppressed,header=T)
dsubgraphtaxa <- merge(dsubgraphtaxa,dsupp,by='ID',all.x=T)

cat('\n Adding infection date information \n')

if(args$infdate==1){
	dinf <- read.csv(infile.infdate,header=T)
}else{
	dinf <- read.csv(infile.diagdate,header=T)
}
dsubgraphtaxa <- merge(dsubgraphtaxa,dinf,by='ID',all.x=T)

cat(paste0('\n Calculating sequence sampling fraction for index cases at start of ',args$start_d, '\n'))

outfile.rho.index <- file.path(args$source_dir,'data','patient_data',paste0('sampling_fraction_index_cases_',args$start_d,'.csv'))
rho_m <- read.csv(outfile.rho.index, header=T)
rho_m <- subset(rho_m, TRANSM==args$trsm)

cat(paste0('\n Load subgraph classifications by ', args$start_d, '\n'))

# label phylogenetic subgraphs by whether they pre-existed by start date, or emerged since start date
if(args$infdate==1){
	dc <- read.csv(infile.chaintype.infd, header=T)
}else{
	dc <- read.csv(infile.chaintype.diagd, header=T)
}

cat(paste('\n Summarising subgraph size distribution \n'))

cases <- subgraph_sizes(dsubgraphtaxa,rho_m,dc,args$start_d,args$end_d,args$trsm,args$index_in_pre) 

if(args$nonB==1){
	cat(paste('\n Grouping non-B subtypes together \n'))
	
	cases[,ST2:=ST]
	dsubgraphtaxa[ST!='B', ST:='nonB']
	cases[ST!='B', ST:='nonB']
}

# save subgraph data
outfile <- file.path(args$outdir,'subgraphs_withmetadata.RDS')
saveRDS(dsubgraphtaxa, file=outfile)

freqs <- cases[, list(N=length(NAME)), by=c('SELECT','ST','icases','jcases','SIZE')]
freqs2 <- cases[, list(N=length(NAME)), by=c('SELECT','ST2','icases','jcases','SIZE')]

# plot for any aggregated subtypes
freqs[, `Index cases`:=icases]
freqs2[, `Index cases`:=icases]
plot_cs <- ggplot(data=freqs) +
	geom_bar(aes(x=jcases,y=N,fill=ST),stat='identity', position = "dodge") +
	scale_x_continuous(expand=c(0,0),breaks=c(seq(0,max(freqs$N),1)),labels=c(seq(0,max(freqs$N),1))) +
  labs(x='\nNumber of new cases',y="Number of chains",fill="Subtype") +
	theme_bw() +
	ggsci::scale_fill_npg() +
	facet_wrap(.~`Index cases`,ncol=1,scales="free_y",labeller = label_both) +
	theme(legend.position="bottom",
				strip.background=element_blank())
ggsave(file=file.path(args$job_dir, paste0('chainsizes',args$trsm,'.png')),plot_cs,w=8, h=12)

# plot for all subtypes
freqs2 <- subset(freqs2,icases>=0)
freqs2[,  `Index cases`:=factor(`Index cases`,levels=seq(0,max(icases),1))]
freqs2$ST2 <- factor(freqs2$ST2,levels=c('B','01AE','02AG','C','A1','G','D','06cpx'))
plot_cs <- ggplot(data=freqs2) +
	geom_bar(aes(x=jcases,y=N,fill=ST2),stat='identity', position = "dodge") +
	scale_x_continuous(expand=c(0,0),breaks=c(seq(0,max(freqs$N),1)),labels=c(seq(0,max(freqs$N),1))) +
	labs(x='Number of new cases',y="Number of chains",fill="Subtype") +
	theme_bw(30) +
	ggsci::scale_fill_npg() +
	facet_grid(ST2~`Index cases`,scales="free",labeller = label_both) +
	theme(legend.position="bottom",
				strip.background=element_blank())
ggsave(file=file.path(args$job_dir, paste0('chainsizes_sbt_',args$trsm,'.png')),plot_cs,w=30, h=30)
freqs[,  `Index cases`:=NULL]

# drop all with icases== -1 (pre-existing subgraphs with no observed index cases and no new cases, assumed to have died out)
freqs <- subset(freqs,icases>=0)

cat("\nSaving Subgraphs data\n")
if(args$savedata) saveRDS(freqs,file.path(args$job_dir, paste0('subgraph_sizes_',args$trsm,'.RDS')))

cat("\nLoading undiagnosed data \n")
file <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('p_undiagnosed_average_',args$undiag_job,'_',args$trsm,'.RDS'))

du <- readRDS(file)
dl <- data.table(trsm=c(rep('MSM',5),rep('HSX',4)),mg=c(1,2,3,4,5,1,2,3,4),mwmb=c('G1','G2','G3','NL','Other','G4','G5','NL','Other'))
dl <- subset(dl,trsm==args$trsm)
p_undiagnosed <- merge(dl,subset(du,select=c('trsm','mg','p0.5')),by=c('trsm','mg'),all.x=T)
setnames(p_undiagnosed,'p0.5','du')

cat("\nLoading number diagnosed by risk group \n")

if(args$infdate==1){
	pr <- data.table(read.csv(infile.diagnoses.i, header=T))
}else{
	pr <- data.table(read.csv(infile.diagnoses.d, header=T))
}
pr <- subset(pr,TRANSM==args$trsm)

cat("\nCalculating sequence sampling fraction \n")

rho <- calculate_sampling_fraction_mwmb(freqs,pr,args$infdate,args$start_d,args$end_d,p_undiagnosed)

seq <- rho[[1]]
inf <- rho[[2]]

cat("\nGenerate Stan data \n")

stan.data = generate_stan_data(freqs, seq, inf, p_undiagnosed, args$start_d, args$end_d, max_icases = args$max_index_cases, max_jcases = 37,args$upper.bound.multiplier,args$index_flag,args$infdate)

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_stanin.RData')) )

## init values
stan_init <- list()
r0 <- rep(c(0.5,0.25,0.75), 3) 
stan_init$r0 <- r0[args$chain]
stan_init$vmr_minus_one <- 0.05
stan_init$rho <- 0.5
stan_init$r0 <- rep(stan_init$r0, stan.data$N_sbts)

# save stan.data object
if(args$cmdstan){
  rstan::stan_rdump( names(stan_init), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstaninit.R')), envir=list2env(stan_init))
  rstan::stan_rdump( names(stan.data), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstanin.R')), envir=list2env(stan.data))  	
} else{
  model <- rstan::stan_model(args$file_stanModel)
  fit <- rstan::sampling(model,data=stan.data,iter=10,warmup=5,chains=1,seed=args$seed,init=list(stan_init),verbose=TRUE)
}



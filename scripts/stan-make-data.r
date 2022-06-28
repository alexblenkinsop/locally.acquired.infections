library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
library(ggsci)

args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	#source_dir= '~/Documents/GitHub/bpm',
	#indir='~/Box\ Sync/Roadmap/RQ1 Estimating introductions',
	#outdir= '~/Box\ Sync/Roadmap/RQ1 Estimating introductions/branching_process_model',
	stanModelFile= 'branching_process_210810b_cmdstan',
	analysis= 'analysis_211101',
	#analysis= 'analysis_200917',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	job_tag= 'ECDC_undiagnosed',
	trsm= 'MSM',
	cmdstan = 1L,
	max_index_cases = 22,
	start_d = 2014,
	end_d = 2019,
	index_flag = 1,
	upper.bound.multiplier = 10,
	last_case = 10,
	infdate=1,
	nonB=1,
	sensitivity_infdate=0,
	undiag_job='2010_2012_notrunc'
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-analysis')
	stopifnot(args_line[[7]]=='-seed', !is.na(as.integer(args_line[[8]])))	
	stopifnot(args_line[[9]]=='-chain', !is.na(as.integer(args_line[[10]])))
	stopifnot(args_line[[11]]=='-indir')
	stopifnot(args_line[[13]]=='-outdir')
	stopifnot(args_line[[15]]=='-jobtag')
	stopifnot(args_line[[17]]=='-trsm')
	stopifnot(args_line[[19]]=='-cmdstan')
	stopifnot(args_line[[21]]=='-max_index_cases')
	stopifnot(args_line[[23]]=='-start_d')
	stopifnot(args_line[[25]]=='-end_d')
	stopifnot(args_line[[27]]=='-index_flag')
	stopifnot(args_line[[29]]=='-infdate')
	stopifnot(args_line[[31]]=='-nonB')
	stopifnot(args_line[[33]]=='-sensitivity_infdate')
	stopifnot(args_line[[35]]=='-undiag_job')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['seed']] <- as.integer(args_line[[8]])
	args[['chain']] <- as.integer(args_line[[10]])
	args[['indir']] <- args_line[[12]]
	args[['outdir']] <- args_line[[14]]
	args[['job_tag']] <- args_line[[16]]
	args[['trsm']] <- args_line[[18]]
	args[['cmdstan']] <- as.integer(args_line[[20]])  
	args[['max_index_cases']] <- as.integer(args_line[[22]])  
	args[['start_d']] <- as.integer(args_line[[24]])  
	args[['end_d']] <- as.integer(args_line[[26]])  
	args[['index_flag']] <- as.integer(args_line[[28]])  
	args[['infdate']] <- as.integer(args_line[[30]])  
	args[['nonB']] <- args_line[[32]]  
	args[['sensitivity_infdate']] <- args_line[[34]]  
	args[['undiag_job']] <- args_line[[36]]  
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
args$with_subtypes = 0
if(grepl("1209|1215|0203|0205|0210|0215|0217|0322|0414|0505|0802|0810", args$file_stanModel)){
  args$with_subtypes = 1
}
args$index_in_pre = 0

str(args)
set.seed(args$seed)

## make job dir
dir.create( args$job_dir )


## save input args
saveRDS( args, file=file.path(args$job_dir, paste0(basename(args$job_dir), '_args.RDS')))


cat("\nRead phylogenetic subgraphs \n")

indir.phsc <- file.path(args$indir, args$analysis, 'subgraphs_dated')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '200917_sequence_labels.rda')
infile.seq <-	file.path(args$indir, 'Data', 'data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
infile.geo <- file.path(args$indir,'misc','NEWGEO.csv')

cat('\n Extract subgraph taxa\n')
dsubgraphtaxa <- extract_subgraphs(indir.phsc)

cat("\nRead patient metadata \n")
load(infile.seq)
load(infile.meta)
setnames(dind, gsub("CENS_D", "RECART_D", names(dind)))

# flag which diagnosed individuals have a sequence
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
dind <- merge(dind,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)

# add meta data from pre-processed persons file
dind <- as.data.table(dind)
setnames(dind, c('PATIENT','BIRTH_Y','BIRTH_CNTRY'), c('ID','BIRTH_YEAR','BIRTH_COUNTRY'))
tmp <- subset(dind, select=c(ID,BIRTH_YEAR,BIRTH_COUNTRY,LOC_BIRTH,CITY,TRANSM,GENDER,RECART_D,HIV1_POS_D,HIV1_POS_D_lower,HIV1_POS_D_upper))
set(tmp, NULL, 'BIRTH_YEAR', tmp[, as.numeric(BIRTH_YEAR)])
dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')

cat('\n Add viral load data \n')
infile.cd4 <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_CD4.csv')
du <- add_cd4_counts(dsubgraphtaxa,infile.cd4,args$start_d,args$trsm)

infile.rna <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_RNA.csv')
infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
dsubgraphtaxa <- add_viral_loads(dsubgraphtaxa,infile.rna,infile.indinfo,args$start_d)
dind <- add_viral_loads(dind,infile.rna,infile.indinfo,args$start_d)

cat('\n Correct misrecorded dates of diagnosis/ART start \n')
dsubgraphtaxa <- correct_misrecorded_dates(dsubgraphtaxa,infile.indinfo)


cat('\n Add birth place region \n')
dsubgraphtaxa <- add_georeg_bplace_sg_data(dsubgraphtaxa,dind,infile.geo)
dsubgraphtaxa <- map_mwmb_regions(dsubgraphtaxa)
dind <- add_georeg_bplace_ind_data(dind,infile.geo)
dind <- map_mwmb_regions(dind)

cat('\n Update HIV infection date using VL data \n')
dt_m <- file.path(args[['indir']],'branching_process_model',paste0('undiagnosed_211102-',args[['undiag_job']]),paste0('median_timetodiagnosis_',args[['undiag_job']],'_MSM.RDS'))
dt_h <- file.path(args[['indir']],'branching_process_model',paste0('undiagnosed_211102-',args[['undiag_job']]),paste0('median_timetodiagnosis_',args[['undiag_job']],'_HSX.RDS'))
dt_m <- readRDS(dt_m)
dt_h <- readRDS(dt_h)
dt_inf <- rbind(dt_m,dt_h)
setnames(dt_inf,'trsm','TRANSM')

infile.inftime <- file.path(args$indir,'Data','infection_time_estimates','roadmap_cd4_vl_est.csv')
if(args$sensitivity_infdate==1){
	dsubgraphtaxa <- add_infection_time_midpoint_seroconv_diagnosis(dsubgraphtaxa,infile.inftime)
}else{
	dsubgraphtaxa <- add_infection_time(dsubgraphtaxa,infile.inftime,dt_inf)
}
dind <- add_infection_time(dind,infile.inftime,dt_inf)

cat(paste('\n Summarise subgraph size distribution \n'))
subgraphs_inf <- subgraph_sizes_infdate(dsubgraphtaxa,dind,args$start_d,args$end_d,args$trsm,args$index_in_pre) 
subgraphs_diag <- subgraph_sizes_diagdate(dsubgraphtaxa,dind,args$start_d,args$end_d,args$trsm,args$index_in_pre)

cases_inf <- subgraphs_inf[[1]]
cases_diag <- subgraphs_diag[[1]]

dind_inf <- subgraphs_inf[[2]]
dind_diag <- subgraphs_diag[[2]]
	
if(args$infdate==1){
	cases <- cases_inf
	dind <- dind_inf
}else{
	cases <- cases_diag
	dind <- dind_diag
}

if(args$nonB==1){
	cases[,ST2:=ST]
	dsubgraphtaxa[ST!='B', ST:='nonB']
	dind[ST!='B', ST:='nonB']
	cases[ST!='B', ST:='nonB']
}

# save subgraph data
outfile <- file.path(args$outdir,'subgraphs_withmetadata.RDS')
saveRDS(dsubgraphtaxa, file=outfile)

if(args$with_subtypes){
  freqs <- cases[, list(N=length(NAME)), by=c('SELECT','ST','icases','jcases','SIZE')]
  freqs2 <- cases[, list(N=length(NAME)), by=c('SELECT','ST2','icases','jcases','SIZE')]
  
  freqs[, `Index cases`:=icases]
  freqs2[, `Index cases`:=icases]
  plot_cs <- ggplot(data=freqs) +
  	geom_bar(aes(x=SIZE,y=N,fill=ST),stat='identity', position = "dodge") +
  	scale_x_continuous(expand=c(0,0),breaks=c(seq(0,max(freqs$N),1)),labels=c(seq(0,max(freqs$N),1))) +
  	#scale_y_continuous(breaks = integer_breaks(),expand = expansion(mult = c(0, .1))) +
    labs(x='\nNumber of new cases',y="Number of chains",fill="Subtype") +
  	theme_bw() +
  	ggsci::scale_fill_npg() +
  	facet_wrap(.~`Index cases`,ncol=1,scales="free_y",labeller = label_both) +
  	theme(legend.position="bottom",
  				strip.background=element_blank())
  ggsave(file=file.path(args$job_dir, paste0('chainsizes',args$trsm,'.png')),plot_cs,w=8, h=12)
  freqs2[,  `Index cases`:=factor(`Index cases`,levels=seq(0,max(icases),1))]
  freqs2$ST2 <- factor(freqs2$ST2,levels=c('B','01AE','02AG','C','A1','G','D','06cpx'))
  plot_cs <- ggplot(data=freqs2) +
  	geom_bar(aes(x=jcases,y=N,fill=ST2),stat='identity', position = "dodge") +
  	scale_x_continuous(expand=c(0,0),breaks=c(seq(0,max(freqs$N),1)),labels=c(seq(0,max(freqs$N),1))) +
  	#scale_y_continuous(breaks = integer_breaks(),expand = expansion(mult = c(0, .1))) +
  	labs(x='Number of new cases',y="Number of chains",fill="Subtype") +
  	theme_bw(30) +
  	ggsci::scale_fill_npg() +
  	facet_grid(ST2~`Index cases`,scales="free",labeller = label_both) +
  	theme(legend.position="bottom",
  				strip.background=element_blank())
  ggsave(file=file.path(args$job_dir, paste0('chainsizes_sbt_',args$trsm,'.png')),plot_cs,w=30, h=30)
  freqs[,  `Index cases`:=NULL]
  
  # drop all with icases== -1
  	freqs <- subset(freqs,icases>=0)
  # set m=1 index case and combine with other m=1,j=0 subgraph counts
}else{
  freqs <- cases[, list(N=length(NAME)), by=c('SELECT','jcases','icases','SIZE')]
  plot_cs  <- ggplot(data=freqs) +
  	geom_bar(aes(x=SIZE,y=N),stat='identity', position = "dodge") +
  	scale_x_continuous(expand=c(0,0),breaks=c(seq(0,max(freqs$N),1)),labels=c(seq(0,max(freqs$N),1))) +
  	labs(x='\nSize of likely transmission chains \nattributed to Amsterdam',y="Number of chains") +
  	facet_wrap(.~icases,ncol=1,scales="free_y",labeller = label_both) +
  	theme_bw() +
  	theme(legend.position="bottom")
  ggsave(file=file.path(args$job_dir, paste0('chainsizes',args$trsm,'.png')),plot_cs,w=8, h=12)
}
cat("\nSave Subgraphs data\n")
if(args$savedata) saveRDS(freqs,file.path(args$job_dir, paste0('subgraph_sizes_',args$trsm,'.RDS')))

cat("\nLoad undiagnosed data \n")
file <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('p_undiagnosed_average_',args$undiag_job,'_',args$trsm,'.RDS'))

du <- readRDS(file)
dl <- data.table(trsm=c(rep('MSM',5),rep('HSX',4)),mg=c(1,2,3,4,5,1,2,3,4),mwmb=c('G1','G2','G3','NL','Other','G4','G5','NL','Other'))
dl <- subset(dl,trsm==args$trsm)
p_undiagnosed <- merge(dl,subset(du,select=c('trsm','mg','p0.5')),by=c('trsm','mg'),all.x=T)
setnames(p_undiagnosed,'p0.5','du')

tmp <- copy(freqs)
rho <- calculate_sampling_fraction_mwmb(infile.geo,dind,tmp,args$infdate,args$start_d,args$end_d,args$trsm,p_undiagnosed)

seq <- rho[[1]]
inf <- rho[[2]]

#
# Generate stan data
cat("\nGenerate Stan data \n")

if(args$with_subtypes){
  stan.data = generate_stan_data_sbt(args$trsm,tmp, dind, seq, inf, p_undiagnosed, args$start_d, args$end_d, max_icases = args$max_index_cases, max_jcases = 37,args$upper.bound.multiplier,args$index_flag,args$infdate)
} else{
  stan.data = generate_stan_data(args$trsm,tmp, dind, args$start_d, args$end_d, max_icases = args$max_index_cases, max_jcases = 37,args$upper.bound.multiplier)
}

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
if(args$with_subtypes){
  stan_init$r0 <- rep(stan_init$r0, stan.data$N_sbts)
} 

# save stan.data object
if(args$cmdstan){
  rstan::stan_rdump( names(stan_init), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstaninit.R')), envir=list2env(stan_init))
  rstan::stan_rdump( names(stan.data), file=file.path(args$job_dir, paste0(basename(args$job_dir), '_cmdstanin.R')), envir=list2env(stan.data))  	
} else{
  model <- rstan::stan_model(args$file_stanModel)
  fit <- rstan::sampling(model,data=stan.data,iter=10,warmup=5,chains=1,seed=args$seed,init=list(stan_init),verbose=TRUE)
}



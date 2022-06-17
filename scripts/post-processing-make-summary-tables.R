# make-tables.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Make tables.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810a_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810a_cmdstan-bugfix_cases_bplace_2014-2019_HSX'
	args_dir[['job_tag']] <- 'bugfix_cases_bplace_2014-2019_HSX'
	args_dir[['trsm']] <- 'HSX'
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
source(file.path(args_dir$source_dir, 'R', 'post-processing-summary-functions.R'))

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
pars.basic <- readRDS(file)
file <- paste0(outfile.base,'-stanout-subgraphs-gqs.RDS')
cat("\n read RDS:", file)
fit.subgraphs <- readRDS(file)
file <- paste0(outfile.base,'-stanout-origins-gqs.RDS')
cat("\n read RDS:", file)
fit.origins <- readRDS(file)
file <- paste0(outfile.base,'-stanout-transmission_pars.RDS')
cat("\n read RDS:", file)
plot.pars.trmspars <- readRDS(file)

cat(" \n -------------------------------- Reading data -------------------------------- \n")

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)



cat('\nSummarise average size of subgraph...')

## Average size of subgraph
cs_obs <- as.data.table(reshape2:: melt(stan.data$cs_obs))
setnames(cs_obs, 1:4, c('index_cases','new_cases','subtype','freq'))
# shift index cases for m=0 row
cs_obs$index_cases <- cs_obs$index_cases - 1
# if index case assumed part of subgraph make 1 new case the index case
#if(stan.data$index_flag==1){
#	cs_obs[index_cases==0,new_cases:= as.integer(new_cases - 1)]
#}
# correct generated cases as col 1 = 0
cs_obs[,new_cases:= new_cases - 1]
# only keep non-negative cases (shouldn't have any)
cs_obs <- subset(cs_obs,new_cases>=0)
# correct index cases for subgraphs starting since 2010 to 1
#cs_obs[index_cases==0,index_cases:=1]
# aggregate subgraphs post-2010 with the rest
cs_obs <- cs_obs[, list(freq=sum(freq)),by=c('subtype','index_cases','new_cases')]
cs_obs[, N:=new_cases*freq]

cat('\nSummarise predicted chains...')
cs_actual <- readRDS(file=paste0(outfile.base,'-','predicted_chains.rds'))
unobs_N <- readRDS(file=paste0(outfile.base,'-','unobserved_chains.rds'))

#dall <- rbind(subset(cs_actual,select=c('iteration','subtype','size')),subset(unobs_N,select=c('iteration','subtype','size')))
##du <- unobs_N$size[unobs_N$subtype]
#mean_act <- dall[, list(size = mean(size)),by=c('subtype')]
dall <- rbind(subset(cs_actual,select=c('iteration','subtype','new_cases')),subset(unobs_N,select=c('iteration','subtype','new_cases')))
mean_act <- dall[, list(new_cases = mean(new_cases)),by=c('subtype')]
#mean_a <- dall[, list(mean_cases= quantile(new_cases, prob=ps,na.rm=T),
#													q_label=p_labs),by=c('subtype')]		

dc <- cs_actual[, list(N=length(subgraph)),by=c('iteration','subtype')]
du <- unobs_N[, list(N_u=length(subgraph)),by=c('iteration','subtype')]
dc <- merge(dc,du,by=c('iteration','subtype'),all=T)
dc[is.na(N), N:=0]
dc[is.na(N_u), N_u:=0]
dc[, chains := N+N_u]
dc <- merge(dc,pars.basic$ds,by.x='subtype',by.y='subtypes')

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

# 	summarise
par <- dc[, list(qs= quantile(chains, prob=ps), qlab=p_labs), by=c('subtypes_name')]
par <- dcast.data.table(par, subtypes_name~qlab, value.var='qs')
par[, L:= paste0( round(M, d=0), ' [',  round(CL, d=0),'-', round(CU, d=0),']')]

cat('\nSummarise local infections by subtype...')

ex <- readRDS(file=paste0(outfile.base,'-','external_importations_sbt_samples.RDS'))
ex <- ex[, list(ext_imp= quantile(ext_imp, prob=ps,na.rm=T),inf_Ams = quantile(inf_Ams, prob=ps,na.rm=T),
								q_label=p_labs),by=c('TRM','subtype')]		

ans <- dcast.data.table(ex, TRM+subtype~q_label, value.var='inf_Ams')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,pars.basic$ds,by.x='subtype',by.y='subtypes')

cat('\nSummarise Rt by subtype...')

## Rt
Rt <- exp(plot.pars.trmspars$log_r0_sbts)
Rt <- make_posterior_summaries(Rt,pars.basic$ds)
Rt[, L:= paste0( round(M, d=2), ' [',  round(CL, d=2),'-', round(CU, d=2),']')]

vmr <- plot.pars.trmspars$vmr_minus_one + 1
vmr <- make_posterior_summaries(vmr,pars.basic$ds)
vmr[, L:= paste0( round(M, d=2), ' [',  round(CL, d=2),'-', round(CU, d=2),']')]

cat(" \n -------------------------------- Create table -------------------------------- \n")

sbts <- unique(pars.basic$ds$subtypes_name)
	
sbt <- list()
for(st in sbts){
	stid <- pars.basic$ds$subtypes[pars.basic$ds$subtypes_name==st]
	
	# mean size of obs subgraphs
	ds <- subset(cs_obs,subtype==pars.basic$ds$subtypes[pars.basic$ds$subtypes_name==st])
	#ds[, size:=index_cases+new_cases]
	#d <- rep(ds$size, ds$freq)
	d <- rep(ds$new_cases, ds$freq)
	
	# mean size of actual chains
	#mean_a <- mean_act$size[mean_act$subtype==stid]
	mean_a <- mean_act$new_cases[mean_act$subtype==stid]
	
	sbt[[st]] <- data.table(subtype=st,
													obs_cases=sum(cs_obs$N[cs_obs$subtype==stid]),
									 sgs_n=sum(stan.data$cs_obs[,,stid]),
									 #sgs_size=round(mean(d),2),
									 sgs_newcases=round(mean(d),2),
									 chains_n=par[subtypes_name==st,L],
									 #chains_size=round(mean_a,2),
									 chains_newcases=round(mean_a,2),
									 Rt=Rt[subtypes_name==st,L],
									 vmr=vmr[subtypes_name==st,L],
									 local=ans[subtypes_name==st,L]
										)
}

dt <- do.call('rbind',sbt)

setnames(dt, c('subtype','obs_cases','sgs_n','sgs_newcases','chains_n','chains_newcases','local'),
				 c('Subtype','Observed new cases','Observed subgraphs','Average observed new cases','Actual transmission chains','Average actual new cases','Infections acquired in Amsterdam (%)'))
dt$Subtype <- factor(dt$Subtype,levels=c('B','nonB'),labels=c('B','Non-B'))

saveRDS(dt,file=paste0(outfile.base,'-','model_summary_table.RDS'))
write.csv(dt,file=paste0(outfile.base,'-','model_summary_table.csv'))

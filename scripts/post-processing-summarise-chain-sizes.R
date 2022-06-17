
cat(" \n -------------------------------- \n \n Summarise-chain-sizes.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'undiagnosed_notrunc'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
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

file <- paste0(outfile.base,'-stanout-subgraphs-gqs.RDS')
cat("\n read RDS:", file)
fit.subgraphs <- readRDS(file)

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

#### obs
cs_obs <- as.data.table(reshape2:: melt(stan.data$cs_obs))
setnames(cs_obs, 1:4, c('index_cases','new_cases','subtype','freq'))
# shift index cases for m=0 row
cs_obs$index_cases <- cs_obs$index_cases - 1
# distinguish the emergent chains
cs_obs$chains <- 'pre-existing'
cs_obs[index_cases==0,chains:='emergent']
# correct generated cases as col 1 = 0
cs_obs[,new_cases:= new_cases - 1]
# only keep non-negative cases (shouldn't have any)
cs_obs <- subset(cs_obs,new_cases>=0)
# aggregate subgraphs post-2010 with the rest
cs_obs <- cs_obs[, list(freq=sum(freq)),by=c('subtype','index_cases','new_cases','chains')]

obs <- cs_obs[, list(QN=sum(freq)),by=c('new_cases','chains')]
obs_t <- cs_obs[, list(tot=sum(freq),tot_grew=sum(freq[new_cases>0])),by=c('chains')]
obs <- merge(obs,obs_t,by=c('chains'),all.x=T)
obs[, p:=QN/tot]
obs[,analysis:= 'observed']
obs[,P:= 'p0.5']

obs_N <- obs[, list(QN_c=sum(QN)), by=c('chains')]
obs_N[, analysis:='observed']
obs_N[, P:='p0.5']
obs_i <- obs[, list(QN_i=sum(QN*new_cases)), by=c('chains')]
obs_i[, analysis:='observed']
obs_i[, P:='p0.5']

cat('\nLoad predicted chain sizes...')
cs_ex <- readRDS(file=paste0(outfile.base,'-','predicted_chains.rds'))
cs_em <- readRDS(file=paste0(outfile.base,'-','unobserved_chains.rds'))
dt <- rbind(cs_ex,cs_em)
dt <- data.table(dt)
dt[, analysis:="predicted"]
N_ind <- dt[, list(N=sum(new_cases)),by=c('iteration','chains')]

cat('\nSummarise...')

cs_act <- dt[, list(freq=length(subgraph)),by=c('iteration','new_cases','chains')]
cs_act_tot <- dt[, list(freq_tot=length(subgraph),freq_tot_grew=length(subgraph[new_cases>0])),by=c('iteration','chains')]
cs_act <- merge(cs_act,cs_act_tot,by=c('iteration','chains'),all.x=T)
cs_act[, p:=freq/freq_tot]

total.replicates	<- cs_act[, length(unique(iteration))]
ans	<- cs_act[, list(P=ps, 
										 QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ps),
										 tot=quantile(c(freq_tot, rep(0, total.replicates-length(freq_tot))), p=ps),
										 tot_grew=quantile(c(freq_tot_grew, rep(0, total.replicates-length(freq_tot_grew))), p=ps),
										 QN_p=quantile(c(p, rep(0, total.replicates-length(p))), p=ps)),
										 by=c('chains','new_cases')]
set(ans, NULL, 'P', ans[, paste0('p',P)])
setnames(ans,'QN_p','p')
ans[,analysis := 'predicted']

# summarise number of chains
N_c <- dt[, list(N=length(subgraph)),by=c('iteration','chains')]
N_c <- N_c[, list(P=ps, 
									QN_c=quantile(c(N, rep(0, total.replicates-length(N))), p=ps)), by=c('chains')]
set(N_c, NULL, 'P', N_c[, paste0('p',P)])

N_ind <- N_ind[, list(P=ps, 
											QN_i=quantile(c(N, rep(0, total.replicates-length(N))), p=ps)), by=c('chains')]
set(N_ind, NULL, 'P', N_ind[, paste0('p',P)])
N_ind[, analysis:='predicted']
N_ind <- merge(N_ind,obs_i,by=c('chains','analysis','P','QN_i'),all=T)
N_c[, analysis:='predicted']
N_c <- merge(N_c,obs_N,by=c('chains','analysis','P','QN_c'),all=T)

### merge
all <- merge(ans,subset(obs,select=c('chains','analysis','new_cases','QN','tot','tot_grew','p','P')),by=c('chains','analysis','new_cases','tot','tot_grew','p','P','QN'),all=T)
all <- merge(all,N_ind,by=c('chains','analysis','P'),all=T)

saveRDS(all,file=paste0(outfile.base,'-','obs_actual_cs_distribution_long_all.rds'))

cat('\nGroup chains of 7+ new cases together...')

cs_obs[, newcases:=factor(new_cases)]
cs_obs[new_cases>=7,newcases:='7+']
cs_obs[, new_cases:=NULL]
setnames(cs_obs,'newcases','new_cases')

obs <- cs_obs[, list(QN=sum(freq)),by=c('new_cases','chains')]
obs_t <- cs_obs[, list(tot=sum(freq),tot_grew=sum(freq[new_cases!='0'])),by=c('chains')]
obs <- merge(obs,obs_t,by=c('chains'),all.x=T)
obs[, p:=QN/tot]
obs[,analysis:= 'observed']
obs[,P:= 'p0.5']
obs_N <- obs[, list(QN_c=sum(QN)), by=c('chains')]
obs_N[, analysis:='observed']
obs_N[, P:='p0.5']

dt[, newcases:=factor(new_cases)]
dt[new_cases>=7,newcases:='7+']
dt[, new_cases:=NULL]
setnames(dt,'newcases','new_cases')

cs_act <- dt[, list(freq=length(subgraph)),by=c('iteration','new_cases','chains')]
cs_act_tot <- dt[, list(freq_tot=length(subgraph),freq_tot_grew=length(subgraph[new_cases!='0'])),by=c('iteration','chains')]
cs_act <- merge(cs_act,cs_act_tot,by=c('iteration','chains'),all.x=T)
cs_act[, p:=freq/freq_tot]

total.replicates	<- cs_act[, length(unique(iteration))]
ans	<- cs_act[, list(P=ps, 
										 QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ps),
										 tot=quantile(c(freq_tot, rep(0, total.replicates-length(freq_tot))), p=ps),
										 tot_grew=quantile(c(freq_tot_grew, rep(0, total.replicates-length(freq_tot_grew))), p=ps),
										 QN_p=quantile(c(p, rep(0, total.replicates-length(p))), p=ps)),
							by=c('chains','new_cases')]
set(ans, NULL, 'P', ans[, paste0('p',P)])
setnames(ans,'QN_p','p')

ans[,analysis := 'predicted']

# summarise number of chains
N_c <- dt[, list(N=length(subgraph)),by=c('iteration','chains')]
N_c <- N_c[, list(P=ps, 
									QN_c=quantile(c(N, rep(0, total.replicates-length(N))), p=ps)), by=c('chains')]
set(N_c, NULL, 'P', N_c[, paste0('p',P)])
N_c[, analysis:='predicted']
N_c <- merge(N_c,obs_N,by=c('chains','analysis','P','QN_c'),all=T)

### merge
all <- merge(ans,subset(obs,select=c('chains','analysis','new_cases','QN','tot','tot_grew','p','P')),by=c('chains','analysis','new_cases','tot','tot_grew','p','P','QN'),all=T)
all <- merge(all,N_ind,by=c('chains','analysis','P'),all=T)

saveRDS(all,file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))


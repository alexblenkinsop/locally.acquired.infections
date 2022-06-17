cat(" \n -------------------------------- \n \n Posterior predictive check.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810f_cmdstan'
	args_dir[['job_name']] <- 'lower_r0_prior'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['period']] <- '2014-2019'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['with_subtypes']] <- 1

	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_',	args_dir[['trsm']])
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
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

file <- paste0(outfile.base,'-stanout-posteriorpredictivecheck-gqs.RDS')
cat("\n read RDS:", file)
fit.csobs <- readRDS(file)

cat(" \n -------------------------------- Reading data -------------------------------- \n")

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)


cat('\nLoading predicted observed chain sizes...')
existing <- fit.csobs$obs_cs_pre 
emergent <- fit.csobs$obs_cs_post 

existing <- as.data.table( reshape2::melt( existing ) )
setnames(existing, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))

emergent <- as.data.table( reshape2::melt( emergent ) )
setnames(emergent, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))

existing[,chains:= 'existing']
emergent[,chains:= 'emergent']


existing <- subset(existing,is.finite(new_cases))
existing <- subset(existing,new_cases>0)
existing[chains=='existing', new_cases:=  new_cases - 1]
existing[, size:= index_cases + new_cases]

emergent <- subset(emergent,is.finite(new_cases))
emergent <- subset(emergent,new_cases>0)
emergent[, size:= index_cases + new_cases]


cs_obs_pr <- rbind(existing,emergent)

cat('\nLoading predicted chain sizes...')
saveRDS(existing,file=paste0(outfile.base,'-','predicted_chains_existing.rds'))
saveRDS(emergent,file=paste0(outfile.base,'-','predicted_chains_emergent.rds'))

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

cat(" \n -------------------------------- Posterior predictive check: new cases -------------------------------- \n")

cat('\nSummarising observed subgraphs...')

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


cat('\nPlot by subtype...')

obs <- cs_obs[, list(QN=sum(freq)),by=c('subtype','new_cases','chains')]
obs[,analysis:= 'observed']
obs[,P:= 'p0.5']
obs[,N:= new_cases*QN]

N_obs <- obs[, list(N=sum(QN)),by=c('subtype','chains')]

# count frequencies
cs_act <- cs_obs_pr[, list(freq=length(subgraph)),by=c('iteration','subtype','new_cases','chains')]
cs_act[,N:= new_cases*freq]
cs_act[chains=='existing',chains := 'pre-existing']

# actual chain size distribution
total.replicates	<- cs_act[, length(unique(iteration))]
act	<- cs_act[, list(P=ps, 
										 QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ps)), by=c('subtype','new_cases','chains')]
set(act, NULL, 'P', act[, paste0('p',P)])
act[,analysis := 'predicted']

obs[, N:= NULL]
ans <- merge(obs,act,by=c('analysis','subtype','new_cases','chains','P','QN'),all=T)

tmp			<- dcast.data.table(ans, analysis+subtype+chains+new_cases~P, value.var='QN')

tmp <- merge(tmp,pars.basic$ds,by.x="subtype",by.y="subtypes")
tmp <- merge(tmp,N_obs,by=c('subtype','chains'))

tmp$subtypes_name <- as.factor(tmp$subtypes_name)
tmp$subtypes_name <- factor(tmp$subtypes_name,levels=c('B','nonB','02AG','A1','01AE','C','G','D','06cpx'))

ob <- subset(obs,select=c('subtype','chains','new_cases','QN'))
tmp <- merge(tmp,ob,by=c('subtype','chains','new_cases'),all=T)

tmp[is.na(QN),QN:=0]

tmp[p0.025<=QN & QN<=p0.975, in_range:=1]
tmp[p0.025>QN | QN>p0.975, in_range:=0]

pr <- tmp[, list(in_range=sum(in_range,na.rm=T),all=length(subtype[!is.na(in_range)])),by=c('subtype','chains')]
pr[, p:=in_range/all]
tmp <- merge(tmp,subset(pr,select=c('subtype','chains','p')),by=c('subtype','chains'),all.x=T)

#dat <- subset(tmp,subtypes_name!='06cpx')
#dat <- subset(dat,subtypes_name!='D')
#if(args$trsm=='MSM'){
#	dat <- subset(dat,subtypes_name!='G')
#}
#if(args$trsm=='HSX'){
#	dat <- subset(dat,subtypes_name!='01AE')
#}

dat <- tmp
dat$chains <- factor(dat$chains,levels=c('pre-existing','emergent'))

plot <- ggplot(data=dat) +
	geom_bar(aes(x=new_cases,y=p0.5,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=p0.025, ymax=p0.975,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(aes(12, N,label=paste("N= ", N, "p=", round(p,2)))) +
	xlim(-1,15) +
	scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
	facet_grid(subtypes_name~chains,scales="free_y") +
	labs(x='\nNumber of cases since 2015',y="Number of chains",fill="") +
	theme_bw(base_size=14) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_newcases_bysubtype_chaintype.png'),plot,w=10, h=12)


plot <- ggplot(data=subset(dat,analysis!='predicted')) +
	geom_bar(aes(x=new_cases,y=p0.5,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=p0.025, ymax=p0.975,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	geom_text(aes(12, 25,label=paste("N= ", N))) +
	xlim(-1,15) +
	scale_y_continuous(trans='log10',expand = expand_scale(mult = c(0, .1))) +
	ylim(0,max(dat$p0.5)) +
	facet_wrap(subtypes_name~chains,ncol=2,scales="free_y") +
	labs(x='\nNumber of cases since 2015',y="Number of chains",fill="") +
	theme_bw() +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_newcases_bysubtype_chaintype_natscale.png'),plot,w=8, h=12)

saveRDS(dat,paste0(outfile.base,'-posteriorpredictivecheck.RDS'))


cat(" \n -------------------------------- Re-do aggregating over subtypes -------------------------------- \n")

existing <- as.data.table( reshape2::melt( existing ) )
setnames(existing, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))

emergent <- as.data.table( reshape2::melt( emergent ) )
setnames(emergent, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))

existing[,chains:= 'existing']
emergent[,chains:= 'emergent']


existing <- subset(existing,is.finite(new_cases))
existing <- subset(existing,new_cases>0)
existing[chains=='existing', new_cases:=  new_cases - 1]
existing[, size:= index_cases + new_cases]

emergent <- subset(emergent,is.finite(new_cases))
emergent <- subset(emergent,new_cases>0)
emergent[, size:= index_cases + new_cases]


cs_obs_pr <- rbind(existing,emergent)


cat('\nSummarising observed subgraphs...')

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


cat('\nPlot by subtype...')

obs <- cs_obs[, list(QN=sum(freq)),by=c('new_cases','chains')]
obs[,analysis:= 'observed']
obs[,P:= 'p0.5']
obs[,N:= new_cases*QN]

N_obs <- obs[, list(N=sum(QN)),by=c('chains')]

# count frequencies
cs_act <- cs_obs_pr[, list(freq=length(subgraph)),by=c('iteration','new_cases','chains')]
cs_act[,N:= new_cases*freq]
cs_act[chains=='existing',chains := 'pre-existing']

# actual chain size distribution
total.replicates	<- cs_act[, length(unique(iteration))]
act	<- cs_act[, list(P=ps, 
										 QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ps)), by=c('new_cases','chains')]
set(act, NULL, 'P', act[, paste0('p',P)])
act[,analysis := 'predicted']

obs[, N:= NULL]
ans <- merge(obs,act,by=c('analysis','new_cases','chains','P','QN'),all=T)

tmp			<- dcast.data.table(ans, analysis+chains+new_cases~P, value.var='QN')

tmp <- merge(tmp,N_obs,by=c('chains'))

ob <- subset(obs,select=c('chains','new_cases','QN'))
tmp <- merge(tmp,ob,by=c('chains','new_cases'),all=T)

tmp[is.na(QN),QN:=0]

tmp[p0.025<=QN & QN<=p0.975, in_range:=1]
tmp[p0.025>QN | QN>p0.975, in_range:=0]

dat <- tmp
dat$chains <- factor(dat$chains,levels=c('pre-existing','emergent'))

saveRDS(dat,paste0(outfile.base,'-posteriorpredictivecheck_allsbts.RDS'))

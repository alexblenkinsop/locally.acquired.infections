
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
dsubgraphtaxa <- subset(dsubgraphtaxa,SELECT %in% c('AmsMSM','AmsHSX'))

dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
do <- dsubgraphtaxa
# count unique subgraphs with origin in each region (by clade too)
do <- do[,list(N=length(unique(FULL_NAME))),by=c('REP','ORIGIN','SELECT','ST','ST_CLADE')]
# aggregate over clades for B
do <- do[,list(N=sum(N)),by=c('REP','ORIGIN','SELECT','ST')]
do <- subset(do, ORIGIN!='Unknown')
do_0 <- do[REP=='000', list(ORIGIN=ORIGIN,N=N,pct=round(N/sum(N)*100,1)),by=c('REP','SELECT','ST')]
do <- do[REP!='000', list(ORIGIN=ORIGIN,N=N,pct=N/sum(N)),by=c('REP','SELECT','ST')]
da <- do[, list(N= round(quantile(N, prob=ps,na.rm=T),0),
								pct= round(quantile(pct, prob=ps,na.rm=T)*100,1),
								q_label=p_labs),by=c('SELECT','ST','ORIGIN')]		
dn <- dcast(da,SELECT+ST+ORIGIN~q_label,value.var="N")
dp <- dcast(da,SELECT+ST+ORIGIN~q_label,value.var="pct")

# merge central analysis with BS
dn <- merge(dn,subset(do_0,select=c('SELECT','ST','ORIGIN','N'),by=c('SELECT','ST','ORIGIN'),all=T))
dp <- merge(dp,subset(do_0,select=c('SELECT','ST','ORIGIN','pct'),by=c('SELECT','ST','ORIGIN'),all=T))

dn[, N_CI:= paste0(N, " [",CL,"-",CU,"]")]
dp[, pct_CI:= paste0(pct, "% [",CL,"-",CU,"%]")]

dt <- merge(subset(dn,select=c('SELECT','ST','ORIGIN','N_CI')),
						subset(dp,select=c('SELECT','ST','ORIGIN','pct_CI')),
						by=c('SELECT','ST','ORIGIN'),all=T)
dt[is.na(ORIGIN), ORIGIN:='Unknown']

dt$SELECT <- factor(dt$SELECT,levels=c('AmsMSM','AmsHSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
dt$ORIGIN <- factor(dt$ORIGIN,levels=c('AmsnonMSM','AmsnonHSX','NL','Africa','WEurope','EEuropeCentralAsia','NorthAm','LaAmCar','MENA','Asia','Oceania'),
										labels=c('Amsterdam - other risk group','Amsterdam - other risk group','Netherlands','Sub-Saharan Africa','Western Europe','Eastern Europe & Central Asia',
														 'North America','South America & Caribbean','Middle East & North Africa','South and South-East Asia','Oceania'))
dt$ST <- factor(dt$ST,levels=c('B','01AE','02AG','C','A1','G','D','06cpx'))
dt <- dt[order(SELECT,ST,ORIGIN),]

df <- dcast(dt,ST+ORIGIN~SELECT,value.var='N_CI')
df[is.na(`Amsterdam MSM`), `Amsterdam MSM`:='-']
df[is.na(`Amsterdam heterosexuals`), `Amsterdam heterosexuals`:='-']

write.csv(df,file=paste0(outfile.base,'-','observed_phylo_origins_bysubtype_alltime.csv'))
saveRDS(df,file=paste0(outfile.base,'-','observed_phylo_origins_bysubtype_alltime.RDS'))

dt[, N:= paste0(N_CI,' (',pct_CI,')')]
set(dt,NULL,c('N_CI','pct_CI'),NULL)
dt <- dcast(dt,ST+ORIGIN~SELECT,value.var='N')

dt[is.na(`Amsterdam MSM`), `Amsterdam MSM`:='-']
dt[is.na(`Amsterdam heterosexuals`), `Amsterdam heterosexuals`:='-']

d_origins <- copy(dt)

write.csv(dt,file=paste0(outfile.base,'-','observed_phylo_origins_bysubtype_alltime_table.csv'))
saveRDS(dt,file=paste0(outfile.base,'-','observed_phylo_origins_bysubtype_alltime_table.RDS'))


cat(" \n -------------------------------- summarise growing chains -------------------------------- \n")

# Remove individuals with unknown HIV positive date
#dsubgraphtaxa <- subset(dsubgraphtaxa, INF_D<end_d & INF_D>1912)
#dind <- subset(dind, INF_D<=end_d)

# include those infected in 2019 here]
args$end_d <- 2020

cat(paste('\n Count index cases and new cases \n'))
# summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
dsubgraphsize <- dsubgraphtaxa[, list(icases=length(ID[INF_D<args$start_d & supp==0]),jcases=length(ID[INF_D>=args$start_d & INF_D<args$end_d])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
cases <- subset(dsubgraphsize,select=c('REP','NAME','SELECT','ST','ST_CLADE','icases','jcases'))
cases <- subset(cases,SELECT %in% c('AmsMSM','AmsHSX'))

mind <- subset(dsubgraphtaxa)[,list(MIND=min(INF_D,na.rm=T),MAXD=max(INF_D,na.rm=T)),by=c('REP','NAME','ST','SELECT','ST_CLADE')]

cases <- merge(cases,mind,by=c('REP','NAME','ST','SELECT','ST_CLADE'))
cases$period <- "pre-existing"
cases$period[cases$MIND>=args$start_d] <- "emergent"

# final subgraph sizes from individuals not in ART by start date
cases$SIZE <- cases$icases + cases$jcases

# summarise by all suppressed, and growth/no further growth
cases[, allsupp:=0]
cases[icases==0 & jcases==0, allsupp:=1]
cases[, grew:= 0]
cases[jcases>0, grew:= 1]

ds <- cases[, list(N_pre=length(NAME[period=='pre-existing']),
									 N_allsupp=length(NAME[allsupp==1]),
									 N_nogrowth=length(NAME[period=='pre-existing' & grew==0]),
									 N_grew=length(NAME[period=='pre-existing' & grew==1]),
									 N_emergent=length(NAME[period=='emergent']),
									 median_size_pre=median(as.numeric(SIZE)[period=='pre-existing' & grew==1]),
									 median_size_em=median(as.numeric(SIZE)[period=='emergent'])),
						by=c('REP','SELECT','ST')]

da <- ds[, list(N_pre= round(quantile(N_pre, prob=ps,na.rm=T),0),
								N_allsupp = round(quantile(N_allsupp, prob=ps,na.rm=T),0),
								N_nogrowth = round(quantile(N_nogrowth, prob=ps,na.rm=T),0),
								N_grew = round(quantile(N_grew, prob=ps,na.rm=T),0),
								N_emergent = round(quantile(N_emergent, prob=ps,na.rm=T),0),
								median_size_pre = round(quantile(median_size_pre, prob=ps,na.rm=T),0),
								median_size_em = round(quantile(median_size_em, prob=ps,na.rm=T),0),
								q_label=p_labs),by=c('SELECT','ST')]	
#dn <- da[, list(L=paste0())]

dn <- dcast(da,SELECT+ST~q_label,value.var=c('N_pre','N_allsupp','N_nogrowth','N_grew','N_emergent','median_size_pre','median_size_em'))

dn[, N_pre:= paste0(N_pre_M,' [',N_pre_CL,'-',N_pre_CU,']')]
dn[, N_allsupp:= paste0(N_allsupp_M,' [',N_allsupp_CL,'-',N_allsupp_CU,']')]
dn[, N_nogrowth:= paste0(N_nogrowth_M,' [',N_nogrowth_CL,'-',N_nogrowth_CU,']')]
dn[, N_grew:= paste0(N_grew_M,' [',N_grew_CL,'-',N_grew_CU,']')]
dn[, N_emergent:= paste0(N_emergent_M,' [',N_emergent_CL,'-',N_emergent_CU,']')]
dn[, median_size_pre:= paste0(median_size_pre_M,' [',median_size_pre_CL,'-',median_size_pre_CU,']')]
dn[, median_size_em:= paste0(median_size_em_M,' [',median_size_em_CL,'-',median_size_em_CU,']')]
dn <- subset(dn,select=c('SELECT','ST','N_pre','N_allsupp','N_nogrowth','N_grew','N_emergent','median_size_pre','median_size_em'))
dn[median_size_pre=='NA [NA-NA]',median_size_pre:='-']
dn[median_size_em=='NA [NA-NA]',median_size_em:='-']

dt <- copy(dn)
dt$SELECT <- factor(dt$SELECT,levels=c('AmsMSM','AmsHSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
dt$ST <- factor(dt$ST,levels=c('B','01AE','02AG','C','A1','G','D','06cpx'))
dt <- dt[order(SELECT,ST),]


#dp <- dcast(da,SELECT+ST+ORIGIN~q_label,value.var="pct")
#dn[, N_CI:= paste0(M, " [",CL,"-",CU,"]")]
#dp[, pct_CI:= paste0(M, "% [",CL,"-",CU,"%]")]

#dt <- merge(subset(dn,select=c('SELECT','ST','ORIGIN','N_CI')),
#						subset(dp,select=c('SELECT','ST','ORIGIN','pct_CI')),
#						by=c('SELECT','ST','ORIGIN'),all=T)

write.csv(dt,file=paste0(outfile.base,'-','observed_phylo_growth_bysubtype.csv'))
saveRDS(dt,file=paste0(outfile.base,'-','observed_phylo_growth_bysubtype.RDS'))


cat(" \n -------------------------------- summarise chains -------------------------------- \n")

# include those infected in 2019 here]
args$end_d <- 2020

dsubgraphsize <- dsubgraphtaxa[, list(chain_size=length(ID)), by=c('ST','ST_CLADE','REP','SELECT','NAME')]

dsubgraphsize_0 <- dsubgraphsize[REP=='000', list(chains=length(NAME),
																			chain_size_1=length(NAME[chain_size==1]),
																			chain_size_25=length(NAME[chain_size>=2 & chain_size<5]),
																			chain_size_510=length(NAME[chain_size>=5 & chain_size<10]),
																			chain_size_10=length(NAME[chain_size>=10])),by=c('ST','REP','SELECT')]

dsubgraphsize <- dsubgraphsize[REP!='000', list(chains=length(NAME),
																			chain_size_1=length(NAME[chain_size==1]),
																			chain_size_25=length(NAME[chain_size>=2 & chain_size<5]),
																			chain_size_510=length(NAME[chain_size>=5 & chain_size<10]),
																			chain_size_10=length(NAME[chain_size>=10])),by=c('ST','REP','SELECT')]

da <- dsubgraphsize[, list(chains= round(quantile(chains, prob=ps,na.rm=T),0),
													 chain_size_1 = round(quantile(chain_size_1, prob=ps,na.rm=T),0),
													 chain_size_25 = round(quantile(chain_size_25, prob=ps,na.rm=T),0),
													 chain_size_510 = round(quantile(chain_size_510, prob=ps,na.rm=T),0),
													 chain_size_10 = round(quantile(chain_size_10, prob=ps,na.rm=T),0),
								q_label=p_labs),by=c('SELECT','ST')]	

dn <- dcast(da,SELECT+ST~q_label,value.var=c('chains','chain_size_1','chain_size_25','chain_size_510','chain_size_10'))

dn <- merge(dn,dsubgraphsize_0,by=c('SELECT','ST'),all=T)

dn[, chains:= paste0(chains,' [',chains_CL,'-',chains_CU,']')]
dn[, chain_size_1:= paste0(chain_size_1,' [',chain_size_1_CL,'-',chain_size_1_CU,']')]
dn[, chain_size_25:= paste0(chain_size_25,' [',chain_size_25_CL,'-',chain_size_25_CU,']')]
dn[, chain_size_510:= paste0(chain_size_510,' [',chain_size_510_CL,'-',chain_size_510_CU,']')]
dn[, chain_size_10:= paste0(chain_size_10,' [',chain_size_10_CL,'-',chain_size_10_CU,']')]
dn <- subset(dn,select=c('SELECT','ST','chains','chain_size_1','chain_size_25','chain_size_510','chain_size_10'))

dt <- copy(dn)

dt$SELECT <- factor(dt$SELECT,levels=c('AmsMSM','AmsHSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
dt$ST <- factor(dt$ST,levels=c('B','01AE','02AG','C','A1','G','D','06cpx'))
dt <- dt[order(SELECT,ST),]

write.csv(dt,file=paste0(outfile.base,'-','observed_phylo_chainsize_bysubtype.csv'))
saveRDS(dt,file=paste0(outfile.base,'-','observed_phylo_chainsize_bysubtype.RDS'))


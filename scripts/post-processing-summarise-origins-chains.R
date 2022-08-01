
cat(" \n -------------------------------- \n \n Summarise-origins-chains.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_211101'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['job_tag']] <- paste0('test_refactor_gqs_2014-2018_',args_dir[['trsm']])
#	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
#	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810b_cmdstan-',args_dir[['job_tag']])   
	args_dir[['in_dir']] <- '/Users/alexb/Documents/Roadmap/refactor_code'
	args_dir[['out_dir']] <- paste0('/Users/alexb/Documents/Roadmap/refactor_code/branching_process_model/branching_process_210810b_cmdstan-',args_dir[['job_tag']])   
	args_dir[['source_dir']] <- '~/git/locally.acquired.infections-private'
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

cat(paste('\n Subset data to individuals diagnosed between ',args_dir$start_d,'-',args_dir$end_d,' \n'))
# select patients in study window
dsubgraphtaxa[inf_after_startd==1 & inf_after_endd==0, keep:=1]

dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',args$trsm) & keep==1)

dsubgraphtaxa[, TRANSM:= gsub('Ams','',SELECT)]
dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
do <- dsubgraphtaxa
# count unique subgraphs with origin in each region (by clade too)
do <- do[,list(N=length(unique(FULL_NAME))),by=c('REP','ORIGIN','TRANSM','ST','ST_CLADE')]
# aggregate over clades for B
do <- do[,list(N=sum(N)),by=c('REP','ORIGIN','TRANSM','ST')]
do <- subset(do, ORIGIN!='Unknown')
do <- do[, list(ORIGIN=ORIGIN,N=N,pct=N/sum(N)),by=c('REP','TRANSM','ST')]
da <- do[, list(N= round(quantile(N, prob=ps,na.rm=T),0),
								pct= round(quantile(pct, prob=ps,na.rm=T)*100,1),
													q_label=p_labs),by=c('TRANSM','ST','ORIGIN')]		
dn <- dcast(da,TRANSM+ST+ORIGIN~q_label,value.var="N")
dp <- dcast(da,TRANSM+ST+ORIGIN~q_label,value.var="pct")
dn[, N_CI:= paste0(M, " [",CL,"-",CU,"]")]
dp[, pct_CI:= paste0(M, "% [",CL,"-",CU,"%]")]

dt <- merge(subset(dn,select=c('TRANSM','ST','ORIGIN','N_CI')),
						subset(dp,select=c('TRANSM','ST','ORIGIN','pct_CI')),
						by=c('TRANSM','ST','ORIGIN'),all=T)

write.csv(da,file=paste0(outfile.base,'-','observed_phylo_origins.csv'))

cat(" \n -------------------------------- summarise predicted origins -------------------------------- \n")

cat('\nSummarise predicted origins...')

do <- dsubgraphtaxa
do <- do[,list(N=length(FULL_NAME)),by=c('ORIGIN','TRANSM','ST','ST_CLADE')]
setnames(do,c('TRANSM'),c('TRM_GROUP'))

do <- subset(do, ORIGIN!='Unknown')

da <- do[, list(N=sum(N)), by=c('ORIGIN','TRM_GROUP','ST')]
da <- merge(da, da[, list(TOTAL=sum(N)), by=c('TRM_GROUP','ST')], by=c('TRM_GROUP','ST'))
tmp <- dcast.data.table(da,ORIGIN~ST,value.var='N')

loc <- data.table(location=seq(1:length(unique(tmp$ORIGIN))),
									loc_label=unique(tmp$ORIGIN))

cat('\nLoad predicted origins of subgraphs...')
origins <- fit.origins$origins_subgraphs_x
origins_em <- fit.origins$origins_subgraphs_e

cat('\nPre-existing subgraphs...')
do <- as.data.table( reshape2::melt( origins ) )
setnames(do, 1:5, c('iteration','subtype','subgraph','index_cases','location'))
do <- subset(do,is.finite(location))
do <- subset(do,location>0)

do <- merge(do,loc,by='location')
do[,chains := 'pre-existing']

cat('\nEmergent subgraphs...')
du <- as.data.table( reshape2::melt( origins_em ) )
setnames(du, 1:5, c('iteration','subtype','subgraph','index_cases','location'))
du <- subset(du,is.finite(location))
du <- subset(du,location>0)
du <- subset(du,index_cases==1)

# just keep the predicted origins for the number of individuals in unobserved subgraphs
du <- merge(du,loc,by='location')
du[,chains := 'emergent']

cat('\nCombine subgraph origins...')
do <- merge(do,du,by=c('iteration','subgraph','subtype','index_cases','location','loc_label','chains'),all=T)
do <- merge(do, pars.basic$ds,by.x='subtype',by.y='subtypes')

#	summarise subgraphs by origins
cat('\nSummarise origins...')
ans <- do[, list(N=length(subgraph)),by=c('iteration','loc_label','subtype')]
total <- do[, list(total=length(subgraph)),by=c('iteration','subtype')]
ans <- merge(ans,total,by=c('iteration','subtype'))
ans[,origin:= N/total]
ans <- ans[, list(N = quantile(N, prob=ps),
									pct = round(quantile(origin, prob=ps)*100,1),
									q_label=p_labs), 
					 by=c('loc_label','subtype')]		

ans_n <- dcast.data.table(ans, loc_label+subtype~q_label, value.var='N')
ans_p <- dcast.data.table(ans, loc_label+subtype~q_label, value.var='pct')

ans_n[, N_CI:= paste0(M, " [",round(CL,0),"-",round(CU,0),"]")]
ans_p[, pct_CI:= paste0(M, "% [",CL,"-",CU,"%]")]

ans <- merge(subset(ans_n,select=c('subtype','loc_label','N_CI')),
						subset(ans_p,select=c('subtype','loc_label','pct_CI')),
						by=c('subtype','loc_label'),all=T)

ans <- merge(ans,pars.basic$ds,by.x='subtype',by.y='subtypes')
setnames(ans,c('subtypes_name','loc_label','N_CI','pct_CI'),c('ST','ORIGIN','N_e_CI','pct_e_CI'))

ans <- merge(dt,ans,by=c('ST','ORIGIN'),all=T)

ans <- subset(ans,select=c('TRANSM','ST','ORIGIN','N_CI','pct_CI','N_e_CI','pct_e_CI'))

ans[ORIGIN=='AmsnonHSX', ORIGIN:='Amsterdam - other risk group']
ans[ORIGIN=='AmsnonMSM', ORIGIN:='Amsterdam - other risk group']
ans[ORIGIN=='Africa', ORIGIN:='Sub-Saharan Africa']
ans[ORIGIN=='CEurope', ORIGIN:='Central Europe']
ans[ORIGIN=='EEurope', ORIGIN:='Eastern Europe']
ans[ORIGIN=='EEuropeCentralAsia', ORIGIN:='Eastern Europe & Central Asia']
ans[ORIGIN=='FormerCurrDutchColonies', ORIGIN:='Suriname, Curacao & Aruba']
ans[ORIGIN=='LaAmCar', ORIGIN:='South America & Caribbean']
ans[ORIGIN=='MENA', ORIGIN:='Middle East & North Africa']
ans[ORIGIN=='NorthAm', ORIGIN:='North America']
ans[ORIGIN=='WEurope', ORIGIN:='Western Europe']
ans[ORIGIN=='NL', ORIGIN:='Netherlands']

ans[TRANSM=='MSM', TRANSM:='Amsterdam MSM']
ans[TRANSM=='HSX', TRANSM:='Amsterdam heterosexual']

write.csv(ans,file=paste0(outfile.base,'-','observed_pred_phylo_origins.csv'))
saveRDS(ans,file=paste0(outfile.base,'-','observed_pred_phylo_origins.RDS'))



require(data.table)
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))

args_dir <- list()
args_dir[['stanModelFile']] <- 'branching_process_210414r_cmdstan'
#args_dir[['stanModelFile']] <- 'branching_process_210414q_cmdstan'
args_dir[['analysis']] <- 'analysis_200917'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414r_cmdstan-inftime_undiagnosed_until2018_2010-2014_MSM'
#args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414r_cmdstan-inftime_undiagnosed_until2018_2015-2019_MSM'
#args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414q_cmdstan-inftime_undiagnosed_until2018_2010-2014_HSX'
#args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414q_cmdstan-inftime_undiagnosed_until2018_2015-2019_HSX'
args_dir[['job_tag']] <- 'inftime_undiagnosed_until2018_2010-2014_MSM'
#args_dir[['job_tag']] <- 'inftime_undiagnosed_until2018_2015-2019_MSM'
#args_dir[['job_tag']] <- 'inftime_undiagnosed_until2018_2010-2014_HSX'
#args_dir[['job_tag']] <- 'inftime_undiagnosed_until2018_2015-2019_HSX'
#args_dir[['trsm']] <- 'MSM'
args_dir[['trsm']] <- 'HSX'
args_dir[['with_subtypes']] <- 1
args_dir[['source_dir']] <- '~/git/bpm'

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

cs_obs[, newcases:=factor(new_cases)]
cs_obs[new_cases>=7,newcases:='7+']
cs_obs[, new_cases:=NULL]
setnames(cs_obs,'newcases','new_cases')

obs <- cs_obs[, list(QN=sum(freq)),by=c('new_cases','chains')]
obs[,analysis:= 'observed']
obs[,P:= 'p0.5']

################################################
cat('\nSaving predicted chain sizes...')
cs_ex <- readRDS(file=paste0(outfile.base,'-','preexisting_chain_sizes.rds'))
cs_em <- readRDS(file=paste0(outfile.base,'-','emergent_chain_sizes.rds'))
dt <- rbind(cs_ex,cs_em)
dt <- data.table(dt)
dt[, analysis:="predicted"]
N_ind <- dt[, list(N=sum(new_cases)),by=c('iteration','chains')]
dt[, newcases:=factor(new_cases)]
dt[new_cases>=7,newcases:='7+']
dt[, new_cases:=NULL]
setnames(dt,'newcases','new_cases')

cs_act <- dt[, list(freq=length(subgraph)),by=c('iteration','new_cases','chains')]

total.replicates	<- cs_act[, length(unique(iteration))]
ans	<- cs_act[, list(P=ps, 
										 QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ps)), by=c('chains','new_cases')]
set(ans, NULL, 'P', ans[, paste0('p',P)])
ans[,analysis := 'predicted']

# summarise number of chains
N_c <- dt[, list(N=length(subgraph)),by=c('iteration','chains')]
N_c <- N_c[, list(P=ps, 
									QN_c=quantile(c(N, rep(0, total.replicates-length(N))), p=ps)), by=c('chains')]
set(N_c, NULL, 'P', N_c[, paste0('p',P)])

N_ind <- N_ind[, list(P=ps, 
											QN_i=quantile(c(N, rep(0, total.replicates-length(N))), p=ps)), by=c('chains')]
set(N_ind, NULL, 'P', N_ind[, paste0('p',P)])

### merge
all <- merge(ans,subset(obs,select=c('chains','analysis','new_cases','QN','P')),by=c('chains','analysis','new_cases','P','QN'),all=T)
tot <- all[, list(tot=sum(QN)),by=c('chains','analysis','P')]
all <- merge(all,tot,by=c('chains','analysis','P'),all.x=T)
all[,p:=QN/tot]
all <- merge(all,N_c,by=c('chains','P'),all=T)
all <- merge(all,N_ind,by=c('chains','P'),all=T)

saveRDS(all,file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))


cat(" \n -------------------------------- Combine analyses -------------------------------- \n")

outfile.base <- "/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414r_cmdstan-fix_cobs_gqs_2010-2014_MSM/branching_process_210414r_cmdstan-fix_cobs_gqs_2010-2014_MSM"
msm1 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
msm1[, lab:='MSM 2010-2014']
outfile.base <- "/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414r_cmdstan-fix_cobs_gqs_2015-2019_MSM/branching_process_210414r_cmdstan-fix_cobs_gqs_2015-2019_MSM"
msm2 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
msm2[, lab:='MSM 2015-2019']
outfile.base <- "/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414q_cmdstan-fix_cobs_gqs_2010-2014_HSX/branching_process_210414q_cmdstan-fix_cobs_gqs_2010-2014_HSX"
hsx1 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
hsx1[, lab:='HSX 2010-2014']
outfile.base <- "/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414q_cmdstan-fix_cobs_gqs_2015-2019_HSX/branching_process_210414q_cmdstan-fix_cobs_gqs_2015-2019_HSX"
hsx2 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
hsx2[, lab:='HSX 2015-2019']

all <- rbind(msm1,msm2,hsx1,hsx2)

all[,N:=paste0( QN, ' (',round(p*100, d=1), '%)')]
QN_L <- subset(all,P=='p0.025',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_L,c('QN','p'),c('QN_L','p_L'))
QN_U <- subset(all,P=='p0.975',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_U,c('QN','p'),c('QN_U','p_U'))

all <- merge(subset(all,P=='p0.5'),QN_L,by=c('lab','analysis','chains','new_cases'),all=T)
all <- merge(all,QN_U,by=c('lab','analysis','chains','new_cases'),all=T)
all[analysis=='predicted',N:=paste0( QN, ' [', QN_L,'-',QN_U,']',' (',round(p*100, d=1), '%' , ' [', round(p_L*100, d=1),'-',round(p_U*100, d=1),'%',']' ,')')]

cat(" \n -------------------------------- Make table -------------------------------- \n")

tmp <- dcast.data.table(subset(all,P=="p0.5"), lab + analysis + new_cases ~ chains, value.var='N')
tmp <- subset(tmp,new_cases %in% c('0','1','2','3','4','5','6','7+'))
# print table
tmp

# summarise number of chains
chains <- all[, list(N=sum(QN)),by=c('lab', 'analysis','chains')]

cat(" \n -------------------------------- Make figure -------------------------------- \n")

# make figures
all[chains=='existing', chains:='pre-existing']

tmp <- dcast.data.table(subset(all,P=="p0.5"), lab + analysis + new_cases ~ chains, value.var='N')
tmp <- subset(tmp,new_cases %in% c('0','1','2','3','4','5','6','7+'))

dat <- subset(all,P=='p0.5' & new_cases %in% c('0','1','2','3','4','5','6','7+'))
dat[,trsm:='Amsterdam MSM']
dat[lab %in% c('HSX 2010-2014','HSX 2015-2019'),trsm:='Amsterdam heterosexual']

dat[chains=='pre-existing', chains:='pre-existing*']
dat$chains <- factor(dat$chains,levels=c('pre-existing*','emergent'))

dat$trsm <- factor(dat$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexual'))

dat[, time:= '2010-2014']
dat[lab %in% c('HSX 2015-2019','MSM 2015-2019'), time:= '2015-2019']

plot <- ggplot(data=subset(dat,time=='2010-2014')) +
	geom_bar(aes(x=new_cases,y=QN,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=QN_L, ymax=QN_U,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~chains,scales="free_y") +
	labs(x='\nNumber of cases between 2010-2014',y="Number of transmission chains",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_fig3_2010-14.png'),plot,w=15, h=12)

plot <- ggplot(data=subset(dat,time=='2015-2019')) +
	geom_bar(aes(x=new_cases,y=QN,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=QN_L, ymax=QN_U,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~chains,scales="free_y") +
	labs(x='\nNumber of cases between 2015-2019',y="Number of transmission chains",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_fig3_2015-19.png'),plot,w=15, h=12)

plot <- ggplot(data=subset(dat)) +
	geom_bar(aes(x=new_cases,y=QN,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=QN_L, ymax=QN_U,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~time+chains,scales="free_y") +
	labs(x='\nNumber of cases',y="Number of transmission chains",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_fig3_all.png'),plot,w=15, h=12)



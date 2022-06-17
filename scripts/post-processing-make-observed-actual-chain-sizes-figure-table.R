
cat(" \n -------------------------------- \n \n Summarise-chain-sizes.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

args_dir <- list()
args_dir[['stanModelFileMSM']] <- 'branching_process_210810b_cmdstan'
args_dir[['stanModelFileHSX']] <- 'branching_process_210810m_cmdstan'
args_dir[['period']] <- '2014-2018'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['job_name']] <- 'undiagnosed_untilmay'


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFileMSM')	
	stopifnot(args_line[[3]]=='-stanModelFileHSX')	
	stopifnot(args_line[[5]]=='-in_dir')
	stopifnot(args_line[[7]]=='-job_name')
	stopifnot(args_line[[9]]=='-period')
	stopifnot(args_line[[11]]=='-source_dir')
	
	args_dir <- list()
	args_dir[['stanModelFileMSM']] <- args_line[[2]]
	args_dir[['stanModelFileHSX']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['job_name']] <- args_line[[8]]
	args_dir[['period']] <- args_line[[10]]
	args_dir[['source_dir']] <- args_line[[12]]
	args_dir[['infdate']] <- 1
	
} 


cat(" \n -------------------------------- Reading data -------------------------------- \n")

#### MSM
args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileMSM']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_MSM')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_MSM')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileMSM , "-", args_dir$job_tag)
msm1 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
msm1[, lab:='MSM 2014-2019']
msm1_a <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long_all.rds'))
msm1_a[, lab:='MSM 2014-2019']

#### HSX
args_dir[['out_dir']] <- paste0(args_dir$in_dir,'/branching_process_model/',args_dir[['stanModelFileHSX']],'-',args_dir[['job_name']],'_',args_dir[['period']],'_HSX')
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_HSX')

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFileHSX , "-", args_dir$job_tag)

hsx1 <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long.rds'))
hsx1[, lab:='HSX 2014-2019']
hsx1_a <- readRDS(file=paste0(outfile.base,'-','obs_actual_cs_distribution_long_all.rds'))
hsx1_a[, lab:='HSX 2014-2019']

cat(" \n -------------------------------- Combine analyses -------------------------------- \n")

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

all <- rbind(msm1,hsx1)

all[,N:=paste0( QN, ' (',round(p*100, d=1), '%)')]
QN_L <- subset(all,P=='p0.025',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_L,c('QN','p'),c('QN_L','p_L'))
QN_U <- subset(all,P=='p0.975',select=c('lab','analysis','chains','new_cases','QN','p'))
setnames(QN_U,c('QN','p'),c('QN_U','p_U'))

#tot <- all[, list(new_cases='Total',QN=QN_c[1]),by=c('lab','analysis','chains','P')]
tot <- all[, list(new_cases='Total',QN=tot[1]),by=c('lab','analysis','chains','P')]
tot_grew <- all[, list(new_cases='Total_grew',QN=tot_grew[1]),by=c('lab','analysis','chains','P')]
tot <- rbind(tot,tot_grew)
tot <- dcast.data.table(tot, lab + new_cases + analysis ~ chains+P, value.var='QN')
tot[analysis=='predicted',`pre-existing`:=paste0( `pre-existing_p0.5`, ' [', `pre-existing_p0.025`,'-',`pre-existing_p0.975`,']')]
tot[analysis=='predicted',emergent:=paste0( emergent_p0.5, ' [', emergent_p0.025,'-',emergent_p0.975,']')]
tot[analysis=='observed',`pre-existing`:=`pre-existing_p0.5`]
tot[analysis=='observed',emergent:=emergent_p0.5]
set(tot,NULL,c('emergent_p0.025','emergent_p0.5','emergent_p0.975','pre-existing_p0.025','pre-existing_p0.5','pre-existing_p0.975'),NULL)
tot1 <- dcast.data.table(tot,lab+new_cases ~ analysis, value.var='pre-existing')
tot2 <- dcast.data.table(tot,lab+new_cases ~ analysis, value.var='emergent')
setnames(tot1,c('observed','predicted'),c('pre-existing_observed','pre-existing_predicted'))
setnames(tot2,c('observed','predicted'),c('emergent_observed','emergent_predicted'))
tot <- merge(tot1,tot2,by=c('lab','new_cases'),all=T)

all <- merge(subset(all,P=='p0.5'),QN_L,by=c('lab','analysis','chains','new_cases'),all=T)
all <- merge(all,QN_U,by=c('lab','analysis','chains','new_cases'),all=T)
all[analysis=='predicted',N:=paste0( QN, ' [', round(QN_L, d=0),'-',round(QN_U, d=0),']',' (',format(round(p*100, d=1),nsmall=1), '%' , ' [', format(round(p_L*100, d=1),nsmall=1),'-',format(round(p_U*100, d=1),nsmall=1),'%',']' ,')')]

cat(" \n -------------------------------- Figure for chains and those which grew -------------------------------- \n")
dg <- all[, list(tot_grew_em=sum(QN[new_cases!="0" & new_cases!="1"])),by=c('lab','analysis','chains','P')]
all <- merge(all,dg,by=c('lab','analysis','chains','P'),all.x=T)
all[chains=='emergent',tot_grew:=tot_grew_em]
all[, pct:=paste0(round(tot_grew/tot*100,d=0),"%")]
all$chaintype <- factor(all$chains,levels=c('pre-existing','emergent'),labels=c('pre-existing','emergent'))
all$trsmlab <- factor(all$lab,levels=c('MSM 2014-2019','HSX 2014-2019'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))

col_pal = pal_npg("nrc")(2)

plot <- ggplot(data=subset(all,analysis=='predicted' & P=='p0.5')) +
	geom_bar(aes(x=chaintype,y=tot, fill=col_pal[2]), stat="identity", position ="identity", alpha=1) +
	geom_bar(aes(x=chaintype,y=tot_grew, fill=col_pal[1]), stat="identity", position="identity", alpha=1) +
	geom_text(aes(x=chaintype,y=tot_grew,label=pct), position=position_dodge(width=0.9), vjust=-0.75,size=8) +
	facet_grid(.~trsmlab,scales="free") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections \n 2014-2019",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_text(size=24)) +
	scale_fill_manual(values = c( col_pal[2], col_pal[1]),
										labels = c("Total chains", "Chain which grew")) +
	guides(fill=guide_legend(override.aes=list(fill=c("Total chains"=col_pal[2],"Chain which grew"=col_pal[1]))))
ggsave(file=paste0(outfile.base,'-chains_grew_2014-2019.png'),plot,w=15, h=8)
cat(" \n -------------------------------- Make table -------------------------------- \n")

tmp <- dcast.data.table(subset(all,P=="p0.5"), lab + new_cases ~ chains+analysis, value.var='N')
tmp <- subset(tmp,new_cases %in% c('0','1','2','3','4','5','6','7+'))
tmp <- tmp[,c('lab', 'new_cases','pre-existing_observed','pre-existing_predicted','emergent_observed','emergent_predicted')]
tmp <- merge(tmp,tot,by=c('lab','new_cases','pre-existing_observed','pre-existing_predicted','emergent_observed','emergent_predicted'),all=T)
tmp$lab <- factor(tmp$lab,levels=c('MSM 2014-2019','HSX 2014-2019'))
tmp$new_cases <- factor(tmp$new_cases,levels=c('0','1','2','3','4','5','6','7+','Total_grew','Total'))
tmp <- tmp[order(lab,new_cases),]

dt <-dcast.data.table(subset(all,P=="p0.5" & new_cases==1), lab ~ chains+analysis, value.var='QN_i')
dt <- dt[order(lab),]
# save table
write.csv(tmp,file=paste0(outfile.base,'-obs-actual-chain-sizes-table.csv'))
saveRDS(tmp,file=paste0(outfile.base,'-obs-actual-chain-sizes-table.RDS'))

# summarise number of chains
chains <- all[, list(N=sum(QN)),by=c('lab', 'analysis','chains')]

cat(" \n -------------------------------- Make figure -------------------------------- \n")

# make figures
all[chains=='existing', chains:='pre-existing']

tmp <- dcast.data.table(subset(all,P=="p0.5"), lab + analysis + new_cases ~ chains, value.var='N')
tmp <- subset(tmp,new_cases %in% c('0','1','2','3','4','5','6','7+'))

dat <- subset(all,P=='p0.5' & new_cases %in% c('0','1','2','3','4','5','6','7+'))
dat[,trsm:='Amsterdam MSM']
dat[lab %in% c('HSX 2014-2019'),trsm:='Amsterdam heterosexual']

dat[chains=='pre-existing', chains:='pre-existing*']
dat$chains <- factor(dat$chains,levels=c('pre-existing*','emergent'))

dat$trsm <- factor(dat$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexual'))

dat[, time:= '2014-2019']

# pad data with 0s to keep bars same width
dat2 <- as.data.table(tidyr::crossing(time=unique(dat$time),analysis=unique(dat$analysis),
																			chains=unique(dat$chains),
																			trsm=unique(dat$trsm),new_cases=unique(droplevels(dat$new_cases))))
dat <- merge(dat,dat2,by=c('time','analysis','chains','trsm','new_cases'),all=T)
dat[is.na(QN),QN:=0]

plot <- ggplot(data=subset(dat,time=='2014-2019')) +
	geom_bar(aes(x=new_cases,y=QN,fill=analysis),stat='identity', position = "dodge") +
	geom_errorbar(aes(x=new_cases,ymin=QN_L, ymax=QN_U,fill=analysis),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~chains,scales="free_y") +
	labs(x='\nNumber of cases between 2014-2018',y="Number of transmission chains",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-posteriorpredictivecheck_fig3_2014-2019.png'),plot,w=15, h=12)


cat(" \n -------------------------------- Get number of new cases -------------------------------- \n")

all <- rbind(msm1_a,hsx1_a)
all <- all[, list(N_i=QN_i[1]),by=c('lab','analysis','chains','P')]
all <- all[, list(N_i=paste0(N_i[P=='p0.5'],' [',round(N_i[P=='p0.025'],d=0),'-',round(N_i[P=='p0.975'],d=0),']')),by=c('lab','analysis','chains')]

N_i <- dcast.data.table(all, lab + analysis ~ chains, value.var='N_i')

N1 <- dcast.data.table(N_i, lab ~ analysis, value.var='pre-existing')
N2 <- dcast.data.table(N_i, lab ~ analysis, value.var='emergent')
setnames(N1,c('observed','predicted'),c('pre-existing_observed','pre-existing_predicted'))
setnames(N2,c('observed','predicted'),c('emergent_observed','emergent_predicted'))
N_i <- merge(N1,N2,by=c('lab'),all=T)

N_i$lab <- factor(N_i$lab,levels=c('MSM 2014-2019','HSX 2014-2019'))
N_i <- N_i[order(lab),]

# save table
write.csv(N_i,file=paste0(outfile.base,'-obs-actual-new_cases_table.csv'))

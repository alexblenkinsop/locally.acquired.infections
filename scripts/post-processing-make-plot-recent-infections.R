
cat(" \n -------------------------------- \n \n Running post-processing-plot-recent-hiv-infections.R\n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(bayesplot, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggpubr, quietly = TRUE))
suppressMessages(library(viridis, quietly = TRUE))
suppressMessages(library(forcats, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

args_dir <- list()
args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
args_dir[['period']] <- '2014-2018'
args_dir[['trsm']] <- 'MSM'
args_dir[['source_dir']] <- '~/git/bpm'
args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810m_cmdstan-undiagnosed_untilmay_2014-2018_HSX'
args_dir[['job_name']] <- 'undiagnosed_untilmay'
args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFile')	
	stopifnot(args_line[[3]]=='-in_dir')
	stopifnot(args_line[[5]]=='-out_dir')
	stopifnot(args_line[[7]]=='-job_tag')
	stopifnot(args_line[[9]]=='-trsm')

	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['in_dir']] <- args_line[[4]]
	args_dir[['out_dir']] <- args_line[[6]]
	args_dir[['job_tag']] <- args_line[[8]]
	args_dir[['trsm']] <- args_line[[10]]
} 

args_dir[['out_dir']] <- paste0(args_dir[['in_dir']],'/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)


dat <- data.table(study=c('BLMM','BLMM','Avidity Assay','Avidity Assay'),
									trsm=c('Amsterdam MSM','Amsterdam Heterosexual','Amsterdam MSM','Amsterdam Heterosexual'),
									value=c(0.256,0.0804,0.2,0.1))

dat$trsm <- factor(dat$trsm,levels=c('Amsterdam MSM','Amsterdam Heterosexual'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
dat$study <- factor(dat$study,levels=c('BLMM','Avidity Assay'))

plot <- ggplot(data=dat) +
	geom_bar(aes(x=trsm,y=value,fill=study),stat='identity', position = "dodge") +
	scale_y_continuous(expand=expansion(c(0,0)), breaks=seq(0,1,0.1), minor_breaks=seq(0,0.1,0.05),
										 limits = c(0, 1),labels=scales::percent_format(accuracy = 1L)) +
	labs(x='',y=paste0("% of infections estimated to have been acquired\n",intToUtf8(8804)," 6 months before diagnosis"),fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-recent_hiv_infections.png'),plot,w=8, h=8)



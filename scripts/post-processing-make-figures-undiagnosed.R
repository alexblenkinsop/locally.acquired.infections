
cat(" \n -------------------------------- \n \n Running post-processing-make-figures-undiagnosed.R \n \n -------------------------------- \n")

require(ggplot2)
require(ggpubr)
require(ggsci)
require(patchwork)
require(grid)
require(data.table)

# setup
if(0){
	args <- list( 
		#source_dir= '/rds/general/user/ablenkin/home/git/bpm',
		#indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
		#outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
		source_dir= '~/Documents/GitHub/locally.acquired.infections-private',
		indir='/Users/alexb/Documents/Roadmap/refactor_code',
		outdir= '/Users/alexb/Documents/Roadmap/refactor_code/branching_process_model',
		stanModelFile= 'undiagnosed_211102',
		analysis= 'analysis_211101',
		#analysis= 'analysis_200917',
		hmc_stepsize= 0.02,
		hmc_num_samples= 15,
		hmc_num_warmup= 10,			
		seed= 42,
		chain= 1,
		job_tag= 'undiag_untilmay2019_weights',
		sens=F
	)
}

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-analysis')
	stopifnot(args_line[[7]]=='-indir')
	stopifnot(args_line[[9]]=='-outdir')
	stopifnot(args_line[[11]]=='-jobtag')
	stopifnot(args_line[[13]]=='-sens')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['indir']] <- args_line[[8]]
	args[['outdir']] <- args_line[[10]]
	args[['job_tag']] <- args_line[[12]]
	args[['sens']] <- args_line[[14]]
} 

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag)) 
job_tag <- args[['job_tag']]
outdir <- args[['outdir']]

infile.inftimes.m <- file.path(args$source_dir,'data','infection_times',paste0('time_to_diagnosis_birthplace-','MSM','.rds'))
infile.inftimes.h <- file.path(args$source_dir,'data','infection_times',paste0('time_to_diagnosis_birthplace-','HSX','.rds'))

# load data ----
cat(" \n -------------------------------- \n Load data \n -------------------------------- \n")

do <- readRDS(infile.inftimes.m)
tmp <- readRDS(infile.inftimes.h)
do <- rbind(do,tmp)

# plot of undiagnosed by year ----
cat(" \n -------------------------------- \n Make plot of undiagnosed by year \n -------------------------------- \n")

do$mlab <- factor(do$migrant_group,levels=c('NL','W.Europe,N.America,Oceania','E. & C. Europe','S. America & Caribbean','Sub-Saharan Africa','Other'),
									labels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))
do$mlab2 <- factor(do$mlab, levels=c('Other','Sub-Saharan\nAfrica','S. America &\n Caribbean','E. & C. Europe','W.Europe,\nN.America,\nOceania','NL'))
do$trsm <- factor(do$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

col_msm = pal_npg("nrc")(9)[c(1:5)]
col_hsx = pal_npg("nrc")(9)[c(8,6,7,9)]

## bar plot 
msm_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_MSM",'.RDS')))
hsx_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_HSX",'.RDS')))
du <- rbind(msm_s,hsx_s)
du <- reshape2::dcast(du,trsm+mg+year~qlabel,value.var='p')
du <- data.table(du)
du$year <- as.integer(as.character(du$year))

du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du <- du[!is.na(du$trsm),]

du$trsm <- factor(du$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexuals'))

plot_msm <- ggplot(data=subset(du,trsm=='Amsterdam MSM' & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by end of 2018",fill="Ethnicity") +
	ggtitle('Amsterdam MSM') +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=30),
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5)) +
	scale_fill_manual(values=col_msm) 
plot_hsx <- ggplot(data=subset(du,trsm=='Amsterdam heterosexuals' & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by end of 2018",fill="") +
	ggtitle('Amsterdam heterosexuals') +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=30),
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5)) +
	scale_fill_manual(values=col_hsx)
plot <- plot_msm | plot_hsx
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb.png'),plot,w=34, h=16)
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb.pdf'),plot,w=34, h=16)

saveRDS(do,file=file.path(outdir,'metadata_with_inftime_mgroup.RDS'))

# plot empirical/fitted CDFs ----
cat(" \n -------------------------------- \n Plot empirical/fitted CDFs \n -------------------------------- \n")

msm_s <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_MSM",'.rds')))
hsx_s <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_HSX",'.rds')))

shape_msm <- data.table(reshape::melt(msm_s$wb_shape_grp))
setnames(shape_msm,c('iterations','Var.2'),c('iter','mg'))
shape_hsx <- data.table(reshape::melt(hsx_s$wb_shape_grp))
setnames(shape_hsx,c('iterations','Var.2'),c('iter','mg'))
shape_msm[, trsm:='MSM']
shape_hsx[, trsm:='HSX']
shape_msm[, par:='shape']
shape_hsx[, par:='shape']

scale_msm <- data.table(reshape::melt(msm_s$wb_scale_grp))
setnames(scale_msm,c('iterations','Var.2'),c('iter','mg'))
scale_hsx <- data.table(reshape::melt(hsx_s$wb_scale_grp))
setnames(scale_hsx,c('iterations','Var.2'),c('iter','mg'))
scale_msm[, trsm:='MSM']
scale_hsx[, trsm:='HSX']
scale_msm[, par:='scale']
scale_hsx[, par:='scale']

ds <- rbind(shape_msm,shape_hsx,scale_msm,scale_hsx)
ds <- ds[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg','par')]
ds <- dcast(ds,trsm+mg+par~qlabel,value.var="p")

ds <- ds[, list(x=seq(0,10,0.1)),by=c('trsm','mg','par','p0.025','p0.5','p0.975')]
ds <- dcast(ds,trsm+mg+x~par,value.var=c('p0.025','p0.5','p0.975'))
ds[, cdf_L:= pweibull(x,shape=`p0.025_shape`,scale=`p0.025_scale`)]
ds[, cdf_M:= pweibull(x,shape=`p0.5_shape`,scale=`p0.5_scale`)]
ds[, cdf_U:= pweibull(x,shape=`p0.975_shape`,scale=`p0.975_scale`)]
ds[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
ds[trsm=='MSM', trsm := 'Amsterdam MSM']
ds <- merge(ds,unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all.x=T)
do[, mg:= mgid]

dat_msm <- subset(ds,trsm=='Amsterdam MSM')
dat_msm$mlab <- factor(dat_msm$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Other'))
dat_msm$mlab <- factor(dat_msm$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))

g_msm <- ggplot(subset(do,trsm %in% c('Amsterdam MSM')), aes(time_to_diagnosis)) +
	stat_ecdf(geom = "step",size=1.5,colour="black") +
	geom_ribbon(data=dat_msm,aes(x=x,ymin = cdf_L,ymax = cdf_U,fill=mlab),alpha = 0.4) +
	geom_line(data=dat_msm,aes(x=x,y = cdf_M,colour=mlab), show.legend = F) +
	scale_fill_manual(values=col_msm) +
	scale_colour_manual(values=col_msm) +
	labs(x='Time to diagnosis (years)',y='Probability of being diagnosed',fill='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(trsm~mlab) +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5)) +
	ggtitle('Amsterdam MSM')

dat_hsx <- subset(ds,trsm=='Amsterdam heterosexuals')
dat_hsx$mlab <- factor(dat_hsx$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))
g_hsx <- ggplot(subset(do,trsm %in% c('Amsterdam heterosexuals')), aes(time_to_diagnosis)) +
	stat_ecdf(geom = "step",size=1.5,colour="black") +
	geom_ribbon(data=dat_hsx,aes(x=x,ymin = cdf_L,ymax = cdf_U,fill=mlab),alpha = 0.4) +
	geom_line(data=dat_hsx,aes(x=x,y = cdf_M,colour=mlab), show.legend = F) +
	scale_fill_manual(values=col_hsx) +
	scale_colour_manual(values=col_hsx) +
	labs(x='Time to diagnosis (years)',y='Probability of being diagnosed',fill='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(trsm~mlab) +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5)) +
	ggtitle('Amsterdam heterosexual')

g_cdf <- ggarrange(g_msm,g_hsx,ncol=1,align="hv")
ggsave(file=file.path(outdir,'timetodiag_cdfs_facets.png'),g_cdf,w=35, h=25)

# make violin plot ----
cat(" \n -------------------------------- \n Make violin plot \n -------------------------------- \n")
mypal = pal_npg("nrc")(6)[1:6]

g <- ggplot(subset(do,trsm %in% c('Amsterdam MSM','Amsterdam heterosexuals')), aes(x=mlab2, y=time_to_diagnosis,colour=mlab2)) + 
	geom_violin() +
	coord_flip() +
	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5,binwidth=0.2) +
	scale_y_continuous(breaks=c(seq(0,15,1))) +
	facet_wrap(~trsm,scales="free") +
	labs(y='Time to diagnosis (years)', x='Ethnicity') +
	#scale_color_npg() +
	scale_color_manual(values=rev(mypal)) +
	theme_bw(base_size=40) +
	theme(strip.background=element_blank(),
				#strip.text = element_blank(),
				legend.position="none",
				#axis.text.y=element_blank(),
				axis.ticks.y=element_blank())
ggsave(file.path(outdir,'timetodiagnosis_violinplot.png'),g,w=30, h=20)

# make time to diagnosis plot ----
cat(" \n -------------------------------- \n Make time to diagnosis plot \n -------------------------------- \n")

msm_s <- readRDS(file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_MSM",'.RDS')))
hsx_s <- readRDS(file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_HSX",'.RDS')))

du <- rbind(msm_s,hsx_s)
du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du[, mg:=as.integer(mg)]

du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du <- du[!is.na(du$trsm),]
du[is.na(mg),mlab:='All']

du$trsm <- factor(du$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexuals'))

du$mlab <- factor(du$mlab,levels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe',
																	 'S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'),
									labels=c('All','NL',' W.Europe,\nN.America,Oceania','E. & C. Europe',
													 'S. America &\n Caribbean','Sub-Saharan Africa','Other'))
du[, mwmb_lab_trsm:=paste0(trsm,' ',mlab)]
order <- unique(subset(du,select=c('trsm','mlab','mwmb_lab_trsm')))
order <- order[order(trsm,mlab),]
du$mwmb_lab_trsm <- factor(du$mwmb_lab_trsm,levels=unique(order$mwmb_lab_trsm))

col_msm = pal_npg("nrc")(9)[c(1:5)]
col_hsx = pal_npg("nrc")(9)[c(8,6,7,9)]

plot_msm <- ggplot(data=subset(du,trsm=='Amsterdam MSM' & !is.na(mwmb_lab_trsm) & mlab!='All' ),aes(x=mlab,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
	#										 limits = c(0, 1)) +
	coord_flip() +
	labs(x='',y="Posterior median time to diagnosis (years)",fill="MSM, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_msm) + guides(fill = guide_legend(title.position = "top"))
ggsave(file.path(outdir,'timetodiag_MSM_CI.png'),plot_msm,w=34, h=16)

plot_msm <- ggplot(data=subset(du,trsm=='Amsterdam heterosexuals' & !is.na(mwmb_lab_trsm) & mlab!='All' ),aes(x=mlab,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_wrap(. ~ trsm,scales="free_y") +
	coord_cartesian(ylim=c(0,1)) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1)) +
	coord_flip() +
	labs(x='',y="Posterior median time to diagnosis (years)",fill="MSM, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_msm) + guides(fill = guide_legend(title.position = "top"))
ggsave(file.path(outdir,'timetodiag_HSX_CI.png'),plot_msm,w=34, h=16)

col_pal_c = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

plot <- ggplot(data=subset(du,!is.na(mwmb_lab_trsm)),aes(x=mlab,y=p0.5,fill=mwmb_lab_trsm)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mwmb_lab_trsm),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_wrap(. ~ trsm,scales="free_y") +
	coord_flip() +
	labs(x='',y="Posterior median time to diagnosis (years)",fill="MSM, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="none",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_pal_c) + guides(fill = guide_legend(title.position = "top"))
ggsave(file.path(outdir,'timetodiag_CI.png'),plot,w=34, h=16)

write.csv(du,file=file.path(outdir,'timetodiag_CI.csv'))

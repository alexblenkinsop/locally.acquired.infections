cat(" \n -------------------------------- \n \n make figure-phylogenetic-subgraphs.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))
suppressMessages(library(ggpubr, quietly = TRUE))

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810m_cmdstan'
	args_dir[['analysis']] <- 'analysis_211101'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810m_cmdstan-undiagnosed_weighted_inf_rate_2014-2018_HSX'   
	args_dir[['job_tag']] <- 'undiagnosed_weighted_inf_rate_2014-2018_HSX'
	args_dir[['trsm']] <- 'HSX'
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
	#args_dir[['source_dir']] <- '~/Documents/GitHub/bpm'
	args_dir[['start_d']] <- 2014
	args_dir[['end_d']] <- 2019
}
args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	#source_dir= '~/Documents/GitHub/bpm',
	#indir='~/Box\ Sync/Roadmap/RQ1 Estimating introductions',
	#outdir= '~/Box\ Sync/Roadmap/RQ1 Estimating introductions/branching_process_model',
	stanModelFile= 'branching_process_210810b_cmdstan',
	analysis= 'analysis_200917',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	job_tag= 'undiagnosed_untilmay',
	trsm= 'HSX',
	cmdstan = 1L,
	max_index_cases = 12,
	start_d = 2014,
	end_d = 2020,
	index_flag = 1,
	upper.bound.multiplier = 10,
	last_case = 10,
	keep_prop_dead_chains = 0,
	rho1=0,
	infdate=1,
	p_undiag=1,
	p_undiag_b=0.13,
	p_undiag_nb=0.18,
	rho="infdate",
	pre="infdate",
	em="infdate",
	nonB=1
)

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
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# get loc labels
## read subgraphs with metadata
cat('\nReading subgraph data...')
infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')

dsubgraphtaxa <- readRDS(infile.subgraphs)

cat(paste('\n Subset data to individuals diagnosed up until ',args_dir$end_d,' \n'))
# Remove individuals with unknown HIV positive date
dsubgraphtaxa <- subset(dsubgraphtaxa, INF_D<=2020 & INF_D>1912)

regex.tip.label <- '^([A-Za-z0-9]+)_+([0-9]+)_([0-9]+)_([a-zA-Z0-9-]+)$'
ylab <- gsub(regex.tip.label,'\\1',dsubgraphtaxa$FULL_NAME)

data <- subset(dsubgraphtaxa,SELECT!='Ams' & REP=='000')
data <- data[order(SEQ_YEAR),]
# find min/max seq date per subgraph
mind <- data[,list(MIND=min(INF_D,na.rm=T),MAXD=max(INF_D,na.rm=T)),by=c('FULL_NAME','ST_CLADE')]
data <- merge(data,mind,by=c('FULL_NAME','ST_CLADE'))


# get index cases
dsubgraphsize <- dsubgraphtaxa[, list(icasesart=length(ID[INF_D<args_dir$start_d & RECART_D>=args_dir$start_d]),
																			icases=length(ID[INF_D<args_dir$start_d & supp==0]),
																			jcases=length(ID[INF_D>=args_dir$start_d & INF_D<args_dir$end_d]),
																			mindicase=min(INF_D[INF_D<args_dir$start_d & supp==0],na.rm=T)),
															 by=c('ST','ST_CLADE','REP','SELECT','NAME')]
cases <- subset(dsubgraphsize,REP=='000',select=c('NAME','SELECT','ST','ST_CLADE','icases','jcases','mindicase'))

data2 <- merge(data,cases,by=c('NAME','SELECT','ST','ST_CLADE'),all.x=T)
data2[mindicase==Inf,mindicase:=MIND]
data2$sdate <- as.numeric(as.Date(0, origin = "2014-01-01",format="%Y-%m-%d"))

data2[,SG:= FULL_NAME]
data2$SG <- as.factor(data2$SG)
data2[!is.na(ST_CLADE),SG:= paste0(SG,'_',ST_CLADE)]
sgs <- data2$SG[order(data2$MIND)]
data2$SG <- factor(data2$SG,levels=sgs,labels=sgs)
data2 <- data2[order(SELECT,data2$MIND),]
data2[SELECT=='AmsHSX', id:=seq(1:nrow(data2[data$SELECT=='AmsHSX']))]
data2[SELECT=='AmsMSM', id:=seq(1:nrow(data2[data$SELECT=='AmsMSM']))]

hsx <- unique(subset(data2,SELECT=='AmsHSX', select=c('SELECT','SG')))
msm <- unique(subset(data2,SELECT=='AmsMSM', select=c('SELECT','SG')))
hsx$SG_id <- factor(hsx$SG,labels=seq(1:nrow(hsx)))
msm$SG_id <- factor(msm$SG,labels=seq(1:nrow(msm)))
sg_ids <- rbind(hsx,msm)

data2 <- merge(data2,sg_ids,by=c('SELECT','SG'),all.x=T)

# add variables to indicate whether chains continued or had no new cases for colourings
data2[, chaintype:='Pre-existed by 2014, continued to grow']
data2[MIND>=2014, chaintype:='Emergent since 2014']
data2[MAXD<2014 & icases==0, chaintype:='All members virally suppressed by 2014, no further growth']
data2[MAXD<2014 & icases>0, chaintype:='Some members not virally suppressed by 2014, no further growth']
data2[, chaintype:=factor(chaintype,levels=c('All members virally suppressed by 2014, no further growth',
																						 'Some members not virally suppressed by 2014, no further growth',
																						 'Pre-existed by 2014, continued to grow','Emergent since 2014'))]

# labels
data2[, trsm:= 'Amsterdam heterosexuals']
data2[SELECT=='AmsMSM', trsm:= 'Amsterdam MSM']
data2$trsm <- factor(data2$trsm, levels=c('Amsterdam MSM','Amsterdam heterosexuals'))

data2[, SG_id:=NULL]
data3 <- data2[with(data2, order(SELECT,MAXD)),]
hsx <- unique(subset(data3,SELECT=='AmsHSX', select=c('SELECT','SG')))
msm <- unique(subset(data3,SELECT=='AmsMSM', select=c('SELECT','SG')))
hsx[,SG_id:=seq(1:nrow(hsx))]
msm[,SG_id:=seq(1:nrow(msm))]
sg_ids <- rbind(hsx,msm)
sg_ids[, SG_id:=factor(SG_id)]

data3 <- merge(data3,sg_ids,by=c('SELECT','SG'),all.x=T)

#mypal = pal_npg("nrc")(5)[1:4]
mypal = pal_npg("nrc")(6)[c(1,2,5,4)]

# panel : 1980-2020
p1 <- ggplot(data=subset(data3)) +
	geom_line(col="black",size=0.6,aes(x=INF_D,y=SG_id)) +
	geom_point(size = 6,aes(x=INF_D,y=SG_id,col=chaintype)) +
	scale_x_continuous(expand=c(0,0),labels=seq(1980,2020,5),breaks = seq(1980, 2020, 5)) +
	scale_y_discrete(expand=c(0,0), breaks=levels(data3$SG_id)[c(rep(F, 199),T)]) +
	coord_cartesian(xlim=c(1979.9,2020)) +#,ylim=c(0,1350)) +
	facet_grid(trsm~., space='free_y',scales='free_y') +
	labs(x="Estimated date of infection",y="Phylogenetically observed transmission chains",col="") +
	scale_colour_manual(values=rev(mypal)) +
	theme_bw(base_size=40) +
	theme(axis.text.x=element_text(angle=60, vjust = 0, hjust=0),
		strip.background=element_blank(),
		strip.text = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		axis.ticks.y = element_blank(),
		panel.background = element_blank(),
		legend.position="bottom",
		legend.text=element_text(size=30),
		panel.spacing = unit(2, "lines"),
		plot.margin = unit(c(1,1.2,1,1), "cm")) +
	guides(col=guide_legend(nrow=4,byrow=TRUE))
ggsave(paste0(outfile.base,'-fig2_panel1_sidebyside_order_date.png'),p1,w=15, h=27)
ggsave(paste0(outfile.base,'-fig2_panel1_sidebyside_order_date_wide_legright.png'),p1 + theme(legend.position="right",),w=32, h=24)


### panel 2: 2014-2019, coloured by emergent, and continued since 2010 with new cases
data2[, chaintype:='Pre-existed by 2014']
data2[MIND>=2014, chaintype:='Emergent since 2014']
data2[, chaintype:=factor(chaintype,levels=c('Pre-existed by 2014','Emergent since 2014'))]
data2[, SG_id:=NULL]

data_trim <- subset(data2,MAXD>=2014)
data_trim[, SG_id:=NULL]

data4 <- data_trim[with(data_trim, order(SELECT,MAXD)),]

hsx <- unique(subset(data4,SELECT=='AmsHSX', select=c('SELECT','SG')))
msm <- unique(subset(data4,SELECT=='AmsMSM', select=c('SELECT','SG')))
hsx[,SG_id:=seq(1:nrow(hsx))]
msm[,SG_id:=seq(1:nrow(msm))]

sg_ids <- rbind(hsx,msm)
sg_ids[, SG_id:=factor(SG_id)]

data4 <- merge(data4,sg_ids,by=c('SELECT','SG'),all.x=T)

mypal = pal_npg("nrc")(5)[1:2]
# panel : 1980-2020
p2 <- ggplot(data=subset(data4,MAXD>=2014)) +
	geom_line(col="black",size=0.6,aes(x=INF_D,y=SG_id)) +
	geom_point(size = 6,aes(x=INF_D,y=SG_id,col=chaintype)) +
	scale_x_continuous(expand=c(0,0),labels=seq(2014,2020,1),breaks = seq(2014, 2020, 1)) +
	coord_cartesian(xlim=c(2013.9,2020)) + #,ylim=c(0,200)) +
	scale_y_discrete(expand=c(0,0), breaks=levels(data3$SG_id)[c(rep(F, 19),T)]) +
	facet_grid(trsm~., space='free_y',scales='free_y') +
	labs(x="Estimated date of infection",y="",col="") +
	scale_colour_manual(values=rev(mypal)) +
	theme_bw(base_size=40) +
	theme(axis.text.x=element_text(angle=60, vjust = 0, hjust=0),
		strip.background=element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		axis.ticks.y = element_blank(),
		panel.background = element_blank(),
		legend.position="bottom",
		panel.spacing = unit(2, "lines"),
		plot.margin = unit(c(1,1.2,1,1), "cm")) +
	guides(col=guide_legend(nrow=2,byrow=TRUE)) #+ 
ggsave(paste0(outfile.base,'-fig2_rightpanel_sidebyside_order_date.png'),p2,w=15, h=27)


g <- ggarrange(p1,p2,ncol=2,widths=c(1,1),labels=c('A','B'),font.label=list(size=40),align="hv")
ggsave(paste0(outfile.base,'-fig2_2panels_2014-2019_order_date.png'),w=30, h=30)

# replot just right panel with some tweaks
data4[is.na(trsm), trsm:= "Amsterdam heterosexuals"]
data4[, trsm := factor(trsm,levels=c("Amsterdam MSM","Amsterdam heterosexuals"),labels=c("Amsterdam MSM","Amsterdam\nheterosexuals"))]
p2 <- ggplot(data=subset(data4,MAXD>=2014)) +
	geom_line(col="black",size=0.6,aes(x=INF_D,y=SG_id)) +
	geom_point(size = 4,aes(x=INF_D,y=SG_id,col=chaintype)) +
	scale_x_continuous(expand=c(0,0),labels=seq(2014,2020,1),breaks = seq(2014, 2020, 1)) +
	coord_cartesian(xlim=c(2013.9,2020)) + #,ylim=c(0,200)) +
	scale_y_discrete(expand=c(0,0), breaks=levels(data3$SG_id)[c(rep(F, 19),T)]) +
	facet_grid(trsm~., space='free_y',scales='free_y') +
	labs(x="Estimated date of infection",y="Phylogenetically observed transmission chains",col="") +
	scale_colour_manual(values=rev(mypal)) +
	theme_bw(base_size=40) +
	theme(axis.text.x=element_text(angle=60, vjust = 0, hjust=0),
				strip.background=element_blank(),
				panel.grid.minor = element_blank(),
				panel.grid.major = element_blank(),
				axis.ticks.y = element_blank(),
				panel.background = element_blank(),
				legend.position="bottom",
				panel.spacing = unit(2, "lines"),
				plot.margin = unit(c(1,1.2,1,1), "cm")) +
	guides(col=guide_legend(nrow=2,byrow=TRUE)) #+ 
ggsave(paste0(outfile.base,'-fig2_rightpanel_order_date_labs.png'),p2,w=15, h=27)
ggsave(paste0(outfile.base,'-fig2_rightpanel_order_date_labs.pdf'),p2,w=15, h=27)



require(ggplot2)
require(ggpubr)
require(ggsci)
require(patchwork)
require(grid)
require(data.table)

# setup
if(0){
		home <- '/Users/alexb/Box Sync/Roadmap'
	
	dir.name <- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
	outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
	outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model')
	
	args <- list()
	args$trsm <- 'MSM'
}
args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	#source_dir= '~/Documents/GitHub/bpm',
	#indir='~/Box\ Sync/Roadmap/RQ1 Estimating introductions',
	#outdir= '~/Box\ Sync/Roadmap/RQ1 Estimating introductions/branching_process_model',
	stanModelFile= 'undiagnosed_211102',
	analysis= 'analysis_211101',
	#analysis= 'analysis_200917',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	job_tag= 'undiag_untilmay2019',
	sens=F
)

## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-source_dir')
	stopifnot(args_line[[3]]=='-stanModelFile')
	stopifnot(args_line[[5]]=='-analysis')
	stopifnot(args_line[[7]]=='-seed', !is.na(as.integer(args_line[[8]])))	
	stopifnot(args_line[[9]]=='-indir')
	stopifnot(args_line[[11]]=='-outdir')
	stopifnot(args_line[[13]]=='-jobtag')
	stopifnot(args_line[[15]]=='-sens')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['seed']] <- as.integer(args_line[[8]])
	args[['indir']] <- args_line[[10]]
	args[['outdir']] <- args_line[[12]]
	args[['job_tag']] <- args_line[[14]]
	args[['sens']] <- args_line[[16]]
} 

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag)) 
#home <- file.path(args$indir,args$analysis,'Data')
#outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outpath <- file.path(args$outdir)
#outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model',job_tag)
outdir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag))

# load data
file.seqlabels <- file.path(args$indir,args$analysis,'misc/200917_sequence_labels.rda')
infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','roadmap_cd4_vl_est.csv')
geo.file <- file.path(args$indir,'misc/NEWGEO.csv')

load(file.seqlabels)
dinf <- read.csv(infile.inftime,header=T)
geo <- data.table(read.csv(geo.file))

geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

dind <- data.table(dind)
dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
dinf <- unique(dinf)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D')),by.x='id',by.y='PATIENT',all.x=T)
do <- data.table(dinf)
do[, time:=estsctodiagMedian]

do[, INF_D:=HIV1_POS_D - time]
do <- merge(do,subset(dind,select=c('PATIENT','ORIGIN')),by.x='id', by.y='PATIENT',all.x=T)
do <- merge(do,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)

# reclassify by main migrant groups
## msm
do[TRANSM=='MSM', mwmb:="Other"]
do[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
# western countires (non-NL)
do[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania") & ORIGIN!='NL', mwmb:="G1"]
# eastern and central europe
do[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="G2"]
# caribbean and south america
do[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]

## hsx
do[TRANSM=='HSX', mwmb:="Other"]
do[TRANSM=='HSX' & ORIGIN %in% c("NL"), mwmb:="NL"]
# sub-saharan africa
do[TRANSM=='HSX' & WRLD_born %in% c("Africa"), mwmb:="G4"]
# caribbean and south america
do[TRANSM=='HSX' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]

# exclude individuals without an infection time estimate
do <- subset(do,!is.na(time))

# summarise number diagnosed to adjust for undiagnosed
do[, infdate:=INF_D]
do[is.na(INF_D) & HIV1_POS_D>=2014, infdate:=HIV1_POS_D] ## in model we use HIV pos date where there is no infection date - check this
n_diag <- do[, list(N_diag=length(unique(id[infdate>=2014]))),by=c('TRANSM','mwmb')]
n_diag <- n_diag[order(mwmb),]

## exclude patients who were infected outside the NL
do <- subset(do,INF_CNTRY=='Netherlands' | is.na(INF_CNTRY) | INF_CNTRY=='Unknown')
do <- subset(do,!(INF_D<MIG_D & INF_CNTRY!='Netherlands'))


do[mwmb=='G1' & TRANSM=='MSM',mgid:=1]
do[mwmb=='G2' & TRANSM=='MSM',mgid:=2]
do[mwmb=='G3' & TRANSM=='MSM',mgid:=3]
do[mwmb=='NL' & TRANSM=='MSM',mgid:=4]
do[mwmb=='Other' & TRANSM=='MSM',mgid:=5]
do[mwmb=='G4' & TRANSM=='HSX',mgid:=1]
do[mwmb=='G5'& TRANSM=='HSX',mgid:=2]
do[mwmb=='NL'& TRANSM=='HSX',mgid:=3]
do[mwmb=='Other'& TRANSM=='HSX',mgid:=4]

## relabel migrant groups
do[mwmb=='G1' & TRANSM=='MSM', mlab:='W.Europe,\nN.America,\nOceania']
do[mwmb=='G2' & TRANSM=='MSM', mlab:='E. & C. Europe']
do[mwmb=='G3' & TRANSM=='MSM', mlab:='S. America &\n Caribbean']
do[mwmb=='NL' & TRANSM=='MSM', mlab:='NL']
do[mwmb=='Other' & TRANSM=='MSM', mlab:='Other']

do[mwmb=='G4' & TRANSM=='HSX', mlab:='Sub-Saharan\nAfrica']
do[mwmb=='G5' & TRANSM=='HSX', mlab:='S. America &\n Caribbean']
do[mwmb=='NL' & TRANSM=='HSX', mlab:='NL']
do[mwmb=='Other' & TRANSM=='HSX', mlab:='Other']

do$mlab <- factor(do$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))

####### MAKE VIOLIN PLOTS
do$mlab2 <- factor(do$mlab, levels=c('Other','Sub-Saharan\nAfrica','S. America &\n Caribbean','E. & C. Europe','W.Europe,\nN.America,\nOceania','NL'))
do$trsm <- factor(do$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

require(ggsci)
mypal = pal_npg("nrc")(6)[1:6]

# load stan results
#sd_prior <- 'q80-q50_lnorm(0,1)'
#sd_prior <- '2010-2012_diagd_inf_NL'
#sd_prior <- '2010-2012_update'
sd_prior <- 'undiag_untilmay2019'

#outdir <- file.path(outdir,sd_prior)
msm_s <- readRDS(file=file.path(outdir, paste0('samples_',sd_prior,"_MSM",'.rds')))
hsx_s <- readRDS(file=file.path(outdir, paste0('samples_',sd_prior,"_HSX",'.rds')))

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

##Create normal density points
#dd <- crossing(trsm=unique(do$trsm),mlab=unique(do$mlab[!is.na(do$mlab)]),time=seq(0,10,0.1))
#dd$norm<-with(dd,dnorm(vals,means[as.numeric(State_CD)],
#											 sds[as.numeric(State_CD)]))

plot <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & !is.na(mlab))) +
	geom_histogram(aes(x=time,fill=mlab), colour="black") +
	facet_grid(trsm~mlab,scales="free_y") +
	#scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	#coord_cartesian(ylim=c(0,1)) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
				#strip.text = element_blank()) +
	scale_fill_npg()
ggsave(file.path(outdir,'timetodiagnosis.png'),plot,w=30, h=20)

p1 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='NL' & trsm=='Amsterdam MSM')) +
	geom_histogram(bins = 20,aes(x=time,y =..density..),fill= col_msm[1],colour = "black") +
stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("Netherlands") 
p2 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='W.Europe,\nN.America,\nOceania' & trsm=='Amsterdam MSM')) +
	geom_histogram(bins = 20,aes(x=time,y =..density..),fill= col_msm[2],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("W.Europe,\nN.America,\nOceania") +
	scale_fill_npg()
p3 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='E. & C. Europe' & trsm=='Amsterdam MSM')) +
	geom_histogram(bins = 20,aes(x=time,y =..density..),fill= col_msm[3],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("E. & C. Europe")

p4 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='S. America &\n Caribbean' & trsm=='Amsterdam MSM')) +
	geom_histogram(bins = 50,aes(x=time,y =..density..),fill= col_msm[4],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("S. America &\n Caribbean") 

p5 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='Other' & trsm=='Amsterdam MSM')) +
	geom_histogram(bins = 20,aes(x=time,y =..density..),fill= col_msm[5],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==5], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==5])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("Other") 

p_msm <- ggarrange(p1,p2,p3,p4,p5,align="hv",nrow=1)
ggsave(file.path(outdir,'timetodiagnosis_pdf_MSM.png'),p_msm,w=35, h=10)

p6 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='NL' & trsm=='Amsterdam heterosexuals')) +
	geom_histogram(bins=20,aes(x=time,y =..density..),fill= col_hsx[1],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("Netherlands") 
p7 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='S. America &\n Caribbean' & trsm=='Amsterdam heterosexuals')) +
	geom_histogram(bins=20,aes(x=time,y =..density..),fill= col_hsx[2],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("S. America &\n Caribbean") +
	scale_fill_npg()

p8 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='Sub-Saharan\nAfrica' & trsm=='Amsterdam heterosexuals')) +
	geom_histogram(bins=20,aes(x=time,y =..density..),fill= col_hsx[3],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("Sub-Saharan\nAfrica") 

p9 <- ggplot(data=subset(do,INF_D>=2010 & INF_D<2013 & mlab=='Other' & trsm=='Amsterdam heterosexuals')) +
	geom_histogram(bins=20,aes(x=time,y =..density..),fill= col_hsx[4],colour = "black") +
	stat_function(data=data.frame(x=seq(1,10,1)),fun = dweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = "black", size = 1.5) +
	labs(x='Time to diagnosis (years)',fill="Ethnicity") +
	#facet_grid(trsm~mlab,scales="free_y") +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	ggtitle("Other") 

p_hsx <- ggarrange(p6,p7,p8,p9,align="hv",nrow=1)
ggsave(file.path(outdir,'timetodiagnosis_pdf_HSX.png'),p_hsx,w=35, h=10)

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Amsterdam MSM", angle = 90) + theme_bw(base_size=40) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Amsterdam heterosexuals", angle = 90) + theme_bw(base_size=40) + theme_void()
blank <- ggplot() + annotate(geom = 'text', x=1, y=1, label="", angle = 270) + theme_void()

p_r1 <- p_r1 + annotation_custom(
	grob = textGrob(label = "Amsterdam MSM", hjust = 0, gp = gpar(cex = 1.5)),
	ymin = 0,      # Vertical position of the textGrob
	ymax = 0.3,
	xmin = 11,         # Note: The grobs are positioned outside the plot area
	xmax = 11)

p_r2 <- p_r2 + annotation_custom(
	grob = textGrob(label = "Amsterdam heterosexuals", hjust = 0, gp = gpar(cex = 1.5)),
	ymin = 0,      # Vertical position of the textGrob
	ymax = 0.2,
	xmin = 11,         # Note: The grobs are positioned outside the plot area
	xmax = 11)

p_r1 <- ggarrange(p1,p2,p3,p4,p5,align="hv",ncol=6)
p_r2 <- ggarrange(p6,p7,p8,p9,blank,align="hv",ncol=6)
p_all <- ggarrange(p_r1,p_r2,align="hv",ncol=1)
ggsave(file.path(outdir,'timetodiagnosis_pdf_all.pdf'),p_all,w=40, h=20)

col_mg = pal_npg("nrc")(6)[1:6]
#col_msm = pal_npg("nrc")(6)[c(1:4,6)]
#col_hsx = pal_npg("nrc")(6)[c(1,4:6)]
col_msm = pal_npg("nrc")(9)[c(1:5)]
col_hsx = pal_npg("nrc")(9)[c(8,6,7,9)]

#col_msm = pal_lancet("lanonc")(7)[c(2:5,7)]
#col_hsx = pal_lancet("lanonc")(7)[c(2,5:7)]

g_msm <- ggplot(subset(do,TRANSM %in% c('MSM') & INF_D>=2010 & INF_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step",size=1.5) +
	#scale_color_npg() +
	scale_colour_manual(values=col_msm) +
	#stat_function(fun=pweibull, geom="ribbon", 
	#							mapping = aes(ymin=0,ymax=..y..)) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_msm[1],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_msm[2],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_msm[3],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_msm[4],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==5], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==5])
								, aes(x = x),show.legend=FALSE,colour = col_msm[5],size=1.5) +
	labs(x='Time to diagnosis (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5)) +
	ggtitle('Amsterdam MSM')
ggsave(file.path(outdir,'msm_timetodiag_cdfs.png'),g_msm,w=15, h=15)

g_hsx <- ggplot(subset(do,TRANSM %in% c('HSX') & INF_D>=2010 & INF_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step",size=1.5) +
	#scale_color_npg() +
	scale_colour_manual(values=col_hsx) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[1],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[2],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[3],size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[4],size=1.5) +
	labs(x='Time to diagnosis (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5)) +
	ggtitle('Amsterdam heterosexuals')
ggsave(file.path(outdir,'hsx_timetodiag_cdfs.png'),g_hsx,w=15, h=30)

g_cdf <- ggarrange(g_msm,g_hsx,ncol=2,align="hv")
ggsave(file=file.path(outdir,'timetodiag_cdfs.png'),g_cdf,w=35, h=20)

#stat_function(fun = pweibull, args = list(shape = 1, scale = 1)) + 
	
g_dens <- ggplot(subset(do,TRANSM %in% c('MSM','HSX') & INF_D>=2010 & INF_D<2013), aes(x=time, color=mlab,size=0.1)) +
	geom_density()+
	#geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
	#					 linetype="dashed") +
	scale_x_continuous(breaks=c(seq(0,15,1))) +
	labs(x='Time to diagnosis (years)',y='Density', color='Ethnicity') +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank()) +
	scale_color_npg()
#ggsave(file.path(home,'timetodiagnosis_density_h.png'),g_dens,w=30, h=8)

g <- ggplot(subset(do,TRANSM %in% c('MSM','HSX') & INF_D>=2010 & INF_D<2013), aes(x=mlab2, y=time,colour=mlab2)) + 
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
ggsave(file.path(outdir,'timetodiagnosis_violinplot_h_wlabs.png'),g,w=30, h=20)

# do cdf plot without axis labels on HSX
g_hsx_facet <- ggplot(subset(do,TRANSM %in% c('HSX') & INF_D>=2010 & INF_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step") +
	#scale_color_npg() +
	scale_colour_manual(values=col_hsx) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[1]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[2]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[3]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[4]) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank(),
				axis.title.y=element_blank(),
				axis.text.y=element_blank(),
				axis.ticks.y =element_blank(), 
				plot.title = element_text(hjust=0.5))

g_v_cdf <- ggarrange(g_msm,g_hsx_facet,ncol=2,align="hv")
g_v_cdf <- ggarrange(g,g_v_cdf,ncol=1,align="hv")
ggsave(g_v_cdf,file=file.path(outdir,'timetodiagnosis_violinplot_cdfs.png'),w=20, h=20)


msm_s <- data.table(reshape2::melt(msm_s$p_undiag_year_trunc))
setnames(msm_s,c('iterations','Var2','Var3'),c('iter','mg','year'))
hsx_s <- data.table(reshape2::melt(hsx_s$p_undiag_year_trunc))
setnames(hsx_s,c('iterations','Var2','Var3'),c('iter','mg','year'))
msm_s[, trsm:='MSM']
hsx_s[, trsm:='HSX']

#### bar plot
du <- rbind(msm_s,hsx_s)
du <- du[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg','year')]
du <- dcast(du,trsm+mg+year~qlabel,value.var="p")
du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du[, year:=2019-year]

du$trsm <- factor(du$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexuals'))
plot <- ggplot(data=subset(du,year>=2014 & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by 2019",fill="Ethnicity") +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank()) +
	scale_fill_npg()
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb.png'),plot,w=30, h=16)

###########
## bar plot using truncated probabilities
msm_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_byyear_',sd_prior,"_MSM",'.RDS')))
hsx_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_byyear_',sd_prior,"_HSX",'.RDS')))
du <- rbind(msm_s,hsx_s)
du <- reshape2::melt(du,id.vars=c('trsm','mg','qlabel'))
du <- subset(du,variable!='p') # drop average probability for plot
du$year <- factor(du$variable,levels=c('p_2014','p_2015','p_2016','p_2017','p_2018'),labels=c(2014,2015,2016,2017,2018))
du$year <- as.integer(as.character(du$year))
du <- reshape2::dcast(du,trsm+mg+year~qlabel,value.var="value")
du <- data.table(du)
#### bar plot
du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du <- du[!is.na(du$trsm),]

du$trsm <- factor(du$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexuals'))
plot <- ggplot(data=subset(du,year>=2014 & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by 2019",fill="Ethnicity") +
	theme_bw(base_size=40) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank()) +
	scale_fill_npg()
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb_trunc.png'),plot,w=30, h=16)

## separate legends
plot_msm <- ggplot(data=subset(du,trsm=='Amsterdam MSM' & year>=2014 & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by 2019",fill="Ethnicity") +
	ggtitle('Amsterdam MSM') +
  theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=30),
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5)) +
	scale_fill_manual(values=col_msm) 
plot_hsx <- ggplot(data=subset(du,trsm=='Amsterdam heterosexuals' & year>=2014 & !is.na(mlab)),aes(x=year,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=year,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(.~trsm) +
	scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	labs(x='Year of infection',y="% infections undiagnosed by 2019",fill="") +
	ggtitle('Amsterdam heterosexuals') +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=30),
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5)) +
	scale_fill_manual(values=col_hsx)
plot <- plot_msm | plot_hsx
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb_truncateddist_lancet.png'),plot,w=34, h=16)

#g_und <- ggarrange(g,g_dens,plot,nrow=3,labels = c('A','B','C'),font.label = list(size=50),align="hv",axis="l",heights=c(1,1,1.5))
#g_und <- g / g_dens / plot
#g_und <- g_und + plot_annotation(tag_levels = 'A')
#ggsave(file.path(outdir,'figure_undiagnosed_3panel_byethnicity.png'),g_und,w=30, h=35)


g_msm <- ggplot(subset(do,TRANSM %in% c('MSM') & INF_D>=2010 & INF_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step") +
	#scale_color_npg() +
	scale_colour_manual(values=col_msm) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_msm[1]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_msm[2]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_msm[3]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_msm[4]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='MSM' & ds$mg==5], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='MSM' & ds$mg==5])
								, aes(x = x),show.legend=FALSE,colour = col_msm[5]) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				plot.margin = margin(0, 0, 0, 0, "cm")) +
	ggtitle('Amsterdam MSM')

g_hsx <- ggplot(subset(do,TRANSM %in% c('HSX') & INF_D>=2010 & INF_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step") +
	#scale_color_npg() +
	scale_colour_manual(values=col_hsx) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==3], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[1]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==2], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[2]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==1], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[3]) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$p0.5[ds$par=='shape' & ds$trsm=='HSX' & ds$mg==4], scale = ds$p0.5[ds$par=='scale' & ds$trsm=='HSX' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[4]) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin = margin(-1, -1, -1, -1, "cm")) +
	ggtitle('Amsterdam HSX')

g_cdf <- ggarrange(g_msm+ rremove("xlab"),g_hsx+ rremove("xlab"),ncol=2,align="hv")
annotate_figure(g_cdf, bottom = textGrob("Time since infection (years)", gp = gpar(cex = 1.3)))

g_und <- (g_msm+ rremove("xlab") | g_hsx+ rremove("xlab"))
annotate_figure(g_und, bottom = textGrob("Time since infection (years)", gp = gpar(cex = 1.3)))

#g_und <- (g_und) / plot
g_und <- ((g_msm | g_hsx)  + plot_layout(tag_level = 'new')) / plot  + plot_annotation(tag_levels = c('A','i'))
#g_und <- g_und + plot_annotation(tag_levels = c('A','i'))
ggsave(file.path(outdir,'figure_undiagnosed_2panel_byethnicity.png'),g_und,w=30, h=25)

g_und <- ggarrange(g_msm + rremove("xlab"), g_hsx+ rremove("xlab"),ncol=2,nrow=1)
annotate_figure(g_und, bottom = textGrob("Time since infection (years)", gp = gpar(cex = 1.3)))
g_n <- ggarrange(g_und,plot,labels = c("A", "B"),ncol=1,align="hv",heights=c(0.8,1))
ggsave(file.path(outdir,'figure_undiagnosed_2panel_byethnicity_v2.png'),g_n,w=30, h=25)

leg <- get_legend(plot)
############
## compare model with fistdistrplus using 2010-2012 diagnoses
dp <- readRDS("/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/branching_process_model/wb_estimates_fitdistrplus.RDS")
ds <- readRDS("/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/hierarchical_model/model_211102/wb_estimates_stan.RDS")
ds <- subset(ds,select=c('trsm','mg','par','p0.5'))
ds <- dcast(ds,trsm+mg~par,value.var='p0.5')

g_msm <- ggplot(subset(do,TRANSM %in% c('MSM')), aes(time,colour=mlab)) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='MSM' & ds$mg==4], scale = ds$scale[ds$trsm=='MSM' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_msm[1],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='MSM' & ds$mg==1], scale = ds$scale[ds$trsm=='MSM' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_msm[2],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='MSM' & ds$mg==2], scale = ds$scale[ds$trsm=='MSM' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_msm[3],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='MSM' & ds$mg==3], scale = ds$scale[ds$trsm=='MSM' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_msm[4],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='MSM' & ds$mg==5], scale = ds$scale[ds$trsm=='MSM' & ds$mg==5])
								, aes(x = x),show.legend=FALSE,colour = col_msm[5],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==4], scale = dp$scale[dp$trsm=='MSM' & dp$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_msm[1],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==1], scale = dp$scale[dp$trsm=='MSM' & dp$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_msm[2],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==2], scale = dp$scale[dp$trsm=='MSM' & dp$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_msm[3],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==3], scale = dp$scale[dp$trsm=='MSM' & dp$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_msm[4],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==5], scale = dp$scale[dp$trsm=='MSM' & dp$mg==5])
								, aes(x = x),show.legend=FALSE,colour = col_msm[5],linetype=2,size=1.5) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
	ggtitle('Amsterdam MSM')
ggsave(file.path(outdir,'weibull_msm_fistdistrplus_stanmodel.png'),g_msm,w=20, h=15)

g_hsx <- ggplot() + 	
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='HSX' & ds$mg==3], scale = ds$scale[ds$trsm=='HSX' & ds$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[1],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='HSX' & ds$mg==2], scale = ds$scale[ds$trsm=='HSX' & ds$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[2],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='HSX' & ds$mg==1], scale = ds$scale[ds$trsm=='HSX' & ds$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[3],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = ds$shape[ds$trsm=='HSX' & ds$mg==4], scale = ds$scale[ds$trsm=='HSX' & ds$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[4],linetype=1,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[ds$trsm=='HSX' & dp$mg==3], scale = dp$scale[dp$trsm=='HSX' & dp$mg==3])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[1],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='HSX' & dp$mg==2], scale = dp$scale[dp$trsm=='HSX' & dp$mg==2])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[2],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='HSX' & dp$mg==1], scale = dp$scale[dp$trsm=='HSX' & dp$mg==1])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[3],linetype=2,size=1.5) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='HSX' & dp$mg==4], scale = dp$scale[dp$trsm=='HSX' & dp$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_hsx[4],linetype=2,size=1.5) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
	ggtitle('Amsterdam heterosexual')
ggsave(file.path(outdir,'weibull_hsx_fistdistrplus_stanmodel.png'),g_hsx,w=20, h=15)

g_compare <- ggarrange(g_msm,g_hsx,ncol=2,align="hv",common.legend=T)
g_compare <- ggarrange(g_compare,leg,ncol=1,heights=c(1,0.2))
ggsave(file.path(outdir,'weibull_fistdistrplus_stanmodel_leg.png'),g_compare,w=30, h=15)
ggsave(file.path(outdir,'weibull_fistdistrplus_stanmodel_leg.pdf'),g_compare,w=30, h=15)




###### plot data vs. fitdistrplus
g_msm <- ggplot(subset(do,TRANSM %in% c('MSM') & mlab=='NL' & HIV1_POS_D>=2010 & HIV1_POS_D<2013), aes(time,colour=mlab)) +
	stat_ecdf(geom = "step") +
	scale_colour_manual(values=col_msm) +
	stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==4], scale = dp$scale[dp$trsm=='MSM' & dp$mg==4])
								, aes(x = x),show.legend=FALSE,colour = col_msm[1],linetype=1,size=1.5) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
	ggtitle('Amsterdam MSM')
ggsave(file.path(outdir,'weibull_msm_fistdistrplus_stanmodel.png'),g_msm,w=20, h=15)


## time to diagnosis bar plot

msm_s <- readRDS(file=file.path(outdir,paste0('median_timetodiagnosis_',args$job_tag,"_MSM",'.RDS')))
hsx_s <- readRDS(file=file.path(outdir,paste0('median_timetodiagnosis_',args$job_tag,"_HSX",'.RDS')))
du <- rbind(msm_s,hsx_s)
#du <- subset(du,mg!='All')
du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du[, mg:=as.integer(mg)]

du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du <- du[!is.na(du$trsm),]
du[is.na(mg),mlab:='All']

#du <- merge(du,da,by=c('trsm','mlab','p0.025','p0.5','p0.975'),all=T)

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
	#facet_wrap(. ~ trsm,scales="free_y") +
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
	coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
	#										 limits = c(0, 1)) +
	coord_flip() +
	labs(x='',y="Posterior median time to diagnosis (years)",fill="MSM, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="none",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_pal_c) + guides(fill = guide_legend(title.position = "top"))
ggsave(file.path(outdir,'timetodiag_CI.png'),plot,w=34, h=16)


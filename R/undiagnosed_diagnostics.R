require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
library(cmdstanr)
library(posterior)

# setup
home <- '/Users/alexb/Box Sync/Roadmap'

dir.name <- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model')

args <- list()
args$trsm <- 'MSM'

# load data
file.seqlabels <- file.path(home,'analysis_200917/misc/200917_sequence_labels.rda')
infile.inftime <- file.path(home,'RQ1 Estimating introductions/Data/infection_time_estimates','roadmap_cd4_vl_allpts_est.csv')
geo.file <- file.path(home,'RQ1 Estimating introductions/analysis_200917/misc/NEWGEO.csv')

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

# make stan data
sd_prior <- "2010-2012_lnorm01"
outdir <- file.path(outdir,sd_prior)

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2013)
#dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2003 & do$INF_D<2013)
dt[mwmb=='G1',mgid:=1]
dt[mwmb=='G2',mgid:=2]
dt[mwmb=='G3',mgid:=3]
dt[mwmb=='NL',mgid:=4]
dt[mwmb=='Other',mgid:=5]
N_diag <- subset(n_diag,TRANSM==args$trsm)

q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))

#exp_wb_shape <- (log(5) - log(2)) / (q80-q50)
#wb_shape <- log(exp_wb_shape)
wb_shape <- (log(log(5)) - log(log(2))) / (log(q80)-log(q50))
#wb_scale <- q50 / (log(2)^(1/shape))
wb_scale <- exp(log(q50) / (log(log(2))/(wb_shape)))

x <- seq(0,15,0.1)
y <- dweibull(x,shape=wb_shape,scale=wb_scale)
plot(x,y)

ds <- rweibull(463,shape=wb_shape,scale=wb_scale)
ggplot() + geom_histogram(aes(x=ds),binwidth=0.25) +
	theme_bw() +
	labs(x='Time to diagnosis (years)',y='simulated counts')
ggsave(file.path(outdir,'simulated_times_MSM.png'),h=5,w=5)


# calculate empirical cdf myself
dt <- subset(do,TRANSM=='MSM' & do$BIRTH_CNTRY=='Netherlands' & do$INF_D>=2010 & do$INF_D<2013)
dat <- data.table(time=dt$time[dt$time>0])
dat <- dat[order(time),]
dat[, id:=seq(1,nrow(dat),1)]
dat <- dat[, list(id=id,time=time,p=id/nrow(dat))]

	d_blmsm <- subset(do,TRANSM %in% c('MSM') & mlab=='NL' & INF_D>=2010 & INF_D<2013)
g_msm <- ggplot(dat) +
	stat_ecdf(data=subset(do,TRANSM %in% c('MSM') & mlab=='NL' & INF_D>=2010 & INF_D<2013),geom = "step",aes(x=time,colour=mlab)) +
	geom_step(mapping=aes(x=time, y=p), direction="vh", linetype=3) +
	#scale_colour_manual(values=col_msm) +
	#stat_function(data=data.frame(x=seq(0,10,1)),fun = pweibull, args = list(shape = dp$shape[dp$trsm=='MSM' & dp$mg==4], scale = dp$scale[dp$trsm=='MSM' & dp$mg==4])
	#							, aes(x = x),show.legend=FALSE,colour = col_msm[1],linetype=1,size=1.5) +
	labs(x='Time since infection (years)',y='Probability of being diagnosed',colour='Ethnicity') +
	scale_y_continuous(expand = c(0,0))  +
	#facet_grid(.~trsm) +
	theme_bw(base_size=40) +
	theme(legend.position="none",
				strip.background=element_blank(),
				strip.text = element_blank(),
				plot.title = element_text(hjust=0.5),
				plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
	ggtitle('Amsterdam MSM')
ggsave(file.path(outdir,'weibull_msm_ecdf.png'),g_msm,w=20, h=15)



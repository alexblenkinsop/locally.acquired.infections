require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
library(cmdstanr)
library(posterior)

# setup
home <- '/Users/alexb/Box Sync/Roadmap'
sd_prior <- "2010-2012_inf_NL_sens_midpoint_SC"

dir.name <- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model',sd_prior)

args <- list()
args$trsm <- 'HSX'

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
dinf <- subset(dinf,select=c('id','u','estsctodiagMedian','hiv_pos_d','hiv1_neg_d'))
dinf <- unique(dinf)
dinf <- subset(dinf,!is.na(hiv1_neg_d)) # just keep patients with a last negative test for sensitivity analysis
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D','INF_CNTRY','MIG_D')),by.x='id',by.y='PATIENT',all.x=T)
do <- data.table(dinf)
do[, estsctodiagMedian:=u/2] # use midpoint between time at risk (last neg test and first pos test)
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

# make stan data

args <- list()

args$trsm <- 'HSX'

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2013)
#dt <- subset(do,TRANSM==args$trsm & do$HIV1_POS_D>=2010 & do$HIV1_POS_D<2013)
#dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2003 & do$INF_D<2013)
dt[mwmb=='G4',mgid:=1]
dt[mwmb=='G5',mgid:=2]
dt[mwmb=='NL',mgid:=3]
dt[mwmb=='Other',mgid:=4]
N_diag <- subset(n_diag,TRANSM==args$trsm)

stan_data <- list()
q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mwmb)), idx_to_obs_array=dt$mgid, y = dt$time, 
									 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$N_diag)

options(mc.cores=parallel::detectCores())
warmup <- 1000
# model using rstan
#model <- rstan::stan_model(file.path('stan-models','undiagnosed_211102.stan'))
#fit <- rstan::sampling(model,data=stan_data,iter=2000,warmup=warmup,chains=3,verbose=TRUE)

# decrease step size
#fit <- rstan::sampling(model,data=stan_data,iter=2000,warmup=warmup,chains=3,verbose=TRUE, control=list(adapt_delta=0.99))

# model using cmdstan
model = cmdstan_model(file.path('stan-models','undiagnosed_211102.stan'))
cmdstanfit = model$sample(chains=3,data=stan_data,
													iter_warmup=500,iter_sampling=2000,
													adapt_delta=.99)
fit <- rstan::read_stan_csv(cmdstanfit$output_files())

fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp','wb_log_q80_q50_overall','wb_log_q80_q50_grp','wb_log_q50_sd','wb_log_q80_q50_sd',
										 'wb_shape_grp[1]','wb_scale_grp[1]')
gqs.pars <- c('p_undiag_av')
summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars,gqs.pars), permuted = FALSE, inc_warmup = TRUE))
print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# model summary
mixture_su <- summary(fit)$summary
saveRDS(mixture_su,file=file.path(outdir, paste0('model_summary_fit_',sd_prior,"_",args$trsm,'.rds')))

# get samples from the chains
samples <- extract(fit, inc_warmup = FALSE)
saveRDS(samples,file=file.path(outdir, paste0('samples_',sd_prior,"_",args$trsm,'.rds')))

#	traces	
color_scheme_set("mix-blue-red")
model.pars <- c('q50_overall','q50_r','q80_overall','q80_r','q50_sd','q80_sd')
p <- rstan::traceplot(fit, pars=c(fit.target.pars,'lp__'),inc_warmup=FALSE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior_",sd_prior,"_",args$trsm,".pdf")), w=10, h=20)
print(p)
dev.off()

#	pair plots	
model.pars <- c('q50_overall','q50_r[1]','q50_r[2]','q50_r[3]','q50_r[4]','q50_r[5]',
								'q80_overall','q80_r[1]','q80_r[2]','q80_r[3]','q80_r[4]','q80_r[5]',
								'q50_sd','q80_sd',"lp__")

p <- mcmc_pairs(rstan::extract(fit, pars=c(fit.target.pars,'lp__'), permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_sdprior_",sd_prior,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()

color_scheme_set("darkgray")
p <- mcmc_pairs(as.array(fit), np = nuts_params(fit), pars =c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q80_q50_overall',
																															'wb_log_q80_q50_grp[1]','wb_log_q50_sd','wb_log_q80_q50_sd',
																															'wb_shape_grp[1]','wb_scale_grp[1]',"lp__" ),
								diag_fun = "dens",off_diag_args = list(size = 0.75))
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_divergences_sdprior_",sd_prior,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()


fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q50_grp[2]','wb_log_q50_grp[3]','wb_log_q50_grp[4]','p_undiag_av[1]',
										 'p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','wb_log_q50_sd','wb_log_q80_q50_sd')
color_scheme_set("blue")
po <- rstan::extract(fit,inc_warmup=TRUE,permuted=FALSE)
p <- mcmc_intervals(po, pars = fit.target.pars)
p
ggsave(file=file.path(outdir,paste0('marginal_pairs_',sd_prior,"_",args$trsm,'.pdf')), p, w=6, h=6)

### summarise the prop undiagnosed

p <- mcmc_intervals(po, pars = ,c('p_undiag_av[1]','p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]'))
ggsave(file=file.path(outdir,paste0('marginal_pairs_undiagnosed_',sd_prior,"_",args$trsm,'.pdf')), p, w=6, h=6)

du <- data.table(reshape2::melt(samples$p_undiag_av))
setnames(du,c('iterations','Var2'),c('iter','mg'))
du <- du[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('mg')]
du <- dcast(du,mg~qlabel,value.var="p")

### calculate undiagnosed using truncated distribution
#samples <- readRDS(file=file.path(outdir, paste0('samples_',sd_prior,"_",args$trsm,'.rds')))

shape_hsx <- data.table(reshape::melt(samples$wb_shape_grp))
setnames(shape_hsx,c('iterations','Var.2'),c('iter','mg'))
shape_hsx[, trsm:='HSX']
shape_hsx[, par:='shape']

scale_hsx <- data.table(reshape::melt(samples$wb_scale_grp))
setnames(scale_hsx,c('iterations','Var.2'),c('iter','mg'))
scale_hsx[, trsm:='HSX']
scale_hsx[, par:='scale']

ds <- rbind(shape_hsx,scale_hsx)
ds <- dcast(ds,trsm+mg+iter~par,value.var="value")

ds <- ds[, list(p_2014=1 - ptrunc(2019-2014,"weibull",shape=shape,scale=scale, a=0, b=8),
								p_2015=1 - ptrunc(2019-2015,"weibull",shape=shape,scale=scale, a=0, b=8),
								p_2016=1 - ptrunc(2019-2016,"weibull",shape=shape,scale=scale, a=0, b=8),
								p_2017=1 - ptrunc(2019-2017,"weibull",shape=shape,scale=scale, a=0, b=8),
								p_2018=1 - ptrunc(2019-2018,"weibull",shape=shape,scale=scale, a=0, b=8)),
				 by=c('trsm','iter','mg')]
#ds[, av_undiagnosed:=mean(p_2014,p_2015,p_2016,p_2017,p_2018)]
ds[, av_undiagnosed:= rowMeans(ds[,4:8])]
saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_samples_',sd_prior,"_",args$trsm,'.RDS')))

ds <- ds[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
								p_2014=quantile(p_2014,prob=c(0.025,0.5,0.975)),
								p_2015=quantile(p_2015,prob=c(0.025,0.5,0.975)),
								p_2016=quantile(p_2016,prob=c(0.025,0.5,0.975)),
								p_2017=quantile(p_2017,prob=c(0.025,0.5,0.975)),
								p_2018=quantile(p_2018,prob=c(0.025,0.5,0.975)),
									qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg')]
saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',sd_prior,"_",args$trsm,'.RDS')))

ds <- dcast(ds,trsm+mg~qlabel,value.var="p")
saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',sd_prior,"_",args$trsm,'.RDS')))
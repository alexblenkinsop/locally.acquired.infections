require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
require(cmdstanr)
require(posterior)
require(truncdist)
require(ggsci)

# setup
home <- '/Users/alexb/Box Sync/Roadmap'
sd_prior <- "2010-2012_inf_NL_sens_midpoint_SC"

dir.name <- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model',sd_prior)

dir.create( outdir )

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

da <- subset(do,do$INF_D>=2010 & do$INF_D<2013 & TRANSM %in% c('HSX','MSM'))
############
# make table

tab <- data.table(var1='TRANSM',var2='TOTAL',table(da$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(da$mwmb,da$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_all <- copy(tab)

###########
## exclude patients who were infected outside the NL
do <- subset(do,INF_CNTRY=='Netherlands' | is.na(INF_CNTRY) | INF_CNTRY=='Unknown')
do <- subset(do,!(INF_D<MIG_D & INF_CNTRY!='Netherlands'))

de <- subset(do,do$INF_D>=2010 & do$INF_D<2013 & TRANSM %in% c('HSX','MSM'))

dexcl <- da[!(da$id %in% de$id),]

tab <- data.table(var1='TRANSM',var2='TOTAL',table(dexcl$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(dexcl$mwmb,dexcl$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]
tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_excl <- copy(tab)

tab <- data.table(var1='TRANSM',var2='TOTAL',table(de$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(de$mwmb,de$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

time <- de[, list(var='INFTIME',
									V1='Estimated time to diagnosis (years)',
									q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
													 " [",
													 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
													 "-",
													 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
													 "]")),
					 by=c('TRANSM','mwmb')]
time_all <- de[, list(mwmb='TOTAL',
											var='INFTIME',
											V1='Estimated time to diagnosis (years)',
											q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
															 " [",
															 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
															 "-",
															 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
															 "]")),
							 by=c('TRANSM')]
time <- rbind(time,time_all)
time <- dcast(time,TRANSM+mwmb~V1,value.var='q')
tab[var2=='TOTAL',V1:=var2]
tab <- merge(tab,time,by.y=c('TRANSM','mwmb'),by.x=c('V2','V1'),all=T)
tab <- subset(tab,!(N==0 & is.na(`Estimated time to diagnosis (years)`)))
tab_incl <- copy(tab)

setnames(tab_all,'N','N_all')
setnames(tab_excl,'N','N_excluded')
setnames(tab_incl,'N','N_included')
tab <- merge(tab_all,tab_excl,by=c('var1','var2','V1','V2'),all=T)
tab <- merge(tab,tab_incl,by=c('var1','var2','V1','V2'),all=T)

tab[V1=='G1', V1:='W.Europe, N.America, Oceania']
tab[V1=='G2', V1:='E. & C. Europe']
tab[V1=='G3', V1:='S. America & Caribbean']
tab[V1=='G4', V1:='Sub-Saharan Africa']
tab[V1=='G5', V1:='S. America & Caribbean']
tab[V1=='NL', V1:='Netherlands']
tab[V1=='Other', V1:='Other']
tab[V1=='TOTAL', V1:='All']

tab$V2 <- factor(tab$V2,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexual'))
tab <- tab[order(V2),]
saveRDS(tab,file=file.path(outdir, paste0("characteristics_patients_undiagnosed.RDS")))

### make labels for birthplaces

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

########################

# make stan data

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2013)
#dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2003 & do$INF_D<2013)
dt[mwmb=='G1',mgid:=1]
dt[mwmb=='G2',mgid:=2]
dt[mwmb=='G3',mgid:=3]
dt[mwmb=='NL',mgid:=4]
dt[mwmb=='Other',mgid:=5]
N_diag <- subset(n_diag,TRANSM==args$trsm)

ps <- c(0.5,0.025,0.975)
p_labs <- c('M','CL','CU')
time_m_ci <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
											 q_label=p_labs),by=c('TRANSM','mwmb')]
time_m_ci <- subset(time_m_ci,TRANSM %in% c('MSM','HSX'))
time_m_ci <- dcast(time_m_ci,TRANSM+mwmb~q_label,value.var='q')
time_m_ci <- subset(time_m_ci,select=c('TRANSM','mwmb','M','CL','CU'))

time_m_ci_all <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
											 q_label=p_labs),by=c('TRANSM')]
time_m_ci_all <- subset(time_m_ci_all,TRANSM %in% c('MSM','HSX'))
time_m_ci_all <- dcast(time_m_ci_all,TRANSM~q_label,value.var='q')
time_m_ci_all <- subset(time_m_ci_all,select=c('TRANSM','M','CL','CU'))


q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mwmb)), idx_to_obs_array=dt$mgid, y = dt$time, 
									 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$N_diag)

## init values
stan_init <- list()
stan_init$q50_r <- c(0.3,0.6,0.9)
stan_init$q80_r <- c(0.3,0.6,0.9)
stan_init$q50_r <- rep(stan_init$q50_r, stan_data$r)
stan_init$q80_r <- rep(stan_init$q80_r, stan_data$r)

options(mc.cores=parallel::detectCores())
#options(mc.cores=1)
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

#	examine neff and rhat
fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp','wb_log_q80_q50_overall','wb_log_q80_q50_grp','wb_log_q50_sd','wb_log_q80_q50_sd',
										 'wb_shape_grp[1]','wb_scale_grp[1]')
gqs.pars <- c('p_undiag_av','undiagnosed[1]')
summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars,gqs.pars), permuted = FALSE, inc_warmup = TRUE))
print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# model summary
mixture_su <- summary(fit)$summary
saveRDS(mixture_su,file=file.path(outdir, paste0('model_summary_fit_',sd_prior,"_",args$trsm,'.rds')))

# get samples from the chains
samples <- rstan::extract(fit, inc_warmup = FALSE)
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


fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q50_grp[2]','wb_log_q50_grp[3]','wb_log_q50_grp[4]','wb_log_q50_grp[5]','p_undiag_av[1]',
										 'p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','p_undiag_av[5]','wb_log_q50_sd','wb_log_q80_q50_sd')
color_scheme_set("blue")
po <- rstan::extract(fit,inc_warmup=TRUE,permuted=FALSE)
p <- mcmc_intervals(po, pars = fit.target.pars)
p
ggsave(file=file.path(outdir,paste0('marginal_pairs_',sd_prior,"_",args$trsm,'.pdf')), p, w=6, h=6)

### summarise the prop undiagnosed

p <- mcmc_intervals(po, pars = ,c('p_undiag_av[1]','p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','p_undiag_av[5]'))
ggsave(file=file.path(outdir,paste0('marginal_pairs_undiagnosed_',sd_prior,"_",args$trsm,'.pdf')), p, w=6, h=6)

du <- data.table(reshape2::melt(samples$p_undiag_av))
setnames(du,c('iterations','Var2'),c('iter','mg'))
du <- du[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('mg')]
du <- dcast(du,mg~qlabel,value.var="p")

### calculate undiagnosed using truncated distribution
#samples <- readRDS(file=file.path(outdir, paste0('samples_',sd_prior,"_",args$trsm,'.rds')))

shape_msm <- data.table(reshape::melt(samples$wb_shape_grp))
setnames(shape_msm,c('iterations','Var.2'),c('iter','mg'))
shape_msm[, trsm:='MSM']
shape_msm[, par:='shape']

scale_msm <- data.table(reshape::melt(samples$wb_scale_grp))
setnames(scale_msm,c('iterations','Var.2'),c('iter','mg'))
scale_msm[, trsm:='MSM']
scale_msm[, par:='scale']

ds <- rbind(shape_msm,scale_msm)
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

## calculate actual # infected with CIs

#setDT(ds)
#ds <- ds[, .(av_undiagnosed= rowMeans(.SD)), by=c('trsm','iter','mg'), .SDcols = c("p_2014", "p_2015", "p_2016", "p_2017", "p_2018")]

##########################################
## HSX model
##########################################
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
save(stan_data, file=file.path(outdir, 'stanin.RData'))
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

###########################################
## PLOT
###########################################
do$mlab2 <- factor(do$mlab, levels=c('Other','Sub-Saharan\nAfrica','S. America &\n Caribbean','E. & C. Europe','W.Europe,\nN.America,\nOceania','NL'))
do$trsm <- factor(do$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

mypal = pal_npg("nrc")(6)[1:6]

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
ggsave(file.path(outdir,'undiagnosedbyyear_mwmb_truncateddist_lancet_sens_midpoint_SC.png'),plot,w=34, h=16)


#############
## plot CDF against main analysis
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

col_msm = pal_lancet("lanonc")(7)[c(2:5,7)]
col_hsx = pal_lancet("lanonc")(7)[c(2,5:7)]

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
ggsave(file=file.path(outdir,'timetodiag_cdfs_lancet.png'),g_cdf,w=35, h=20)


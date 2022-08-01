require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
require(tidyr)
require(reshape)

## set up ----
if(0){
    # setup
    home <- '/Users/alexb/Box Sync/Roadmap'
    home <- file.path(args$in_dir,args$analysis,'Data')
    job_tag <- "undiag_untilmay2019_weights"
    #job_tag <- "2010-2012_inf_NL_sens_midpoint_SC"
    sens <- F # whether to run sensitivity analysis
}

args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/locally.acquired.infections-private',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	stanModelFile= 'undiagnosed_211102',
	analysis= 'analysis_211101',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	job_tag= 'undiag_weighted_ECDC',
	sens=F,
	weights='ECDC'
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
	stopifnot(args_line[[17]]=='-weights')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['seed']] <- as.integer(args_line[[8]])
	args[['indir']] <- args_line[[10]]
	args[['outdir']] <- args_line[[12]]
	args[['job_tag']] <- args_line[[14]]
	args[['sens']] <- args_line[[16]]
	args[['weights']] <- args_line[[18]]
} 


args$trsm <- 'MSM'
job_tag <- args$job_tag

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag)) 
outpath <- file.path(args$outdir)
outdir <- args$outdir

infile.inftimes.m <- file.path(args$source_dir,'data','infection_times',paste0('time_to_diagnosis_birthplace-','MSM','.rds'))
infile.inftimes.h <- file.path(args$source_dir,'data','infection_times',paste0('time_to_diagnosis_birthplace-','HSX','.rds'))
infile.diagnoses.i <- file.path(args$source_dir,'data','patient_data',paste0('N_diagnosed_by_mg_infdate_since_','2014','.csv'))

#dir.create( outdir )

## make stan data ----
cat(" \n -------------------------------- \n Load data \n -------------------------------- \n")

infile.ecdc.msm <- file.path(args$source_dir,'data','infection_times','AMS_MSM_infections_ECDC.csv')
infile.ecdc.nonmsm <- file.path(args$source_dir,'data','infection_times','AMS_NO_MSM_infections_ECDC.csv')

# load ECDC estimates
dw <- data.table(read.csv(infile.ecdc.msm))
dw[, trsm:='MSM']
tmp <- data.table(read.csv(infile.ecdc.nonmsm))
tmp[, trsm:='HSX']
dw <- rbind(dw,tmp)

cat(" \n -------------------------------- \n Generate weights \n -------------------------------- \n")

dw <- subset(dw, select=c(trsm,year,N_Inf_M), year %in% c(2014:2018))
dw <- dw[, list(year=year,w=N_Inf_M/sum(N_Inf_M)),by=c('trsm')]

cat(" \n -------------------------------- \n Load infection times \n -------------------------------- \n")

do <- readRDS(infile.inftimes.m)
tmp <- readRDS(infile.inftimes.h)
do <- rbind(do,tmp)

n_diag <- data.table(read.csv(infile.diagnoses.i, header=T))

# summarise number diagnosed to adjust for undiagnosed
n_diag <- n_diag[order(mwmb),]

## MSM model ----
cat(" \n -------------------------------- \n MSM model: Make stan data \n -------------------------------- \n")

args$trsm <- 'MSM'

dt <- subset(do,trsm==args$trsm)
N_diag <- subset(n_diag,TRANSM==args$trsm)

data_mg <- data.table(trsm='MSM',mgid=dt$mgid,time_to_diagnosis=dt$time)
data_mg[, migrant_group:=factor(mgid,
														 levels=c(1,2,3,4,5),
														 labels=c('W.Europe,N.America,Oceania','E. & C. Europe','S. America & Caribbean',
														 				 'NL','Other'))]

q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mgid)), idx_to_obs_array=dt$mgid, y = dt$time, 
										 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$diag)
save(stan_data, file=file.path(outdir, paste0(job_tag,'-',args$trsm,'-stanin.RData')))

## init values
init_list <- list(
	list(wb_log_q50_overall=-1.2, wb_log_q80_q50_overall=0.7, wb_log_q50_sd=0.05,wb_log_q80_q50_sd=0.01),
	list(wb_log_q50_overall=-0.5, wb_log_q80_q50_overall=1.3, wb_log_q50_sd=0.2,wb_log_q80_q50_sd=0.1),
	list(wb_log_q50_overall=-0.1, wb_log_q80_q50_overall=1.6, wb_log_q50_sd=0.3,wb_log_q80_q50_sd=0.2)
)

cat(" \n -------------------------------- \n MSM model: run model \n -------------------------------- \n")

options(mc.cores=parallel::detectCores())
#options(mc.cores=1)
warmup <- 1000
# model using rstan
#model <- rstan::stan_model(file.path('stan-models','undiagnosed_211102.stan'))
#fit <- rstan::sampling(model,data=stan_data,iter=2000,warmup=warmup,chains=3,verbose=TRUE)

# decrease step size
#fit <- rstan::sampling(model,data=stan_data,iter=2000,warmup=warmup,chains=3,verbose=TRUE, control=list(adapt_delta=0.99))

# model using cmdstan
if(args$sens!=F){
	model = rstan::stan_model(file.path(args$source_dir,'stan-models','undiagnosed_211116.stan'))
}else{
  model = rstan::stan_model(file.path(args$source_dir,'stan-models','undiagnosed_211102.stan'))
}
#fit = rstan::sampling(model,chains=3,data=stan_data,
#													iter_warmup=500,iter_sampling=2000,
#													adapt_delta=.99,
#													init = init_list)
# stan update
fit = rstan::sampling(model,chains=3,data=stan_data,
											warmup=500,iter=2000,
											control=list(adapt_delta=.99),
											init = init_list)
saveRDS(fit,file=file.path(outdir, paste0('stanfit_',job_tag,"_",args$trsm,'.rds')))

#	examine neff and rhat
fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp','wb_log_q80_q50_overall','wb_log_q80_q50_grp','wb_log_q50_sd','wb_log_q80_q50_sd',
										 'wb_shape_grp[1]','wb_scale_grp[1]')
gqs.pars <- c('p_undiag_av','undiagnosed[1]')
summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars,gqs.pars), permuted = FALSE, inc_warmup = TRUE))
print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# model summary
mixture_su <- summary(fit)$summary
saveRDS(mixture_su,file=file.path(outdir, paste0('model_summary_fit_',job_tag,"_",args$trsm,'.rds')))

# get samples from the chains
samples <- rstan::extract(fit, inc_warmup = FALSE)
saveRDS(samples,file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))

#	traces	
color_scheme_set("mix-blue-red")
model.pars <- c('q50_overall','q50_r','q80_overall','q80_r','q50_sd','q80_sd')
p <- rstan::traceplot(fit, pars=c(fit.target.pars,'lp__'),inc_warmup=FALSE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior_",job_tag,"_",args$trsm,".pdf")), w=10, h=20)
print(p)
dev.off()

# worst parameter
lp <- which(grepl("lp__",rownames(mixture_su)))
wrst <- rownames(mixture_su)[which.min(mixture_su[-lp,'n_eff'])]
small.neff <- mixture_su[wrst,]
write.csv(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.csv'))
saveRDS(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.rds'),version = 2)

p <- rstan::traceplot(fit, pars=wrst,inc_warmup=TRUE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior-worst_par-",job_tag,"_",args$trsm,".pdf")), w=8, h=3)
print(p)
dev.off()

#	pair plots	
model.pars <- c('q50_overall','q50_r[1]','q50_r[2]','q50_r[3]','q50_r[4]','q50_r[5]',
								'q80_overall','q80_r[1]','q80_r[2]','q80_r[3]','q80_r[4]','q80_r[5]',
								'q50_sd','q80_sd',"lp__")

p <- mcmc_pairs(rstan::extract(fit, pars=c(fit.target.pars,'lp__'), permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()

color_scheme_set("darkgray")
p <- mcmc_pairs(as.array(fit), np = nuts_params(fit), pars =c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q80_q50_overall',
																															'wb_log_q80_q50_grp[1]','wb_log_q50_sd','wb_log_q80_q50_sd',
																															'wb_shape_grp[1]','wb_scale_grp[1]',"lp__" ),
								diag_fun = "dens",off_diag_args = list(size = 0.75))
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_divergences_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()


fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q50_grp[2]','wb_log_q50_grp[3]','wb_log_q50_grp[4]','wb_log_q50_grp[5]','p_undiag_av[1]',
										 'p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','p_undiag_av[5]','wb_log_q50_sd','wb_log_q80_q50_sd')
color_scheme_set("blue")
po <- rstan::extract(fit,inc_warmup=TRUE,permuted=FALSE)
p <- mcmc_intervals(po, pars = fit.target.pars)
p
ggsave(file=file.path(outdir,paste0('marginal_pairs_',job_tag,"_",args$trsm,'.pdf')), p, w=6, h=6)

### summarise the prop undiagnosed

p <- mcmc_intervals(po, pars = ,c('p_undiag_av[1]','p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','p_undiag_av[5]'))
ggsave(file=file.path(outdir,paste0('marginal_pairs_undiagnosed_',job_tag,"_",args$trsm,'.pdf')), p, w=6, h=6)

du <- data.table(reshape2::melt(samples$p_undiag_av))
setnames(du,c('iterations','Var2'),c('iter','mg'))
du <- du[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('mg')]
du <- dcast(du,mg~qlabel,value.var="p")

cat(" \n -------------------------------- \n MSM model: calculate undiagnosed by month/year from samples \n -------------------------------- \n")
samples <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))

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

# calculate median
ds[, median:= scale*(log(2))^(1/shape)]
dm <- ds[, list(p=quantile(median,prob=c(0.025,0.5,0.975)),
								qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg')]
dm <- dcast(dm,trsm+mg~qlabel,value.var="p")
dm_all <- ds[, list(mg='All',
								p=quantile(median,prob=c(0.025,0.5,0.975)),
								qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm')]
dm_all <- dcast(dm_all,trsm+mg~qlabel,value.var="p")
dm <- rbind(dm,dm_all)
saveRDS(dm,file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_",args$trsm,'.RDS')))

# calculate prob(undiagnosed) for a given month/year of infection
dat <- tidyr::crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- tidyr::crossing(year=seq(2014,2018,1))
}
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
if(args$weights=='ECDC'){
	ds[, time:=(2019 - year)]
}else{
	ds[, time:=(2019+(4/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
}
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]

if(args$weights=='ECDC'){
	ds <- merge(ds,subset(dw,trsm=='MSM'),by=c('year','trsm'),all.x=T)
	# calculate mean prob across all months (per mc sample)
	mean <- ds[, list(av_undiagnosed=sum(p*w)),
						 by=c('trsm','mg','iter')]
	saveRDS(mean,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_",args$trsm,'.RDS')))
	
	# summarise quantiles for each year
	mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
													qlabel=c('p0.025','p0.5','p0.975')),
									 by=c('trsm','mg','year')] # summarise quantiles for each year
	saveRDS(mean_y,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
										qlabel=c('p0.025','p0.5','p0.975')),
						 by=c('trsm','mg')] # summarise quantiles for all months
	ds <- dcast(ds,trsm+mg~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_",args$trsm,'.RDS')))
}else{
	
	# calculate mean prob across all months (per mc sample)
	mean_s <- ds[, list(av_undiagnosed=mean(p)),
									by=c('trsm','mg','iter')]
	saveRDS(mean_s,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_",args$trsm,'.RDS')))
	
	# calculate mean prob for January of each year (per mc sample)
	mean_y <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','mg','iter','year')]
	# summarise quantiles for each year
	ds <- mean_y[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','mg','year')] # summarise quantiles for each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean_s[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
															qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','mg')] # summarise quantiles for all months
	ds <- dcast(ds,trsm+mg~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_",args$trsm,'.RDS')))

}
cat(" \n -------------------------------- \n MSM model: calculate undiagnosed for all MSM using overall parameters \n -------------------------------- \n")

shape_msm <- data.table(reshape::melt(samples$wb_shape_overall))
setnames(shape_msm,c('iterations'),c('iter'))
shape_msm[, trsm:='MSM']
shape_msm[, par:='shape']

scale_msm <- data.table(reshape::melt(samples$wb_scale_overall))
setnames(scale_msm,c('iterations'),c('iter'))
scale_msm[, trsm:='MSM']
scale_msm[, par:='scale']

ds <- rbind(shape_msm,scale_msm)
ds <- dcast(ds,trsm+iter~par,value.var="value")

# calculate median
ds[, median:= scale*(log(2))^(1/shape)]
dm <- ds[, list(p=quantile(median,prob=c(0.025,0.5,0.975)),
								qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm')]
dm <- dcast(dm,trsm~qlabel,value.var="p")
saveRDS(ds,file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_all_",args$trsm,'.RDS')))

# calculate prob(undiagnosed) for a given month/year of infection
dat <- tidyr::crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- tidyr::crossing(year=seq(2014,2018,1))
}
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
#ds[, time:=(2019-year)-(month/12)] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
if(args$weights=='ECDC'){
	ds[, time:=(2019 - year)]
}else{
	ds[, time:=(2019+(4/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
}
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]
ds_c <- copy(ds)

if(args$weights=='ECDC'){
	ds <- merge(ds,subset(dw,trsm=='MSM'),by=c('year','trsm'),all.x=T)
	# calculate mean prob across all months (per mc sample)
	mean <- ds[, list(av_undiagnosed=sum(p*w)),
						 by=c('trsm','iter')]
	saveRDS(mean,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_all_",args$trsm,'.RDS')))
	
	ds <- copy(ds_c)
	ds <- merge(ds,subset(dw,trsm=='MSM'),by=c('year','trsm'),all.x=T)
	
	# summarise quantiles for each year
	mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','year')] # summarise quantiles for Jan of each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm')] # summarise quantiles for all months
	ds <- dcast(ds,trsm~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",args$trsm,'.RDS')))
	
}else{
	mean_s <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','iter')]
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# calculate mean prob for each year (per mc sample)
	mean_y <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','iter','year')]
	# summarise quantiles for each year
	ds <- mean_y[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','year')] # summarise quantiles for Jan of each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean_s[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm')] # summarise quantiles for all months
	ds <- dcast(ds,trsm~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",args$trsm,'.RDS')))
	
}

## HSX model ----
cat(" \n -------------------------------- \n HSX model: Make stan data \n -------------------------------- \n")

args$trsm <- 'HSX'

dt <- subset(do,trsm==args$trsm)
N_diag <- subset(n_diag,TRANSM==args$trsm)

data_mg <- data.table(trsm='HSX',mgid=dt$mgid,time_to_diagnosis=dt$time)
data_mg[, migrant_group:=factor(mgid,
																levels=c(1,2,3,4),
																labels=c('Sub-Saharan Africa','S. America & Caribbean',
																				 'NL','Other'))]
stan_data <- list()
q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mgid)), idx_to_obs_array=dt$mgid, y = dt$time, 
									 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$diag)
save(stan_data, file=file.path(outdir, paste0(job_tag,'-',args$trsm,'-stanin.RData')))

options(mc.cores=parallel::detectCores())
warmup <- 1000

# model using cmdstan
#if(args$sens!=F){
#	model = cmdstan_model(file.path('stan-models','undiagnosed_211116.stan'))
#}else{
#	model = cmdstan_model(file.path('stan-models','undiagnosed_211102.stan'))
#}
#cmdstanfit = model$sample(chains=3,data=stan_data,
#													iter_warmup=500,iter_sampling=2000,
#													adapt_delta=.99)
#fit <- rstan::read_stan_csv(cmdstanfit$output_files())

if(args$sens!=F){
	model = rstan::stan_model(file.path(args$source_dir,'stan-models','undiagnosed_211116.stan'))
}else{
	model = rstan::stan_model(file.path(args$source_dir,'stan-models','undiagnosed_211102.stan'))
}

fit = rstan::sampling(model,chains=3,data=stan_data,
											warmup=500,iter=2000,
											control=list(adapt_delta=.99),
											init = init_list)
saveRDS(fit,file=file.path(outdir, paste0('stanfit_',job_tag,"_",args$trsm,'.rds')))

fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp','wb_log_q80_q50_overall','wb_log_q80_q50_grp','wb_log_q50_sd','wb_log_q80_q50_sd',
										 'wb_shape_grp[1]','wb_scale_grp[1]')
gqs.pars <- c('p_undiag_av')
summary <- rstan::monitor(rstan::extract(fit, pars=c(fit.target.pars,gqs.pars), permuted = FALSE, inc_warmup = TRUE))
print(summary,probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# model summary
mixture_su <- summary(fit)$summary
saveRDS(mixture_su,file=file.path(outdir, paste0('model_summary_fit_',job_tag,"_",args$trsm,'.rds')))

# get samples from the chains
samples <- rstan::extract(fit, inc_warmup = FALSE)
saveRDS(samples,file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))

#	traces	
color_scheme_set("mix-blue-red")
model.pars <- c('q50_overall','q50_r','q80_overall','q80_r','q50_sd','q80_sd')
p <- rstan::traceplot(fit, pars=c(fit.target.pars,'lp__'),inc_warmup=TRUE, ncol = 1)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior_",job_tag,"_",args$trsm,".pdf")), w=10, h=20)
print(p)
dev.off()

# worst parameter
lp <- which(grepl("lp__",rownames(mixture_su)))
wrst <- rownames(mixture_su)[which.min(mixture_su[-lp,'n_eff'])]
small.neff <- mixture_su[wrst,]
write.csv(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.csv'))
saveRDS(small.neff,file=paste0(outdir,'/undiag_untilmay2019',"_",args$trsm,'-pars-with-smallest-neff.rds'),version = 2)

p <- rstan::traceplot(fit, pars=wrst,inc_warmup=TRUE, ncol = 1,)
pdf(file=file.path(outdir, paste0("undiagnosed_traces_sdprior-worst_par-",job_tag,"_",args$trsm,".pdf")), w=8, h=3)
print(p)
dev.off()

#	pair plots	
model.pars <- c('q50_overall','q50_r[1]','q50_r[2]','q50_r[3]','q50_r[4]','q50_r[5]',
								'q80_overall','q80_r[1]','q80_r[2]','q80_r[3]','q80_r[4]','q80_r[5]',
								'q50_sd','q80_sd',"lp__")

p <- mcmc_pairs(rstan::extract(fit, pars=c(fit.target.pars,'lp__'), permuted=FALSE, inc_warmup=FALSE), diag_fun = "dens", off_diag_fun = "hex")
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()

color_scheme_set("darkgray")
p <- mcmc_pairs(as.array(fit), np = nuts_params(fit), pars =c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q80_q50_overall',
																															'wb_log_q80_q50_grp[1]','wb_log_q50_sd','wb_log_q80_q50_sd',
																															'wb_shape_grp[1]','wb_scale_grp[1]',"lp__" ),
								diag_fun = "dens",off_diag_args = list(size = 0.75))
pdf(file=file.path(outdir, paste0("undiagnosed_mcmc_pairs_divergences_sdprior_",job_tag,"_",args$trsm,".pdf")), w=40, h=40)
print(p)
dev.off()


fit.target.pars <- c('wb_log_q50_overall','wb_log_q50_grp[1]','wb_log_q50_grp[2]','wb_log_q50_grp[3]','wb_log_q50_grp[4]','p_undiag_av[1]',
										 'p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]','wb_log_q50_sd','wb_log_q80_q50_sd')
color_scheme_set("blue")
po <- rstan::extract(fit,inc_warmup=TRUE,permuted=FALSE)
p <- mcmc_intervals(po, pars = fit.target.pars)
p
ggsave(file=file.path(outdir,paste0('marginal_pairs_',job_tag,"_",args$trsm,'.pdf')), p, w=6, h=6)

p <- mcmc_intervals(po, pars = ,c('p_undiag_av[1]','p_undiag_av[2]','p_undiag_av[3]','p_undiag_av[4]'))
ggsave(file=file.path(outdir,paste0('marginal_pairs_undiagnosed_',job_tag,"_",args$trsm,'.pdf')), p, w=6, h=6)

du <- data.table(reshape2::melt(samples$p_undiag_av))
setnames(du,c('iterations','Var2'),c('iter','mg'))
du <- du[, list(p=quantile(value,prob=c(0.025,0.5,0.975)),qlabel=c('p0.025','p0.5','p0.975')),by=c('mg')]
du <- dcast(du,mg~qlabel,value.var="p")
du

cat(" \n -------------------------------- \n HSX model: calculate prop. undiagnosed \n -------------------------------- \n")

samples <- readRDS(file=file.path(outdir, paste0('samples_',job_tag,"_",args$trsm,'.rds')))

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

# calculate median
ds[, median:= scale*(log(2))^(1/shape)]
dm <- ds[, list(p=quantile(median,prob=c(0.025,0.5,0.975)),
								qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm','mg')]
dm <- dcast(dm,trsm+mg~qlabel,value.var="p")
dm_all <- ds[, list(mg='All',
										p=quantile(median,prob=c(0.025,0.5,0.975)),
										qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm')]
dm_all <- dcast(dm_all,trsm+mg~qlabel,value.var="p")
dm <- rbind(dm,dm_all)
saveRDS(dm,file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_",args$trsm,'.RDS')))

# calculate prob(undiagnosed) for a given month/year of infection
dat <- tidyr::crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- tidyr::crossing(year=seq(2014,2018,1))
}
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
#ds[, time:=(2019-year)-(month/12)] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
if(args$weights=='ECDC'){
	ds[, time:=(2019 - year)]
}else{
	ds[, time:=(2019+(4/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
}
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]

if(args$weights=='ECDC'){
	ds <- merge(ds,subset(dw,trsm=='HSX'),by=c('year','trsm'),all.x=T)
	# calculate mean prob across all months (per mc sample)
	mean <- ds[, list(av_undiagnosed=sum(p*w)),
						 by=c('trsm','mg','iter')]
	saveRDS(mean,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_",args$trsm,'.RDS')))
	
	# calculate mean prob for each year (per mc sample)
	mean_y <- ds[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
													qlabel=c('p0.025','p0.5','p0.975')),
									 by=c('trsm','mg','year')] # summarise quantiles for each year
	saveRDS(mean_y,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
										qlabel=c('p0.025','p0.5','p0.975')),
						 by=c('trsm','mg')] # summarise quantiles for all months
	ds <- dcast(ds,trsm+mg~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_",args$trsm,'.RDS')))
}else{
	
	# calculate mean prob across all months (per mc sample)
	mean_s <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','mg','iter')]
	saveRDS(mean_s,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_",args$trsm,'.RDS')))
	
	# calculate mean prob for January of each year (per mc sample)
	mean_y <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','mg','iter','year')]
	# summarise quantiles for each year
	ds <- mean_y[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','mg','year')] # summarise quantiles for each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean_s[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','mg')] # summarise quantiles for all months
	ds <- dcast(ds,trsm+mg~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_",args$trsm,'.RDS')))
	
}

cat(" \n -------------------------------- \n HSX model: prop undiagnosed using overall parameters \n -------------------------------- \n")

shape_hsx <- data.table(reshape::melt(samples$wb_shape_overall))
setnames(shape_hsx,c('iterations'),c('iter'))
shape_hsx[, trsm:='HSX']
shape_hsx[, par:='shape']

scale_hsx <- data.table(reshape::melt(samples$wb_scale_overall))
setnames(scale_hsx,c('iterations'),c('iter'))
scale_hsx[, trsm:='HSX']
scale_hsx[, par:='scale']

ds <- rbind(shape_hsx,scale_hsx)
ds <- dcast(ds,trsm+iter~par,value.var="value")

# calculate median
ds[, median:= scale*(log(2))^(1/shape)]
dm <- ds[, list(p=quantile(median,prob=c(0.025,0.5,0.975)),
								qlabel=c('p0.025','p0.5','p0.975')),by=c('trsm')]
dm <- dcast(dm,trsm~qlabel,value.var="p")
saveRDS(ds,file=file.path(outdir,paste0('median_timetodiagnosis_',job_tag,"_all_",args$trsm,'.RDS')))

# calculate prob(undiagnosed) for a given month/year of infection
dat <- tidyr::crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- tidyr::crossing(year=seq(2014,2018,1))
}
ds <- merge(dat,ds,all=T)
ds <- data.table(ds)
#ds[, time:=(2019-year)-(month/12)] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
if(args$weights=='ECDC'){
	ds[, time:=(2019 - year)]
}else{
	ds[, time:=(2019+(4/12))-(year + (month/12))] # +(1/12) to add 1 month so we count until end Dec 2018/start of Jan 2019. but I think it makes sense someone infected Dec 2018 has 0 prob of being diagnosed by end of 2018 due to time taken to detect virus
}
ds[, p:=1 - pweibull(time,shape=shape,scale=scale)]
ds_c <- copy(ds)

if(args$weights=='ECDC'){
	ds <- merge(ds,subset(dw,trsm=='HSX'),by=c('year','trsm'),all.x=T)
	# calculate mean prob across all months (per mc sample)
	mean <- ds[, list(av_undiagnosed=sum(p*w)),
						 by=c('trsm','iter')]
	saveRDS(mean,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
										qlabel=c('p0.025','p0.5','p0.975')),
						 by=c('trsm')] # summarise quantiles for all months
	ds <- dcast(ds,trsm~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",args$trsm,'.RDS')))
	
	ds <- copy(ds_c)
	ds <- merge(ds,subset(dw,trsm=='HSX'),by=c('year','trsm'),all.x=T)
	
	# summarise quantiles for each year
	ds <- mean_y[, list(p=quantile(p,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','year')] # summarise quantiles for Jan of each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
										qlabel=c('p0.025','p0.5','p0.975')),
						 by=c('trsm')] # summarise quantiles for all months
	ds <- dcast(ds,trsm~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",args$trsm,'.RDS')))
	
}else{
	mean_s <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','iter')]
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_samples_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# calculate mean prob for each year (per mc sample)
	mean_y <- ds[, list(av_undiagnosed=mean(p)),
							 by=c('trsm','iter','year')]
	# summarise quantiles for each year
	ds <- mean_y[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm','year')] # summarise quantiles for Jan of each year
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_byyear_',job_tag,"_all_",args$trsm,'.RDS')))
	
	# summarise quantiles for average across all months/years
	ds <- mean_s[, list(p=quantile(av_undiagnosed,prob=c(0.025,0.5,0.975)),
											qlabel=c('p0.025','p0.5','p0.975')),
							 by=c('trsm')] # summarise quantiles for all months
	ds <- dcast(ds,trsm~qlabel,value.var="p")
	saveRDS(ds,file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",args$trsm,'.RDS')))
	
}

require(data.table)
require(rstan)
require(bayesplot)
require(hexbin)
require(cmdstanr)
require(posterior)
require(truncdist)
require(scales)
require(ggsci)
require(tidyr)
require(reshape)
require(patchwork)
require(ggpubr)

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
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	dir.ecdc='/rds/general/project/ratmann_roadmap_data_analysis/live/Data/Undiagnosed/ECDC_model',
	#source_dir= '~/Documents/GitHub/bpm',
	#indir='~/Box Sync/Roadmap',
	#outdir= '~/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/hierarchical_model',
	dir.ecdc='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/ECDC_model',
	stanModelFile= 'undiagnosed_211102',
	analysis= 'analysis_211101',
	#analysis= 'analysis_200917',
	hmc_stepsize= 0.02,
	hmc_num_samples= 15,
	hmc_num_warmup= 10,			
	seed= 42,
	chain= 1,
	#job_tag= 'undiag_untilmay2019',
	job_tag= 'undiag_untilmay2019_weights',
	#job_tag= '2010-2012_inf_NL_sens_midpoint_SC',
	#job_tag= 'sens_q40',
	#sens='mp',
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
	stopifnot(args_line[[17]]=='-dir.ecdc')
	stopifnot(args_line[[19]]=='-weights')
	
	args <- list()
	args[['source_dir']] <- args_line[[2]]
	args[['stanModelFile']] <- args_line[[4]]
	args[['analysis']] <- args_line[[6]]
	args[['seed']] <- as.integer(args_line[[8]])
	args[['indir']] <- args_line[[10]]
	args[['outdir']] <- args_line[[12]]
	args[['job_tag']] <- args_line[[14]]
	args[['sens']] <- args_line[[16]]
	args[['dir.ecdc']] <- args_line[[18]]
	args[['weights']] <- args_line[[18]]
} 

args$job_dir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag)) 
#home <- file.path(args$indir,args$analysis,'Data')
#outpath <- file.path(home,'RQ1 Estimating introductions','Manuscript')
outpath <- file.path(args$outdir)
#outdir <- file.path(home,'RQ1 Estimating introductions','undiagnosed','hierarchical_model',job_tag)
outdir <- file.path(args$outdir,paste0(args$stanModelFile,'-',args$job_tag))

dir.create( outdir )

args$trsm <- 'MSM'
job_tag <- args$job_tag

## make stan data ----
cat(" \n -------------------------------- \n Load data \n -------------------------------- \n")

# load data
file.seqlabels <- file.path(args$indir,args$analysis,'misc/200917_sequence_labels.rda')
infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','roadmap_cd4_vl_est.csv')
if(args$sens!=F){
	infile.inftime <- file.path(args$indir,'Data/infection_time_estimates','roadmap_cd4_vl_est-quantiles.csv')
}
geo.file <- file.path(args$indir,'misc/NEWGEO.csv')

infile.ecdc.msm <- file.path(args$dir.ecdc,'AMS_MSM_LOCAL_MR_Result_main.csv')
infile.ecdc.nonmsm <- file.path(args$dir.ecdc,'AMS_NO_MSM_LOCAL_MR_Result_main.csv')

load(file.seqlabels)
dinf <- read.csv(infile.inftime,header=T)
geo <- data.table(read.csv(geo.file))

# load ECDC estimates
dw <- data.table(read.csv(infile.ecdc.msm))
dw[, trsm:='MSM']
tmp <- data.table(read.csv(infile.ecdc.nonmsm))
tmp[, trsm:='HSX']
dw <- rbind(dw,tmp)

cat(" \n -------------------------------- \n Generate weights \n -------------------------------- \n")

dw <- subset(dw, select=c(trsm,year,N_Inf_M), year %in% c(2014:2018))
dw <- dw[, list(year=year,w=N_Inf_M/sum(N_Inf_M)),by=c('trsm')]

cat(" \n -------------------------------- \n Define migrant groups \n -------------------------------- \n")

geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

dind <- data.table(dind)
if(args$sens==F){
	dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
}else if(args$sens=='mp'){
	dinf <- subset(dinf,select=c('id','u','estsctodiagMedian','hiv_pos_d','hiv1_neg_d'))
}else{
  dinf <- subset(dinf,select=c('id','estsctodiagMedian','estsctodiag30','estsctodiag40','hiv_pos_d'))
}
dinf <- unique(dinf)
if(args$sens=='mp'){
	dinf <- subset(dinf,!is.na(hiv1_neg_d)) # just keep patients with a last negative test for sensitivity analysis
}
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D','INF_CNTRY','MIG_D')),by.x='id',by.y='PATIENT',all.x=T)
do <- data.table(dinf)
if(args$sens==30){
	do[, time:=estsctodiag30]
}else if(args$sens==40){
	do[, time:=estsctodiag40]
}else if(args$sens=='mp'){
	do[, estsctodiagMedian:=u/2] # use midpoint between time at risk (last neg test and first pos test)	
	do[, time:=estsctodiagMedian]
}else{
	do[, time:=estsctodiagMedian]
}

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

cat(" \n -------------------------------- \n Make summary table of characteristics \n -------------------------------- \n")

tab <- data.table(var1='TRANSM',var2='TOTAL',table(da$TRANSM))
tab[, V2:=V1]

tab <- rbind(tab,data.table(var1='TRANSM',var2='BIRTH_PLACE',table(da$mwmb,da$TRANSM)))
tab <- tab[, c('var1','var2','V1','V2','N')]

tab[var2=='TOTAL',V1:=var2]
tab <- subset(tab,!(N==0))
tab_all <- copy(tab)

# make synthetic subset of data
de <- subset(do,do$INF_D>=2010 & do$INF_D<2013 & TRANSM %in% c('HSX','MSM'))
dexcl <- da

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
tab <- subset(tab,is.na(N_excluded)) # drop the pts infected outside NL
saveRDS(tab,file=file.path(outdir, paste0("characteristics_patients_undiagnosed.RDS")))

## MSM model ----
cat(" \n -------------------------------- \n MSM model: Make stan data \n -------------------------------- \n")

args$trsm <- 'MSM'

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2013)
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

time_m_ci_total <- do[, list(q=quantile(time[INF_D>=2010 & INF_D<2013],prob=ps,na.rm=T),
													 q_label=p_labs)]
time_m_ci_total <- dcast(time_m_ci_total,.~q_label,value.var="q")

data_mg <- data.table(trsm='MSM',mgid=dt$mgid,time_to_diagnosis=dt$time)
data_mg[, migrant_group:=factor(mgid,
														 levels=c(1,2,3,4,5),
														 labels=c('W.Europe,N.America,Oceania','E. & C. Europe','S. America & Caribbean',
														 				 'NL','Other'))]
saveRDS(data_mg,file=file.path(outdir, paste0('time_to_diagnosis_birthplace-',job_tag,"-",args$trsm,'.rds')))
	
q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mwmb)), idx_to_obs_array=dt$mgid, y = dt$time, 
										 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$N_diag)
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
dat <- crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- crossing(year=seq(2014,2018,1))
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
dat <- crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- crossing(year=seq(2014,2018,1))
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

dt <- subset(do,TRANSM==args$trsm & do$INF_D>=2010 & do$INF_D<2013)
dt[mwmb=='G4',mgid:=1]
dt[mwmb=='G5',mgid:=2]
dt[mwmb=='NL',mgid:=3]
dt[mwmb=='Other',mgid:=4]
N_diag <- subset(n_diag,TRANSM==args$trsm)

data_mg <- data.table(trsm='HSX',mgid=dt$mgid,time_to_diagnosis=dt$time)
data_mg[, migrant_group:=factor(mgid,
																levels=c(1,2,3,4),
																labels=c('Sub-Saharan Africa','S. America & Caribbean',
																				 'NL','Other'))]
saveRDS(data_mg,file=file.path(outdir, paste0('time_to_diagnosis_birthplace-',job_tag,"-",args$trsm,'.rds')))

stan_data <- list()
q50 <- quantile(dt$time,probs=c(0.5))
q80 <- quantile(dt$time,probs=c(0.8))
stan_data <- list( n=nrow(dt), r=length(unique(dt$mwmb)), idx_to_obs_array=dt$mgid, y = dt$time, 
									 log_wb_quantiles_priormean=c(log(q50),log(q80)-log(q50)),N_diag=N_diag$N_diag)
#save(stan_data, file=file.path(outdir, 'stanin.RData'))
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
dat <- crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- crossing(year=seq(2014,2018,1))
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
dat <- crossing(year=seq(2014,2018,1),month=seq(1,12,1))
if(args$weights=='ECDC'){
	dat <- crossing(year=seq(2014,2018,1))
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

## make figures ----
cat(" \n -------------------------------- \n Make figures \n -------------------------------- \n")

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

cat(" \n -------------------------------- \n Make plot of undiagnosed by year \n -------------------------------- \n")

do$mlab <- factor(do$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))
do$mlab2 <- factor(do$mlab, levels=c('Other','Sub-Saharan\nAfrica','S. America &\n Caribbean','E. & C. Europe','W.Europe,\nN.America,\nOceania','NL'))
do$trsm <- factor(do$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

col_msm = pal_npg("nrc")(9)[c(1:5)]
col_hsx = pal_npg("nrc")(9)[c(8,6,7,9)]

#col_msm = pal_lancet("lanonc")(9)[c(1:5)]
#col_hsx = pal_lancet("lanonc")(9)[c(8,6,7,9)]

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

g_msm <- ggplot(subset(do,TRANSM %in% c('MSM') & INF_D>=2010 & INF_D<2013), aes(time)) +
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
ggsave(file.path(outdir,'msm_timetodiag_cdfs_facets.png'),g_msm,w=25, h=12)

dat_hsx <- subset(ds,trsm=='Amsterdam heterosexuals')
dat_hsx$mlab <- factor(dat_hsx$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))
g_hsx <- ggplot(subset(do,TRANSM %in% c('HSX') & INF_D>=2010 & INF_D<2013), aes(time)) +
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
ggsave(file.path(outdir,'hsx_timetodiag_cdfs_facets.png'),g_hsx,w=25, h=12)

g_cdf <- ggarrange(g_msm,g_hsx,ncol=1,align="hv")
ggsave(file=file.path(outdir,'timetodiag_cdfs_facets.png'),g_cdf,w=35, h=25)


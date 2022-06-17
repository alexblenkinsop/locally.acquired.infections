cat(" \n -------------------------------- \n \n Running post-processing.r \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(abind, quietly = TRUE))


`%notin%` <- Negate(`%in%`)

# for testing
if(0){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210414r_cmdstan'
	args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210414r_cmdstan-rho_infd_pre_infd_em_infd_1000sg_2015-2019_MSM'
	args_dir[['job_tag']] <- 'rho_infd_pre_infd_em_infd_1000sg_2015-2019_MSM'
	args_dir[['numb_chains']] <- 3
	args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['trsm']] <- 'MSM'
}

# save args for report before loading those from running session 
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFile')	
	stopifnot(args_line[[3]]=='-out_dir')
	stopifnot(args_line[[5]]=='-job_tag')
	stopifnot(args_line[[7]]=='-numb_chains')	
	stopifnot(args_line[[9]]=='-source_dir')	
	stopifnot(args_line[[11]]=='-trsm')	
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['out_dir']] <- args_line[[4]]
	args_dir[['job_tag']] <- args_line[[6]]
	args_dir[['numb_chains']] <- args_line[[8]]
	args_dir[['source_dir']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
} 

# start script
args_dir[['work_dir']] <- getwd()
source(file.path(args_dir$source_dir, 'R', 'post-processing-summary-functions.R'))

cat(" \n --------------------------------  with arguments -------------------------------- \n")
str(args_dir)

do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
cat(paste("\n", nrow(do),"/",args_dir$numb_chains, "chains are finished \n"))

outfile.base <- unique( do[, file.path(dirname(dirname(F)), paste0(args_dir$stanModelFile,'-',args_dir$job_tag))] )
stopifnot(length(outfile.base)==1 )

#	load all input variables for this analysis run
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes
str(args)

cat(" \n -------------------------------- load: fit -------------------------------- \n")
#	reading job output, merge separate stanfits into one consolidated stanfit object
rf <- vector('list', nrow(do))
median_lp_ <- vector('numeric', nrow(do))
for(i in seq_len(nrow(do)))
{
	cat('Loading output in ',do[i,F],'\n')
	z <- load(do[i,F])
	stopifnot('fit' %in% z)
	median_lp_[i] = median(rstan:::extract(fit)$lp__)
	if(all(rstan::summary(fit)$summary[,1] == 0) & all(is.na(rstan::summary(fit)$summary[,2]))) next
	rf[[i]] <- fit
}
fit <- rstan:::sflist2stanfit(rf[lapply(rf,length)>0])
re <- rstan::extract(fit)

cat(" \n -------------------------------- save: fit -------------------------------- \n")
file = paste0(outfile.base,'-stanout-fit.RDS')
cat("\n save file:", file)
io_saveRDS(fit, args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
fit <- NULL
rf <- NULL
gc()

cat(" \n -------------------------------- save: transmission pars -------------------------------- \n")
file =paste0(outfile.base,'-stanout-transmission_pars.RDS')
cat("\n save file:", file)
cat("\n saving objects:", names(re))
io_saveRDS(re, args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
re <- NULL
gc()

cat(" \n -------------------------------- load: generated quantities -------------------------------- \n")
do <- data.table(F=list.files(args_dir$out_dir, pattern='_stangqs.RDS', recursive=TRUE, full.name=TRUE))
do[, STANMF:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\1',basename(F))]
do[, JOB_TAG:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\2',basename(F))]
do[, JOB_ID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\3',basename(F))]

if(nrow(do)!=0){
  rf.gqs <- list()
  re.gqs = list()
  
  chains = as.numeric(unique(gsub(".*\\[(.+)\\].*", "\\1", do[,JOB_ID])))
  icases = as.numeric(unique(gsub(".*_indexcase(.+)_.*", "\\1", do[,JOB_ID])))
  
  for(i in chains)
  {
    
    rf.gqs[[i]] <- list()
    re.gqs[[i]] = list()
    
    for(icase in icases){
      index = which( grepl(paste0("_indexcase", icase, "_"), do[,JOB_ID]) & grepl(paste0("\\[", i, "\\]"), do[,JOB_ID]) )
      cat('Loading output in ', do[index,F],'\n')
      rf.gqs[[i]][[icase]] <- readRDS(do[index,F])
      
      for (var in names(rf.gqs[[i]][[icase]])){
      	cat('Adding ', var,'\n')
      	re.gqs[[i]][[var]] =  array(unlist(lapply(rf.gqs[[i]], "[[", var)), dim = c(dim(rf.gqs[[i]][[1]][[var]]), length(icases)))
      }
      #rf.gqs[[i]][[icase]] <- NULL
    }
  }
  
  #vars = names(re.gqs[[chains[1]]])[names(re.gqs[[chains[1]]]) %notin% names(re)]
  #vars = names(re.gqs[[chains[1]]])
  vars = names(re.gqs[[chains[1]]])[names(re.gqs[[chains[1]]]) %in% c("actual_cs" , "actual_cs_unobs", "N_sgs_unobs" ,"N_sgs_e","pr_unsampled","unsampled_ch",
  																																		"origins_subgraphs_x","origins_subgraphs_e","subt_x","subt_e","obs_cs_pre" , "obs_cs_post",
  																																		"origins_ind_x","origins_ind_e")]
  
  for (var in vars){
    listvar = list( do.call("abind",list(lapply(re.gqs, "[[", var), along = 1)) )
    names(listvar) = var
    #re <- c(re, listvar)
    if(var==vars[1]){
    	re <- listvar
    }else{
    	re <- c(re, listvar)
    }
  }
  
  gc()
}

if(length(unique(lapply(re, function(x) dim(x)[1])))!=1)
{
  stop('number of iterations in Stan and gqs don\'t match')
}

cat(" \n -------------------------------- save: basic quantities -------------------------------- \n")

#	subtypes ID and subtypes names
file <- list.files(args_dir$out_dir, pattern=paste0('subgraph_sizes_',args_dir$trsm,'.RDS'), recursive=TRUE, full.name=TRUE)
freqs <- readRDS(file[1])
ds <- data.table(subtypes=seq_along(unique(freqs$ST)), subtypes_name=unique(freqs$ST))
icases_map <- data.frame(id=seq(1,length(icases),1),index_cases=icases)

basic <- list(subtypes=ds$subtypes,
              ds=ds,
							icases_map=icases_map,
              stan.data=stan.data,
              freqs=freqs,
              dind=dind,
              cases=cases,
              with_subtypes=with_subtypes,
              outfile.base=outfile.base,
							trsm = args_dir$trsm,
              args=args)
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n save file:", file)
io_saveRDS(basic, args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
basic <- NULL
gc()

cat(" \n -------------------------------- save: pmfsum gqs -------------------------------- \n")
pmfsum_gqs = c("pmfsum_pr","pmfsum","pmfqsum")
if(all(pmfsum_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-pmfsum-gqs.RDS')
	io_saveRDS(re[pmfsum_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[pmfsum_gqs] <- NULL
}

cat(" \n -------------------------------- save: pmf gqs -------------------------------- \n")
pmf_gqs = c("cs_actual_pmf","cs_obs_pre_pmf","cs_obs_post_pmf")
if(all(pmf_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-pmf-gqs.RDS')
	io_saveRDS(re[pmf_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[pmf_gqs] <- NULL
}

cat(" \n -------------------------------- save: subgraphs gqs -------------------------------- \n")
subgraph_gqs = c("actual_cs" , "actual_cs_unobs", "N_sgs_unobs" ,"N_sgs_e")
if(all(subgraph_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-subgraphs-gqs.RDS')
	io_saveRDS(re[subgraph_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[subgraph_gqs] <- NULL
}

cat(" \n -------------------------------- save: origins gqs -------------------------------- \n")
origins_gqs = c("origins_subgraphs_x","origins_subgraphs_e")
if(all(origins_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-origins-gqs.RDS')
	io_saveRDS(re[origins_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[origins_gqs] <- NULL
}

cat(" \n -------------------------------- save: posterior predictive check gqs -------------------------------- \n")
csobs_gqs = c("obs_cs_pre" , "obs_cs_post")
if(all(csobs_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-posteriorpredictivecheck-gqs.RDS')
	io_saveRDS(re[csobs_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[csobs_gqs] <- NULL
}

cat(" \n -------------------------------- save: subtypes gqs -------------------------------- \n")
subt_gqs = c("subt_x","subt_e")
if(all(subt_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-subtypes-gqs.RDS')
	io_saveRDS(re[subt_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[subt_gqs] <- NULL
}


cat(" \n -------------------------------- save: pr(unsampled) gqs -------------------------------- \n")
unsamp_gqs = c("pr_unsampled","unsampled_ch")
if(all(unsamp_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-pr_chain_unsampled-gqs.RDS')
	io_saveRDS(re[unsamp_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[unsamp_gqs] <- NULL
}

cat(" \n -------------------------------- save: birthplace gqs -------------------------------- \n")
bplace_gqs = c("origins_ind_x","origins_ind_e")
if(all(bplace_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-bplace-gqs.RDS')
	io_saveRDS(re[bplace_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[bplace_gqs] <- NULL
}

cat(" \n -------------------------------- save: newcases gqs -------------------------------- \n")
cases_gqs = c("N_inf_x","N_inf_e")
if(all(cases_gqs %in% names(re))){
	file = paste0(outfile.base,'-stanout-cases-gqs.RDS')
	io_saveRDS(re[cases_gqs], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[cases_gqs] <- NULL
}

cat(" \n -------------------------------- save: iterations gqs -------------------------------- \n")
iter = c("iterations")
if(all(iter %in% names(re))){
	cat("\n ---saving iterations in which pmf does not add to 1--- \n")
	file = paste0(outfile.base,'-stanout-iterationsflag-gqs.RDS')
	io_saveRDS(re[iter], args_dir[['work_dir']], dirname(file), basename(file), check_if_saved_n=10)
	re[iter] <- NULL
}


cat(" \n -------------------------------- \n \n End post-processing.r \n \n -------------------------------- \n")


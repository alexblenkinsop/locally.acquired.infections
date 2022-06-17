# post-processing-transmission-pars.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running post-processing-transmission-pars.R\n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(bayesplot, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggpubr, quietly = TRUE))

# for dev purposes
if(0){
  args_dir <- list()
  args_dir[['stanModelFile']] <- 'branching_process_210217a_cmdstan'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210217a_cmdstan-bysubtype_randomeffects_R0_indexcaseout_2015-2019_MSM'
  args_dir[['job_tag']] <- 'bysubtype_randomeffects_R0_indexcaseout_2015-2019_MSM'
  args_dir[['source_dir']] <- '~/git/bpm'
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-stanModelFile')	
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-job_tag')
  stopifnot(args_line[[7]]=='-source_dir')
  args_dir <- list()
  args_dir[['stanModelFile']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['job_tag']] <- args_line[[6]]
  args_dir[['source_dir']] <- args_line[[8]]
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'post-processing-summary-functions.R'))

outfile.base <- paste0(args_dir$out_dir, "/",
                       args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
plot.pars.basic <- readRDS(file)
file <- paste0(outfile.base,'-stanout-transmission_pars.RDS')
cat("\n read RDS:", file)
plot.pars.trmspars <- readRDS(file)

tryCatch(
  if("rho" %in% names(plot.pars.trmspars)) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["rho"]])
    g_rho <- make_posterior_intervals_R0( tmp, 
                                         "rho", 
                                         xtick = NULL, 
                                         xintercept=NULL, 
                                         xmin=NULL, 
                                         xlab=expression(rho),
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)


tryCatch(
  if("r0" %in% names(plot.pars.trmspars) & !plot.pars.basic$with_subtypes) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["r0"]])
    g_R0 <- make_posterior_intervals_R0( tmp, 
                                         "r0", 
                                         xtick = NULL, 
                                         xintercept=1, 
                                         xmin=NULL, 
                                         xlab=expression(R[0]),
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)


tryCatch(
  if(any(c("r0","log_r0_sbts","logit_r0_sbts") %in% names(plot.pars.trmspars)) & plot.pars.basic$with_subtypes) 
  {			
    if("r0" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["r0"]])
    if("log_r0_sbts" %in% names(plot.pars.trmspars)) tmp <- as.matrix(exp(plot.pars.trmspars[["log_r0_sbts"]]))
    if("logit_r0_sbts" %in% names(plot.pars.trmspars)) tmp <- as.matrix(gtools::inv.logit(plot.pars.trmspars[["logit_r0_sbts"]]))
    g_R0 <- make_posterior_intervals_R0( tmp, 
                                         "r0", 
                                         xtick = plot.pars.basic$ds$subtypes_name, 
                                         xintercept=1, 
                                         xmin=NULL, 
                                         xlab=expression(R[0]),
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch(
  if(any(c("sd_r0","log_r0_sd","logit_r0_sd") %in% names(plot.pars.trmspars))) 
  {			
    if("sd_r0" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["sd_r0"]])
    if("log_r0_sd" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["log_r0_sd"]])
    if("logit_r0_sd" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["logit_r0_sd"]])
    g_sd_R0_rnde <- make_posterior_intervals_R0(tmp, 
                                                "sd_r0", 
                                                xtick = NULL, 
                                                xintercept=NULL, 
                                                xmin=NULL, 
                                                xlab="sd_r0",
                                                outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)


tryCatch(
  if("vmr_minus_one" %in% names(plot.pars.trmspars) & !plot.pars.basic$with_subtypes) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["vmr_minus_one"]])
    g_vmr_minus_one <- make_posterior_intervals_R0( tmp, 
                                         "vmr_minus_one", 
                                         xtick = NULL, 
                                         xintercept=NULL, 
                                         xmin=NULL, 
                                         xlab='vmr_minus_one',
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch(
  if(any(c("vmr_minus_one","vmr_minus_one_sbts") %in% names(plot.pars.trmspars)) & plot.pars.basic$with_subtypes) 
  {			
    if("vmr_minus_one" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["vmr_minus_one"]])
    if("vmr_minus_one_sbts" %in% names(plot.pars.trmspars)) tmp <- as.matrix(plot.pars.trmspars[["vmr_minus_one_sbts"]])
    g_vmr_minus_one <- make_posterior_intervals_R0( tmp, 
                                         "vmr_minus_one", 
                                         xtick = plot.pars.basic$ds$subtypes_name, 
                                         xintercept=NULL, 
                                         xmin=NULL, 
                                         xlab='vmr_minus_one',
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch(
  if("sd_vmr_minus_one" %in% names(plot.pars.trmspars)) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["sd_vmr_minus_one"]])
    g_sd_vmr_rnde <- make_posterior_intervals_R0(tmp, 
                                             "sd_vmr_minus_one", 
                                             xtick = NULL, 
                                             xintercept=1, 
                                             xmin=NULL, 
                                             xlab="sd_vmr_minus_one",
                                             outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch(
  if("kappa" %in% names(plot.pars.trmspars) & !plot.pars.basic$with_subtypes) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["kappa"]])
    g_kappa <- make_posterior_intervals_R0( tmp, 
                                         "kappa", 
                                         xtick = NULL, 
                                         xintercept=1, 
                                         xmin=NULL, 
                                         xlab=expression(kappa),
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch(
  if("kappa" %in% names(plot.pars.trmspars) & plot.pars.basic$with_subtypes) 
  {			
    tmp <- as.matrix(plot.pars.trmspars[["kappa"]])
    g_kappa <- make_posterior_intervals_R0( tmp, 
                                         "kappa", 
                                         xtick = plot.pars.basic$ds$subtypes_name, 
                                         xintercept=1, 
                                         xmin=NULL, 
                                         xlab=expression(kappa),
                                         outfile.base = outfile.base)
  }
  , error = function(err) { warning(err) }
)

tryCatch({
  gpl <- list()	
  if(exists('g_rho'))
    gpl[[length(gpl)+1]] <- g_rho
  if(exists('g_R0'))
    gpl[[length(gpl)+1]] <- g_R0
  if(exists('g_sd_R0_rnde'))
    gpl[[length(gpl)+1]] <- g_sd_R0_rnde
  if(exists('g_vmr_minus_one'))
    gpl[[length(gpl)+1]] <- g_vmr_minus_one
  if(exists('g_sd_vmr_rnde'))
    gpl[[length(gpl)+1]] <- g_sd_vmr_rnde
  if(exists('g_kappa'))
    gpl[[length(gpl)+1]] <- g_kappa

  make_transmission_parameter_summary_plot(gpl, outfile.base)
})


cat(" \n -------------------------------- \n \n End post-processing-transmission-pars.R\n \n -------------------------------- \n")


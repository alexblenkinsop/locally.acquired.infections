# post-processing-knit-report.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n post-processing-knit-report.R \n \n -------------------------------- \n")

suppressMessages(library(rmarkdown, quietly = TRUE))

#	for dev purposes: olli
if(1)
{
  args_dir <- list()
  args_dir[['stanModelFile']] <- 'branching_process_201209_cmdstan'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_201209_cmdstan-heterogeneity_sbt'
  args_dir[['job_tag']] <- 'heterogeneity_sbt'
  args_dir[['rmd_file']] <- "scripts/post-processing-make-report.Rmd"
  args_dir[['source_dir']] <- '~/git/bpm'
  args_dir[['report_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/reports'
}

#	for runtime
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-rmd_file')
  stopifnot(args_line[[3]]=='-stanModelFile')	
  stopifnot(args_line[[5]]=='-out_dir')
  stopifnot(args_line[[7]]=='-job_tag')
  stopifnot(args_line[[9]]=='-report_dir')
  stopifnot(args_line[[11]]=='-source_dir')
  args_dir <- list()
  args_dir[['rmd_file']] <- args_line[[2]]
  args_dir[['stanModelFile']] <- args_line[[4]]
  args_dir[['out_dir']] <- args_line[[6]]
  args_dir[['job_tag']] <- args_line[[8]]
  args_dir[['report_dir']] <- args_line[[10]]
  args_dir[['source_dir']] <- args_line[[12]]
} 

## start script
cat(" \n -------------------------------- \n with post-processing arguments \n -------------------------------- \n")
args_dir[['report_path_to_file']] <- file.path(args_dir[['report_dir']], paste0("report_", args_dir[['stanModelFile']], "-", args_dir[['job_tag']], ".html") )
args_dir[['rmd_path_to_file']] <- file.path(args_dir$source_dir, args_dir[['rmd_file']])

str(args_dir)

outfile.base <- paste0(args_dir$out_dir, "/", args_dir$stanModelFile , "-", args_dir$job_tag)

##	make report
cat(paste("\n ----------- create report ----------- \n"))


rmarkdown::render( args_dir[['rmd_path_to_file']], 
                   output_file= args_dir[['report_path_to_file']], 
                   params = list(
                     stanModelFile = args_dir$stanModelFile,
                     job_dir= args_dir$out_dir,
                     job_tag= args_dir$job_tag
                   ))

cat(paste("\n -------------------------------- \n \n post-processing-knit-report.R \n \n -------------------------------- \n"))

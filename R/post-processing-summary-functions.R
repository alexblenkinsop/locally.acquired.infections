io_saveRDS <- function(obj, work_dir, out_dir, base_name, check_if_saved_n = 0)
{
  cat('\nSave to file ',file.path(out_dir, base_name),'...')	
  tmp <- check_if_saved_n
  repeat
  {
    #   comp_9 = xzfile(tmp, compression = 9)
    # 	saveRDS(fit.gqs, comp_9)	
    tryCatch(
      {
        saveRDS(obj, file=file.path(work_dir, base_name))
        if(work_dir!=out_dir)
        {
          file.copy(file.path(work_dir, base_name), 
                    file.path(out_dir, base_name), 
                    overwrite = TRUE, 
                    recursive = FALSE,
                    copy.mode = TRUE, 
                    copy.date = TRUE
          )							
        }
      }, error = function(err) { warning(err) } )
    if(check_if_saved_n<1)
      break
    check_if_saved <- try(readRDS(file=file.path(out_dir, base_name)))		
    if(!'try-error'%in%class(check_if_saved))
      break	
    tmp <- tmp-1
    if(tmp<=0)
    {
      stop('Failed to save ',check_if_saved_n,' times')
    }
  }	
}

make_posterior_summaries <- function(par,subtypes)
{
  cat("\n ----------- summarise posterior ----------- \n")
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  #	do calculations per subtype to minimise memory footprint
  
  stopifnot( dim(par)[2]==length(subtypes$subtypes) )
  par <- as.data.table( reshape2::melt( par ) )
  setnames(par, 1:3, c('iteration','subtypes','value'))
  
  # 	summarise
  par <- par[, list(qs= quantile(value, prob=ps), qlab=p_labs), by=c('subtypes')]
  par <- dcast.data.table(par, subtypes~qlab, value.var='qs')
  
  # make human readable labels
  par <- merge(par, subtypes, by='subtypes')
  setnames(par, c('subtypes'), c('subtype'))
  
  return(par)
}

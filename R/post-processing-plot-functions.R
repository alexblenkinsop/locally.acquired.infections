make_posterior_intervals_R0 <- function(tmp,pars,xtick=NULL,xintercept=NULL,xmin=NULL,
                                        xlab=NULL,outfile.base){
  cat("\n ----------- make ", pars ," figure ----------- \n")
  if(is.null(xlab)){
    xlab <- pars
  }
  
  if(is.null(xtick)){
    xtick <- 1:ncol(tmp)
  }
  
  colnames(tmp) <- xtick
  
  gv <- bayesplot::mcmc_intervals(tmp,prob = .9) +
    geom_vline(xintercept=xintercept) +
    labs(x=pars) +
    theme_bw() +
    scale_x_log10()
  
  ggsave(paste0(outfile.base,'-',pars,'.png'),gv,width=4,height=(length(xtick)/2 + 0.3))
  
  ###
  
  gh <- bayesplot::mcmc_intervals(tmp,prob = .9) +
    geom_vline(xintercept=xintercept) +
    labs(x=xlab) +
    coord_flip() +
    scale_y_discrete(position = "right",labels=colnames(tmp)) +
    scale_x_log10() +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=18,angle=90),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14,angle=90, vjust = 0.5, hjust=0))
  
  if(!is.null(xmin)){
    gh <- gh + scale_x_continuous(expand=c(0,0),limits = c(0, NA))
  }
  ggsave(paste0(outfile.base,'-',pars,'.png'),gh,width=4,height=(length(xtick)/2 + 0.3))
  
  return(gh)
}

make_transmission_parameter_summary_plot <- function(ggplot.list, outfile.base){
  cat("\n ----------- make model_transmission_parameters figure ----------- \n")
  
  g <- ggarrange(plotlist=ggplot.list,ncol=1,align="v",heights=c(1.3,rep(1,length(ggplot.list)-1)))
  ggsave(paste0(outfile.base,'-model_transmission_parameters','.png'), g, w = 8, h = 2*(length(ggplot.list)))
}

make_gqs_plots <- function(tmp,xpar,ylab=NULL,outfile.base,filename){
  cat("\n ----------- make ", filename ," figure ----------- \n")
  g <- ggplot(tmp, aes(x=xpar)) +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=0.3) +
    geom_point(aes(y=M)) +			
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() +
    theme(axis.text.x=element_text(size=12,angle=90, vjust = 0.5, hjust=0)) +
    labs(y=ylab, x='')
  ggsave(paste0(outfile.base,'-',filename,'.png'), w=6, h=6)
  
  return(g)
}

make_gqs_plots_st <- function(tmp,xpar,ylab=NULL,outfile.base,filename){
  cat("\n ----------- make ", filename ," figure ----------- \n")
  g <- ggplot(ans, aes(x=xpar)) +
    geom_errorbar(aes(ymin=CL, ymax=CU,col=subtypes_name), width=0.3,position=position_dodge(width=0.6), width=0.3) +
    geom_point(aes(y=M,col=subtypes_name), position=position_dodge(width=0.6)) +			
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() +
    theme(axis.text.x=element_text(size=12,angle=90, vjust = 0.5, hjust=0),
          legend.position="bottom") +
    labs(y=ylab, x='',col='')
  ggsave(paste0(outfile.base,'-',filename,'.png'), w=6, h=6)
  
  return(g)
}

make_plot_ratio <- function(tmp,xpar,ylab=NULL,sqrt_y=NULL,ylim=NULL,outfile.base,filename){
  cat("\n ----------- make ", filename ," figure ----------- \n")
  g <- ggplot(tmp, aes(x=xpar)) +
    geom_errorbar(aes(ymin=CL, ymax=CU), width=0.3) +
    geom_point(aes(y=M)) +	
    theme_bw() +
    theme(axis.text.x=element_text(size=12,angle=90, vjust = 0.5, hjust=0)) +
    labs(y=ylab, x='')
  if(!is.null(sqrt_y)){
    g <- g + scale_y_sqrt()
  }
  if(!is.null(ylim)){
    g <- g + coord_cartesian(ylim=c(ylim[1], ylim[2]))
  }
  ggsave(paste0(outfile.base,'-',filename,'.png'), w=6, h=6)
  
  return(g)
}


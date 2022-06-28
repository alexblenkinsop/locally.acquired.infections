extract_subgraphs = function(indir.phsc){
  infiles <- data.table(F=list.files(indir.phsc, pattern='rda$', full.names=TRUE, recursive=TRUE))
  infiles <- subset(infiles, grepl('subgraphs_',F))
  #    label by transmission group
  infiles[, SELECT:= gsub('^(.*)subgraphs_([A-Za-z]+).rda','\\2',basename(F))]
  #    ST stores the subtype, and ST_CLADE the large subtree (for ST B only)
  infiles[, ST:= gsub('^.*subtype_([^_]+)_.*\\.rda','\\1',basename(F))]
  infiles[, ST_CLADE:= as.integer(gsub('[^c]*c([0-9]+)','\\1',ST))]
  infiles[, ST:= gsub('([^c]*)c([0-9]+)','\\1',ST)]
  #    ID of bootstrap replicate, 000 for analysis on real alignment
  infiles[, REP:= gsub('^.*wOutgroup_([0-9]+)_.*\\.rda','\\1',basename(F))]
  dsubgraphtaxa <- infiles[, {
    infile <- F
    cat('Process',infile,'\n')
    load(infile)
    if(length(subgraphs)==1){
      subgraph.names <- rep(subgraphs[[1]]$subgraph.name, length(subgraphs[[1]]$tip.label))
      subgraph.taxa <- subgraphs[[1]]$tip.label
      subgraph.parent.state <- subgraphs[[1]]$subgraph.parent.state
    }  else{
      subgraph.names <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.name, length(subgraph$tip.label)) ))
      subgraph.taxa <- unlist(sapply( subgraphs, function(subgraph)  subgraph$tip.label))
      subgraph.parent.state <- unlist(sapply( subgraphs, function(subgraph)  rep(subgraph$subgraph.parent.state, length(subgraph$tip.label))))
    }
    list(NAME=subgraph.names,
         TAXA= subgraph.taxa,
         ORIGINHOST= subgraph.parent.state
    )
  }, by=c('ST','ST_CLADE','REP','SELECT')]
  #    add meta data from taxa names
  regex.tip.label <- '^([A-Za-z]+)_+(T[0-9]+)_([0-9]+)_([a-zA-Z]+)_([a-zA-Z]+)_([a-zA-Z]+)-([a-zA-Z]+)_([a-zA-Z_]+)_([0-9]+)-([0-9-]+)-([0-9-]+)_([0-9]+)$'
  dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
  dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\9',TAXA))]
  dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]
  
  return(dsubgraphtaxa)
}

hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x	
}

add_cd4_counts = function(dsubgraphtaxa,infile.cd4,startd,trsm){
  
  dat <- read.csv(infile.cd4,header=T)
  dat <- subset(dat,CD4_U==1)
  dat$CD4_D <- as.Date(dat$CD4_D,format=c("%Y-%m-%d"))
  
  # keep first CD4 count (and first obs if multiple from same day)
  tmp <- dat %>%
    group_by(PATIENT) %>%
    filter(CD4_D == min(CD4_D)) %>% 
    filter(1:n() == 1)
  setnames(tmp,c('CD4_D','CD4_V'),c('CD4_first_D','CD4_first_V'))

  dsubgraphtaxa <- merge(dsubgraphtaxa,subset(tmp,select=c('PATIENT','CD4_first_D','CD4_first_V')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)
  
  # just look at diagnoses from end 2018
  #ds <- subset(dsubgraphtaxa,REP=='000' & SELECT==paste0('Ams',trsm) & HIV1_POS_D>=2018.5)
  
  # nb there is one person with a CD4 count from 2016 though diagnosed since 2017
  
  # get proportion of these diagnoses with CD4 <350
  #ds[, lowcd4:=0]
  #ds[CD4_first_V<350, lowcd4:=1]
  
  #du <- ds[, list(undiag=sum(lowcd4)/nrow(ds))]
}

add_viral_loads = function(dsubgraphtaxa,infile.rna,infile.indinfo,startd){
  
  dat <- read.csv(infile.rna,header=T)
  dat$RNA_D <- as.Date(dat$RNA_D,format=c("%Y-%m-%d"))
  
  # get last viral load before start of analysis
  startd <- as.Date(0, origin = paste0(startd,"-01-01"),format="%Y-%m-%d")
  tmp <- subset(dat, RNA_D<as.Date(startd))
  
  tmp <- tmp %>%
    group_by(PATIENT) %>%
    filter(RNA_D == max(RNA_D))

  dsubgraphtaxa <- merge(dsubgraphtaxa,subset(tmp,select=c('PATIENT','RNA_D','RNA_V')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)
  
  tmp <- dat %>%
    group_by(PATIENT) %>%
    filter(RNA_D == max(RNA_D))
    setnames(tmp,c('RNA_D','RNA_V'),c('RNA_last_D','RNA_last_V'))
  
  dsubgraphtaxa <- merge(dsubgraphtaxa,subset(tmp,select=c('PATIENT','RNA_last_D','RNA_last_V')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)

  # get death data in case some individuals died before analysis period (no potential to transmit)
  dbas <- read.csv(infile.indinfo,header=T)
  dbas <- subset(dbas,select=c(PATIENT,DEATH_D))
  dbas[,'DEATH_D'] <- as.Date(dbas[,'DEATH_D'],format="%Y-%m-%d")
  
  dsubgraphtaxa <- merge(dsubgraphtaxa,dbas,by.x=c('ID'),by.y=c('PATIENT'),all.x=T)
  
  # flag patients virally suppressed by start date
  dsubgraphtaxa[, supp:=1]
  dsubgraphtaxa[RNA_V>100 & HIV1_POS_D<as.Date(startd), supp:=0]
  # add patients who died before start date to virally suppressed
  dsubgraphtaxa[DEATH_D<startd, supp:=1]
  
  return(dsubgraphtaxa)
  
}

correct_misrecorded_dates = function(dsubgraphtaxa,infile.indinfo){
  
  dbas <- read.csv(infile.indinfo,header=T)
  dbas <- subset(dbas,select=c(PATIENT,RECART_D,HIV1_POS_D,HIV1_POS_D_A))
  dbas[,'HIV1_POS_D'] <- as.Date(dbas[,'HIV1_POS_D'],format="%d/%m/%Y")
  dbas[,'RECART_D'] <- as.Date(dbas[,'RECART_D'],format="%Y-%m-%d")
  setnames(dbas,c('HIV1_POS_D','RECART_D'),c('HIV1_POS_DATE','RECART_DATE'))
  dbas <- merge(dsubgraphtaxa, dbas, by.x=c('ID'),by.y=('PATIENT'),all.x=T)
  # Replace all RECART_D from 1911 with HIV1_POS_D
  dbas[RECART_D<HIV1_POS_D & floor(RECART_D)==1911,RECART_D:=HIV1_POS_D + (1/365)]
  # Replace all HIV_POS_D which are uncertain with lower bound
  dbas[HIV1_POS_D_A==1 & RECART_D<HIV1_POS_D,HIV1_POS_D:=HIV1_POS_D_lower]
  # Replace all RECART_D which have been rounded to after HIV_POS_D
  dbas[format(as.Date(RECART_DATE, format="%Y-%m-%d"),"%m")=='07' & 
                  format(as.Date(RECART_DATE, format="%Y-%m-%d"),"%d")=='01' & RECART_D<HIV1_POS_D,RECART_D:=HIV1_POS_D + (1/365)]
  # For remaining individuals, in which not clear which date is correct, assume ART date started after diagnosis date
  dbas[RECART_D<HIV1_POS_D,RECART_D:=HIV1_POS_D + (1/365)]
  return(dbas)
  
}

add_georeg_bplace_ind_data = function(dat,infile.geo){
  geo <- data.table(read.csv(infile.geo))
  geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
  geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
  setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
  
  dat <- merge(dat,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
  
  return(dat)
}

add_georeg_bplace_sg_data = function(dat,dind,infile.geo){
  geo <- data.table(read.csv(infile.geo))
  geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
  geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
  setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
  
  dat <- merge(dat,subset(dind,select=c('ID','ORIGIN'),by=c('ID'),all.x=T))
  dat <- merge(dat,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
  
  return(dat)
}

add_infection_time = function(dat,infile.inftime,dt_inf){

  dinf <- read.csv(infile.inftime,header=T)
  dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d')) #rem pos d
  dinf <- unique(dinf)
  dat <- merge(dat,dinf,by.x='ID',by.y='id',all.x=T)
  
  # add the migrant group IDs from the undiagnosed model
  dat[TRANSM=='MSM' & mwmb=='G1',mg:=1]
  dat[TRANSM=='MSM' & mwmb=='G2',mg:=2]
  dat[TRANSM=='MSM' & mwmb=='G3',mg:=3]
  dat[TRANSM=='MSM' & mwmb=='NL',mg:=4]
  dat[TRANSM=='MSM' & mwmb=='Other',mg:=5]
  dat[TRANSM=='HSX' & mwmb=='G4',mg:=1]
  dat[TRANSM=='HSX' & mwmb=='G5',mg:=2]
  dat[TRANSM=='HSX' & mwmb=='NL',mg:=3]
  dat[TRANSM=='HSX' & mwmb=='Other',mg:=4]
  dat[, mg:=as.character(mg)]
  # add median times to diagnosis
  dat <- merge(dat,dt_inf,by=c('TRANSM','mg'),all.x=T)
  dat$hiv_pos_d <- as.Date(dat$hiv_pos_d,format=c("%Y-%m-%d")) # del
  # estimate infection time
  dat[, INF_D:= HIV1_POS_D - estsctodiagMedian]
  #dat[is.na(INF_D),INF_D:=HIV1_POS_D]
  # where no infection time estimate, subtract median time to diagnosis from diagnosis date
  dat[is.na(INF_D),INF_D:=HIV1_POS_D - `p0.5`]
  return(dat)
  
}

add_infection_time_midpoint_seroconv_diagnosis = function(dat,infile.inftime){
  
  dinf <- data.table(read.csv(infile.inftime,header=T))
  dinf <- subset(dinf,select=c('id','u')) ## u = time between time at risk (last neg test or 15 y/o)
  dinf <- unique(dinf)
  dinf[, u:=u/2] # take mid-point of time at risk and subtract from diagnosis date
  dat <- merge(dat,dinf,by.x='ID',by.y='id',all.x=T)
  dat[, INF_D:= HIV1_POS_D - u]
  dat[is.na(INF_D),INF_D:=HIV1_POS_D]
  return(dat)
  
}

est_prop_undiagnosed = function(infile.undiag,start_d,end_d){
  
  undiag <- data.table(read.csv(infile.undiag))
  
  last_d <- end_d - 1
  if(end_d==2020) last_d <- end_d - 2
  
  # estimate number infected during analysis time period and undiagnosed by 2019 (end of study)
  
  ud <- melt(subset(undiag,select=c('year','N_Inf_M','t_diag_1','t_diag_2','t_diag_3','t_diag_4','t_diag_5','t_diag_6','t_diag_7',
                                    't_diag_7','t_diag_8','t_diag_9','t_diag_10','t_diag_11','t_diag_12','t_diag_13','t_diag_14','t_diag_15')),
             id.vars=c('year','N_Inf_M'))
  
  dat <- data.table(year=seq(start_d,last_d))
  
  dat[,yrs:= 2019 - year]
  for(i in start_d:last_d){
    dat[year==i, N_Inf:= ud$N_Inf_M[ud$year==i][1]]
    dat[year==i,p_undiag:= 1 - sum(unique(ud$value[ud$year==i][1:(2019-i)]))/100]
    dat[year==i, N_undiag:= ud$N_Inf_M[ud$year==i][1]*p_undiag]
  }
  
  # weight probabilities by number infected in each year
  undiagnosed <- dat
  undiagnosed[, infected:=N_Inf]
  undiagnosed[, weight:= N_Inf/sum(N_Inf)]
  undiagnosed[, p_undiag_y:=p_undiag*weight]
  p_undiagnosed <- sum(undiagnosed$p_undiag_y)
  
  return(p_undiagnosed)
}

make_standata_prop_undiagnosed = function(infile.seqdat,infile.seqlabs,infile.inftime,infile.geo,trsm){
  
  load(infile.seqdat)
  load(infile.seqlabs)
  dinf <- read.csv(infile.inftime,header=T)
  geo <- data.table(read.csv(infile.geo))
  
  geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
  geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
  setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
  
  dind <- data.table(dind)
  dind$SEQ <- dind$PATIENT %in% ds$PATIENT
  dinf <- read.csv(infile.inftime,header=T)
  dinf$SEQ <- dinf$id %in% ds$PATIENT
  dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
  dinf <- unique(dinf)
  dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D')),by.x='id',by.y='PATIENT',all.x=T)
  do <- data.table(dinf)
  do[, time:=estsctodiagMedian]
  ## Estimate undiagnosed by migrant groups
  do[, INF_D:=HIV1_POS_D - time]
  do <- merge(do,subset(dind,select=c('PATIENT','ORIGIN')),by.x='id', by.y='PATIENT',all.x=T)
  do <- merge(do,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
  
  do <- map_mwmb_regions(do)
  
  dt <- subset(do,TRANSM==trsm & do$INF_D>=2010 & do$INF_D<2013)
  
  if(trsm=='HSX'){
    n <- max(length(dt$time[dt$mwmb=='G4']), length(dt$time[dt$mwmb=='G5']),
             length(dt$time[dt$mwmb=='NL']),
             length(dt$time[dt$mwmb=='Other']))
    n_r <- c(length(dt$time[dt$mwmb=='G4']), length(dt$time[dt$mwmb=='G5']),
             length(dt$time[dt$mwmb=='NL']),
             length(dt$time[dt$mwmb=='Other']))
    time <- data.table(G4=dt$time[dt$mwmb=='G4'][1:n],
                       G5=dt$time[dt$mwmb=='G5'][1:n],
                       NL=dt$time[dt$mwmb=='NL'][1:n],
                       Other=dt$time[dt$mwmb=='Other'][1:n])
  }else if(trsm=='MSM'){
    n <- max(length(dt$time[dt$mwmb=='G1']), length(dt$time[dt$mwmb=='G2']),
             length(dt$time[dt$mwmb=='G3']),length(dt$time[dt$mwmb=='NL']),
             length(dt$time[dt$mwmb=='Other']))
    n_r <- c(length(dt$time[dt$mwmb=='G1']), length(dt$time[dt$mwmb=='G2']),
             length(dt$time[dt$mwmb=='G3']),length(dt$time[dt$mwmb=='NL']),
             length(dt$time[dt$mwmb=='Other']))
    time <- data.table(G1=dt$time[dt$mwmb=='G1'][1:n],
                       G2=dt$time[dt$mwmb=='G2'][1:n],
                       G3=dt$time[dt$mwmb=='G3'][1:n],
                       NL=dt$time[dt$mwmb=='NL'][1:n],
                       Other=dt$time[dt$mwmb=='Other'][1:n])
  }
  time[is.na(time)] <- Inf
  stan.data <- list( n=n, n_r = n_r, r=length(unique(dt$mwmb)), time = time)
  
  return(stan.data)
}

subgraph_sizes_infdate = function(dsubgraphtaxa,dind,start_d,end_d,trsm,index_in_pre){
  cat(paste('\n Subset data to individuals diagnosed up until ',end_d,' \n'))
  # Remove individuals with unknown HIV positive date
  dsubgraphtaxa <- subset(dsubgraphtaxa, INF_D<end_d & INF_D>1912)
  dind <- subset(dind, INF_D<=end_d)
  
  cat(paste('\n Count index cases and new cases \n'))
  # summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
  dsubgraphsize <- dsubgraphtaxa[, list(icasesart=length(ID[INF_D<start_d & RECART_D>=start_d]),icases=length(ID[INF_D<start_d & supp==0]),jcases=length(ID[INF_D>=start_d & INF_D<end_d])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
  cases <- subset(dsubgraphsize,REP=='000',select=c('NAME','SELECT','ST','ST_CLADE','icases','jcases'))
  cases <- subset(cases,SELECT==paste0('Ams',trsm))
  
  cat('\n Calculate expected index cases from observed using sampling fraction at start of analysis \n')
  dat <- dind[dind$CITY=='Amsterdam' & dind$TRANSM==trsm & dind$INF_D<start_d & dind$RECART_D>=start_d,]
  
  inf.startd <- length(unique(dat$ID)) 
  seq.startd <- length(unique(dat$ID[dat$SEQ==T])) 
  sampling.prob <- seq.startd/inf.startd
  
  cases$m <- round(cases$icases/sampling.prob,0)
  
  mind <- subset(dsubgraphtaxa,REP=='000')[,list(MIND=min(INF_D,na.rm=T),MAXD=max(INF_D,na.rm=T)),by=c('NAME','ST','SELECT','ST_CLADE')]
  
  cases <- merge(cases,mind,by=c('NAME','ST','SELECT','ST_CLADE'))
  cases$period <- "pre-startd"
  cases$period[cases$MIND>=start_d] <- "post-startd"
  
  # distinguish the pre-existing subgraphs with 0 index, 0 new cases
  cases[m==0 & jcases==0 & period=='pre-startd',m:= -1]
  # for subgraphs starting prior to cut-date with no observed index case, assume one unobserved index case
  if(index_in_pre==1){
    cases[m==0 & period=='pre-startd',jcases:= jcases-1]
  }
  cases[m==0 & cases$period=='pre-startd',m:= 1]
  
  # final subgraph sizes from individuals not in ART by start date
  cases$SIZE <- cases$m + cases$jcases
  cases$icases <- NULL
  setnames(cases,'m','icases')
  
  return(list(cases,dind))
  
}

subgraph_sizes_diagdate = function(dsubgraphtaxa,dind,start_d,end_d,trsm,index_in_pre){
  cat(paste('\n Subset data to individuals diagnosed up until ',end_d,' \n'))
  # Remove individuals with unknown HIV positive date
  dsubgraphtaxa <- subset(dsubgraphtaxa, HIV1_POS_D<end_d & HIV1_POS_D>1912)
  dind <- subset(dind, HIV1_POS_D<=end_d)
  
  cat(paste('\n Count index cases and new cases \n'))
  # summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
  dsubgraphsize <- dsubgraphtaxa[, list(icasesart=length(ID[HIV1_POS_D<start_d & RECART_D>=start_d]),icases=length(ID[HIV1_POS_D<start_d & supp==0]),jcases=length(ID[HIV1_POS_D>=start_d & HIV1_POS_D<end_d])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
  
  cases <- subset(dsubgraphsize,REP=='000',select=c('NAME','SELECT','ST','ST_CLADE','icases','jcases'))
  cases <- subset(cases,SELECT==paste0('Ams',trsm))
  
  cat('\n Calculate expected index cases from observed using sampling fraction at start of analysis \n')
  dat <- dind[dind$CITY=='Amsterdam' & dind$TRANSM==trsm & dind$HIV1_POS_D<start_d & dind$RECART_D>=start_d,]
  
  inf.startd <- length(unique(dat$ID)) 
  seq.startd <- length(unique(dat$ID[dat$SEQ==T])) 
  sampling.prob <- seq.startd/inf.startd
  
  cases$m <- round(cases$icases/sampling.prob,0)
  
  mind <- subset(dsubgraphtaxa,REP=='000')[,list(MIND=min(HIV1_POS_D,na.rm=T),MAXD=max(HIV1_POS_D,na.rm=T)),by=c('NAME','ST','SELECT','ST_CLADE')]
  
  cases <- merge(cases,mind,by=c('NAME','ST','SELECT','ST_CLADE'))
  cases$period <- "pre-startd"
  cases$period[cases$MIND>start_d] <- "post-startd"
  
  # distinguish the pre-existing subgraphs with 0 index, 0 new cases
  cases[m==0 & jcases==0 & period=='pre-startd',m:= -1]
  # for subgraphs starting prior to cut-date with no observed index case, assume one unobserved index case
  if(index_in_pre==1){
    cases[m==0 & period=='pre-startd',jcases:= jcases-1]
  }
  cases[m==0 & cases$period=='pre-startd',m:= 1]
  
  # final subgraph sizes from individuals not in ART by start date
  cases$SIZE <- cases$m + cases$jcases
  cases$icases <- NULL
  setnames(cases,'m','icases')
  
  return(list(cases,dind))
  }

generate_stan_data = function(trsm, tmp, dind, start_d, end_d, max_icases, max_jcases,upper.bound.multiplier){
  inf <- list()
  seq <- list()

  # max icases 
  max_icases = max(tmp$icases[which(tmp$icases <= max_icases)]) # trim index cases 
  # max jcases 
  max_jcases =  max(tmp$jcases[which(tmp$jcases <= max_jcases)]) # trim secondary cases 
  cat("\nmaximum index cases ", max_icases, ' and maximum secondary cases ', max_jcases, '\n')
  
  # exclude those indiviuals on ART by start date (also excluded from subgraph sizes)
  dind$prestartd <- 0
  dind$prestartd[dind$HIV1_POS_D<start_d & dind$RECART_D<start_d & !is.na(dind$RECART_D)] <-	1
  dind <- subset(dind,prestartd==0)
  
  MSM <- subset(dind,CITY=='Amsterdam' & TRANSM=='MSM' & prestartd==0 & HIV1_POS_D<end_d)
  HSX <- subset(dind,CITY=='Amsterdam' & TRANSM=='HSX' & prestartd==0 & HIV1_POS_D<end_d)
  inf[['MSM']] <- length(unique(MSM$ID)) 
  seq[['MSM']] <- length(unique(MSM$ID[MSM$SEQ==T])) 
  inf[['HSX']] <- 	length(unique(HSX$ID)) 
  seq[['HSX']] <- length(unique(HSX$ID[HSX$SEQ==T])) 
  
  # max obs chain size
  N_cs_obs <- max(tmp$SIZE)
  # fill counts of 0
  list <- as.data.table(tidyr::crossing(jcases=seq(0:max(tmp$jcases)),icases=seq(1:max(tmp$icases))))
  list$jcases <- list$jcases - 1
  tmp <- merge(list,tmp,by=c("jcases","icases"),all=TRUE)
  # cs_obs = matrix of index cases by generated cases
  tmp <- dcast.data.table(tmp,icases~jcases,value.var='N')
  tmp[is.na(tmp)] <- 0
  # trim/set max dimensions for observed subgraphs
  tmp <- tmp[1:max_icases,1:(max_jcases+2)]

  stan.data <- list()
  stan.data$index_cases <- sort(unique(tmp$icases))
  tmp <- tmp[,!1]
  stan.data$cs_obs <- tmp
  stan.data$N_cs_obs <- N_cs_obs # max obs chain size = max intial cases + max generated cases (minus first col 0 generated)
  stan.data$N_icases <- length(stan.data$index_cases)
  stan.data$N_jcases <- ncol(tmp)
  stan.data$sampling_n <- inf[[trsm]]				# actual number infected in Amsterdam not on ART by start date (incl index cases)
  stan.data$sampling_k <- seq[[trsm]]		# actual number sequenced in Amsterdam not on ART by start date
  stan.data$N_cs_actual <- ceiling(stan.data$N_cs_obs / (stan.data$sampling_k /stan.data$sampling_n) * upper.bound.multiplier)	#	set upper bound for infinite sum approximation
  
  return(stan.data)
}

map_mwmb_regions= function(dat){
  ## msm
  dat[TRANSM=='MSM', mwmb:="Other"]
  dat[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
  # western countires (non-NL)
  dat[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania") & ORIGIN!='NL', mwmb:="G1"]
  # eastern and central europe
  dat[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="G2"]
  # caribbean and south america
  dat[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]
  
  ## hsx
  dat[TRANSM=='HSX', mwmb:="Other"]
  dat[TRANSM=='HSX' & ORIGIN %in% c("NL"), mwmb:="NL"]
  # sub-saharan africa
  dat[TRANSM=='HSX' & WRLD_born %in% c("Africa"), mwmb:="G4"]
  # caribbean and south america
  dat[TRANSM=='HSX' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]
  
  return(dat)
}

calculate_sampling_fraction_mwmb = function(infile.geo,dind,tmp,infdate,start_d,end_d,trsm,p_undiagnosed){
  
  # exclude those individuals virally suppressed by start date (also excluded from subgraph sizes)
  dind$prestartd <- 0
  if(infdate==1){
    dind$prestartd[dind$INF_D<start_d & dind$supp==1 & !is.na(dind$RECART_D)] <-	1
  }else{
    dind$prestartd[dind$HIV1_POS_D<start_d & dind$supp==1 & !is.na(dind$RECART_D)] <-	1
  }
  dat <- subset(dind,prestartd==0)
  
  # sampling fraction for new cases since start_d
  if(infdate==1){
    dat <- subset(dat,CITY=='Amsterdam' & TRANSM==trsm & prestartd==0 & INF_D>=start_d & INF_D<=end_d)
  }else{
    dat <- subset(dat,CITY=='Amsterdam' & TRANSM==trsm & prestartd==0 & HIV1_POS_D>=start_d & HIV1_POS_D<=end_d)
  }
  
  # # sequenced = number of new cases in subgraph data
  tmp[,f:=jcases * N]
  seq <- sum(tmp$f)
  # 
  #inf <- length(unique(dat$ID)) 
  
  pr <- dat[, list(diag=length(ID)),by=c('TRANSM','mwmb')]
  pr <- merge(pr,p_undiagnosed,by=c('mwmb'),all=T)
  # update estimated number infected using prop undiagnosed
  pr[, inf:= round(diag/(1-du),digits=0) ]
  inf <- sum(pr$inf)

  return(list(seq,inf))
}

generate_stan_data_sbt = function(trsm, tmp, dind, seq, inf, p_undiagnosed, start_d, end_d, max_icases, max_jcases, upper.bound.multiplier, index_flag, infdate){
  
  # max icases 
  max_icases = max(tmp$icases[which(tmp$icases <= max_icases & tmp$jcases <= max_jcases)]) # trim index cases  
  # max jcases 
  max_jcases =  max(tmp$jcases[which(tmp$jcases <= max_jcases & tmp$icases <= max_icases)]) # trim secondary cases 
  cat("\nmaximum index cases ", max_icases, ' and maximum secondary cases ', max_jcases, '\n')
  
  # subtypes
  subtypes = unique(tmp$ST)
  
  stan.data <- list()
  stan.data$sampling_n <- inf				# actual number infected in Amsterdam
  stan.data$sampling_k <- seq		# actual number sequenced infected in Amsterdam	
  #stan.data$inf_m <- inf_m
  
  # max obs chain size
  stan.data$N_cs_obs <- max_icases+max_jcases
  # number of subtype
  stan.data$N_sbts = length(unique(tmp$ST))
  # max index cases
  stan.data$N_icases = max_icases
  
  # flag indicating whether to assume index case is part of observed subgraph for chains starting since cutoff date
  stan.data$index_flag = index_flag
  
  # max number of generated cases in data by subtypes
  stan.data$N_cs_actual = vector(mode = 'integer', length = stan.data$N_sbts)
  # secondary cases by subtypes
  stan.data$N_jcases = vector(mode = 'integer', length = stan.data$N_sbts)
  # index cases by subtypes
  stan.data$index_cases = array(dim = c(max_icases, stan.data$N_sbts), 0)
  # matrix of index cases by generated cases
  stan.data$cs_obs = array(dim = c(max_icases+1, max_jcases+1, stan.data$N_sbts), 0)

  for(s in 1:stan.data$N_sbts){
    #s = 6
    
    tmp1 = subset(tmp, ST == subtypes[s] & icases <= max_icases & jcases <= max_jcases)
    
    # index cases
    tmp2 = subset(tmp1,icases!=0)
    index_cases = sapply(1:max_icases, function(i) i %in% sort(unique(tmp2$icases)))
    index_cases[index_cases == 0] = -1
    index_cases[index_cases == 1] = sort(unique(tmp2$icases))
    stan.data$index_cases[,s] = index_cases
    
    # number of generated cases in data
    stan.data$N_jcases[s] <- max(tmp1$jcases) + 1
    stan.data$N_cs_actual[s] <- ceiling(max(tmp$jcases) / (stan.data$sampling_k /stan.data$sampling_n) * upper.bound.multiplier)
    
    # cs_obs = matrix of index cases by generated cases
    list <- as.data.table(tidyr::crossing(jcases=seq(0:max_jcases),icases=seq(0:max_icases)))
    list$icases <- list$icases - 1
    list$jcases <- list$jcases - 1
    tmp1 <- merge(list,tmp1,by=c("jcases","icases"),all=TRUE)
    tmp1 <- dcast.data.table(tmp1,icases~jcases,value.var='N')
    tmp1[is.na(tmp1)] <- 0
    tmp1 = tmp1[, -'icases']
    # row 1 (m=0) are the newly emergent subgraphs
    stan.data$cs_obs[,,s] = as.matrix(tmp1)
  }
  
  return(stan.data)
}

stan_data_add_origins = function(trsm,stan.data,outdir){
  
  dsubgraphtaxa <- readRDS(file.path(outdir,'subgraphs_withmetadata.RDS'))
  
  dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
  dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
  do <- dsubgraphtaxa
  do <- do[,list(N=length(TAXA)),by=c('ORIGINHOST','TRANSM')]
  setnames(do,c('ORIGINHOST','TRANSM'),c('PARENT_STATE','TRM_GROUP'))
  
  do <- subset(do, PARENT_STATE!='NA')
  da <- do[, list(N=sum(N)), by=c('PARENT_STATE','TRM_GROUP')]
  da <- merge(da, da[, list(TOTAL=sum(N)), by='TRM_GROUP'], by='TRM_GROUP')
  da[, pct:=N/TOTAL]
  
  dt <- list()
  trm <- c('HSX','MSM')
  for(i in trm)
  {
    dt[[i]] <- subset(da, TRM_GROUP==i)
  }
  stan.data$N_origins <- length(unique(dt[[trsm]]$PARENT_STATE))
  stan.data$N_subgraphs <- unique(dt[[trsm]]$TOTAL)
  stan.data$obs_origins <- dt[[trsm]]$N
  stan.data$pr_origins <- dt[[trsm]]$pct
  
  return(stan.data)
}

stan_data_add_origins_sbt = function(trsm,stan.data,outdir,indir.results,start_d,end_d,infdate){
  `%notin%` <- Negate(`%in%`)
  
  dsubgraphtaxa <- readRDS(file.path(outdir,'subgraphs_withmetadata.RDS'))
  tmp <- readRDS(file.path(indir.results,paste0('subgraph_sizes_',trsm,'.RDS')))
  
  cat(paste('\n Subset data to individuals diagnosed between ',start_d,'-',end_d,' \n'))
  # look at origins in time window of analysis
  # keep only at new cases since start date
  if(infdate==1){
    dsubgraphtaxa[INF_D>=start_d & INF_D<end_d, keep:=1]
  }else{
    dsubgraphtaxa[HIV1_POS_D>=start_d & HIV1_POS_D<end_d, keep:=1]
  }

  dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',trsm) & keep==1)

  # ensure subtype indices match rest of stan.data
  subtypes = unique(tmp$ST)
  
  dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
  dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
  do <- dsubgraphtaxa
  # count unique subgraphs with origin in each region (by clade too)
  do <- do[,list(N=length(unique(FULL_NAME))),by=c('REP','ORIGIN','TRANSM','ST','ST_CLADE')]
  # aggregate over clades for B
  do <- do[,list(N=sum(N)),by=c('REP','ORIGIN','TRANSM','ST')]
  setnames(do,c('TRANSM'),c('TRM_GROUP'))
  
  do <- subset(do, ORIGIN!='Unknown')
  
  da <- do[, list(N=sum(N)), by=c('ORIGIN','TRM_GROUP','ST')]
  da <- merge(da, da[, list(TOTAL=sum(N)), by=c('TRM_GROUP','ST')], by=c('TRM_GROUP','ST'))
  da[, pct:=N/TOTAL]
  
  tmp <- dcast.data.table(da,ORIGIN~ST,value.var='N')
  tmp[is.na(tmp)] <- 0
  # reorder columns so subtypes match rest of stan data
  colorder <- c("ORIGIN",subtypes)
  # add 0 cols for subtypes with no new cases
  mis <- length(colorder[colorder %notin% colnames(tmp)])
  tmp2 <- cbind(tmp,matrix(0, nrow = nrow(tmp), ncol = mis))
  colnames(tmp2) <- c(colnames(tmp),colorder[colorder %notin% colnames(tmp)])
  tmp <- tmp2
  setcolorder(tmp,colorder)
  # only keep subtypes in stan data (with >0 subgraphs)
  tmp <- tmp[, colnames(tmp) %in% subtypes, with=FALSE]
  
  stan.data$N_subgraphs <- unique(da$TOTAL[da$ST %in% subtypes])
  stan.data$N_origins <- length(unique(da$ORIGIN[da$ST %in% subtypes]))
  stan.data$obs_origins <- t(as.matrix(tmp))
  
  tmp <- dcast.data.table(da,ORIGIN~ST,value.var='pct')
  tmp[is.na(tmp)] <- 0
  # reorder columns so subtypes match rest of stan data
  colorder <- c("ORIGIN",subtypes)
  # add 0 cols for subtypes with no new cases
  mis <- length(colorder[colorder %notin% colnames(tmp)])
  tmp2 <- cbind(tmp,matrix(1/(nrow(tmp)), nrow = nrow(tmp), ncol = mis))
  colnames(tmp2) <- c(colnames(tmp),colorder[colorder %notin% colnames(tmp)])
  tmp <- tmp2
  setcolorder(tmp,colorder)
  # only keep subtypes in stan data (with >0 subgraphs)
  tmp <- tmp[, colnames(tmp) %in% subtypes, with=FALSE]
  
  stan.data$pr_origins <- t(as.matrix(tmp))
  
  dsubgraphtaxa[, loc_Ams:= as.character(factor(grepl('Ams',ORIGIN), levels=c(TRUE,FALSE), labels=c('Ams','External')))]
  do <- dsubgraphtaxa[,list(N=length(FULL_NAME)),by=c('loc_Ams','ST')]
  tmp <- dcast.data.table(do,loc_Ams~ST,value.var='N')
  tmp[is.na(tmp)] <- 0
 
  colorder <- c("loc_Ams",subtypes)
  # add 0 cols for subtypes with no new cases
  mis <- length(colorder[colorder %notin% colnames(tmp)])
  tmp2 <- cbind(tmp,matrix(0, nrow = nrow(tmp), ncol = mis))
  colnames(tmp2) <- c(colnames(tmp),colorder[colorder %notin% colnames(tmp)])
  tmp <- tmp2
  setcolorder(tmp,colorder)
  
  tmp <- tmp[, colnames(tmp) %in% subtypes, with=FALSE]
  tmp <- tmp + 0.5
  
  stan.data$origin_Ams <- t(as.matrix(tmp))
  
  return(stan.data)
}


stan_data_add_sgs_m_sbt <- function(stan.data,m,outdir){

  stan.data$N_sgs_m_s <- array(dim = c(stan.data$N_sbts), 0)
  for(s in 1:stan.data$N_sbts){
   stan.data$N_sgs_m_s[s] <- sum( stan.data$cs_obs[m,,s])
  }
  stan.data$N_sgs_m <- max(stan.data$N_sgs_m_s)

  return(stan.data)
}

stan_data_add_sgs_m <- function(stan.data,m,outdir){
  
  stan.data$N_sgs_m <- max(rowSums(stan.data$cs_obs))
  
  return(stan.data)
}

stan_data_add_B_nonB = function(trsm,stan.data,outdir,indir.results,start_d,end_d,infdate){
  `%notin%` <- Negate(`%in%`)
  
  dsubgraphtaxa <- readRDS(file.path(outdir,'subgraphs_withmetadata.RDS'))
  tmp <- readRDS(file.path(indir.results,paste0('subgraph_sizes_',trsm,'.RDS')))
  
  cat(paste('\n Subset data to individuals diagnosed between ',start_d,'-',end_d,' \n'))
  # look at origins in time window of analysis
  # keep only at new cases since start date
  if(infdate==1){
    dsubgraphtaxa[INF_D>=start_d & INF_D<end_d, keep:=1]
  }else{
    dsubgraphtaxa[HIV1_POS_D>=start_d & HIV1_POS_D<end_d, keep:=1]
  }
  
  dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',trsm) & keep==1)
  
  # ensure subtype indices match rest of stan.data
  subtypes = unique(tmp$ST)
  
  do <- dsubgraphtaxa
  do <- do[,list(N=length(TAXA)),by=c('TRANSM','ST')]
  setnames(do,c('TRANSM'),c('TRM_GROUP'))
  
  da <- do[, list(N=sum(N)), by=c('TRM_GROUP','ST')]
  da <- merge(da, da[, list(TOTAL=sum(N)), by=c('TRM_GROUP')], by=c('TRM_GROUP'))
  da[, pct:=N/TOTAL]
  
  tmp <- dcast.data.table(da,TRM_GROUP~ST,value.var='pct')
  tmp[is.na(tmp)] <- 0
  # reorder columns so subtypes match rest of stan data
  colorder <- c("TRM_GROUP",subtypes)
  setcolorder(tmp,colorder)

  stan.data$st_p <- t(as.matrix(tmp[,-1]))
  return(stan.data)
}

stan_data_add_bplace_sbt = function(trsm,stan.data,outdir,indir.results,start_d,end_d,infdate){
  `%notin%` <- Negate(`%in%`)
  
  dsubgraphtaxa <- readRDS(file.path(outdir,'subgraphs_withmetadata.RDS'))
  tmp <- readRDS(file.path(indir.results,paste0('subgraph_sizes_',trsm,'.RDS')))
  
  cat(paste('\n Subset data to individuals diagnosed between ',start_d,'-',end_d,' \n'))
  
  # look at origins in time window of analysis
  # keep only at new cases since start date
  if(infdate==1){
    dsubgraphtaxa[INF_D>=start_d & INF_D<end_d, keep:=1]
  }else{
    dsubgraphtaxa[HIV1_POS_D>=start_d & HIV1_POS_D<end_d, keep:=1]
  }
  
  dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',trsm) & keep==1)
  
  geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
  geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
  
  load(file='/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200917/misc/200917_sequence_labels.rda')
  dseq <- merge(dseq,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
  setnames(dseq,'ORIGIN','BIRTH_COUNTRY_ISO')
  dsubgraphtaxa <- merge(dsubgraphtaxa,unique(subset(dseq,select=c('PATIENT','BIRTH_COUNTRY_ISO','WRLD'))),by.x='ID',by.y='PATIENT',all.x=T)
  
  
  # ensure subtype indices match rest of stan.data
  subtypes = unique(tmp$ST)
  
  dsubgraphtaxa[TRANSM=='MSM', mwmb:="Other"]
  dsubgraphtaxa[TRANSM=='MSM' & BIRTH_COUNTRY %in% c("Netherlands"), mwmb:="NL"]
  # western countires (non-NL)
  dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("WEurope","NorthAm","Oceania") & BIRTH_COUNTRY!='Netherlands', mwmb:="G1"]
  # eastern and central europe
  dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("EEurope", "CEurope"), mwmb:="G2"]
  # caribbean and south america
  dsubgraphtaxa[TRANSM=='MSM' & WRLD %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]
  
  ## hsx
  dsubgraphtaxa[TRANSM=='HSX', mwmb:="Other"]
  dsubgraphtaxa[TRANSM=='HSX' & BIRTH_COUNTRY %in% c("Netherlands"), mwmb:="NL"]
  # sub-saharan africa
  dsubgraphtaxa[TRANSM=='HSX' & WRLD %in% c("Africa"), mwmb:="G4"]
  # caribbean and south america
  dsubgraphtaxa[TRANSM=='HSX' & WRLD %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]
  
  do <- dsubgraphtaxa
  do <- do[,list(N=length(FULL_NAME)),by=c('mwmb','TRANSM','ST','ST_CLADE')]
  setnames(do,c('TRANSM'),c('TRM_GROUP'))

  da <- do[, list(N=sum(N)), by=c('mwmb','TRM_GROUP','ST')]
  da <- da[, list(mwmb=mwmb,pct=N/sum(N)),by=c('ST')]

  tmp <- dcast.data.table(da,mwmb~ST,value.var='pct')
  tmp[is.na(tmp)] <- 0
  # reorder columns so subtypes match rest of stan data
  colorder <- c("mwmb",subtypes)
  # add 0 cols for subtypes with no new cases
  mis <- length(colorder[colorder %notin% colnames(tmp)])
  tmp2 <- cbind(tmp,matrix(0, nrow = nrow(tmp), ncol = mis))
  colnames(tmp2) <- c(colnames(tmp),colorder[colorder %notin% colnames(tmp)])
  tmp <- tmp2
  setcolorder(tmp,colorder)
  # only keep subtypes in stan data (with >0 subgraphs)
  tmp <- tmp[, colnames(tmp) %in% subtypes, with=FALSE]
  
  stan.data$N_bplace <- length(unique(da$mwmb))
  stan.data$pr_bplace <- t(as.matrix(tmp))
  
  return(stan.data)
}



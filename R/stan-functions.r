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
	regex.tip.label <- '^([A-Za-z]+)_+([0-9]+)_([0-9]+)$'
	dsubgraphtaxa[, ID:= as.numeric(gsub(regex.tip.label,'\\2',TAXA))]
	dsubgraphtaxa[, SEQ_YEAR:= as.numeric(gsub(regex.tip.label,'\\3',TAXA))]
	dsubgraphtaxa[, FULL_NAME:= paste0(ST,'_',REP,'_',SELECT,'_',NAME)]
	
	dsubgraphtaxa <- unique(dsubgraphtaxa)
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
  dinf <- subset(dinf,select=c('id','estsctodiagMedian')) #rem pos d
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
  # estimate infection time
  dat[, INF_D:= HIV1_POS_D - estsctodiagMedian]
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

subgraph_sizes = function(dsubgraphtaxa,rho_m,dc,start_d,end_d,trsm,index_in_pre){
	cat(paste('\n Subset data to individuals diagnosed up until ',end_d,' \n'))
	# Remove individuals diagnosed after end date and with unknown HIV positive date
	dsubgraphtaxa <- subset(dsubgraphtaxa, inf_after_endd==0 & diagd_unknown==0)

	cat(paste('\n Count index cases and new cases \n'))
	# summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
	cases <- dsubgraphtaxa[, list(icases=length(ID[inf_after_startd==0 & supp==0]),jcases=length(ID[inf_after_startd==1])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]
	cases <- subset(cases,SELECT==paste0('Ams',trsm))
	
	cat('\n Calculate expected index cases from observed using sampling fraction at start of analysis \n')

	inf.startd <- length(unique(rho_m$ID)) 
	seq.startd <- length(unique(rho_m$ID[rho_m$SEQ==T])) 
	sampling.prob <- seq.startd/inf.startd
	
	cases$m <- round(cases$icases/sampling.prob,0)
	
	cases <- merge(cases,dc,by=c('SELECT','NAME','ST','ST_CLADE'),all.x=T)

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
	
	return(cases)
	
}

map_mwmb_regions= function(dat){
  ## msm
  dat[TRANSM=='MSM', mwmb:="Other"]
  dat[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
  # western countries (non-NL)
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

calculate_sampling_fraction_mwmb = function(tmp,pr,infdate,start_d,end_d,p_undiagnosed){
  
  # sequenced = number of new cases since start date in subgraph data
  tmp[,f:=jcases * N]
  seq <- sum(tmp$f)
  
  # infected = number of new diagnosed cases infected since start date
  pr <- merge(pr,p_undiagnosed,by=c('mwmb'),all=T)
  # update estimated number infected using prop undiagnosed
  pr[, inf:= round(diag/(1-du),digits=0) ]
  inf <- sum(pr$inf)

  return(list(seq,inf))
}

generate_stan_data = function(tmp, seq, inf, p_undiagnosed, start_d, end_d, max_icases, max_jcases, upper.bound.multiplier, index_flag, infdate){
  
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

stan_data_add_origins_sbt = function(trsm,stan.data,outdir,indir.results,start_d,end_d,infdate){
  `%notin%` <- Negate(`%in%`)
  
  dsubgraphtaxa <- readRDS(file.path(outdir,'subgraphs_withmetadata.RDS'))
  tmp <- readRDS(file.path(indir.results,paste0('subgraph_sizes_',trsm,'.RDS')))
  
  # keep subgraphs with new cases since start date
  dsubgraphtaxa[inf_after_startd==1 & inf_after_endd==0, keep:=1]
  
  dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',trsm) & keep==1)

  # ensure subtype indices match rest of stan.data
  subtypes = unique(tmp$ST)
  
  dsubgraphtaxa[, ORIGIN := ORIGINHOST]
  dsubgraphtaxa[is.na(ORIGIN), ORIGIN:= 'Unknown']
  do <- dsubgraphtaxa
  do[, TRANSM:= gsub('Ams','',SELECT)]
  
  # count unique subgraphs with origin in each region (by clade too)
  do <- do[,list(N=length(unique(FULL_NAME))),by=c('REP','ORIGIN','TRANSM','ST','ST_CLADE')]
  # aggregate over clades for B
  do <- do[,list(N=sum(N)),by=c('REP','ORIGIN','TRANSM','ST')]
  setnames(do,c('TRANSM'),c('TRM_GROUP'))
  # exclude the unknowns
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
  
  return(stan.data)
}

stan_data_add_bplace_sbt = function(trsm,stan.data,outdir,indir.results,infile.bplaces,start_d,end_d,infdate){
  `%notin%` <- Negate(`%in%`)
  
  tmp <- readRDS(file.path(indir.results,paste0('subgraph_sizes_',trsm,'.RDS')))
  da <- data.table(read.csv(infile.bplaces))
  da <- subset(da,TRM_GROUP==trsm)
  # ensure subtype indices match rest of stan.data
  subtypes = unique(tmp$ST)
  
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



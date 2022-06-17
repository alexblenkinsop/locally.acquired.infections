## ---- rmd.chunk.roadmap.200203.alignment.make.bootstraps ----
roadmap.200303.alignment.make.bootstraps <- function(analysis)
{
  require(ape)
  require(data.table)
  
	analysis <- 'analysis_200917'
  
  home <- '/Users/alexb/Box Sync/Roadmap'
  #home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  indir <- file.path(home,analysis,'alignments_subtype')
  outdir <- file.path(home,analysis,'alignments_100bs')
  bsn <- 100
  seed <- 42L
  #
  df <- data.table(F=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',F))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,F]
    cat('\nprocess ',infile.fasta)
    sq                  <- read.dna(infile.fasta, format='fa')
    set.seed(seed)
    for(b in 0:bsn)
    {
      outfile <- gsub('\\.fasta',paste0('_',sprintf("%03d",b),'.fasta'),infile.fasta)
      outfile <- file.path(outdir,basename(outfile))
      bs <- 1:ncol(sq)
      if(b>0)
      {
        bs<- sample(1:ncol(sq), ncol(sq), replace=TRUE)
      }
      write.dna(sq[,bs], file=outfile, format='fa', colsep='', nbcol=-1)
    }
  }
}


## ---- rmd.chunk.roadmap.200203.fastree ----
roadmap.200303.fastree<- function(analysis)
{
  require(big.phylo)
  require(data.table)
  #source("/rdsgpfs/general/project/ratmann_roadmap_data_analysis/live/bigphylo.cmd.2.R")
	
	analysis <- 'analysis_200917'
	
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

  indir <- file.path(home,analysis,'alignments_bs')
  outdir <- file.path(home,analysis,'fasttree')
  
  #   get UNIX commands for all trees to be processed
  df <- data.table(FIN=list.files(indir, pattern='fasta$',full.names=TRUE))
  df <- subset(df, grepl('noDRM',FIN))
  df[, ST:= gsub('.*_subtype_([A-Z0-9a-z]+)_.*','\\1',basename(FIN))]
  df[, FOUT:= file.path(outdir,gsub('\\.fasta','_ft.newick',basename(FIN)))]
  cmds <- vector('list', nrow(df))
  for(i in seq_len(nrow(df)))
  {
    infile.fasta <- df[i,FIN]
    outfile <- df[i,FOUT]
    tmp <- cmd.fasttree(infile.fasta, outfile=outfile, pr.args='-nt -gtr -gamma', check.binary=TRUE)
    cmds[[i]] <- tmp
  }
  df[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(df[1,CMD])
  
  #   run on HPC as array job
  #df <- subset(df, ST=='B')
  df[, CASE_ID:= 1:nrow(df)]
  
  #   make PBS header
  #hpc.load    <- "module load R/3.4.0"
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate Renv3"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 72                      # walltime
  hpc.q       <- NA                       # PBS queue
  hpc.mem     <- "6gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":00:00,pcput=", hpc.walltime, ":00:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("trs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.roadmap.200203.sequence.labels ----
roadmap.200303.sequence.labels<- function(analysis)
{
  require(data.table)
  require(ape)
  require(tidyverse)
	#require(big.phylo)
  
	analysis <- 'analysis_200917'
  source('/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo functions.R')
  home <- '/Users/alexb/Box Sync/Roadmap'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

    infile.indinfo <- file.path(home,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
    infile.nlinfo <- file.path(home,'Data','data_200316_Netherlands','SHM_1902_ROADMAP_200316_tblBAS.csv')
  	infile.subtypes <- file.path(home,'Data','data_200821','ROADMAP_200821_All_Taxa_withsubtype.rda')
    infile.seqinfo <- file.path(home,'Data','data_200821','SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
    infile.georeg <- file.path(home,'misc','UN_geographic_locations.csv')
    infile.countrycodes <- file.path(home,'misc','ISO_country_codes.csv')
    infiles.lanl <- file.path(home,'alignments','data_200305_LANL_alignment.fasta')
    outfile.noorigin <- file.path(home,analysis,'misc','No_birth_country.csv')
    outfile.demodata <- file.path(home,analysis,'misc','Baseline_data_Ams_NL.csv')
    outfile.base <- file.path(home,analysis,'misc','200917_')
    
    #
  # collect country codes - from UN (geographic regions) and ISO.org (for country codes to match)
  #
  headers = read.csv(infile.georeg, header = F, nrows = 1, as.is = T)
  dgeo <- read.csv(infile.georeg,header=F,skip=26)
  colnames(dgeo)= headers
  names(dgeo)[5] <- 'Alpha.3.code'
  dgeo <- dgeo[,c('Location','Alpha.3.code','GeoRegName')]
    iso <- read.csv(infile.countrycodes,header=T)
    iso <- merge(iso,dgeo,by='Alpha.3.code',all=T)
    iso <- iso[!is.na(iso[,'Alpha.2.code']),]
    
    # Eastern Europe from United Nations regional groups
    # Central Asia from UNAids
    iso  <- as_tibble(iso) %>%
      mutate(Alpha.2.code:=as.character(Alpha.2.code),
             WRLD:= as.character(GeoRegName),
             CNTRY:= as.character(Location)) %>%
      select(-English.short.name,-Numeric,-Location,-GeoRegName) %>%
      mutate(WRLD:= case_when(Alpha.2.code=='SX'~'FormerCurrDutchColonies',
                              Alpha.2.code=='BQ'~'FormerCurrDutchColonies',
                              Alpha.2.code=='AW'~'FormerCurrDutchColonies',
                              Alpha.2.code=='CW'~'FormerCurrDutchColonies',
                              CNTRY=='Suriname'~'FormerCurrDutchColonies',
                              CNTRY=='Albania'~'EEuropeCentralAsia',
                              CNTRY=='Armenia'~'EEuropeCentralAsia',
                              CNTRY=='Azerbijan'~'EEuropeCentralAsia',
                              CNTRY=='Belarus'~'EEuropeCentralAsia',
                              CNTRY=='Bosnia and Herzegovina'~'EEuropeCentralAsia',
                              CNTRY=='Bulgaria'~'EEuropeCentralAsia',
                              CNTRY=='Croatia'~'EEuropeCentralAsia',
                              CNTRY=='Czechia'~'EEuropeCentralAsia',
                              CNTRY=='Estonia'~'EEuropeCentralAsia',
                              CNTRY=='Georgia'~'EEuropeCentralAsia',
                              CNTRY=='Hungary'~'EEuropeCentralAsia',
                              CNTRY=='Kazakhstan'~'EEuropeCentralAsia',
                              CNTRY=='Kyrgyzstan'~'EEuropeCentralAsia',
                              CNTRY=='Latvia'~'EEuropeCentralAsia',
                              CNTRY=='Lithuania'~'EEuropeCentralAsia',
                              CNTRY=='North Macedonia'~'EEuropeCentralAsia',
                              CNTRY=='Republic of Moldova'~'EEuropeCentralAsia',
                              CNTRY=='Montenegro'~'EEuropeCentralAsia',
                              #CNTRY=='Morocco'~'Morocco',
                              CNTRY=='Poland'~'EEuropeCentralAsia',
                              CNTRY=='Romania'~'EEuropeCentralAsia',
                              CNTRY=='Russian Federation'~'EEuropeCentralAsia',
                              CNTRY=='Serbia'~'EEuropeCentralAsia',
                              CNTRY=='Slovakia'~'EEuropeCentralAsia',
                              CNTRY=='Slovenia'~'EEuropeCentralAsia',
                              CNTRY=='Tajikistan'~'EEuropeCentralAsia',
                              #CNTRY=='Turkey'~'Turkey',
                              #CNTRY=='Thailand'~'Thailand',
                              CNTRY=='Turkmenistan'~'EEuropeCentralAsia',
                              CNTRY=='Ukraine'~'EEuropeCentralAsia',
                              CNTRY=='Uzbekistan'~'EEuropeCentralAsia',
                              CNTRY=='Algeria'~'MENA',
                              CNTRY=='Bahrain'~'MENA',
                              CNTRY=='Djibouti'~'MENA',
                              CNTRY=='Egypt'~'MENA',
                              CNTRY=='Iraq'~'MENA',
                              CNTRY=='Iran (Islamic Republic of)'~'MENA',
                              CNTRY=='Jordan'~'MENA',
                              CNTRY=='Kuwait'~'MENA',
                              CNTRY=='Lebanon'~'MENA',
                              CNTRY=='Libya'~'MENA',
                              CNTRY=='Morocco'~'MENA',
                              CNTRY=='Oman'~'MENA',
                              CNTRY=='Qatar'~'MENA',
                              CNTRY=='Saudi Arabia'~'MENA',
                              CNTRY=='Somalia'~'MENA',
                              CNTRY=='Sudan'~'MENA',
                              CNTRY=='Syrian Arab Republic'~'MENA',
                              CNTRY=='Tunisia'~'MENA',
                              CNTRY=='Turkey'~'MENA',
                              CNTRY=='United Arab Emirates'~'MENA',
                              CNTRY=='Yemen'~'MENA',
                              WRLD=='Latin America and the Caribbean'~'LaAmCar',
                              WRLD=='Northern America'~'NorthAm',
                              WRLD=='Europe'~'WEurope',
                              Alpha.3.code=='UMI'~'NorthAm',
                              Alpha.3.code=='SJM'~'WEurope',
                              Alpha.3.code=='SGS'~'LaAmCar',
                              Alpha.3.code=='PCN'~'Oceania',
                              Alpha.3.code=='NFK'~'Oceania',
                              Alpha.3.code=='JEY'~'WEurope',
                              Alpha.3.code=='IOT'~'Asia',
                              Alpha.3.code=='HMD'~'Oceania',
                              Alpha.3.code=='GGY'~'WEurope',
                              Alpha.3.code=='CXR'~'Asia',
                              Alpha.3.code=='CCK'~'Asia',
                              Alpha.3.code=='ATF'~'Africa',
                              Alpha.3.code=='ALA'~'WEurope',
                              TRUE ~ WRLD),
             CNTRY:= case_when(CNTRY=='United States of America'~'USA',
                               CNTRY=='United Kingdom'~'UK',
                               Alpha.3.code=='UMI'~'United States Minor Outlying Islands (the)',
                               Alpha.3.code=='SJM'~'Svalbard and Jan Mayen',
                               Alpha.3.code=='SGS'~'South Georgia and the South Sandwich Islands',
                               Alpha.3.code=='PCN'~'Pitcairn',
                               Alpha.3.code=='NFK'~'Norfolk Island',
                               Alpha.3.code=='JEY'~'Jersey',
                               Alpha.3.code=='IOT'~'British Indian Ocean Territory (the)',
                               Alpha.3.code=='HMD'~'Heard Island and McDonald Islands',
                               Alpha.3.code=='GGY'~'Guernsey',
                               Alpha.3.code=='CXR'~'Christmas Island',
                               Alpha.3.code=='CCK'~'Cocos (Keeling) Islands (the)',
                               Alpha.3.code=='ATF'~'French Southern Territories (the)',
                               Alpha.3.code=='ALA'~'Åland Islands',
                               TRUE~CNTRY))
    
  #
  # read Amsterdam individual data and merge with Dutch data
  #
  dind <- read.csv(infile.indinfo,header=T)
  dnl <- read.csv(infile.nlinfo,header=T)
  dind[,'CITY'] <- 'Amsterdam'
  dnl[,'CITY'] <- 'Non-Amsterdam'
  dind <- merge(dind,dnl,all=T)
  
# Format dates to numeric and add ranges for ones which may be inaccurate
  dind[,'BIRTH_D'] <- as.Date(dind[,'BIRTH_D'])
  dind[,'MIG_D'] <- as.Date(dind[,'MIG_D'],format="%d/%m/%Y")
  dind[,'MIG_D_pq'] <- as.Date(dind[,'MIG_D_pq'],format="%d/%m/%Y")
  dind[,'MIG_D_aq'] <- as.Date(dind[,'MIG_D_aq'],format="%d/%m/%Y")
  dind$MIG_D_lower <- hivc.db.Date.lower(dind$MIG_D)
  dind$MIG_D_upper <- hivc.db.Date.upper(dind$MIG_D)
  dind$MIG_D_lower[!is.na(as.POSIXlt(dind$MIG_D_pq))] <- dind$MIG_D_pq[!is.na(as.POSIXlt(dind$MIG_D_pq))]
  dind$MIG_D_upper[!is.na(as.POSIXlt(dind$MIG_D_aq))] <- dind$MIG_D_aq[!is.na(as.POSIXlt(dind$MIG_D_aq))]
  dind$MIG_D_lower <- as.Date(dind$MIG_D_lower)
  dind$MIG_D_upper <- as.Date(dind$MIG_D_upper)
  
  dind[,'HIV1_NEG_D'] <- as.Date(dind[,'HIV1_NEG_D'],format="%d/%m/%Y")
  dind[,'HIV1_NEG_D_pq'] <- as.Date(dind[,'HIV1_NEG_D_pq'],format="%d/%m/%Y")
  dind[,'HIV1_NEG_D_aq'] <- as.Date(dind[,'HIV1_NEG_D_aq'],format="%d/%m/%Y")
  dind$HIV1_NEG_D_lower <- hivc.db.Date.lower(dind$HIV1_NEG_D)
  dind$HIV1_NEG_D_upper <- hivc.db.Date.upper(dind$HIV1_NEG_D)
  dind$HIV1_NEG_D_lower[!is.na(as.POSIXlt(dind$HIV1_NEG_D_pq))] <- dind$HIV1_NEG_D_pq[!is.na(as.POSIXlt(dind$HIV1_NEG_D_pq))]
  dind$HIV1_NEG_D_upper[!is.na(as.POSIXlt(dind$HIV1_NEG_D_aq))] <- dind$HIV1_NEG_D_aq[!is.na(as.POSIXlt(dind$HIV1_NEG_D_aq))]
  dind$HIV1_NEG_D_lower <- as.Date(dind$HIV1_NEG_D_lower)
  dind$HIV1_NEG_D_upper <- as.Date(dind$HIV1_NEG_D_upper)
  
  dind[,'HIV1_POS_D'] <- as.Date(dind[,'HIV1_POS_D'],format="%d/%m/%Y")
  dind[,'HIV1_POS_D_pq'] <- as.Date(dind[,'HIV1_POS_D_pq'],format="%d/%m/%Y")
  dind[,'HIV1_POS_D_aq'] <- as.Date(dind[,'HIV1_POS_D_aq'],format="%d/%m/%Y")
  dind$HIV1_POS_D_lower <- hivc.db.Date.lower(dind$HIV1_POS_D)
  dind$HIV1_POS_D_upper <- hivc.db.Date.upper(dind$HIV1_POS_D)
  dind$HIV1_POS_D_lower[!is.na(as.POSIXlt(dind$HIV1_POS_D_pq))] <- dind$HIV1_POS_D_pq[!is.na(as.POSIXlt(dind$HIV1_POS_D_pq))]
  dind$HIV1_POS_D_upper[!is.na(dind$HIV1_POS_D_aq)] <- dind$HIV1_POS_D_aq[!is.na(dind$HIV1_POS_D_aq)]
  dind$HIV1_POS_D_lower <- as.Date(dind$HIV1_POS_D_lower)
  dind$HIV1_POS_D_upper <- as.Date(dind$HIV1_POS_D_upper)
  
  dind[,'T0'] <- as.Date(dind[,'T0'],format="%d/%m/%Y")
  dind$T0_lower <- as.Date(hivc.db.Date.lower(dind$T0))
  dind$T0_upper <- as.Date(hivc.db.Date.upper(dind$T0))
  
  dind[,'RECART_D'] <- as.Date(dind[,'RECART_D'],format="%Y-%m-%d")
  dind$RECART_D_lower <- as.Date(hivc.db.Date.lower(dind$RECART_D))
  dind$RECART_D_upper <- as.Date(hivc.db.Date.upper(dind$RECART_D))
  
  dind[,'CARE_FRS_D'] <- as.Date(dind[,'CARE_FRS_D'],format="%d/%m/%Y")
  dind$CARE_FRS_D_lower <- as.Date(hivc.db.Date.lower(dind$CARE_FRS_D))
  dind$CARE_FRS_D_upper <- as.Date(hivc.db.Date.upper(dind$CARE_FRS_D))
  dind[,'DROP_D'] <- as.Date(dind[,'DROP_D'],format="%Y-%m-%d")
  dind$DROP_D_lower <- as.Date(hivc.db.Date.lower(dind$DROP_D))
  dind$DROP_D_upper <- as.Date(hivc.db.Date.upper(dind$DROP_D))
  
  dind[,'BIRTH_D'] <- hivc.db.Date2numeric(dind[,'BIRTH_D'])
  dind[,'BIRTH_Y'] <- floor(dind[,'BIRTH_D'])
  dind[,'MIG_D'] <- hivc.db.Date2numeric(dind$MIG_D)
  dind[,'MIG_D_lower'] <- hivc.db.Date2numeric(dind$MIG_D_lower)
  dind[,'MIG_D_upper'] <- hivc.db.Date2numeric(dind$MIG_D_upper)
  dind[,'HIV1_NEG_D'] <- hivc.db.Date2numeric(dind[,'HIV1_NEG_D'])
  dind[,'HIV1_NEG_D_lower'] <- hivc.db.Date2numeric(dind$HIV1_NEG_D_lower)
  dind[,'HIV1_NEG_D_upper'] <- hivc.db.Date2numeric(dind$HIV1_NEG_D_upper)
  dind[,'HIV1_POS_D'] <- hivc.db.Date2numeric(dind[,'HIV1_POS_D'])
  dind[,'HIV1_POS_D_lower'] <- hivc.db.Date2numeric(dind$HIV1_POS_D_lower)
  dind[,'HIV1_POS_D_upper'] <- hivc.db.Date2numeric(dind$HIV1_POS_D_upper)
  dind[,'T0'] <- hivc.db.Date2numeric(dind[,'T0'])
  dind[,'T0_lower'] <- hivc.db.Date2numeric(dind[,'T0_lower'])
  dind[,'T0_upper'] <- hivc.db.Date2numeric(dind[,'T0_upper'])
  dind[,'RECART_D'] <- hivc.db.Date2numeric(dind[,'RECART_D'])
  dind[,'RECART_D_lower'] <- hivc.db.Date2numeric(dind[,'RECART_D_lower'])
  dind[,'RECART_D_upper'] <- hivc.db.Date2numeric(dind[,'RECART_D_upper'])
  dind[,'CARE_FRS_D_lower'] <- hivc.db.Date2numeric(dind[,'CARE_FRS_D_lower'])
  dind[,'CARE_FRS_D_upper'] <- hivc.db.Date2numeric(dind[,'CARE_FRS_D_upper'])
  dind[,'DROP_D'] <-  hivc.db.Date2numeric(dind[,'DROP_D'])
  dind[,'DROP_D_lower'] <-  hivc.db.Date2numeric(dind[,'DROP_D_lower'])
  dind[,'DROP_D_upper'] <-  hivc.db.Date2numeric(dind[,'DROP_D_upper'])
  
  # Recode baseline variables
    dind <- as_tibble(dind) %>%
      select(PATIENT,BIRTH_Y,GENDER,MODE,MODE_OTH,ORIGIN,MIG_D,MIG_D_lower,MIG_D_upper,HIV1_NEG_D,HIV1_NEG_D_lower,HIV1_NEG_D_upper,
             HIV1_NEG_D_A,HIV1_POS_D,HIV1_POS_D_lower,HIV1_POS_D_upper,HIV1_POS_D_A,INF_NL,INF_COUNTRY_1,
             REG_FRS_GGD,REG_LST_GGD,CITY,RECART_Y,T0,T0_lower,T0_upper,RECART_D,RECART_D_lower,RECART_D_upper,NAIVE_Y,
             CARE_FRS_D,CARE_FRS_D_lower,CARE_FRS_D_upper,DROP_D,DROP_D_lower,DROP_D_upper,DROP_RS) %>%
      mutate( GENDER:= case_when(GENDER==1~'Male',
                                 GENDER==2~'Female'),
              TRANSM:= case_when(MODE==1~'MSM',
                                MODE==2~'IDU',
                                MODE==4~'Other',
                                MODE==5~'Other',
                                MODE==6~'HSX',
                                MODE==8~'Other',
                                MODE==9~'Other',
                                MODE==90~'Other',
                                MODE==99~'Unknown'),
              DUTCH_BORN:= case_when(ORIGIN=='Unknown'~'Unknown',
                                     ORIGIN==''~'Unknown',
                                     ORIGIN=="NL"~'Dutch',
                                     TRUE~'Other'),
              INF_NL:= case_when(INF_COUNTRY_1=='NL'~'Dutch',
                                 INF_COUNTRY_1=='-1'~'Unknown',
                                 INF_COUNTRY_1==""~'Unknown',
                                 INF_COUNTRY_1!='NL'&INF_COUNTRY_1!='-1'&INF_COUNTRY_1!=""~'Non-NL'),
              GGD_FIRST:= case_when(REG_FRS_GGD==111~'Groningen',
                                    REG_FRS_GGD==706~'Drenthe',
                                    REG_FRS_GGD==1009~'IJsselland',
                                    REG_FRS_GGD==1106~'Twente',
                                    REG_FRS_GGD==1406~'Gelre_IJssel',
                                    REG_FRS_GGD==1906~'Hulpverlening_Gelderland_Midden',
                                    REG_FRS_GGD==2006~'Rivierenland',
                                    REG_FRS_GGD==2106~'Nijmegen',
                                    REG_FRS_GGD==2209~'Flevoland',
                                    REG_FRS_GGD==2406~'Utrecht',
                                    REG_FRS_GGD==2506~'Midden_Nederland',
                                    REG_FRS_GGD==2707~'Hollands_Noorden',
                                    REG_FRS_GGD==3109~'Kennemerland',
                                    REG_FRS_GGD==3406~'Amsterdam',
                                    REG_FRS_GGD==3606~'Gooi_Vechtstreek',
                                    REG_FRS_GGD==3906~'Den_Haag',
                                    REG_FRS_GGD==4106~'Zuid_Holland_West',
                                    REG_FRS_GGD==4506~'Hollands_Midden',
                                    REG_FRS_GGD==4607~'Rotterdam_Rijnmond',
                                    REG_FRS_GGD==4810~'Zuid_Holland_Zuid',
                                    REG_FRS_GGD==5006~'Zeeland',
                                    REG_FRS_GGD==5206~'West_Brabant',
                                    REG_FRS_GGD==5406~'Hart_voor_Brabant',
                                    REG_FRS_GGD==5608~'Brabant_Zuidoost',
                                    REG_FRS_GGD==6011~'Limburg-Noord',
                                    REG_FRS_GGD==6106~'Zuid_Limburg',
                                    REG_FRS_GGD==7206~'Fryslan',
                                    REG_FRS_GGD==7306~'Zaanstreek_Waterland'),
              GGD_LAST:= case_when(REG_LST_GGD==111~'Groningen',
                                    REG_LST_GGD==706~'Drenthe',
                                    REG_LST_GGD==1009~'IJsselland',
                                    REG_LST_GGD==1106~'Twente',
                                    REG_LST_GGD==1406~'Gelre_IJssel',
                                    REG_LST_GGD==1906~'Hulpverlening_Gelderland_Midden',
                                    REG_LST_GGD==2006~'Rivierenland',
                                    REG_LST_GGD==2106~'Nijmegen',
                                    REG_LST_GGD==2209~'Flevoland',
                                    REG_LST_GGD==2406~'Utrecht',
                                    REG_LST_GGD==2506~'Midden_Nederland',
                                    REG_LST_GGD==2707~'Hollands_Noorden',
                                    REG_LST_GGD==3109~'Kennemerland',
                                    REG_LST_GGD==3406~'Amsterdam',
                                    REG_LST_GGD==3606~'Gooi_Vechtstreek',
                                    REG_LST_GGD==3906~'Den_Haag',
                                    REG_LST_GGD==4106~'Zuid_Holland_West',
                                    REG_LST_GGD==4506~'Hollands_Midden',
                                    REG_LST_GGD==4607~'Rotterdam_Rijnmond',
                                    REG_LST_GGD==4810~'Zuid_Holland_Zuid',
                                    REG_LST_GGD==5006~'Zeeland',
                                    REG_LST_GGD==5206~'West_Brabant',
                                    REG_LST_GGD==5406~'Hart_voor_Brabant',
                                    REG_LST_GGD==5608~'Brabant_Zuidoost',
                                    REG_LST_GGD==6011~'Limburg-Noord',
                                    REG_LST_GGD==6106~'Zuid_Limburg',
                                    REG_LST_GGD==7206~'Fryslan',
                                    REG_LST_GGD==7306~'Zaanstreek_Waterland'
                                 )) %>%
      mutate(REGION_FIRST:= case_when(GGD_FIRST=='Amsterdam'~'Amsterdam',
                                      GGD_FIRST=='Rotterdam'~'Rotterdam',
                                      GGD_FIRST=='The Hague'~'The Hague',
                                      GGD_FIRST=='Utrecht'~'Utrecht',
                                      TRUE~'Other'),
            REGION_LAST:= case_when(GGD_LAST=='Amsterdam'~'Amsterdam',
                                      GGD_LAST=='Rotterdam'~'Rotterdam',
                                      GGD_LAST=='The Hague'~'The Hague',
                                      GGD_LAST=='Utrecht'~'Utrecht',
                                      TRUE~'Other'
                                      )) %>%
      mutate( Alpha.2.code:= as.character(ORIGIN),
              ORIGIN:= as.character(ORIGIN)) %>%
      mutate(ORIGIN:= case_when(ORIGIN==""~'Unknown',
                                is.na(ORIGIN)~'Unknown',
                         TRUE ~ ORIGIN))
    
    # Obtain world region for country of origin and recode a few manually not in geog regions data
    dind <- dind %>% left_join(iso,by='Alpha.2.code') %>%
            mutate(BIRTH_CNTRY:= case_when(Alpha.2.code=='AN'~'Netherlands Antilles',
                                           Alpha.2.code=='YU'~'Yugoslavia',
                                           Alpha.2.code=='CS'~'Serbia and Montenegro',
                                           CNTRY=='Unknown'~'Unknown',
                                           is.na(CNTRY)~'Unknown',
                                           TRUE ~ CNTRY),
                   LOC_BIRTH:= case_when(Alpha.2.code=='AN'~'FormerCurrDutchColonies',
                                         Alpha.2.code=='YU'~'EEuropeCentralAsia',
                                         Alpha.2.code=='CS'~'EEuropeCentralAsia',
                                         CNTRY=='Unknown'~'Unknown',
                                         is.na(CNTRY)~'Unknown',
                                         TRUE ~ WRLD))  %>%
            select(-WRLD,-CNTRY,-Alpha.2.code,-Alpha.3.code) %>%
            mutate(Alpha.2.code:= as.character(INF_COUNTRY_1)) %>%
            left_join(iso,by='Alpha.2.code') %>%
            mutate(INF_CNTRY:= case_when(Alpha.2.code=='AN'~'Netherlands Antilles',
                                     Alpha.2.code=='YU'~'Yugoslavia',
                                     Alpha.2.code=='CS'~'Serbia and Montenegro',
                                     CNTRY=='Unknown'~'Unknown',
                                     is.na(CNTRY)~'Unknown',
                                     TRUE ~ CNTRY),
                   LOC_INF:= case_when(Alpha.2.code=='AN'~'FormerCurrDutchColonies',
                                   Alpha.2.code=='YU'~'EEuropeCentralAsia',
                                   Alpha.2.code=='CS'~'EEuropeCentralAsia',
                                   CNTRY=='Unknown'~'Unknown',
                                   is.na(CNTRY)~'Unknown',
                                   TRUE ~ WRLD))

  # Export individuals without a birth country recorded for investigation
  noorigin <- dind[is.na(dind$ORIGIN) | dind$ORIGIN=="" | dind$ORIGIN=="Unknown",]
  write.csv(noorigin,file=outfile.noorigin)
  
  # Save demographic data
  write.csv(dind,file=outfile.demodata)
  
  # Read in subtype file
  st <- load(infile.subtypes)
  st.l <- as_tibble(ds) %>% select(TAXA_L,SUBTYPE_L) %>% distinct()
  st.a <- as_tibble(ds) %>% select(FASTASampleCode,SUBTYPE) %>% distinct()
  colnames(st.a)[1] <- 'SEQ_LABEL'
  
  #
  # read sequence labels for Amsterdam and NL and add in subtype and geographical data
  #
  load(infile.seqinfo)
  dseq <- ds
  dseq <- as_tibble(dseq) %>% inner_join(st.a, by='SEQ_LABEL')
  
  dseq <- dind %>% inner_join(dseq, by='PATIENT')
  dseq <- dseq %>% select(PATIENT, SEQ_LABEL, SEQ_ID, SEQ_D, SEQ_L, SUBTYPE, ORIGIN, GENDER, TRANSM, CITY, CNTRY, LOC_BIRTH) %>%
          rename(SEQ_DATE=SEQ_D, SEQ_LENGTH=SEQ_L, BIRTHCOUNTRY=CNTRY)

  #
  # collect LANL labels
  #
  dl <- read.dna(infiles.lanl,format='fa')
  dl <- tibble(TAXA=unlist(rownames(dl))) %>%
    filter(!grepl('-.-.-.',TAXA))  %>%
    mutate( SUBTYPE:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',1),
            Alpha.2.code:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',2),
            YEAR:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',3),
            GENBANK:= sapply(strsplit(TAXA,'.',fixed=TRUE),'[[',5)) %>%
    distinct()
  dl$TAXA[grepl('K03455',dl$TAXA)] <- 'HXB2'
  
  outfile <- paste0(outfile.base,'sequence_labels.rda')
  save(dseq, dind, dl, iso, file= outfile)
  
  #
  #   make sequence labels for AmsMSM
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'&TRANSM=='MSM'~'AmsMSM',
                           CITY=='Amsterdam'&TRANSM!='MSM'~'AmsnonMSM',
                           CITY=='Amsterdam'&is.na(TRANSM)~'AmsnonMSM',
                           CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
                            TRUE~WRLD)) %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW)  %>%
    rename(GRP=WRLD) %>%
    rbind(tmp)  %>%
    arrange(GRP)
  outfile <- paste0(outfile.base,'sequence_labels_AmsMSM.csv')
  write.csv(tmp, file=outfile)
  #
  #   make sequence labels for AmsHSX
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'&TRANSM=='HSX'~'AmsHSX',
                            CITY=='Amsterdam'&TRANSM!='HSX'~'AmsnonHSX',
                           CITY=='Amsterdam'&is.na(TRANSM)~'AmsnonHSX',
                           CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
                            TRUE~WRLD)) %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW) %>%
    rename(GRP=WRLD) %>%
    rbind(tmp) %>%
    arrange(GRP)
  outfile <- paste0(outfile.base,'sequence_labels_AmsHSX.csv')
  write.csv(tmp, file=outfile)
  #
  #   make sequence labels for Amsterdam overall
  #
  tmp <- dseq %>% rename(TAXA:= SEQ_LABEL) %>%
    mutate(GRP:= case_when(CITY=='Amsterdam'~'Ams',
                            CITY!='Amsterdam'~'NL')) %>%
    mutate(TAXA_NEW:= paste0(GRP,'___',SEQ_ID,'_',PATIENT,'_',CITY,'_',ORIGIN,'_',LOC_BIRTH,'-',GENDER,'_',TRANSM,'_',SEQ_DATE,'_',SEQ_LENGTH)) %>%
    select(SUBTYPE,GRP,TAXA,TAXA_NEW)
  tmp <- dl %>% left_join(iso, by='Alpha.2.code') %>%
    mutate(WRLD:= case_when(is.na(WRLD)~'Unknown',
                            TRUE~WRLD)) %>%
    mutate(TAXA_NEW:= paste0(WRLD,'___',TAXA)) %>%
    select(SUBTYPE,WRLD,TAXA,TAXA_NEW) %>%
    rename(GRP=WRLD) %>%
    rbind(tmp) %>%
    arrange(GRP)
  outfile <- paste0(outfile.base,'sequence_labels_Ams.csv')
  write.csv(tmp, file=outfile)
}


## ---- rmd.chunk.roadmap.200203.phyloscanner.nonB ----
roadmap.200303.phyloscanner.nonB <- function(analysis)
{
  require(dplyr)
  require(data.table)
  require(ape)
  require(adephylo)
  require(phytools)
  require(phangorn)
	
	analysis <- 'analysis_200917'
	
  
  plot.phylogenies <- 1
  max.Ntip <- 5e3
  if(0)
  {
    prog.phyloscanner_analyse_trees <- '/Users/alexb/Documents/software/phyloscanner/phyloscanner/phyloscanner_analyse_trees.R'
    home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  }
  if(1)
  {
    prog.phyloscanner_analyse_trees <- '/rds/general/project/ratmann_roadmap_data_analysis/live/R/phyloscanner_analyse_trees.R'
    home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  }
  
  
  indir.trees <- file.path(home,analysis,'fasttree')
   infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_Ams.csv'),
                                           file.path(home,analysis,'misc','200917_sequence_labels_AmsMSM.csv'),
                                           file.path(home,analysis,'misc','200917_sequence_labels_AmsHSX.csv'))
   )
  outdir <- file.path(home,analysis,'phyloscanner')
  infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
  infiles[, DUMMY:= 1L]
  infile.labels[, DUMMY:= 1L]
  infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
  infiles[, DUMMY:= NULL]
  infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
  infiles[, SELECT:= gsub('^.*_labels_([A-Z]+)\\.csv$','\\1', basename(FLABEL))]
  infiles[, DUMMY:= paste0(gsub('\\.newick$','_rerooted_',basename(infiles$FIN)), infiles$SELECT)]
  tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
  tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
  infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
  infiles <- subset(infiles, is.na(FOUT) & ST!='B', c(ST, FIN, FLABEL))
  cmds <- vector('list',nrow(infiles))
  for(i in seq_len(nrow(infiles)))
  {
    #i<- 1
    cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
    infile <- infiles[i,FIN]
    ph <- read.tree(infile)
    ph$node.label <- NULL
    #   update tip labels: add world region to start of label
    local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
    dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
    dl[, X:=NULL]
    stopifnot(!any(ph$tip.label==''))
    dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
    dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
    dp$TAXA_NEW[grepl('K03455',dp$TAXA)] <- 'HXB2'
    dp$SUBTYPE[grepl('K03455',dp$TAXA)] <- 'B'
    dp$WRLD[grepl('K03455',dp$TAXA)] <- 'Europe'
    dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <-  gsub('Outgroup.([0-9_A-Z]+).*','\\1',basename(dp$TAXA)[grepl('Outgroup',dp$TAXA)])
    dp$TAXA_NEW[grepl('Outgroup',dp$TAXA)] <- dp$TAXA[grepl('Outgroup',dp$TAXA)]
    stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
    stopifnot( !any(duplicated(ph$tip.label)) )
    ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
    #   re-root
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    tmp <- tmp[order(-NST),][2,SUBTYPE]
    tmp <- subset(dp, SUBTYPE==tmp,)[,TAXA_NEW]
    root <- getMRCA(ph, tmp)
    ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
    #   drop other subtypes
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    # Drop all subtypes which are not the most common
    tmp <- tmp[order(-NST),][-1,SUBTYPE]
    tmp <- subset(dp, SUBTYPE%in%tmp,)[,TAXA_NEW]
    ph <- drop.tip(ph, tmp)
    stopifnot(is.binary(ph))
    #   write to file
    intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
    write.tree(ph, file=intree.phsc)
    if(plot.phylogenies)
    {
      pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(ph)/10)
      plot(ph, show.node.label=TRUE, cex=0.3)
      dev.off()
    }
    #   make phyloscanner UNIX command
    infile <- intree.phsc
    outputString <- paste0(gsub('\\.newick','_',intree.phsc))
    tip.regex <- "^([A-Za-z]+)___.*$"
    cmd <- paste("CWD=$(pwd)\n",sep='')
    cmd <- paste(cmd,"echo $CWD\n",sep='')
    tmpdir.prefix <- paste('phsc_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
    tmpdir <- paste("$CWD/",tmpdir.prefix,sep='')
    tmp.in <- file.path(tmpdir, basename(infile))
    cmd <- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
    cmd <- paste(cmd,'cp "',infile,'" ',tmp.in,'\n', sep='')
    cmd <- paste(cmd,'cd ', tmpdir,'\n', sep='')
    cmd <- paste0(cmd,'Rscript ',prog.phyloscanner_analyse_trees,' ',basename(infile),' ',basename(outputString))
    #   don t use -m multifurcation threshold to ensure that output tree is still binary
    #cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda -m 1e-5\n')
    cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
    cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')
    cmd <- paste(cmd, "cd $CWD\n",sep='')
    cmd <- paste(cmd, "rm ", "-r ", tmpdir,'\n',sep='')
    cmds[[i]] <- cmd
  }
  infiles[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(infiles[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles)
  df[, CASE_ID:= 1:nrow(df)]
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate phylo"
  export.path <- 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH'
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "4gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, export.path, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.roadmap.200203.phyloscanner.B ----
roadmap.200303.phyloscanner.B <- function(analysis)
{
  require(dplyr)
  require(data.table)
  require(ape)
  require(adephylo)
  require(phytools)
  require(phangorn)
 
	analysis <- 'analysis_200917'
	
  plot.phylogenies <- 1
  max.Ntip <- 5e3
  if(0)
  {
    prog.phyloscanner_analyse_trees <- '/Users/alexb/Documents/software/phyloscanner/phyloscanner/phyloscanner_analyse_trees.R'
    home <- '/Users/alexb/Documents/Roadmap/Data/SHM_1902_ROADMAP_191223'
  }
  if(1)
  {
    prog.phyloscanner_analyse_trees <- '/rds/general/project/ratmann_roadmap_data_analysis/live/R/phyloscanner_analyse_trees.R'
    home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  }
  
  
  indir.trees <- file.path(home,analysis,'fasttree')
   infile.labels <- data.table(FLABEL= c(  file.path(home,analysis,'misc','200917_sequence_labels_Ams.csv'),
                                           file.path(home,analysis,'misc','200917_sequence_labels_AmsMSM.csv'),
                                           file.path(home,analysis,'misc','200917_sequence_labels_AmsHSX.csv'))
   )
  outdir <- file.path(home,analysis,'phyloscanner')
  infiles <- data.table(FIN=list.files(indir.trees, pattern='\\.newick$',full.names=TRUE))
  infiles[, DUMMY:= 1L]
  infile.labels[, DUMMY:= 1L]
  infiles <- merge(infiles, infile.labels, by='DUMMY', allow.cartesian=TRUE)
  infiles[, DUMMY:= NULL]
  infiles[, ST:= gsub('.*_subtype_([A-Za-z0-9]+)_.*','\\1',basename(FIN))]
  infiles <- subset(infiles, ST=='B')
  set.seed(42L)
  for(i in seq_len(nrow(infiles)))
  {
    #i<- 1
    cat('\nprocess ',infiles[i,FIN],'\nprocess ',infiles[i,FLABEL])
    infile <- infiles[i,FIN]
    ph <- read.tree(infile)
    ph$node.label <- NULL
    #   update tip labels: add world region to start of label
    local.world.region <- gsub('^.*_labels_([A-Za-z]+)\\.csv$','\\1', basename(infiles[i,FLABEL]))
    dl <- as.data.table(read.csv(infiles[i,FLABEL], stringsAsFactors=FALSE))
    dl[, X:=NULL]
    stopifnot(!any(ph$tip.label==''))
    dp <- data.table(IDX=1:Ntip(ph), TAXA=ph$tip.label)
    dp <- merge(dp,dl,by='TAXA',all.x=TRUE)
    dp$TAXA_NEW[grepl('K03455',dp$TAXA)] <- 'HXB2'
    dp$SUBTYPE[grepl('K03455',dp$TAXA)] <- 'B'
    dp$GRP[grepl('K03455',dp$TAXA)] <- 'Europe'
    dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <-  gsub('Outgroup.([0-9_A-Z]+).*','\\1',basename(dp$TAXA)[grepl('Outgroup',dp$TAXA)])
    dp$TAXA_NEW[grepl('Outgroup',dp$TAXA)] <- dp$TAXA[grepl('Outgroup',dp$TAXA)]
    #dp$SUBTYPE[grepl('Outgroup',dp$TAXA)] <- '28_BF'
    stopifnot( nrow(dp[is.na(TAXA_NEW),])==0 )
    stopifnot( !any(duplicated(ph$tip.label)) )
    ph$tip.label <- dp[order(IDX),][, TAXA_NEW]
    #   re-root
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    tmp <- tmp[order(-NST),][2,SUBTYPE]
    tmp <- subset(dp, SUBTYPE==tmp,)[,TAXA_NEW]
    root <- getMRCA(ph, tmp)
    ph <- reroot(ph, root, ph$edge.length[which(ph$edge[,2]==root)]/2)
    #   drop other subtypes
    tmp <- dp[, list(NST=length(TAXA)), by='SUBTYPE']
    # Drop all subtypes which are not the most common
    tmp <- tmp[order(-NST),][-1,SUBTYPE]
    tmp <- subset(dp, SUBTYPE%in%tmp,)[,TAXA_NEW]
    ph <- drop.tip(ph, tmp)
    #   split very large tree into small trees if necessary
    if(is.finite(max.Ntip))
    {
      split.phs <- list()
      repeat({
        local.tips <- which(grepl(paste0('^',local.world.region),ph$tip.label))
        local.tip.ancestors <- Ancestors(ph, sample(local.tips,1))
        cat('\nNumer of local tips left to divide into smaller trees, n=', length(local.tips))
        subtree.mrca <- Ntip(ph)+1L
        for(k in seq_along(local.tip.ancestors))
        {
          tmp <- length(Descendants(ph, local.tip.ancestors[k], type='tips')[[1]])
          if(tmp>max.Ntip)
          {
            subtree.mrca <- local.tip.ancestors[max(1,k-1)]
            break
          }
        }
        stopifnot( subtree.mrca>Ntip(ph) )
        tmp <- extract.clade(ph, subtree.mrca)
        if(Ntip(ph)==Ntip(tmp))
        {
          split.phs[[length(split.phs)+1L]] <- tmp
          break
        }
        if(Ntip(ph)-Ntip(tmp)>50)
        {
          split.phs[[length(split.phs)+1L]] <- tmp
          ph <- drop.tip(ph, split.phs[[length(split.phs)]][['tip.label']])
        }
        if(Ntip(ph)-Ntip(tmp)<=50)
        {
          cat('\nFound very small final tree, retry sampling a local.tip')
        }
      })
    }
    if(!is.finite(max.Ntip))
      split.phs[[1]] <- ph
    #   write to files
    if(length(split.phs)>1)
    {
      for(k in seq_along(split.phs))
      {
        intree.phsc <- gsub('_(subtype_[A-Za-z0-9]+)_',paste0('_\\1c',k,'_'),gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
        intree.phsc <- file.path(outdir,intree.phsc)
        stopifnot(is.binary(split.phs[[k]]))
        write.tree(split.phs[[k]], file=intree.phsc)
        if(plot.phylogenies)
        {
          pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[k]])/10)
          plot(split.phs[[k]], show.node.label=TRUE, cex=0.3)
          dev.off()
        }
      }
    }
    if(length(split.phs)==1)
    {
      intree.phsc <- file.path(outdir,gsub('\\.newick',paste0('_rerooted_',local.world.region,'.newick'),basename(infile)))
      write.tree(split.phs[[1]], file=intree.phsc)
      if(plot.phylogenies)
      {
        pdf(file=gsub('newick','pdf',intree.phsc), w=20, h=10+Ntip(split.phs[[1]])/10)
        plot(split.phs[[1]], show.node.label=TRUE, cex=0.3)
        dev.off()
      }
    }
  }
  
  infiles <- data.table(FIN=list.files(outdir, pattern='\\.newick$',full.names=TRUE))
  infiles[, ST:= gsub('^.*_subtype_([A-Za-z0-9]+)_.*$','\\1',basename(FIN))]
  infiles[, DUMMY:= gsub('\\.newick$','',basename(FIN))]
  tmp <- data.table(FOUT=list.files(outdir, pattern='workspace.rda$',full.names=TRUE))
  tmp[, DUMMY:= gsub('__workspace.rda$','',basename(FOUT))]
  infiles <- merge(infiles, tmp, by='DUMMY', all.x=TRUE)
  infiles <- subset(infiles, is.na(FOUT), c(ST, FIN))
  infiles <- subset(infiles, grepl('Bc[0-9]+',ST))
  cmds <- vector('list',nrow(infiles))
  for(i in seq_len(nrow(infiles)))
  {
    #   make phyloscanner UNIX command
    infile <- infiles[i,FIN]
    outputString <- paste0(gsub('\\.newick','_',infile))
    tip.regex <- "^([A-Za-z]+)___.*$"
    cmd <- paste("CWD=$(pwd)\n",sep='')
    cmd <- paste(cmd,"echo $CWD\n",sep='')
    tmpdir.prefix <- paste('phsc_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
    tmpdir <- paste("$CWD/",tmpdir.prefix,sep='')
    tmp.in <- file.path(tmpdir, basename(infile))
    cmd <- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
    cmd <- paste(cmd,'cp "',infile,'" ',tmp.in,'\n', sep='')
    cmd <- paste(cmd,'cd ', tmpdir,'\n', sep='')
    cmd <- paste0(cmd,'Rscript ',prog.phyloscanner_analyse_trees,' ',basename(infile),' ',basename(outputString))
    #   don t use -m multifurcation threshold to ensure that output tree is still binary
    #cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda -m 1e-5\n')
    cmd <- paste0(cmd,' s,0 -x "',tip.regex,'" -v 1 -ow -rda\n')
    cmd <- paste0(cmd,'mv ',basename(outputString),'* ','"',dirname(outputString),'"','\n')
    cmd <- paste(cmd, "cd $CWD\n",sep='')
    cmd <- paste(cmd, "rm ", "-r ", tmpdir,'\n',sep='')
    cmds[[i]] <- cmd
  }
  infiles[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(infiles[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles)
  df[, CASE_ID:= 1:nrow(df)]
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate phylo"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "4gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}

## ---- rmd.chunk.roadmap.200203.treedater.run ----
# Ancestral state reconstruction
roadmap.200303.treedater.run <- function(analysis)
{
  require(big.phylo)
  # Need to use edited scripts from big.phylo due to update in treedater command
  source('/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo functions.R')
  source("/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo.cmd.2.R")
  require(data.table)
  require(tidyverse)
  require(ape)
  require(phytools)
  require(treedater)
  require(phyloscannerR)
	
	analysis <- 'analysis_200917'
	
  #home <- '~/Box Sync/Roadmap/Data/SHM_1902_ROADMAP_191223'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  infile.seqinfo <- file.path(home,'Data','data_200821','SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
  indir.phsc  <- file.path(home,analysis,'phyloscanner')
  outdir <- file.path(home,analysis,'phyloscanner_dated')
  infiles.phsc <- data.table(F=list.files(indir.phsc, pattern='workspace.rda$', full.names=TRUE, recursive=TRUE))
  infiles.phsc[, BASENAME:= gsub('workspace.rda$','',basename(F))]
  tmp <- data.table(F_TREE=list.files(outdir, pattern='collapsed_dated.newick$', full.names=TRUE, recursive=TRUE))
  tmp[, BASENAME:= gsub('collapsed_dated.newick$','',basename(F_TREE))]
  infiles.phsc <- merge(infiles.phsc, tmp, by='BASENAME', all.x=TRUE)
  infiles.phsc <- subset(infiles.phsc,is.na(infiles.phsc$F_TREE))
  #infiles.phsc <- subset(infiles.phsc,grepl('subtype_B',infiles.phsc$F))
  alignment.length <-  1000 #1302
  
  #
  # read sampling data
  #
  
  load(infile.seqinfo)
  dseq <- ds

  dates <- as.POSIXlt(dseq$SEQ_D)
  dseq$seqy <- 1900 + dates$year
  dseq$seqm <- dates$mon
  dseq <- dseq %>%
    select(SEQ_ID, PATIENT, SEQ_D, seqy, seqm) %>%
    mutate( SEQ_DATE:= seqy + (seqm-1)/12 + 15/365,
            SEQ_DATE_LOWER:= seqy + (seqm-1)/12,
            SEQ_DATE_UPPER:= seqy + (seqm-1)/12 + 30/365,
            SID_PID:= paste(SEQ_ID,PATIENT,sep="_")
    ) %>%
    select(-seqy, -seqm) %>%
    rename(SID= SEQ_ID, PID= PATIENT)
  dseq <- as.data.table(dseq)
  
  #
  #   for each tree:
  #   make data.table of sequence sampling times
  #   remove taxa without data on sampling times
  #
  cmds <- vector('list',nrow(infiles.phsc))
  for(i in seq_len(nrow(infiles.phsc)))
  {
    #   i<- 4
    cat('\nProcess',i)
    infile <- infiles.phsc[i,F]
    load(infile)
    ph <- phyloscanner.trees[[1]][['tree']]
    stopifnot( !any( ph$tip.label=='' ) )
    #stopifnot( is.binary(ph) )
    #
    #   drop tips without sequence dates,
    #   while conserving the ancestral state reconstructions
    #
    
    #   extract taxa names
    dph <- data.table(  TAXA_LABEL=ph$tip.label,
                        TAXA= gsub('^[^_]+___(.*)','\\1',ph$tip.label),
                        SID_PID = gsub('^[^_]+___([^_]*_[^_]*)_(.*)','\\1',ph$tip.label),
                        TAXA_ID= seq_along(ph$tip.label))
    
    #   add Amsterdam sequence dates
    dph <- merge(dph, dseq, by='SID_PID', all.x=TRUE)
    #   extract GenBank sequence dates from taxa names where possible
    dph[, GENBANK_ID:= gsub('^[^_]+___([^\\.]+).([^\\.]+).([^\\.]+).([^\\.]+).(.*)$','\\5',TAXA_LABEL)]
    dph[, GENBANK_SEQDATE:= gsub('^[^_]+___([^\\.]+).([^\\.]+).([^\\.]+).([^\\.]+).(.*)$','\\3',TAXA_LABEL)]
    set(dph, dph[, which(GENBANK_SEQDATE=='-' | !is.na(SEQ_DATE))],'GENBANK_SEQDATE',NA_character_)
    set(dph, NULL, 'GENBANK_SEQDATE_LOWER', dph[, as.numeric(GENBANK_SEQDATE) ])
    set(dph, NULL, 'GENBANK_SEQDATE_UPPER', dph[, as.numeric(GENBANK_SEQDATE) + 364/365 ])
    set(dph, NULL, 'GENBANK_SEQDATE', dph[, as.numeric(GENBANK_SEQDATE) + 1/2 ])
    tmp <- dph[, which(is.na(SEQ_DATE))]
    set(dph, tmp, 'SEQ_DATE', dph[tmp, GENBANK_SEQDATE])
    set(dph, tmp, 'SEQ_DATE_LOWER', dph[tmp, GENBANK_SEQDATE_LOWER])
    set(dph, tmp, 'SEQ_DATE_UPPER', dph[tmp, GENBANK_SEQDATE_UPPER])
    set(dph, NULL, c('GENBANK_SEQDATE','GENBANK_SEQDATE_LOWER','GENBANK_SEQDATE_UPPER'), NULL)
    #   drop tips
    dph.old <- subset(dph, select=c(TAXA, SEQ_DATE, SEQ_DATE_LOWER, SEQ_DATE_UPPER))
    tmp <- subset(dph, is.na(SEQ_DATE))[, TAXA_ID]
    cat('\nDropping tips without sampling date from ', infile,' n=', length(tmp), 'of Ntips=', Ntip(ph), '\n')
    ph <- phyloscanner.to.simmap(ph)
    ph <- phytools:::drop.tip.simmap(ph, ph$tip.label[tmp])
    ph <- simmap.to.phyloscanner(ph)
    
    #
    #   date tree
    #
    
    #   make data.table of sequence sampling times
    dph <- data.table(  TAXA_LABEL=ph$tip.label,
                        TAXA= gsub('^[^_]+___(.*)','\\1',ph$tip.label),
                        TAXA_ID= seq_along(ph$tip.label))
    dph <- merge(dph, dph.old, by= 'TAXA')
    stopifnot( !any(is.na(dph$SEQ_DATE)) )
    dph <- dph[order(TAXA_ID),]
    #   get into format needed for tree.dater
    sampling.times.init <- dph$SEQ_DATE
    names(sampling.times.init) <- dph$TAXA_LABEL
    sampling.times.bounds <- as.data.frame(subset(dph, select=c(SEQ_DATE_LOWER, SEQ_DATE_UPPER)))
    rownames(sampling.times.bounds) <- dph$TAXA_LABEL
    colnames(sampling.times.bounds) <- c('lower','upper')
    #   save files for dating
    outfile.collapsed.phsc <- file.path(outdir, gsub('workspace\\.rda','collapsed_workspace\\.rda',basename(infile)))
    outfile.tree <- file.path(outdir, gsub('workspace\\.rda','collapsed\\.newick',basename(infile)))
    outfile.sampling.times.bounds <- file.path(outdir, gsub('workspace\\.rda','sampling_times_bounds\\.csv',basename(infile)))
    outfile.sampling.times.init <- file.path(outdir, gsub('workspace\\.rda','sampling_times_init\\.csv',basename(infile)))
    save(ph, file=outfile.collapsed.phsc)
    write.tree(ph, file=outfile.tree)
    write.csv(sampling.times.bounds, file=outfile.sampling.times.bounds, row.names=TRUE)
    write.csv(data.table(TAXA=names(sampling.times.init), SEQ_DATE=sampling.times.init), file=outfile.sampling.times.init, row.names=FALSE)
    control <- list( outfile=gsub('.newick$','_dated.newick',outfile.tree),
                     ali.len=alignment.length,
                     root=NA,
                     omega0=NA,
                     temporalConstraints=TRUE,
                     clock='uncorrelated',
                     estimateSampleTimes=outfile.sampling.times.bounds
    )
    cmd <- cmd.treedater.script(outfile.tree, outfile.sampling.times.init, control=control)
    cmds[[i]] <- cmd
  }
  infiles.phsc[, CMD:= unlist(cmds)]
  
  #   submit jobs like this one:
  cat(infiles.phsc[1,CMD])
  
  #   run on HPC as array job
  df <- copy(infiles.phsc)
  df[, CASE_ID:= 1:nrow(df)]
  #   make PBS header
  hpc.load    <- "module load anaconda3/personal"
  r.activate  <- "source activate phylo"
  hpc.select  <- 1                        # number of nodes
  hpc.nproc   <- 1                        # number of processors on node
  hpc.walltime<- 23                       # walltime
  hpc.q       <- NA #"pqeelab"                        # PBS queue
  hpc.mem     <- "4gb"                    # RAM
  hpc.array   <- length(unique(df$CASE_ID))   # number of runs for job array
  pbshead     <- "#!/bin/sh"
  tmp         <- paste("#PBS -l walltime=", hpc.walltime, ":59:00,pcput=", hpc.walltime, ":45:00", sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  tmp         <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead     <- paste(pbshead, tmp, sep = "\n")
  pbshead     <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(!is.na(hpc.array))
    pbshead <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  pbshead     <- paste(pbshead, hpc.load, r.activate, sep = "\n")
  #   make array job
  cmd     <- df[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd     <- cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]
  cmd     <- paste(pbshead,cmd,sep='\n')
  #   submit job
  outfile     <- gsub(':','',paste("phs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
  outfile     <- file.path(outdir, outfile)
  cat(cmd, file=outfile)
  cmd         <- paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}


## ---- rmd.chunk.roadmap.200203.treedater.postprocess ----
#  superimpose the dated branch lengths onto the tree with ancestral state reconstructions, and plot the result
roadmap.200203.treedater.postprocess <- function(analysis)
{
  require(data.table)
  require(ape)
  require(phyloscannerR)
  
	analysis <- 'analysis_200917'
	
  #home <- '~/Box Sync/Roadmap'
  home <- file.path('/rds/general/project/ratmann_roadmap_data_analysis/live',analysis)
  indir <- file.path(home,'phyloscanner_dated')
  infiles <- data.table(F_PHSC=list.files(indir, pattern='collapsed_workspace.rda$', full.names=TRUE, recursive=TRUE))
  infiles[, BASENAME:= gsub('collapsed_workspace.rda$','',basename(F_PHSC))]
  tmp <- data.table(F_TREE=list.files(indir, pattern='collapsed_dated.newick$', full.names=TRUE, recursive=TRUE))
  tmp[, BASENAME:= gsub('collapsed_dated.newick$','',basename(F_TREE))]
  infiles <- merge(infiles, tmp, by='BASENAME', all.x=TRUE)
  #infiles <- subset(infiles,grepl('subtype_B',BASENAME))
  
  stopifnot( nrow(subset(infiles, is.na(F_TREE)))==0 )

  for(i in seq_len(nrow(infiles)))
  {
    #i<- 1
    cat('\nProcess ',infiles[i,F_TREE])
    #   load data
    load(infiles[i,F_PHSC])
    ph.dated <- read.tree(infiles[i,F_TREE])
    stopifnot( all(  ph.dated$tip.label == ph$tip.label ) )
    stopifnot( all( ph.dated$edge == ph$edge ) )
    
    #   since the tree topology is unchanged, we can copy
    #   the branch lenghts in units of time onto the original tree
    #   that has the ancestral state reconstructions
    ph$edge.length <- ph.dated$edge.length
    
    #   plot dated tree to spot obvious errors
    outfile <- gsub('collapsed_workspace\\.rda','annotated_dated_tree.pdf',infiles[i,F_PHSC])
    tmp <- vector('list')
    tmp[['tree']] <- ph
    tmp[['tree']][['node.states']] <- tmp[['tree']][['mapped.edge']] <- tmp[['tree']][['maps']] <- NULL
    attr(tmp[['tree']],'map.order') <- NULL
    attr(tmp[['tree']],'class') <- 'phylo'
    tmp[['read.counts']] <- rep(1, Ntip(ph))
    write.annotated.tree(tmp, outfile, format="pdf", pdf.scale.bar.width = 0.01, pdf.w = 60, pdf.hm = 0.2, verbose = FALSE)
    
    #   save phyloscanner.tree
    outfile <- gsub('collapsed_workspace\\.rda','annotated_dated_tree.rda',infiles[i,F_PHSC])
    save(ph, file=outfile)
  }
}

## ---- rmd.chunk.roadmap.200203.get.dated.subgraphs ----
roadmap.200203.phyloscanner.get.dated.subgraphs<- function(analysis)
{
  require(data.table)
  require(phangorn)
  require(ggplot2)
  require(reshape)
  require(phyloscannerR)
  
  #    working directory with phyloscanner output
  source('/rds/general/project/ratmann_roadmap_data_analysis/live/R/bigphylo functions.R')
	analysis <- 'analysis_200917'
	
  #home <- '~/Box Sync/Roadmap/Data/SHM_1902_ROADMAP_191223'
  home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
  indir.phsc    <- file.path(home,analysis,'phyloscanner_dated')
  outdir <- file.path(home,analysis,'subgraphs_dated')
  infiles        <- data.table(F=list.files(indir.phsc, pattern='_annotated_dated_tree.rda$', full.names=TRUE, recursive=TRUE))
  infiles[, SELECT:= gsub('^.*_rerooted_([A-Za-z0-9]+)_.*$','\\1',basename(F))]

  #    extract subgraphs of dated trees
  for(i in seq_len(nrow(infiles)))
    {
    cat('process', i,'\n')
    infile <- infiles[i, F]
    host <- infiles[i,SELECT]
    load(infile)
    mrcas <- which( attr(ph, 'SUBGRAPH_MRCA') )
    #    some of the tips states are "unknown" and this clashes internally with the NA state, so we need to take some extra care
    #    this is not a problem because the "unknown" and NA state mean the same thing
    attr(ph, 'INDIVIDUAL') <- as.character(attr(ph, 'INDIVIDUAL'))
    attr(ph, 'INDIVIDUAL')[is.na(attr(ph, 'INDIVIDUAL'))] <- 'Unknown'
    mrcas <- mrcas[ attr(ph, 'INDIVIDUAL')[mrcas]==host ]
    stopifnot( !any(is.na(mrcas)) )
    # check that tree is of class simmap
    stopifnot( any(attr(ph,'class')=='simmap') )
    # extract subgraphs
    subgraphs <- lapply(mrcas, function(mrca) extract.subgraph(ph, mrca))
    # save
    outfile <- gsub('_annotated_dated_tree',paste0('_datedsubgraphs_',host),basename(infile))
    outfile <- file.path(outdir,outfile)
    save(subgraphs, file=outfile)
  }
}

### Read phylogenetic subgraphs
require(data.table)
require(phangorn)
require(ggplot2)
require(reshape)

#    working directory with phyloscanner output
home <- '~/Box Sync/Roadmap'
analysis <- 'analysis_200917'
#home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
indir.phsc    <- file.path(home,analysis,'subgraphs_dated')
infile.meta <- file.path(home,analysis,'misc','200917_sequence_labels.rda')
outdir <- file.path(home,analysis,'subgraph_analysis')

#    extract subgraph taxa
infiles <- data.table(F=list.files(indir.phsc, pattern='rda$', full.names=TRUE, recursive=TRUE))
infiles <- subset(infiles, grepl('datedsubgraphs',F))
#    SELECT defines if Ams, AmsHSX, AmsMSM
infiles[, SELECT:= gsub('^(.*)subgraphs_([A-Za-z]+).rda','\\2',basename(F))]
#    for subtype B, I ran separate analyses based on very large subtrees for comp efficiency
#    ST stores the subtype, and ST_CLADE the large subtree
infiles[, ST:= gsub('^.*subtype_([^_]+)_.*\\.rda','\\1',basename(F))]
# For B only
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

# Number of distinct chains
nsubgraph <- dsubgraphtaxa[, list(N_SUBGRAPHS=length(NAME)), by=c('REP','SELECT')]

#   add meta data from pre-processed persons file
load(infile.meta)
dind <- as.data.table(dind)
setnames(dind, c('PATIENT','BIRTH_Y','BIRTH_CNTRY'), c('ID','BIRTH_YEAR','BIRTH_COUNTRY'))
tmp <- subset(dind, select=c(ID,BIRTH_YEAR,BIRTH_COUNTRY,LOC_BIRTH,CITY,TRANSM,GENDER))
set(tmp, NULL, 'BIRTH_YEAR', tmp[, as.numeric(BIRTH_YEAR)])
dsubgraphtaxa <- merge(dsubgraphtaxa,tmp,by='ID')
dsubgraphtaxa[, AGE2020:= cut(2020-BIRTH_YEAR, breaks=c(-Inf,20,25,30,35,40,45,50,55,60,Inf))]

outfile <- file.path(outdir,'subgraphs_withmetadata.rda')
save(dsubgraphtaxa, file=outfile)

# Plot size of subgraphs
dsubgraphsize <- dsubgraphtaxa[, list(SIZE=length(ID)), by=c('ST','REP','SELECT','NAME')]

ggplot(subset(dsubgraphsize, SELECT!='Ams'), aes(x=SIZE, fill=ST)) +
	geom_bar() +
	facet_grid(SELECT~.) +
	scale_y_continuous(expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0)) +
	xlim(0,30) +
	ylim(0,3100) +
	theme_bw() +
	scale_fill_viridis_d(begin=0,end=1,alpha=1,option="viridis") +
	labs(x='\nSize of phylogenetic subgraphs attributed to Amsterdam',y='Count',fill='Subtype/CRF')
outfile    <- file.path(outdir,'dist_subgraph_sizes.pdf')
ggsave(outfile,w=8, h=6)

# plot median prop across reps & CIs
ans.quantiles=c(0.025,0.5,0.975)
ans	<- dsubgraphsize[, list(freq=length(NAME)), by=c('SELECT','REP','SIZE')]
tmp	<- ans[, list(SIZE=SIZE, P=freq/sum(freq)), by=c('SELECT','REP')]
ans	<- merge(ans, tmp, by=c('SELECT','REP','SIZE'))
total.replicates	<- ans[, length(unique(REP))]
ans	<- ans[, list(P=ans.quantiles, 
									QN=quantile(c(freq, rep(0, total.replicates-length(freq))), p=ans.quantiles),
									QF=quantile(c(P, rep(0, total.replicates-length(P))), p=ans.quantiles)
), by=c('SELECT','SIZE')]
set(ans, NULL, 'P', ans[, paste0('p',P)])
save(ans,file=file.path(outdir,'dist_obs_subgraph_sizes.rds'))

tmp			<- dcast.data.table(ans, SELECT+SIZE~P, value.var='QF')
ggplot(subset(tmp,SELECT!='Ams'), aes(x=SIZE, y=p0.5)) + 
	geom_bar(position=position_dodge(), stat='identity') +	
	geom_errorbar(aes(ymin=p0.025, ymax=p0.975),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(~SELECT, scales='free') +
	#xlim(0,30) +
	scale_y_continuous(labels=scales:::percent) + 
	theme_bw() +
	labs(x='\nsubgraph size\n observed in phylogeny', y=paste0('%','\n'), fill='analysis\ngroup')
outfile    <- file.path(outdir,'dist_obs_subgraph_sizes.pdf')
ggsave(outfile,w=8, h=6)

# Proportion singletons
dsubgraphsize[, SINGLETON:= as.character(factor(SIZE==1, levels=c(TRUE,FALSE), labels=c('subgraph size 1','subgraph size >1')))]
singletons <- dsubgraphsize[, list(N=length(NAME)), by=c('SELECT','SINGLETON')]
singletons

# size of subgraphs by origin
originsize <- dsubgraphtaxa[, list(SIZE=length(ID)), by=c('ST','REP','SELECT','NAME','ORIGINHOST')]

originsize <- subset(originsize, ORIGINHOST!='NA' & SELECT!='Ams')
originsize[, PARENT_STATE:= as.character(factor(grepl('Ams',ORIGINHOST), levels=c(TRUE,FALSE), labels=c('Ams','EXT')))]
originsize$PARENT_STATE[grepl('NL',originsize$ORIGINHOST)] <- 'NL'

originsize$PARENT_STATE <- factor(originsize$PARENT_STATE,levels=c("NL","Ams","EXT"))
lm <- lm(SIZE ~ SELECT+PARENT_STATE, data = originsize)
summary(lm)

ggplot(subset(originsize, SELECT!='Ams'), aes(x=SIZE, fill=ORIGINHOST)) +
	geom_bar() +
	facet_grid(SELECT~.) +
	scale_y_continuous(expand=c(0,0)) +
	scale_x_continuous(expand=c(0,0)) +
	xlim(0,30) +
	ylim(0,3100) +
	theme_bw() +
	scale_fill_viridis_d(begin=0,end=1,alpha=1,option="viridis") +
	labs(x='\nOrigin of phylogenetic subgraphs',y='Count',fill='Geographic region \nof origin')
outfile    <- file.path(outdir,'subgraph_sizes_byorigin_georeg.pdf')
ggsave(outfile,w=8, h=6)

os <- originsize[, list(N=mean(SIZE)), by=c('PARENT_STATE','SELECT')]

ggplot(subset(os, SELECT!='Ams'), aes(x=PARENT_STATE, y=N)) +
	geom_bar(stat='identity') +
	facet_grid(SELECT~.) +
	theme_bw() +
	scale_fill_viridis_d(begin=0,end=1,alpha=1,option="viridis") +
	labs(x='\nOrigin of phylogenetic subgraphs attributed to Amsterdam',y='Mean size of subgraphs',fill='Subtype/CRF')
ggsave(file.path(outdir,'subgraph_sizes_byorigin.pdf'),w=8, h=6)

# Origin of introductions
origin <-  dsubgraphtaxa[, list(N=length(NAME)), by=c('SELECT','ORIGINHOST')]
origin$ORIGINHOST[is.na(origin$ORIGINHOST)] <- 'Unknown'
origin$ORIGINHOST <- factor(origin$ORIGINHOST, levels=c('AmsnonHSX','AmsnonMSM','NL','WEurope','EEuropeCentralAsia','Africa','Asia','FormerCurrDutchColonies','LaAmCar','MENA','NorthAm','Oceania','Unknown'),labels=c('Ams (non-HSX)','Ams (non-MSM)','NL','W.Europe','E.Europe & Central Asia','Africa','Asia','Dutch Colonies','Latin America and Caribbean','MENA','N.America','Oceania','Unknown'))
table(dsubgraphtaxa$ORIGIN[dsubgraphtaxa$SELECT=='AmsHSX'],useNA="always")
table(dsubgraphtaxa$ORIGIN[dsubgraphtaxa$SELECT=='AmsMSM'],useNA="always")
ggplot(subset(origin, SELECT!='Ams'), aes(x=SELECT, y=N, fill=ORIGINHOST)) +
	geom_bar(stat='identity') +
	labs(y='#subgraphs',x='',fill='origin of introductions') +
	theme_bw() +
	theme(axis.ticks.x=element_blank()) 
outfile    <- file.path(outdir,'origin_host.png')
ggsave(outfile)

dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
dsubgraphtaxa$ORIGIN <- factor(dsubgraphtaxa$ORIGIN, levels=c('AmsnonHSX','AmsnonMSM','NL','WEurope','EEuropeCentralAsia','Africa','Asia','FormerCurrDutchColonies','LaAmCar','MENA','NorthAm','Oceania','Unknown'),labels=c('Ams (non-HSX)','Ams (non-MSM)','NL','W.Europe','E.Europe & Central Asia','Africa','Asia','Dutch Colonies','Latin America and Caribbean','MENA','N.America','Oceania','Unknown'))
ggplot(subset(dsubgraphtaxa, SELECT!='Ams'), aes(x=SELECT, fill=ORIGIN)) +
	geom_bar(position='fill') +
	labs(y='subgraphs',x='',fill='origin of introductions') +
	theme_bw() +
	scale_y_continuous(labels = scales::percent) +
	theme(axis.ticks.x=element_blank()) 
outfile    <- file.path(outdir,'origin_host_pct.png')
ggsave(outfile)

# save origins for estimating proportion of external introductions
#origins <- dsubgraphtaxa[dsubgraphtaxa$REP=='000']
origins <- dsubgraphtaxa
origins <- origins[,list(N=length(TAXA)),by=c('ORIGINHOST','TRANSM')]
setnames(origins,c('ORIGINHOST','TRANSM'),c('PARENT_STATE','TRM_GROUP'))
outfile    <- file.path(outdir,'200203_subgraph_origins_N.csv')
write.csv(origins,file=outfile)


# Proportion of birth place by world regions
dsubgraphtaxa$LOC_BIRTH[dsubgraphtaxa$BIRTH_COUNTRY=='Netherlands'] <- dsubgraphtaxa$BIRTH_COUNTRY[dsubgraphtaxa$BIRTH_COUNTRY=='Netherlands']
origin <- dsubgraphtaxa[, list(N=length(NAME)), by=c('SELECT','LOC_BIRTH')]
tmp <- dsubgraphtaxa[, list(TOTAL=length(NAME)), by=c('SELECT')]
origin <- merge(origin,tmp)
origin[,'PROP'] <- origin[,N]/origin[,TOTAL]

ggplot(subset(origin, SELECT!='Ams'), aes(x=SELECT, fill=LOC_BIRTH)) +
	geom_bar(position='fill') +
	labs(y='',x='',fill='place of birth') +
	scale_y_continuous(labels = scales::percent) +
	theme_bw() +
	scale_fill_viridis_d(begin=0,end=1,alpha=1,option="viridis") +
	theme(axis.ticks.x=element_blank()) 
outfile    <- file.path(outdir,'birthplace.pdf')
ggsave(outfile)

# Size of subgraphs for all STs/Reps by HSX/MSM (for Stan model)
dsubgraphsize <- dsubgraphtaxa[, list(SIZE=length(ID)), by=c('ST','REP','SELECT','NAME')]
dsubgraphsizefreq <- dsubgraphsize[, list(N=length(NAME)), by=c('SELECT','SIZE')]
dsubgraphsizefreq <- dsubgraphsizefreq[SELECT!='Ams',]
#dsubgraphsize2 <- dsubgraphsizefreq[, list(N=round(N/11,0)), by=c('SELECT','SIZE')]
#dsubgraphsize2 <- dsubgraphsize2[dsubgraphsize2$N!=0,]

freqs <- dsubgraphsizefreq[order(SELECT,SIZE),]
#freqs <- dsubgraphsize2[order(SELECT,SIZE),]
freqshsx <- freqs[SELECT=='AmsHSX',]
freqsmsm <- freqs[SELECT=='AmsMSM',]

listhsx <- data.frame(SELECT=rep('AmsHSX',max(freqshsx$SIZE)),SIZE=1:max(freqshsx$SIZE))
listmsm <- data.frame(SELECT=rep('AmsMSM',max(freqsmsm$SIZE)),SIZE=1:max(freqsmsm$SIZE))

freqshsx <- merge(listhsx,freqshsx,by=c("SELECT","SIZE"),all=TRUE)
freqshsx$N[is.na(freqshsx$N)] <- 0
freqsmsm <- merge(listmsm,freqsmsm,by=c("SELECT","SIZE"),all=TRUE)
freqsmsm$N[is.na(freqsmsm$N)] <- 0

outfilehsx    <- file.path(outdir,'subgraphsizes_hsx_m.csv')
outfilemsm    <- file.path(outdir,'subgraphsizes_msm_m.csv')
write.csv(freqshsx,file=outfilehsx)
write.csv(freqsmsm,file=outfilemsm)

# calculate pr(intros) using partially observed transmission chains
freqshsx$chains <- freqshsx$SIZE*freqshsx$N
sum(freqshsx$N)/sum(freqshsx$chains)
freqsmsm$chains <- freqsmsm$SIZE*freqsmsm$N
sum(freqsmsm$N)/sum(freqsmsm$chains)

# Sample size for phylodynamic analysis of HIV spread in Amsterdam
#   number of individuals in subgraphs of size == 1 or size > 1
dsubgraphind <- dsubgraphsize[, list(N_INDIVIDUALS=sum(SIZE)), by=c('REP','ST','SELECT','SINGLETON')]
dsubgraphind <- dcast.data.table(dsubgraphind, REP+ST+SINGLETON~SELECT,value.var='N_INDIVIDUALS')
dsubgraphind

# Compare subtype B vs. non-B
dsubgraphsize$ST2 <- 'Non-B'
dsubgraphsize$ST2[dsubgraphsize$ST=='B'] <- 'B'
dsubgraphind2 <- dsubgraphsize[, list(N_INDIVIDUALS=sum(SIZE)), by=c('ST2','SELECT','SINGLETON')]
# Proportion of subclades size 1 vs. >1
dsubgraphsin <- dsubgraphsize[, list(N_SINGLETONCLADES=length(NAME)), by=c('ST2','SELECT','SINGLETON')]

# Empirical analysis on group-mixing within subgraphs
#   number of individuals in subgraphs of size == 1 or size > 1
tmp <- subset(dsubgraphsize, SIZE>1, c(ST, REP, SELECT, NAME))
dsubgraphtaxa2 <- merge(dsubgraphtaxa, tmp, by=c('ST','REP','SELECT','NAME'))

# Place of birth
dsubgraphtaxa2[, DUMMY:= as.character(factor(BIRTH_COUNTRY=='Netherlands',levels=c(TRUE,FALSE),labels=c('Netherlands','Other')))]
dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','SELECT','ST','DUMMY')]
tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','SELECT','ST')]
dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','SELECT','ST'))
ggplot(subset(dcnt, SELECT!='Ams'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
  geom_bar(stat='identity') +
  labs(y='#individuals',x='phylo subgraphs',fill='birth place') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(SELECT~ST, scale='free_x', space='free_x')
outfile    <- file.path(indir.phsc,'place_birth.pdf')
ggsave(outfile)

tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]

# Show for all subclades for subtype X, how many/what proportion is >75% NL/>75% other and otherwise mixed?
dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','REP','SELECT','ST','DUMMY')]
dcnt$PROP <- dcnt$N/dcnt$TOTAL
bp <- dcnt[, list(N75= length(FULL_NAME[PROP>=0.75])), by=c('SELECT','REP','ST','DUMMY')]
nsubgraphs <- dsubgraphtaxa2[, list(SIZE=length(unique(NAME))), by=c('REP','ST','SELECT')]
bp <- merge(bp,nsubgraphs, by=c('REP','ST','SELECT'))
mixed <- bp[, list(Mixed= SIZE-sum(N75)), by=c('REP','ST','SELECT')]
bp <- spread(bp, DUMMY, N75)
bp <- merge(bp,mixed,by=c('REP','ST','SELECT'))
bp <- unique(bp)
bp <- melt(bp, id.vars=c('REP', 'ST', 'SELECT', 'SIZE'))
#bp <- bp[order('REP','ST','SELECT'),]

bp$bpp <- factor(bp$variable, levels=c('Netherlands','Other','Mixed'),labels=c(">75% NL", ">75% MW-MB", "Mixed"))

nonbl <- ggplot(subset(bp, SELECT!='Ams' & ST!='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank()) +
  facet_grid(SELECT~., scale='free_x', space='free_x')
legend <- get_legend(nonbl)

nonb <- ggplot(subset(bp, SELECT!='Ams' & ST!='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank(),legend.position="none") +
  facet_grid(SELECT~., scale='free_x', space='free_x')

b <- ggplot(subset(bp, SELECT!='Ams' & ST=='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank(),legend.position="none") +
  facet_grid(SELECT~., scale='free_x', space='free_x')

# Combine the graphs
outfile <- file.path(outpath,"mixing_bp_BnonB.pdf")
bp_bnonb <-   grid.arrange(nonb, b,
                      legend,
                      ncol = 3, widths=c(3.5, 1.5, 1.5))
ggsave(file=outfile, bp_bnonb, w=8, h=6)

# All STs together
ggplot(subset(bp, SELECT!='Ams'), aes(x=ST, fill=bpp)) +
  geom_bar(position = "fill") +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank()) +
facet_grid(SELECT~., scale='free_x', space='free_x')
outfile    <- file.path(indir.phsc,'place_birth_75pct.pdf')
ggsave(outfile)

# Age group
dsubgraphtaxa2[, DUMMY:= AGE2020]
dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','SELECT','ST','DUMMY')]
tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','SELECT','ST')]
dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','SELECT','ST'))
ggplot(subset(dcnt, SELECT!='Ams'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
  geom_bar(stat='identity') +
  labs(y='#individuals',x='phylo subgraphs',fill='age in 2020') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(SELECT~ST, scale='free_x', space='free_x')
outfile    <- file.path(indir.phsc,'age_group.pdf')
ggsave(outfile)

tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]

# Geographic region of birth
dsubgraphtaxa2$LOC_BIRTH[dsubgraphtaxa2$BIRTH_COUNTRY=='Netherlands'] <- dsubgraphtaxa2$BIRTH_COUNTRY[dsubgraphtaxa2$BIRTH_COUNTRY=='Netherlands']
dsubgraphtaxa2[, DUMMY:= factor(LOC_BIRTH, levels=c('Netherlands','WEurope','EEurope','Africa','Asia','LaAmCarr','NorthAm','Oceania','Unknown'),labels=c('Netherlands','WEurope','EEurope','Africa','Asia','LaAmCarr','NorthAm','Oceania','Unknown'))]
dcnt <- dsubgraphtaxa2[, list(N= length(ID)), by=c('FULL_NAME','REP','SELECT','ST','DUMMY')]
tmp <- dsubgraphtaxa2[, list(TOTAL= length(ID)), by=c('FULL_NAME','REP','SELECT','ST')]
dcnt <- merge(dcnt, tmp, by=c('FULL_NAME','REP','SELECT','ST'))
ggplot(subset(dcnt, SELECT!='Ams'), aes(x=FULL_NAME, y=N, fill=DUMMY)) +
  geom_bar(stat='identity') +
  labs(y='#individuals',x='phylo subgraphs',fill='birth place') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  facet_grid(SELECT~ST, scale='free_x', space='free_x')
outfile    <- file.path(indir.phsc,'region_birth.pdf')
ggsave(outfile)

tmp <- dcnt[, list(R= sum( (N/TOTAL)^2), PMAX= max(N/TOTAL) ), by=c('FULL_NAME','SELECT','ST')]
tmp[, list(R_AVG=mean(R), PMAX_AVG=mean(PMAX)), by=c('SELECT')]

# Show for all subclades for subtype X, how many/what proportion is >75% NL/>75% other and otherwise mixed?
dcnt$PROP <- dcnt$N/dcnt$TOTAL
bp <- dcnt[, list(N75= length(FULL_NAME[PROP>=0.75])), by=c('SELECT','REP','ST','DUMMY')]
nsubgraphs <- dsubgraphtaxa2[, list(SIZE=length(unique(NAME))), by=c('REP','ST','SELECT')]
bp <- merge(bp,nsubgraphs, by=c('REP','ST','SELECT'))
mixed <- bp[, list(Mixed= SIZE-sum(N75)), by=c('REP','ST','SELECT')]
bp <- spread(bp, DUMMY, N75)
bp <- merge(bp,mixed,by=c('REP','ST','SELECT'))
bp <- unique(bp)
bp <- melt(bp, id.vars=c('REP', 'ST', 'SELECT', 'SIZE'))

bp$bpp <- factor(bp$variable, levels=c('Netherlands','WEurope','EEurope','Africa','Asia','LaAmCarr','NorthAm','Oceania','Unknown','Mixed'),labels=c('>75% Netherlands','>75% WEurope','>75% EEurope','>75% Africa','>75% Asia','>75% LaAmCarr','>75% NorthAm','>75% Oceania','>75% Unknown','Mixed'))

ggplot(subset(bp, SELECT!='Ams'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank()) +
  facet_grid(SELECT~., scale='free_x', space='free_x')
outfile    <- file.path(indir.phsc,'region_birth_75pct.pdf')
ggsave(outfile)

nonbl <-  ggplot(subset(bp, SELECT!='Ams' & ST!='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank()) +
  facet_grid(SELECT~., scale='free_x', space='free_x')
legend <- get_legend(nonbl)

nonb <- ggplot(subset(bp, SELECT!='Ams' & ST!='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='#subgraphs',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank(),legend.position="none") +
  facet_grid(SELECT~., scale='free_x', space='free_x')

b <-   ggplot(subset(bp, SELECT!='Ams' & ST=='B'), aes(x=ST, y=value, fill=bpp)) +
  geom_bar(stat='identity') +
  labs(y='',x='subtype',fill='birth place') +
  theme(axis.ticks.x=element_blank(),legend.position="none") +
  facet_grid(SELECT~., scale='free_x', space='free_x')

# Combine the graphs
outfile <- file.path(outpath,"mixing_bp_reg_BnonB.pdf")
bp_bnonb <-   grid.arrange(nonb, b,
                           legend,
                           ncol = 3, widths=c(3.5, 1.5, 1.5))
ggsave(file=outfile, bp_bnonb, w=8, h=6)


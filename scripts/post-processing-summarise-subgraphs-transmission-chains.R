
cat(" \n -------------------------------- \n \n summarise subgraphs and chain sizes.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggsci, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_211101'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['period']] <- '2014-2018'
	args_dir[['job_name']] <- 'undiagnosed_weighted_inf_rate'
	args_dir[['job_tag']] <- paste0(args_dir[['job_name']],'_',args_dir[['period']],'_',args_dir[['trsm']])
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/',args_dir[['stanModelFile']],'-',args_dir[['job_tag']])
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
	args_dir[['infdate']] <- 1
	args_dir[['start_y']] <- 2014
	args_dir[['p_undiag_gps']]='G1,G2,G3,NL,Other'
	args_dir[['p_undiag_probs']]='0.1139,0.1729,0.1812,0.13880,0.2268'
	args_dir[['undiag_job']]= 'undiag_untilmay2019_weights'
	
}

# save args for report before loading those from running session 
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
	stopifnot(args_line[[1]]=='-stanModelFile')	
	stopifnot(args_line[[3]]=='-analysis')
	stopifnot(args_line[[5]]=='-in_dir')
	stopifnot(args_line[[7]]=='-out_dir')
	stopifnot(args_line[[9]]=='-job_tag')
	stopifnot(args_line[[11]]=='-trsm')
	stopifnot(args_line[[13]]=='-source_dir')
	stopifnot(args_line[[15]]=='-start_y')
	stopifnot(args_line[[17]]=='-undiag_job')
	
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['analysis']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['out_dir']] <- args_line[[8]]
	args_dir[['job_tag']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
	args_dir[['source_dir']] <- args_line[[14]]
	args_dir[['start_y']] <- args_line[[16]]
	args_dir[['undiag_job']] <- args_line[[18]]
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))
source(file.path(args_dir$source_dir, 'R', 'stan-functions.r' ))

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)

#	load all input variables for this analysis run
do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes
start_d = args$start_d
end_d = args$end_d

p_undiagnosed <- data.table(mwmb=args$p_undiag_gps,du=args$p_undiag_probs)

dsubgraphtaxa <- readRDS(file.path(args_dir[['out_dir']],'subgraphs_withmetadata.RDS'))

infile.meta <- file.path(args$indir, args$analysis, 'misc', '200917_sequence_labels.rda')
infile.seq <-	file.path(args$indir, 'Data', 'data_200821/SHM_1902_ROADMAP_200821_tblLAB_seq.rda')
infile.cd4 <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_CD4.csv')
infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')

cat("\nRead patient metadata \n")
load(infile.seq)
load(infile.meta)
setnames(dind, gsub("CENS_D", "RECART_D", names(dind)))

dind_diag <- copy(dind)
dind_diag$SEQ <- dind_diag$PATIENT %in% ds$PATIENT
dind_diag <- merge(dind_diag,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)
dind_diag <- data.table(dind_diag)

cat('\n Add birth place region \n')
dind <- data.table(dind)
dind <- add_georeg_bplace_ind_data(dind,infile.geo)
dind <- map_mwmb_regions(dind)

cat("\nLoad time to diagnosis estimates \n")
dt_m <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_MSM.RDS'))
dt_h <- file.path(args$indir,'branching_process_model',paste0('undiagnosed_211102-',args$undiag_job),paste0('median_timetodiagnosis_',args$undiag_job,'_HSX.RDS'))
dt_m <- readRDS(dt_m)
dt_h <- readRDS(dt_h)
dt_inf <- rbind(dt_m,dt_h)
setnames(dt_inf,'trsm','TRANSM')

cat("\nUpdate col names if using infection date \n")
if(args$infdate==1){
	infile.inftime <- file.path(args$indir,'Data','infection_time_estimates','roadmap_cd4_vl_est.csv')
	#dsubgraphtaxa <- add_infection_time(dsubgraphtaxa,infile.inftime)
	dind <- data.table(dind)
	dind[, ID:=PATIENT]
	dind <- add_infection_time(dind,infile.inftime,dt_inf)
	dsubgraphtaxa[, HIV1_POS_D:=INF_D]
	dind[, HIV1_POS_D:=INF_D]
}

# flag which diagnosed individuals have a sequence
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
dind <- merge(dind,unique(subset(dsubgraphtaxa,select=c(ID,ST))),by.x='PATIENT',by.y='ID',all.x=T)
#ams <- dind[dind$CITY=='Amsterdam' & dind$HIV1_POS_D>=2010,]
nl <- dind[dind$HIV1_POS_D>=args_dir$start_y,]
ams <- dind[dind$CITY=='Amsterdam' & dind$HIV1_POS_D>=args_dir$start_y & dind$HIV1_POS_D<args$end_d,]

## add CD4 counts at diagnosis
ams <- add_cd4_counts(ams,infile.cd4,args_dir$start_y,NA)
dind_diag[, ID:=PATIENT]
dind_diag <- add_cd4_counts(dind_diag,infile.cd4,args_dir$start_y,NA)
dind[, ID:=PATIENT]
dind <- add_cd4_counts(dind,infile.cd4,args_dir$start_y,NA)

## add georegions/migrant groups

geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))
				 
#ams <- merge(ams,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
dind_diag <- merge(dind_diag,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
ams[BIRTH_CNTRY=='Tajikistan',WRLD_born:='Asia']
ams[WRLD_born=='PacificOceania',WRLD_born:='Asia']
ams[WRLD_born=='Oceania',WRLD_born:='Australia & New Zealand']
ams[BIRTH_CNTRY=='Netherlands',WRLD_born:='Netherlands']

#ams <- map_mwmb_regions(ams)
dind_diag <- map_mwmb_regions(dind_diag)

### summarise subgraphs
chains <- subset(dsubgraphtaxa,REP=="000" & CITY=="Amsterdam" & SELECT!='Ams')
chains[,SG:= FULL_NAME]
chains$SG <- as.factor(chains$SG)
chains[!is.na(ST_CLADE),SG:= paste0(SG,'_',ST_CLADE)]

# get min/max date of diagnosis date in subgraphs
mind <- subset(chains,REP=='000')[,list(MIND=min(HIV1_POS_D,na.rm=T),MAXD=max(HIV1_POS_D,na.rm=T)),by=c('SG')]
chains <- merge(chains,mind,by=c('SG'))

# count chains with no new case since 2010 as dead (not considered in analysis)
chains$dead_p1 <- "ongoing"
chains$dead_p1[chains$MAXD<args_dir$start_y] <- "dead"

# count chains which started before 2015
chains$period_p1 <- "pre-startd"
chains$period_p1[chains$MIND>=args_dir$start_y] <- "post-startd"

dsubgraphtaxa <- subset(dsubgraphtaxa, HIV1_POS_D<end_d & HIV1_POS_D>1912)
dind <- subset(dind, HIV1_POS_D<=end_d)


infile.rna <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_RNA.csv')
infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
#### add last viral load before 2010
dat <- read.csv(infile.rna,header=T)
dat$RNA_D <- as.Date(dat$RNA_D,format=c("%Y-%m-%d"))

# get death date for those who died
dbas <- read.csv(infile.indinfo,header=T)
ddeath <- subset(dbas,select=c('PATIENT','DEATH_D'))
setnames(ddeath,'PATIENT','ID')
ddeath[,'DEATH_D'] <- as.Date(ddeath[,'DEATH_D'],format="%Y-%m-%d")
ddeath[,'DEATH_D'] <- hivc.db.Date2numeric(ddeath[,'DEATH_D'])

# get last viral load before start of analysis
startd <- as.Date(0, origin = paste0(args_dir$start_y,"-01-01"),format="%Y-%m-%d")
tmp <- subset(dat, RNA_D<as.Date(startd))

tmp <- tmp %>%
	group_by(PATIENT) %>%
	filter(RNA_D == max(RNA_D))
setnames(tmp,c('RNA_D','RNA_V'),c('RNA_D_p1','RNA_V_p1'))

#chains <- merge(chains,ddeath,by=c('ID'),all.x=T)
#dsubgraphtaxa <- merge(dsubgraphtaxa,ddeath,by=c('ID'),all.x=T)
chains <- merge(chains,subset(tmp,select=c('PATIENT','RNA_D_p1','RNA_V_p1')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)
dind <- merge(dind,subset(tmp,select=c('PATIENT','RNA_D_p1','RNA_V_p1')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)
dsubgraphtaxa <- merge(dsubgraphtaxa,subset(tmp,select=c('PATIENT','RNA_D_p1','RNA_V_p1')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)

chains[, supp_p1:=1]
chains[RNA_V_p1>100 & HIV1_POS_D<args_dir$start_y, supp_p1:=0]
chains[DEATH_D<args_dir$start_y, supp_p1:=1]

cat(paste('\n Count index cases and new cases \n'))
# summarise the number of individuals in each subgraph diagnosed but not on ART by start date (icases) and individuals diagnosed between start date and end date
cases <- chains[, list(icases_p1=length(ID[HIV1_POS_D<args_dir$start_y & supp_p1==0]),
											 jcases_p1=length(ID[HIV1_POS_D>=args_dir$start_y])),
									 by=c('ST','ST_CLADE','REP','SELECT','FULL_NAME')]

freqs_p1 <- cases[, list(N=length(FULL_NAME)), by=c('SELECT','ST','icases_p1','jcases_p1')]
freqs_i_p1 <- freqs_p1[, list(subgraphs=sum(N)),by=c('icases_p1')]

sizes <- cases[, list(cases_p1=sum(jcases_p1)),by=c('ST','SELECT')]
sizes[SELECT=='AmsHSX',group:='Amsterdam heterosexual']
sizes[SELECT=='AmsMSM',group:='Amsterdam MSM']

cases_pre <- chains[, list(icases_p1=length(ID[HIV1_POS_D<args_dir$start_y & supp_p1==0]),
											 jcases_p1=length(ID[HIV1_POS_D>=args_dir$start_y])),
								by=c('ST','ST_CLADE','REP','SELECT','FULL_NAME','period_p1')]

freqs_pre_p1 <- cases_pre[, list(N=length(FULL_NAME[period_p1=='pre-startd'])), by=c('SELECT','ST','icases_p1','jcases_p1')]
#freqs_pre_i_p1 <- freqs_pre_p1[, list(subgraphs=sum(N)),by=c('icases_p1')]
freqs_pre_i_p1 <- freqs_pre_p1[, list(subgraphs=sum(N)),by=c('SELECT','icases_p1')]

## summarise chains which grew by whether all index cases were suppressed
freqs_pre_p1[, allsupp:=0]
freqs_pre_p1[icases_p1==0, allsupp:=1]
grew <- freqs_pre_p1[, list(pre=sum(N),
														pre_grew=sum(N[jcases_p1>0])),
										 by=c('SELECT','ST','allsupp')]
grew <- grew[order(SELECT,ST,allsupp),]
grew[, grew_pct:=paste0(round((pre_grew/pre)*100,0),"%")]

grew_tsm <- grew[, list(ST='Total',allsupp='Total',
										pre=sum(pre),pre_grew=sum(pre_grew)),
						 by=c('SELECT')]
grew_all <- grew[, list(SELECT='Total',ST='Total',allsupp='Total',
												pre=sum(pre),pre_grew=sum(pre_grew))]
grew_tsm[, grew_pct:=paste0(round((pre_grew/pre)*100,0),"%")]
grew_all[, grew_pct:=paste0(round((pre_grew/pre)*100,0),"%")]

grew[, allsupp:=as.character(allsupp)]
grew <- merge(grew,grew_tsm,by=c('SELECT','ST','allsupp','pre','pre_grew','grew_pct'),all=T)
grew <- merge(grew,grew_all,by=c('SELECT','ST','allsupp','pre','pre_grew','grew_pct'),all=T)
grew$SELECT <- factor(grew$SELECT,levels=c('AmsMSM','AmsHSX','Total'))
grew$ST <- factor(grew$ST,levels=c('B','nonB','Total'))
grew$allsupp <- factor(grew$allsupp,levels=c('1','0','Total'))

write.csv(grew,file=paste0(outfile.base,'-','number_of_subgraphs_table.csv'))

## summarise people in chains
dsubgraphtaxa[, supp_p1:=1]
dsubgraphtaxa[RNA_V_p1>100 & HIV1_POS_D<args_dir$start_y, supp_p1:=0]
dsubgraphtaxa[DEATH_D<args_dir$start_y, supp_p1:=1]
df <- subset(dsubgraphtaxa,REP=='000' & SELECT!='Ams')
cases2 <- chains[, list(icases_p1_notsupp=length(ID[HIV1_POS_D<args_dir$start_y & supp_p1==0]),
												icases_p1_supp=length(ID[HIV1_POS_D<args_dir$start_y & supp_p1!=0]),
												jcases_p1=length(ID[HIV1_POS_D>=args_dir$start_y])),
								 by=c('ST','ST_CLADE','REP','SELECT','FULL_NAME')]

df <- merge(df,cases2,by=c('SELECT','ST','ST_CLADE','FULL_NAME','REP'),all.x=T)
df[, allsupp:=0]
df[icases_p1_notsupp==0, allsupp:=1]
#df[, some_unsupp:=1]
#df[icases_p1_supp>=1, some_unsupp:=1]

ind <- df[, list(N_pre=length(unique(ID[INF_D<args_dir$start_y])),
								 N_pre_unsupp=length(unique(ID[INF_D<args_dir$start_y & supp==0])),
														N_em=length(unique(ID[INF_D>=args_dir$start_y & INF_D<2019]))),
										 by=c('SELECT','ST','allsupp')]
ind[, N_total:=N_pre+N_pre_unsupp+N_em]
ind[, pct_pre:=paste0(round((N_pre/N_total)*100,0),"%")]
ind[, pct_pre_unsupp:=paste0(round((N_pre_unsupp/N_total)*100,0),"%")]
ind[, pct_em:=paste0(round((N_em/N_total)*100,0),"%")]
ind <- ind[, c('SELECT','ST','allsupp','N_total','N_pre','pct_pre','N_pre_unsupp','pct_pre_unsupp','N_em','pct_em')]
tot_i <- ind[, list(ST='Total',allsupp='Total',
										N_total=sum(N_total),N_pre=sum(N_pre),N_pre_unsupp=sum(N_pre_unsupp),N_em=sum(N_em)),
											by=c('SELECT')]
tot_all <- ind[, list(SELECT='Total',ST='Total',allsupp='Total',
											N_total=sum(N_total),N_pre=sum(N_pre),N_pre_unsupp=sum(N_pre_unsupp),N_em=sum(N_em))]
tot_i[, N_total:=N_pre+N_pre_unsupp+N_em]
tot_i[, pct_pre:=paste0(round((N_pre/N_total)*100,0),"%")]
tot_i[, pct_pre_unsupp:=paste0(round((N_pre_unsupp/N_total)*100,0),"%")]
tot_i[, pct_em:=paste0(round((N_em/N_total)*100,0),"%")]
tot_all[, N_total:=N_pre+N_pre_unsupp+N_em]
tot_all[, pct_pre:=paste0(round((N_pre/N_total)*100,0),"%")]
tot_all[, pct_pre_unsupp:=paste0(round((N_pre_unsupp/N_total)*100,0),"%")]
tot_all[, pct_em:=paste0(round((N_em/N_total)*100,0),"%")]

ind[, allsupp:=as.character(allsupp)]
ind <- merge(ind,tot_i,by=c('SELECT','ST','allsupp','N_total','N_pre','pct_pre','N_pre_unsupp','pct_pre_unsupp','N_em','pct_em'),all=T)
ind <- merge(ind,tot_all,by=c('SELECT','ST','allsupp','N_total','N_pre','pct_pre','N_pre_unsupp','pct_pre_unsupp','N_em','pct_em'),all=T)
ind$SELECT <- factor(ind$SELECT,levels=c('AmsMSM','AmsHSX','Total'))
ind$ST <- factor(ind$ST,levels=c('B','nonB','Total'))
ind$allsupp <- factor(ind$allsupp,levels=c('1','0','Total'))

ind <- merge(grew,ind,by=c('SELECT','ST','allsupp'),all=T)

ind <- ind[order(SELECT,ST,allsupp),]
ind[SELECT=='AmsMSM',SELECT:='Amsterdam MSM']
ind[ST=='AmsHSX',SELECT:='Amsterdam heterosexual']
ind[ST=='nonB',ST:='Non-B']
ind[allsupp=='1',allsupp:='Yes']
ind[allsupp=='0',allsupp:='No']

write.csv(ind,file=paste0(outfile.base,'-','number_of_ind_subgraphs_table.csv'))
saveRDS(ind,file=paste0(outfile.base,'-','number_of_ind_subgraphs_table.RDS'))

cat(" \n -------------------------------- table of numbers for text -------------------------------- \n")

summ <- data.table(
	N_diagnoses_ddate=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y]))),
	N_diagnoses_ddate_MSMHSX=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$TRANSM %in% c('MSM','HSX')]))),
	N_diagnoses_ddate_IDU=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$TRANSM %in% c('IDU')]))),
	N_diagnoses_ddate_OTH=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$TRANSM %in% c('Other')]))),
	N_diagnoses_ddate_UKN=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$TRANSM %in% c('Unknown')]))),
	N_diagnoses_ddate_MSM=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y ]))),
	N_diagnoses_ddate_MSM_NL=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='NL']))),
	N_diagnoses_ddate_MSM_G1=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G1']))),
	N_diagnoses_ddate_MSM_G2=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G2']))),
	N_diagnoses_ddate_MSM_G3=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G3']))),
	N_diagnoses_ddate_MSM_Other=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='Other']))),
	N_diagnoses_ddate_HSX=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='HSX' & dind_diag$HIV1_POS_D>=args_dir$start_y ]))),
	N_diagnoses_ddate_HSX_NL=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='HSX' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='NL']))),
	N_diagnoses_ddate_HSX_G4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='HSX' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G4']))),
	N_diagnoses_ddate_HSX_G5=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='HSX' & dind_diag$HIV1_POS_D>=args_dir$start_y  & dind_diag$mwmb=='G5']))),
	N_diagnoses_ddate_HSX_Other=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$TRANSM=='HSX' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='Other']))),
	N_diagnoses_ddate_2009_2013=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=2009 & dind_diag$HIV1_POS_D<2014 & dind_diag$TRANSM %in% c('MSM','HSX')]))),
	N_diagnoses_ddate_2009_2013_late=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=2009 & dind_diag$HIV1_POS_D<2014 & dind_diag$TRANSM %in% c('MSM','HSX') & dind_diag$CD4_first_V<350]))),
	N_diagnoses_infdate_2009_2013=length(unique(na.omit(dind$PATIENT[dind$CITY=='Amsterdam' & dind$INF_D>=2009 & dind$INF_D<2014 & dind$TRANSM %in% c('MSM','HSX')]))),
	N_diagnoses=length(unique(na.omit(ams$PATIENT))),
	N_diagnoses_NL=length(unique(na.omit(nl$PATIENT))),
	rate_diagnoses_Ams=round(length(unique(na.omit(dind_diag$PATIENT[dind_diag$HIV1_POS_D>=args$start_d & dind_diag$CITY=='Amsterdam'])))/854047*100000,d=0),
	rate_diagnoses_NL=round(length(unique(na.omit(dind_diag$PATIENT[dind_diag$HIV1_POS_D>=args$start_d])))/17181084*100000,d=0),
	N_sequences_phylogeny_pre2014=length(unique(na.omit(dind$PATIENT[dind$SEQ==T & dind$CITY=='Amsterdam' & dind$INF_D<args$start_d & dind$INF_D<args$end_d  & !is.na(dind$ST)]))),
	N_sequences_restNL=length(unique(na.omit(dind$PATIENT[dind$SEQ==T & dind$CITY=='Non-Amsterdam']))),
	N_sequences=length(unique(na.omit(ams$PATIENT[ams$SEQ==T]))),
	N_commonsubtypes=length(unique(na.omit(ams$PATIENT[ams$SEQ==T & !is.na(ams$ST)]))),
	N_raresubtypes=length(unique(na.omit(ams$PATIENT[ams$SEQ==T & is.na(ams$ST)]))),
	N_diagnoses_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM']))),
	#N_diagnoses_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CD4_first_V<350]))),
	N_sequences_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T]))),
	N_sequences_MSM_incl=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST)]))),
	N_diagnoses_NL_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN=='NL']))),
	N_diagnoses_nonNL_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL']))),
	N_diagnoses_G1_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$mwmb=='G1']))),
	N_diagnoses_G2_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$mwmb=='G2']))),
	N_diagnoses_G3_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$mwmb=='G3']))),
	N_diagnoses_Other_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$mwmb=='Other']))),
	N_diagnoses_MSM_cd4   =length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$CD4_first_V<350]))),
	N_diagnoses_NL_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$ORIGIN=='NL' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_nonNL_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$ORIGIN!='NL' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_G1_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G1' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_G2_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G2' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_G3_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G3' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_Other_MSM_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='MSM' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='Other' & dind_diag$CD4_first_V<350]))),
	N_sequences_NL_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL']))),
	N_sequences_nonNL_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL']))),
	N_sequences_G1_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='G1']))),
	N_sequences_G2_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='G2']))),
	N_sequences_G3_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='G3']))),
	N_sequences_Other_MSM=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='Other']))),
	N_diagnoses_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX']))),
	#N_diagnoses_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CD4_first_V<350]))),
	N_sequences_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T]))),
	N_sequences_HSX_incl=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST)]))),
	N_diagnoses_NL_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN=='NL']))),
	N_diagnoses_nonNL_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN!='NL']))),
	N_diagnoses_G4_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G4']))),
	N_diagnoses_G5_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G5']))),
	N_diagnoses_Other_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='Other']))),
	N_diagnoses_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$CD4_first_V<350]))),
	N_diagnoses_NL_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$ORIGIN=='NL' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_nonNL_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN!='NL' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$CD4_first_V<350]))),
	N_diagnoses_G4_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G4' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_G5_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='G5' & dind_diag$CD4_first_V<350]))),
	N_diagnoses_Other_HSX_cd4=length(unique(na.omit(dind_diag$PATIENT[dind_diag$TRANSM=='HSX' & dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$mwmb=='Other' & dind_diag$CD4_first_V<350]))),
	N_sequences_G1_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='G4']))),
	N_sequences_G2_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='G5']))),
	N_sequences_Other_HSX=length(unique(na.omit(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$mwmb=='Other']))),
	N_study_p1=length(unique(na.omit(ams$PATIENT[ams$SEQ==T & !is.na(ams$ST) & ams$HIV1_POS_D>=2014]))),
	MSM_chains=length(unique(na.omit(chains$SG[chains$TRANSM=="MSM"]))),
	HSX_chains=length(unique(na.omit(chains$SG[chains$TRANSM=="HSX"]))),
	MSM_N=length(unique(na.omit(chains$ID[chains$TRANSM=="MSM"]))),
	HSX_N=length(unique(na.omit(chains$ID[chains$TRANSM=="HSX"]))),
	MSM_c_pre_p1=length(unique(na.omit(chains$SG[chains$TRANSM=="MSM" & chains$period_p1=="pre-startd"]))),
	HSX_c_pre_p1=length(unique(na.omit(chains$SG[chains$TRANSM=="HSX" & chains$period_p1=="pre-startd"]))),
	MSM_c_pre_p1_newcase=length(unique(na.omit(chains$SG[chains$TRANSM=="MSM" & chains$period_p1=="pre-startd" & chains$HIV1_POS_D>=2014]))),
	HSX_c_pre_p1_newcase=length(unique(na.omit(chains$SG[chains$TRANSM=="HSX" & chains$period_p1=="pre-startd" & chains$HIV1_POS_D>=2014]))),
	MSM_c_em_p1=length(unique(na.omit(chains$SG[chains$TRANSM=="MSM" & chains$MIND>=2014]))),
	HSX_c_em_p1=length(unique(na.omit(chains$SG[chains$TRANSM=="HSX" & chains$MIND>=2014]))),
	MSM_pre_all_supp_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$icases_p1==0 & freqs_pre_i_p1$SELECT=='AmsMSM']),
	MSM_pre_total_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$SELECT=='AmsMSM']),
	MSM_pre_all_supp_p1_new=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1==0 & freqs_pre_p1$jcases_p1>0 & freqs_pre_p1$SELECT=='AmsMSM']),
	MSM_pre_all_supp_p1_nonew=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1==0 & freqs_pre_p1$jcases_p1==0 & freqs_pre_p1$SELECT=='AmsMSM']),
	MSM_pre_someunsupp_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$icases_p1!=0 & freqs_pre_i_p1$SELECT=='AmsMSM']),
	MSM_pre_someunsupp_p1_new=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1!=0 & freqs_pre_p1$jcases_p1!=0 & freqs_pre_p1$SELECT=='AmsMSM']),
	HSX_pre_all_supp_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$icases_p1==0 & freqs_pre_i_p1$SELECT=='AmsHSX']),
	HSX_pre_total_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$SELECT=='AmsHSX']),
	HSX_pre_all_supp_p1_new=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1==0 & freqs_pre_p1$jcases_p1>0 & freqs_pre_p1$SELECT=='AmsHSX']),
	HSX_pre_all_supp_p1_nonew=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1==0 & freqs_pre_p1$jcases_p1==0 & freqs_pre_p1$SELECT=='AmsHSX']),
	HSX_pre_someunsupp_p1=sum(freqs_pre_i_p1$subgraphs[freqs_pre_i_p1$icases_p1!=0 & freqs_pre_i_p1$SELECT=='AmsHSX']),
	HSX_pre_someunsupp_p1_new=sum(freqs_pre_p1$N[freqs_pre_p1$icases_p1!=0 & freqs_pre_p1$jcases_p1!=0 & freqs_pre_p1$SELECT=='AmsHSX'])
)
summ <- reshape2::melt(summ)

write.csv(summ,file=paste0(outfile.base,'-','growth_of_subgraphs.csv'))
saveRDS(summ,file=paste0(outfile.base,'-','growth_of_subgraphs.RDS'))

sizes[,total:=cases_p1]
N_cases <- sizes[order(group,-total)]
#N_cases[, SELECT:=NULL]
N_cases <- N_cases[, c('group','ST','cases_p1')]
setnames(N_cases,c('ST','group','cases_p1'),c('Subtype','Transmission group','Cases (2014-2019)'))

write.csv(N_cases,file=paste0(outfile.base,'-','cases_in_subgraphs.csv'))
saveRDS(N_cases,file=paste0(outfile.base,'-','cases_in_subgraphs.RDS'))

cat(" \n -------------------------------- table of seq/diagnosed/infected -------------------------------- \n")

summ <- data.table(
	N_diagnoses_MSM_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM'])),
	N_sequences_MSM_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_MSM_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_MSM_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN!='NL']))),
	N_diagnoses_HSX_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX'])),
	N_sequences_HSX_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_HSX_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_HSX_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN!='NL'])))
)
#summ[,N_inf_MSM_NL:=round(N_diagnoses_NL_MSM/(1-0.16),digits=0)]
#summ[,N_inf_MSM_nonNL:=round(N_diagnoses_nonNL_MSM/(1-0.18),digits=0)]
#summ[,N_inf_HSX_NL:=round(N_diagnoses_NL_HSX/(1-0.29),digits=0)]
#summ[,N_inf_HSX_nonNL:=round(N_diagnoses_nonNL_HSX/(1-0.34),digits=0)]
summ[,N_inf_MSM_NL:=round(N_diagnoses_NL_MSM/(1-0.17),digits=0)]
summ[,N_inf_MSM_nonNL:=round(N_diagnoses_nonNL_MSM/(1-0.22),digits=0)]
summ[,N_inf_HSX_NL:=round(N_diagnoses_NL_HSX/(1-0.34),digits=0)]
summ[,N_inf_HSX_nonNL:=round(N_diagnoses_nonNL_HSX/(1-0.47),digits=0)]
summ[,N_inf_MSM_NL_ALL:=N_inf_MSM_NL+N_inf_MSM_nonNL]
summ[,N_inf_HSX_NL_ALL:=N_inf_HSX_NL+N_inf_HSX_nonNL]

summ[,N_undiag_MSM_NL:=N_inf_MSM_NL-N_diagnoses_NL_MSM]
summ[,N_undiag_MSM_nonNL:=N_inf_MSM_nonNL-N_diagnoses_nonNL_MSM]
summ[,N_undiag_HSX_NL:=N_inf_HSX_NL-N_diagnoses_NL_HSX]
summ[,N_undiag_HSX_nonNL:=N_inf_HSX_nonNL-N_diagnoses_nonNL_HSX]

summ <- data.table(reshape2::melt(summ))

cat(" \n -------------------------------- plot -------------------------------- \n")

summ[,trsm:='MSM']
summ[grepl('HSX',variable),trsm:='HSX']
summ$trsm <- factor(summ$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

summ[,bplace:='Dutch-born']
summ[grepl('nonNL',variable),bplace:='Foreign-born']
summ[grepl('ALL',variable),bplace:='All']
summ$bplace <- factor(summ$bplace,levels=c('All','Dutch-born','Foreign-born'))

summ[,var:='Sequenced']
summ[grepl('diagnoses',variable),var:='Diagnosed']
summ[grepl('inf',variable),var:='Infected']
summ[grepl('undiag',variable),var:='Undiagnosed']
summ[grepl('diagd',variable),var:='Diagnosed (diag date)']
summ$var <- factor(summ$var,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'))

seqs <- subset(summ,var=='Sequenced')
setnames(seqs,'value','seqs')
summ <- merge(summ,subset(seqs,select=c('trsm','bplace','seqs')),by=c('trsm','bplace'),all.x=T)
summ[, pct:=round(seqs/value*100,digits=0)]
summ[var=='Sequenced',pct:='']
summ[, lab:=paste0(value,'\n(',pct,'%)')]
summ[var=='Sequenced',lab:=value]

plot <- ggplot(data=subset(summ,var!='Undiagnosed' & var!='Diagnosed (diag date)'),aes(x=bplace,y=value,fill=var)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_text(aes(label=lab), position=position_dodge(width=0.9), vjust=-0.25) +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~.,scales="free_y") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-seq_diag_inf_2014-2019.png'),plot,w=8, h=8)

#summ[trsm=='Amsterdam MSM' & bplace=='Dutch-born' & var=='Infected',lab_u:='16%']
#summ[trsm=='Amsterdam MSM' & bplace=='Foreign-born' & var=='Infected',lab_u:='18%']
#summ[trsm=='Amsterdam heterosexuals' & bplace=='Dutch-born' & var=='Infected',lab_u:='29%']
#summ[trsm=='Amsterdam heterosexuals' & bplace=='Foreign-born' & var=='Infected',lab_u:='34%']
summ[trsm=='Amsterdam MSM' & bplace=='Dutch-born' & var=='Infected',lab_u:='17%']
summ[trsm=='Amsterdam MSM' & bplace=='Foreign-born' & var=='Infected',lab_u:='22%']
summ[trsm=='Amsterdam heterosexuals' & bplace=='Dutch-born' & var=='Infected',lab_u:='34%']
summ[trsm=='Amsterdam heterosexuals' & bplace=='Foreign-born' & var=='Infected',lab_u:='47%']

diag <- subset(summ,var=='Diagnosed')
diag[,var2:='Infected']
summ[, var2:=var]
summ <- rbind(summ,diag)
summ$varlab <- factor(summ$var,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'),
																labels=c('Sequenced','Diagnosed \nin 2014-May 2019',
																				 'Estimated \ndiagnosed \nand infected \nin 2014-2018',
																				 'Undiagnosed','Estimated \nundiagnosed, \ndiagnosed and \ninfected \nin 2014-2018'))
summ$varlab2 <- factor(summ$var2,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'),
											labels=c('Sequenced','Diagnosed \nin 2014-May 2019',
															 'Estimated \ndiagnosed \nand infected \nin 2014-2018',
															 'Undiagnosed','Estimated \nundiagnosed, \ndiagnosed and \ninfected \nin 2014-2018'))
col_pal = rev(pal_npg("nrc")(3))

plot <- ggplot(data=subset(summ,var!='Undiagnosed' & var!='Sequenced' & bplace!='All'),aes(x=varlab,y=value,fill=var)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_bar(data=subset(summ,var=='Diagnosed' & var2=='Infected' & bplace!='All'), aes(x=varlab2,y=value),fill=col_pal[2],stat='identity', position = "identity") +
	geom_text(aes(label=lab_u), position=position_dodge(width=0.9), vjust=1.1,size=6,colour="black") +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~bplace,scales="free_y") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	#ggsci::scale_fill_npg()
	scale_fill_manual(values=col_pal)
ggsave(file=paste0(outfile.base,'-diag_inf_2014-2019_oldpct_cols.png'),plot,w=14, h=8)


cat(" \n -------------------------------- table of seq/diagnosed/infected by migrant group -------------------------------- \n")

summ <- data.table(
	N_diagnoses_MSM_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM'])),
	N_sequences_MSM_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_MSM_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_MSM_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_NL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$ORIGIN!='NL']))),
	N_diagnoses_G1_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G1']))),
	N_diagnoses_G2_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G2']))),
	N_diagnoses_G3_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='G3']))),
	N_diagnoses_Oth_MSM_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='MSM' & dind_diag$mwmb=='Other']))),
	N_diagnoses_NL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL'])),
	N_diagnoses_G1_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G1'])),
	N_diagnoses_G2_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G2'])),
	N_diagnoses_G3_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='G3'])),
	N_diagnoses_Oth_MSM=length(unique(ams$PATIENT[ams$TRANSM=='MSM' & ams$ORIGIN!='NL' & ams$mwmb=='Other'])),
	N_diagnoses_HSX_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX'])),
	N_sequences_HSX_incl_ALL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST)])),
	N_sequences_HSX_incl_NL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN=='NL'])),
	N_sequences_HSX_incl_nonNL=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$SEQ==T & !is.na(ams$ST) & ams$ORIGIN!='NL'])),
	N_diagnoses_NL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN=='NL'])),
	N_diagnoses_nonNL_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$ORIGIN!='NL'])),
	N_diagnoses_G4_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G4'])),
	N_diagnoses_G5_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='G5'])),
	N_diagnoses_Oth_HSX=length(unique(ams$PATIENT[ams$TRANSM=='HSX' & ams$mwmb=='Other'])),
	N_diagnoses_NL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN=='NL']))),
	N_diagnoses_nonNL_HSX_diagd=length(unique(na.omit(dind_diag$PATIENT[dind_diag$CITY=='Amsterdam' & dind_diag$HIV1_POS_D>=args_dir$start_y & dind_diag$HIV1_POS_D<args$end_d & dind_diag$TRANSM=='HSX' & dind_diag$ORIGIN!='NL'])))
)
summ[,N_inf_MSM_NL:=round(N_diagnoses_NL_MSM/(1-0.2258),digits=0)]
summ[,N_inf_MSM_G1:=round(N_diagnoses_G1_MSM/(1-0.2102),digits=0)]
summ[,N_inf_MSM_G2:=round(N_diagnoses_G2_MSM/(1-0.2762),digits=0)]
summ[,N_inf_MSM_G3:=round(N_diagnoses_G3_MSM/(1-0.2833),digits=0)]
summ[,N_inf_MSM_Oth:=round(N_diagnoses_Oth_MSM/(1-0.3152),digits=0)]
summ[,N_inf_HSX_NL:=round(N_diagnoses_NL_HSX/(1-0.3570),digits=0)]
summ[,N_inf_HSX_G4:=round(N_diagnoses_G4_HSX/(1-0.6291),digits=0)]
summ[,N_inf_HSX_G5:=round(N_diagnoses_G5_HSX/(1-0.3430),digits=0)]
summ[,N_inf_HSX_Oth:=round(N_diagnoses_Oth_HSX/(1-0.4628),digits=0)]
summ[,N_inf_MSM_nonNL:=N_inf_MSM_G1+N_inf_MSM_G2+N_inf_MSM_G3+N_inf_MSM_Oth]
summ[,N_inf_MSM_NL_ALL:=N_inf_MSM_NL+N_inf_MSM_nonNL]
summ[,N_inf_HSX_nonNL:=N_inf_HSX_G4+N_inf_HSX_G5+N_inf_HSX_Oth]
summ[,N_inf_HSX_NL_ALL:=N_inf_HSX_NL+N_inf_HSX_nonNL]

summ[,N_undiag_MSM_NL:=N_inf_MSM_NL-N_diagnoses_NL_MSM]
summ[,N_undiag_MSM_nonNL:=N_inf_MSM_nonNL-N_diagnoses_nonNL_MSM]
summ[,N_undiag_HSX_NL:=N_inf_HSX_NL-N_diagnoses_NL_HSX]
summ[,N_undiag_HSX_nonNL:=N_inf_HSX_nonNL-N_diagnoses_nonNL_HSX]

summ <- data.table(reshape2::melt(summ))

cat(" \n -------------------------------- plot -------------------------------- \n")

summ[,trsm:='MSM']
summ[grepl('HSX',variable),trsm:='HSX']
summ$trsm <- factor(summ$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))

summ[,bplace:='Dutch-born']
summ[grepl('nonNL',variable),bplace:='Foreign-born']
summ[grepl('ALL',variable),bplace:='All']
summ[grepl('G1',variable),bplace:='G1']
summ[grepl('G2',variable),bplace:='G2']
summ[grepl('G3',variable),bplace:='G3']
summ[grepl('G4',variable),bplace:='G4']
summ[grepl('G5',variable),bplace:='G5']
summ[grepl('Oth',variable),bplace:='Other']
summ <- subset(summ,bplace %in% c('All','Dutch-born','Foreign-born'))
summ$bplace <- factor(summ$bplace,levels=c('All','Dutch-born','Foreign-born'))

summ[,var:='Sequenced']
summ[grepl('diagnoses',variable),var:='Diagnosed']
summ[grepl('inf',variable),var:='Infected']
summ[grepl('undiag',variable),var:='Undiagnosed']
summ[grepl('diagd',variable),var:='Diagnosed (diag date)']
summ$var <- factor(summ$var,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'))

seqs <- subset(summ,var=='Sequenced')
setnames(seqs,'value','seqs')
summ <- merge(summ,subset(seqs,select=c('trsm','bplace','seqs')),by=c('trsm','bplace'),all.x=T)
summ[, pct:=round(seqs/value*100,digits=0)]
summ[var=='Sequenced',pct:='']
summ[, lab:=paste0(value,'\n(',pct,'%)')]
summ[var=='Sequenced',lab:=value]

plot <- ggplot(data=subset(summ,var!='Undiagnosed' & var!='Diagnosed (diag date)'),aes(x=bplace,y=value,fill=var)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_text(aes(label=lab), position=position_dodge(width=0.9), vjust=-0.25) +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~.,scales="free_y") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank()) +
	ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-seq_diag_inf_2014-2019.png'),plot,w=8, h=8)

#summ[trsm=='Amsterdam MSM' & bplace=='Dutch-born' & var=='Infected',lab_u:='16%']
#summ[trsm=='Amsterdam MSM' & bplace=='Foreign-born' & var=='Infected',lab_u:='18%']
#summ[trsm=='Amsterdam heterosexuals' & bplace=='Dutch-born' & var=='Infected',lab_u:='29%']
#summ[trsm=='Amsterdam heterosexuals' & bplace=='Foreign-born' & var=='Infected',lab_u:='34%']
summ[trsm=='Amsterdam MSM' & bplace=='Dutch-born' & var=='Infected',lab_u:='13%']
summ[trsm=='Amsterdam MSM' & bplace=='Foreign-born' & var=='Infected',lab_u:='18%']
summ[trsm=='Amsterdam heterosexuals' & bplace=='Dutch-born' & var=='Infected',lab_u:='27%']
summ[trsm=='Amsterdam heterosexuals' & bplace=='Foreign-born' & var=='Infected',lab_u:='42%']

diag <- subset(summ,var=='Diagnosed')
diag[,var2:='Infected']
summ[, var2:=var]
summ <- rbind(summ,diag)
summ$varlab <- factor(summ$var,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'),
											labels=c('Sequenced','Diagnosed \nin 2014-2019',
															 'Estimated \ndiagnosed \nand infected \nin 2014-2019',
															 'Undiagnosed','Estimated \nundiagnosed, \ndiagnosed and \ninfected \nin 2014-2019'))
summ$varlab2 <- factor(summ$var2,levels=c('Sequenced','Diagnosed (diag date)','Diagnosed','Undiagnosed','Infected'),
											 labels=c('Sequenced','Diagnosed \nin 2014-2019',
											 				 'Estimated \ndiagnosed \nand infected \nin 2014-2019',
											 				 'Undiagnosed','Estimated \nundiagnosed, \ndiagnosed and \ninfected \nin 2014-2019'))
col_pal = rev(pal_npg("nrc")(3))

plot <- ggplot(data=subset(summ,var!='Undiagnosed' & var!='Sequenced' & bplace!='All'),aes(x=varlab,y=value,fill=var)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_bar(data=subset(summ,var=='Diagnosed' & var2=='Infected' & bplace!='All'), aes(x=varlab2,y=value),fill=col_pal[2],stat='identity', position = "identity") +
	geom_text(aes(label=lab_u), position=position_dodge(width=0.9), vjust=1.1,size=6,colour="black") +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(trsm~bplace,scales="free_y") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="none",
				strip.background=element_blank()) +
	#ggsci::scale_fill_npg()
	scale_fill_manual(values=col_pal)
ggsave(file=paste0(outfile.base,'-diag_inf_2014-2019_oldpct_cols.png'),plot,w=14, h=8)

cat(" \n -------------------------------- summarise time to diagnosis -------------------------------- \n")

load(infile.seq)
load(infile.meta)
dind <- data.table(dind)
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
dinf <- read.csv(infile.inftime,header=T)
dinf$SEQ <- dinf$id %in% ds$PATIENT
dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
dinf <- unique(dinf)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D')),by.x='id',by.y='PATIENT',all.x=T)
#setnames(dinf,'hiv_pos_d','HIV1_POS_D')
do <- data.table(dinf)
do[, time:=estsctodiagMedian]
do[, INF_D:=HIV1_POS_D - time]

## HSX
dt <- subset(do,TRANSM=='HSX' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='HSX' & do$BIRTH_CNTRY=='Netherlands' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='HSX' & do$BIRTH_CNTRY!='Netherlands' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)


## MSM
dt <- subset(do,TRANSM=='MSM' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='MSM' & do$BIRTH_CNTRY=='Netherlands' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='MSM' & do$BIRTH_CNTRY!='Netherlands' & do$INF_D>=2010 & do$INF_D<2013)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)


## HSX
dt <- subset(do,TRANSM=='HSX' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='HSX' & do$BIRTH_CNTRY=='Netherlands' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='HSX' & do$BIRTH_CNTRY!='Netherlands' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)


## MSM
dt <- subset(do,TRANSM=='MSM' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='MSM' & do$BIRTH_CNTRY=='Netherlands' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)

dt <- subset(do,TRANSM=='MSM' & do$BIRTH_CNTRY!='Netherlands' & do$INF_D>=2014)
quantile(dt$time,probs=c(0.025,0.5,0.975),na.rm=T)



############

dt <- ams[, list(diag=length(unique(na.omit(PATIENT))),seq=length(unique(na.omit(PATIENT[SEQ==T & !is.na(ST)])))),by=c('TRANSM','mwmb')]
dt <- subset(dt,TRANSM %in% c('MSM','HSX'))
#dt <- data.table(reshape2::melt(dt))
dt$trsm <- factor(dt$TRANSM,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
#dt$mwmb <- factor(dt$mwmb)
dt[mwmb=='G1' & TRANSM=='MSM', mlab:='W.Europe,\nN.America,\nOceania']
dt[mwmb=='G2' & TRANSM=='MSM', mlab:='E. & C. Europe']
dt[mwmb=='G3' & TRANSM=='MSM', mlab:='S. America &\n Caribbean']
dt[mwmb=='NL' & TRANSM=='MSM', mlab:='NL']
dt[mwmb=='Other' & TRANSM=='MSM', mlab:='Other']

dt[mwmb=='G4' & TRANSM=='HSX', mlab:='Sub-Saharan\nAfrica']
dt[mwmb=='G5' & TRANSM=='HSX', mlab:='S. America &\n Caribbean']
dt[mwmb=='NL' & TRANSM=='HSX', mlab:='NL']
dt[mwmb=='Other' & TRANSM=='HSX', mlab:='Other']
dt$mlab <- factor(dt$mlab,levels=c('NL','S. America &\n Caribbean','W.Europe,\nN.America,\nOceania','E. & C. Europe','Sub-Saharan\nAfrica','Other'))

dt[, pct:=round(seq/diag*100,digits=0)]
dt[, lab:=paste0(pct,'%')]

col_pal = pal_npg("nrc")(2)

plot <- ggplot(data=subset(dt)) +
	#geom_bar(stat='identity', position = "dodge") +
	geom_bar(aes(x=mlab,y=diag, fill=col_pal[2]), stat="identity", position ="identity", alpha=1) +
	geom_bar(aes(x=mlab,y=seq, fill=col_pal[1]), stat="identity", position="identity", alpha=1) +
	geom_text(aes(x=mlab,y=diag,label=lab), position=position_dodge(width=0.9), vjust=-0.75,size=8) +
	#geom_text(aes(label=lab), position=position_dodge(width=0.9), vjust=-0.25) +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(.~trsm,scales="free") +
	scale_y_continuous(expand = expansion(mult = c(0, .2)))  +
	labs(x='',y="Number of new infections \n 2014-2019",fill="") +
	theme_bw(base_size=18) +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_text(size=24)) +
	scale_fill_manual(values = c( col_pal[2], col_pal[1]),
										 labels = c("Diagnosed", "Sequenced")) +
	guides(fill=guide_legend(override.aes=list(fill=c("Diagnosed"=col_pal[2],"Sequenced"=col_pal[1]))))
	#ggsci::scale_fill_npg()
ggsave(file=paste0(outfile.base,'-seq_diag_2014-2019.png'),plot,w=15, h=8)


cat(" \n -------------------------------- summary table of characteristics seq/no seq -------------------------------- \n")

ams <- subset(ams, TRANSM %in% c('MSM','HSX'))
ams[is.na(ST),SEQ:=F]

ams[, AGE:=INF_D - BIRTH_Y]
set(ams,NULL,'AGE_GROUP',findInterval(ams$AGE,c(14,25,35,45,60)))
ams$AGE_GROUP <- factor(ams$AGE_GROUP,levels=c(1,2,3,4,5),labels=c('18-24','25-34','35-44','45-59','60+'))

tab <- data.table(var='SEX',table(ams$GENDER,ams$SEQ))

tab <- rbind(tab,data.table(var='TRANSM',table(ams$TRANSM,ams$SEQ)))

tab <- rbind(tab,data.table(var='AGE_GROUP',table(ams$AGE_GROUP,ams$SEQ)))

tab <- rbind(tab,data.table(var='BIRTH_PLACE',table(droplevels(ams$WRLD_born),ams$SEQ)))

time <- ams[, list(var='INFTIME',
									 V1='Estimated time to diagnosis (years)',
									 q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
									 				 " [",
									 				 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
									 				 "-",
									 				 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
									 				 "]")),
									 by=c('SEQ')]
time <- time[, c('var','V1','SEQ','q')]
setnames(time,c('SEQ','q'),c('V2','N'))

tab <- rbind(tab,time)

tab.seq <- dcast(tab,var+V1~V2,value.var='N')

tab.seq$var <- factor(tab.seq$var,levels=c('SEX','TRANSM','AGE_GROUP','BIRTH_PLACE','INFTIME'))
tab.seq <- tab.seq[order(var),]

#############

tab <- data.table(var='SEX',table(ams$GENDER))
tab <- rbind(tab,data.table(var='TRANSM',table(ams$TRANSM)))
tab <- rbind(tab,data.table(var='AGE_GROUP',table(ams$AGE_GROUP)))
tab <- rbind(tab,data.table(var='BIRTH_PLACE',table(droplevels(ams$WRLD_born))))
time <- ams[, list(var='INFTIME',
									 V1='Estimated time to diagnosis (years)',
									 q=paste0(round(quantile(estsctodiagMedian,probs=0.5,na.rm=T),2),
									 				 " [",
									 				 round(quantile(estsctodiagMedian,probs=0.025,na.rm=T),2),
									 				 "-",
									 				 round(quantile(estsctodiagMedian,probs=0.975,na.rm=T),2),
									 				 "]"))]
time <- time[, c('var','V1','q')]
setnames(time,c('q'),c('N'))

tab <- rbind(tab,time)
setnames(tab,'N','All patients')

tab <- merge(tab,subset(tab.seq,select=c('var','V1','TRUE')),by=c('var','V1'),all=T)
setnames(tab,'TRUE','seq')

tab2 <- subset(tab,var!='INFTIME')
tab2$`All patients` <- as.numeric(tab2$`All patients`)
tab2$seq <- as.numeric(tab2$seq)
tab3 <- tab2[, list(V1=V1,
										all=`All patients`,
									all_pct=`All patients`/sum(`All patients`),
									with_seq=seq,
									seq_pct=seq/sum(seq)),
					 by=c('var')]
tab3 <- tab3[,list(var=var,V1=V1,
									 all=paste0(all," (",round(all_pct*100,1),"%)"),
									 seq=paste0(with_seq," (",round(seq_pct*100,1),"%)"))]
inftime <- subset(tab,var=='INFTIME')
setnames(inftime,'All patients','all')
tab3 <- rbind(tab3,inftime)

tab3$var <- factor(tab3$var,levels=c('SEX','TRANSM','AGE_GROUP','BIRTH_PLACE','INFTIME'))
tab3 <- tab3[order(var),]

tab3[V1=='HSX', V1:='Heterosexuals']
tab3[V1=='Africa', V1:='Sub-Saharan Africa']
tab3[V1=='CEurope', V1:='Central Europe']
tab3[V1=='EEurope', V1:='Eastern Europe']
tab3[V1=='FormerCurrDutchColonies', V1:='Suriname, Curacao & Aruba']
tab3[V1=='LaAmCar', V1:='South America & Caribbean']
tab3[V1=='MENA', V1:='Middle East & North Africa']
tab3[V1=='NorthAm', V1:='North America']
tab3[V1=='WEurope', V1:='Western Europe']

write.csv(tab3,file=paste0(outfile.base,'-','table_characteristics.csv'))
saveRDS(tab3,file=paste0(outfile.base,'-','table_characteristics.RDS'))





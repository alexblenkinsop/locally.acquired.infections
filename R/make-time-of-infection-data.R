
library(tidyr)
library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)

args <- list( 
	source_dir= '/rds/general/user/ablenkin/home/git/bpm',
	indir='/rds/general/project/ratmann_roadmap_data_analysis/live',
	outdir= '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model',
	analysis= 'analysis_200917'
)
`%notin%` <- Negate(`%in%`)

home <- '/rds/general/project/ratmann_roadmap_data_analysis/live'

country_continent <- file.path(home,'misc','country_continent.csv')
infile.meta <- file.path(args$indir, args$analysis, 'misc', '200917_sequence_labels.rda')
infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')

load(infile.meta)
setnames(dind, gsub("CENS_D", "RECART_D", names(dind)))
dind <- subset(dind,CITY=='Amsterdam')

# just keep the variables we need
dt <- data.table(subset(dind,select=c('PATIENT', 'GENDER','ORIGIN','BIRTH_CNTRY','TRANSM')))


dbas <- read.csv(infile.indinfo,header=T)
dbas <- subset(dbas,select=c(PATIENT,BIRTH_D,RECART_D,HIV1_POS_D,HIV1_NEG_D,AIDS_D))
dbas[,'HIV1_POS_D'] <- as.Date(dbas[,'HIV1_POS_D'],format="%d/%m/%Y")
dbas[,'HIV1_NEG_D'] <- as.Date(dbas[,'HIV1_NEG_D'],format="%d/%m/%Y")
dbas[,'BIRTH_D'] <- as.Date(dbas[,'BIRTH_D'],format="%Y-%m-%d")
dbas[,'RECART_D'] <- as.Date(dbas[,'RECART_D'],format="%Y-%m-%d")
dbas[,'AIDS_D'] <- as.Date(dbas[,'AIDS_D'],format="%Y-%m-%d")

dt <- merge(dt,dbas,by=c('PATIENT'),all=T)

# generate extra variables
# date of risk onset: since age 15 or since last negative test - 30 days (whichever later)
dt[, fifteen:= BIRTH_D + (365.25*15) ]

dt[, atRiskDate:= apply(subset(dt,select=c('HIV1_NEG_D','fifteen')), 1, max, na.rm=T)]
dt[, atRiskDate:=as.Date(atRiskDate,format="%Y-%m-%d")]
# if person not yet 15, use last hiv neg test
dt[fifteen>as.Date(Sys.Date()),atRiskDate:=HIV1_NEG_D]

dt[, u:= lubridate::time_length(difftime(HIV1_POS_D, atRiskDate), "years")]
# seroconversion date
dt[, SEROCO_D:=NA]
# date first aids symptom
dt[, DIS_D:=NA]
# age at diagnosis
dt[, ageDiag:= lubridate::time_length(difftime(HIV1_POS_D, BIRTH_D), "years")]
# time from 1980-01-01 to diagnosis date
dt[, calendar:= lubridate::time_length(difftime(HIV1_POS_D, "1980-01-01"), "years")]


# generate region_group from georegs
dgeo <- read.csv(country_continent,header=T)

dt <- merge(dt,dgeo,by.x='ORIGIN',by.y='iso.3166.country',all.x=T)

dt[, region_group:='Other/Unknown']
dt[continent.code=='AS', region_group:='Asia']
dt[continent.code=='AF', region_group:='Africa']
dt[continent.code=='EU', region_group:='Europe']

# add cd4/vl data
infile.cd4 <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_CD4.csv')
infile.rna <-	file.path(args$indir, 'Data', 'data_191223_Amsterdam/SHM_1902_ROADMAP_191223_tblLAB_RNA.csv')

dat <- data.table(read.csv(infile.cd4,header=T))
# keep cells/ml (not %)
dat <- subset(dat,CD4_U==1)
dat$CD4_D <- as.Date(dat$CD4_D,format=c("%Y-%m-%d"))
dat[, CD4_V:=sqrt(CD4_V)]
setnames(dat,c('CD4_V','CD4_D'),c('yvar','ExamDate'))
dat[, indi:='cd4']
dat <- subset(dat, select=c('PATIENT','ExamDate','yvar','indi'))
dat[, RNA_L:=NA]
dat[, RNA_UL:=NA]

tmp <- data.table(read.csv(infile.rna,header=T))
tmp$RNA_D <- as.Date(tmp$RNA_D,format=c("%Y-%m-%d"))
tmp <- subset(tmp,RNA_V!=-1)
tmp[, RNA_V:=log10(RNA_V)]
setnames(tmp,c('RNA_V','RNA_D'),c('yvar','ExamDate'))
tmp[, indi:='rna']
tmp <- subset(tmp, select=c('PATIENT','ExamDate','yvar','indi','RNA_L','RNA_UL'))

dat <- rbind(dat,tmp)


# after we remove measurements after AIDS diag/ART start
dat[, consc:=0]
dat[, consr:=0]
dat[indi=='cd4', consc:=1]
dat[indi=='rna', consr:=1]

dat <- merge(dat,dt,by='PATIENT',all.x=T)

# time from diagnosis to current measurement
dat[, dtime:= lubridate::time_length(difftime(ExamDate, HIV1_POS_D), "years")]

# exclude patients with perinatal infections
dat <- subset(dat, ageDiag>=15)

# just keep patients diagnosed since 2010
#dat <- subset(dat,HIV1_POS_D>=as.Date("2010-01-01"))

# calculate time between art start and measurement
dat[, arttime:= lubridate::time_length(difftime(ExamDate, RECART_D), "weeks")]


dat[, keep:= 1]
dat[ExamDate>RECART_D, keep:= 0]
dat[ExamDate>AIDS_D, keep:= 0]
# keep vl measurements within 1 week of ART start 
dat[indi=='vl' & arttime<1, keep:= 1]
# keep cd4 measurements within 1 month of ART start
dat[indi=='cd4' & arttime<4, keep:= 1]

dat <- subset(dat,keep==1)

# replace HIV_NEG_D rounded to 1st July with 1st Jan and recalculate at risk date
dat[,HIV1_NEG_D_L_FLAG:=0]
dat[u<0,HIV1_NEG_D_L_FLAG:=1]

dat[u<0,HIV1_NEG_D:=as.Date(HIV1_NEG_D %m-% months(6))]
dat[u<0, atRiskDate:= HIV1_NEG_D]
dat[u<0, u:= lubridate::time_length(difftime(HIV1_POS_D, atRiskDate), "years")]

mm <- dat[, list(sumc=sum(consc),sumr=sum(consr)),by=c('PATIENT')]
mm[, only:="Both"]
mm[sumc<=0, only:="VL only"]
mm[sumr<=0, only:="CD4 only"]
mm[sumr<=0 & sumc<=0, only:=NA]

da <- merge(dat,subset(mm,select=c('PATIENT','only')),by=c('PATIENT'),all.x=T)
da <- subset(da, !is.na(only))

# deal with pts with no CD4/VL measurements before aids onset/art initiation
# add one record with DIS_D=AIDS_D and yvar=NA
# exclude patients with perinatal infections
dt <- subset(dt, ageDiag>=15)
# just keep patients diagnosed since 2010
#dt <- subset(dt,HIV1_POS_D>=as.Date("2010-01-01"))
dmis <- dt[dt$PATIENT %notin% da$PATIENT,] #1667
dmis <- subset(dmis,AIDS_D<RECART_D) # 707
dmis[,DIS_D:=AIDS_D]

cols <- colnames(dmis)
da[, DIS_D:=as.Date(DIS_D)]
do <- merge(da,dmis,by=c(cols),all=T)

# rename cols
setnames(do,c('PATIENT','TRANSM','HIV1_POS_D','GENDER'),c('id','mode','HIV_POS_D','gender'))

# remove columns not needed
#set(do, NULL, c('ORIGIN','continent.code','BIRTH_CNTRY','RECART_D','HIV1_NEG_D','AIDS_D','fifteen','keep', 'arttime'), NULL)
set(do, NULL, c('ORIGIN','continent.code','BIRTH_CNTRY','fifteen','keep', 'arttime'), NULL)

do <- data.table(do)

setkey(do,id,ExamDate)

do <- subset(do,ExamDate>as.Date('1911-11-11'))

do <- select(do, 'id','u','atRiskDate','yvar','indi','BIRTH_D','HIV1_NEG_D','HIV_POS_D','ExamDate','SEROCO_D','DIS_D','RECART_D','AIDS_D','gender','region_group',
						 'mode','ageDiag','dtime','calendar','consc','consr','only','HIV1_NEG_D_L_FLAG','RNA_L','RNA_UL')

do[is.nan(yvar),yvar:=-1]

# save
write.csv(do,file=file.path(args$indir, args$analysis, 'misc','roadmap_cd4_vl_allpts.csv'))


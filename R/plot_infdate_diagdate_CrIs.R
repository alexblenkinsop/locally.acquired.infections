require(data.table)
require(ggplot2)
require(ggsci)

home <- '/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/branching_process_model'

dir.name <- '/Users/alexb/Box Sync/Roadmap/Data/data_200821'
file			<- file.path(dir.name,"SHM_1902_ROADMAP_200821_tblLAB_seq.rda")
outpath <- '/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript'

dsubgraphtaxa <- readRDS(file.path(home,'subgraphs_withmetadata.RDS'))
do <- subset(dsubgraphtaxa,REP=='000' & SELECT!='Ams')

load(file)
load('/Users/alexb/Box Sync/Roadmap/analysis_200917/misc/200917_sequence_labels.rda')
infile.inftime <- file.path('~/Box Sync/Roadmap/RQ1 Estimating introductions/Data/infection_time_estimates','roadmap_cd4_vl_est.csv')
geo <- data.table(read.csv('/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/analysis_200917/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

dind <- data.table(dind)
dind$SEQ <- dind$PATIENT %in% ds$PATIENT
#infile.inftime <- file.path(args$indir,'Data','infection_time_estimates','roadmap_cd4_vl_est.csv')
dinf <- read.csv(infile.inftime,header=T)
dinf$SEQ <- dinf$id %in% ds$PATIENT
dinf <- subset(dinf,select=c('id','estsctodiagMedian','hiv_pos_d'))
dinf <- unique(dinf)
dinf <- merge(dinf,subset(dind,select=c('PATIENT','TRANSM','BIRTH_CNTRY','HIV1_POS_D')),by.x='id',by.y='PATIENT',all.x=T)
#setnames(dinf,'hiv_pos_d','HIV1_POS_D')
do <- data.table(dinf)
do[, time:=estsctodiagMedian]
## Estimate undiagnosed by migrant groups
do[, INF_D:=HIV1_POS_D - time]
do <- merge(do,subset(dind,select=c('PATIENT','ORIGIN')),by.x='id', by.y='PATIENT',all.x=T)
do <- merge(do,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)

## msm
do[TRANSM=='MSM', mwmb:="Other"]
do[TRANSM=='MSM' & ORIGIN %in% c("NL"), mwmb:="NL"]
# western countires (non-NL)
do[TRANSM=='MSM' & WRLD_born %in% c("WEurope","NorthAm","Oceania") & ORIGIN!='NL', mwmb:="G1"]
# eastern and central europe
do[TRANSM=='MSM' & WRLD_born %in% c("EEurope", "CEurope"), mwmb:="G2"]
# caribbean and south america
do[TRANSM=='MSM' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G3"]

do[mwmb=='G1',mgid:=1]
do[mwmb=='G2',mgid:=2]
do[mwmb=='G3',mgid:=3]
do[mwmb=='NL',mgid:=4]
do[mwmb=='Other',mgid:=5]

# exluce individuals without an infection time estimate
do <- subset(do,!is.na(time))

## hsx
do[TRANSM=='HSX', mwmb:="Other"]
do[TRANSM=='HSX' & ORIGIN %in% c("NL"), mwmb:="NL"]
# sub-saharan africa
do[TRANSM=='HSX' & WRLD_born %in% c("Africa"), mwmb:="G4"]
# caribbean and south america
do[TRANSM=='HSX' & WRLD_born %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]

do[mwmb=='G4',mgid:=1]
do[mwmb=='G5',mgid:=2]
do[mwmb=='NL',mgid:=3]
do[mwmb=='Other',mgid:=4]

do[mwmb=='G1' & TRANSM=='MSM', mlab:='W.Europe,\nN.America,\nOceania']
do[mwmb=='G2' & TRANSM=='MSM', mlab:='E. & C. Europe']
do[mwmb=='G3' & TRANSM=='MSM', mlab:='S. America &\n Caribbean']
do[mwmb=='NL' & TRANSM=='MSM', mlab:='NL']
do[mwmb=='Other' & TRANSM=='MSM', mlab:='Other']

do[mwmb=='G4' & TRANSM=='HSX', mlab:='Sub-Saharan\nAfrica']
do[mwmb=='G5' & TRANSM=='HSX', mlab:='S. America &\n Caribbean']
do[mwmb=='NL' & TRANSM=='HSX', mlab:='NL']
do[mwmb=='Other' & TRANSM=='HSX', mlab:='Other']

do$mlab <- factor(do$mlab,levels=c('NL','W.Europe,\nN.America,\nOceania','E. & C. Europe','S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'))

cat(" \n -------------------------------- load infection times -------------------------------- \n")

home <- '/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed'

dat <- data.table(read.csv('/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/undiagnosed/roadmap_cd4_vl_allpts_est.csv'))
dat[, hiv_pos_d:=as.Date(hiv_pos_d,format="%Y-%m-%d")]

dat[, HIV_POS_D:=hivc.db.Date2numeric(hiv_pos_d)]

dt <- unique(subset(dat,select=c('id','mode','HIV_POS_D','estsctodiagMedian','estsctodiagLL','estsctodiagUL')))

dt[, inf_date:=HIV_POS_D-estsctodiagMedian]
dt[, inf_date_L:=HIV_POS_D-estsctodiagLL]
dt[, inf_date_U:=HIV_POS_D-estsctodiagUL]
dt[is.na(estsctodiagMedian),inf_date:=HIV_POS_D]
length(unique(dt$id[dt$inf_date>=2009 & dt$inf_date<2014]))
length(unique(dt$id[dt$HIV_POS_D>=2009 & dt$HIV_POS_D<2014]))

df <- subset(dt,inf_date>=2014)
df <- merge(df,subset(do,select=c('id','mlab')),by='id',all.x=T)
df$trsm <- factor(df$mode,levels=c('MSM','HSX'),labels=c('Amsterdam MSM', 'Amsterdam heterosexual'))

ggplot(data=subset(df,mode %in% c('MSM','HSX'))) + geom_point(aes(x=HIV_POS_D,y=inf_date)) +
	geom_errorbar(aes(x=HIV_POS_D,ymin=inf_date_L, ymax=inf_date_U),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	scale_y_continuous(expand = c(0,0))  +
	facet_grid(.~trsm) +
	theme_bw() +
	labs(x='Date of diagnosis',y="Estimated date of infection") +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				strip.text = element_blank()) +
	scale_fill_npg()

g <- ggplot(data=subset(df,mode %in% c('MSM','HSX') & !is.na(mlab))) + geom_point(aes(y=HIV_POS_D,x=inf_date,colour=mlab)) +
	geom_errorbar(aes(y=HIV_POS_D,xmin=inf_date_L, xmax=inf_date_U,colour=mlab),position=position_dodge(width=0.9), width=0.1)	+
	scale_y_continuous(expand = c(0,0))  +
	#facet_grid(.~trsm) +
	facet_grid(trsm~.) +
	theme_bw(base_size=26) +
	labs(y='Date of diagnosis',x="Estimated date of infection \n(95% credible intervals)",colour="Ethnicity") +
	theme(legend.position="bottom",
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.margin = margin(1, 1, 0, 1, "cm")) +
	scale_colour_npg()
ggsave(g,file=file.path(home,'infection_date_v_diagnosis_date_tall.png'),w=12,h=17)


dat2 <- data.table(read.csv('/Users/alexb/Documents/GitHub/source.attr.with.infection.time/data_Ams/roadmap_cd4_vl_est.csv'))
dat2[, hiv_pos_d:=as.Date(hiv_pos_d,format="%d/%m/%Y")]

dat2[, HIV_POS_D:=hivc.db.Date2numeric(hiv_pos_d)]

dt2 <- unique(subset(dat2,select=c('id','mode','HIV_POS_D','estsctodiagMedian')))

dt2[, inf_date:=HIV_POS_D-estsctodiagMedian]
#dt2[is.na(estsctodiagMedian) & HIV_POS_D>2014,inf_date:=HIV_POS_D]
length(unique(dt2$id[dt2$inf_date>=2014]))

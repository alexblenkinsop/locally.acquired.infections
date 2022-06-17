# importations by subtype.R
# 
###############################################################################

cat(" \n -------------------------------- \n \n Running importations by subtype.R \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))

# for testing
if(1){
	args_dir <- list()
	args_dir[['stanModelFile']] <- 'branching_process_210810b_cmdstan'
	args_dir[['analysis']] <- 'analysis_200917'
	args_dir[['trsm']] <- 'MSM'
	args_dir[['in_dir']] <- '/rds/general/project/ratmann_roadmap_data_analysis/live'
	args_dir[['job_tag']] <- paste0('undiagnosed_weighted_inf_rate_2014-2018_',args_dir[['trsm']])
	args_dir[['out_dir']] <- paste0('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810b_cmdstan-',args_dir[['job_tag']])   
	args_dir[['with_subtypes']] <- 1
	args_dir[['source_dir']] <- '~/git/bpm'
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
	args_dir <- list()
	args_dir[['stanModelFile']] <- args_line[[2]]
	args_dir[['analysis']] <- args_line[[4]]
	args_dir[['in_dir']] <- args_line[[6]]
	args_dir[['out_dir']] <- args_line[[8]]
	args_dir[['job_tag']] <- args_line[[10]]
	args_dir[['trsm']] <- args_line[[12]]
	args_dir[['source_dir']] <- args_line[[14]]
} 

source(file.path(args_dir$source_dir, 'R', 'post-processing-plot-functions.R'))

#	load all input variables for this analysis run
do <- data.table(F=list.files(args_dir$out_dir, pattern='*_stanout.RData$', recursive=TRUE, full.name=TRUE))
z <- load( gsub('pbs_stanout.RData','pbs_stanin.RData',do[1,F]) )
with_subtypes = 0
if(!is.null(args$with_subtypes)) with_subtypes = args$with_subtypes
start_d = args$start_d
end_d = args$end_d

outfile.base <- paste0(args_dir$out_dir, "/",
											 args_dir$stanModelFile , "-", args_dir$job_tag)

# load inputs for this script
file <- paste0(outfile.base,'-stanout-basic.RDS')
cat("\n read RDS:", file)
pars.basic <- readRDS(file)
file <- paste0(outfile.base,'-stanout-subgraphs-gqs.RDS')
cat("\n read RDS:", file)
fit.subgraphs <- readRDS(file)
file <- paste0(outfile.base,'-stanout-origins-gqs.RDS')
cat("\n read RDS:", file)
fit.origins <- readRDS(file)
file <- paste0(outfile.base,'-stanout-bplace-gqs.RDS')
cat("\n read RDS:", file)
fit.bplace <- readRDS(file)


cat(" \n -------------------------------- Reading data -------------------------------- \n")

### Estimating importations
## read stanin
cat('\nReading Stan input data...')
infile.stanin <- list.files(args_dir$out_dir, pattern='stanin.RData$', recursive=TRUE)[1]
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(args_dir$out_dir, infile.stanin))
stopifnot(c('args','stan.data')%in%tmp)


# get loc labels
## read subgraphs with metadata
cat('\nReading subgraph data...')
infile.subgraphs <- file.path(args_dir$out_dir,'subgraphs_withmetadata.RDS')

dsubgraphtaxa <- readRDS(infile.subgraphs)

if(args$infdate==1){
	dsubgraphtaxa[INF_D>=start_d & INF_D<end_d, keep:=1]
}else{
	dsubgraphtaxa[HIV1_POS_D>=start_d & HIV1_POS_D<end_d, keep:=1]
}
dsubgraphtaxa <- subset(dsubgraphtaxa, SELECT==paste0('Ams',args$trsm) & keep==1)

dsubgraphtaxa$ORIGIN <- dsubgraphtaxa$ORIGINHOST
dsubgraphtaxa$ORIGIN[is.na(dsubgraphtaxa$ORIGIN)] <- 'Unknown'
do <- subset(dsubgraphtaxa,TRANSM==args$trsm)
do <- do[,list(N=length(FULL_NAME)),by=c('ORIGIN','TRANSM','ST','ST_CLADE')]
setnames(do,c('TRANSM'),c('TRM_GROUP'))

do <- subset(do, ORIGIN!='Unknown')

da <- do[, list(N=sum(N)), by=c('ORIGIN','TRM_GROUP','ST')]
da <- merge(da, da[, list(TOTAL=sum(N)), by=c('TRM_GROUP','ST')], by=c('TRM_GROUP','ST'))
tmp <- dcast.data.table(da,ORIGIN~ST,value.var='N')

loc <- data.table(location=seq(1:length(unique(tmp$ORIGIN))),
									loc_label=unique(tmp$ORIGIN))

dsubgraphtaxa <- subset(dsubgraphtaxa, REP=='000')
cat(" \n -------------------------------- Loading predicted subgraphs and origins -------------------------------- \n")

cat('\nPredicting chain sizes for partially unobserved chains...')
# extract predicted subgraphs
cs_ex <- fit.subgraphs$actual_cs # truncate here for subtype-specific
cs_em <- fit.subgraphs$actual_cs_unobs # truncate here for subtype-specific
N_sgs_e <- fit.subgraphs$N_sgs_e
N_sgs_unobs <- fit.subgraphs$N_sgs_unobs

cs_ex <- as.data.table( reshape2::melt( cs_ex ) )
setnames(cs_ex, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))
# remove 0 subgraphs (padding)
cs_ex <- subset(cs_ex,is.finite(new_cases))
cs_ex <- subset(cs_ex,new_cases>0)
# remove 1 from new_cases so they start from 0
cs_ex[, new_cases:=  new_cases - 1]
cs_ex[, size:= index_cases + new_cases]
cs_ex[, chains:= 'pre-existing']


cs_em <- as.data.table( reshape2::melt( cs_em ) )
setnames(cs_em, 1:5, c('iteration','subtype','subgraph','index_cases','new_cases'))
cs_em <- subset(cs_em,is.finite(new_cases))
cs_em <- subset(cs_em,new_cases>0)
# include index case in new cases (don't subtract 1)
cs_em[, size:= index_cases + new_cases]
cs_em[, chains:= 'emergent']

cat('\nSaving predicted chain sizes...')
saveRDS(cs_ex,file=paste0(outfile.base,'-','predicted_chains.rds'))
saveRDS(cs_em,file=paste0(outfile.base,'-','unobserved_chains.rds'))

cat('\nGet number of cases per subtype...')
st <- merge(cs_ex,cs_em,by=c("iteration","subtype","subgraph","index_cases","new_cases",  
															"size","chains"),all=T)
st <- st[, list(new_cases=sum(new_cases,na.rm=T)),by=c('subtype','iteration')]
st <- st[, list(subtype=subtype,prop=new_cases/sum(new_cases,na.rm=T)),by=c('iteration')]

N_sgs_em <- as.data.table( reshape2::melt( N_sgs_e ) )
setnames(N_sgs_em, 1:4, c('iteration','subtype','index_case','NT'))
# just keep index case==1
N_sgs_em <- subset(N_sgs_em,index_case==1,select=c('iteration','subtype','NT'))
# sum over index cases
#N_sgs_em <- N_sgs_em[, list(NT=sum(N)),by=c('iteration','subtype')]

N_sgs_u <- as.data.table( reshape2::melt( N_sgs_unobs ) )
setnames(N_sgs_u, 1:4, c('iteration','subtype','index_case','N_u'))
# just keep index case==1
N_sgs_u <- subset(N_sgs_u,index_case==1,select=c('iteration','subtype','N_u'))
N_sgs_u <- merge(N_sgs_em,N_sgs_u,by=c('iteration','subtype'),all=T)
N_sgs_u[, N_o:= NT - N_u]
saveRDS(N_sgs_u,file=paste0(outfile.base,'-','unobs_emergent_chains.RDS'))

cs_ex <- cs_ex[, list(new_cases=sum(new_cases)),by=c('iteration','subtype','chains')]
cs_em <- cs_em[, list(new_cases=sum(new_cases)),by=c('iteration','subtype','chains')]

# combine pre-exisiting and emergent subgraphs
# add new cases per MC sample
dt <- rbind(cs_ex,cs_em)

dt <- dt[, list(NI=sum(new_cases)),by=c('iteration','subtype')]

# add col for number of emergent subgraphs
dt <- merge(dt,N_sgs_em,by=c('iteration','subtype'),all.x=T)

# sum new cases per MC sample and subtype
#dt <- dt[, list(NI=sum(new_cases)),by=c('iteration','subtype')]

dt <- merge(dt, pars.basic$ds,by.x='subtype',by.y='subtypes')

cat('\nCalculate importations...')

# number of emergent chains/number of new cases
dt[,importations:= NT/NI]

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

dt$TRM <- args_dir$trsm

im <- dt[, list(q= quantile(importations, prob=ps,na.rm=T),
								q_label=p_labs),by=c('TRM','subtypes_name')]		

ans <- dcast.data.table(im, TRM+subtypes_name~q_label, value.var='q')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]

make_gqs_plots_st(ans,
									ans$TRM,
									paste0('Estimated importations between ',args$start_d, ' and ', args$end_d,' \n(posterior median and 95% credibility intervals)\n'),
									outfile.base,
									'importations_sbt_gqs')

saveRDS(ans,file=paste0(outfile.base,'-','importations_sbt_gqs.RDS'))

cat('\nSummarise origins...')

cat('\nLoad predicted origins of subgraphs...')
origins <- fit.origins$origins_subgraphs_x
origins_em <- fit.origins$origins_subgraphs_e

cat('\nPre-existing subgraphs...')
do <- as.data.table( reshape2::melt( origins ) )
setnames(do, 1:5, c('iteration','subtype','subgraph','index_cases','location'))
do <- subset(do,is.finite(location))
do <- subset(do,location>0)

do <- merge(do,loc,by='location')
do[,chains := 'pre-existing']

cat('\nEmergent subgraphs...')
du <- as.data.table( reshape2::melt( origins_em ) )
setnames(du, 1:5, c('iteration','subtype','subgraph','index_cases','location'))
du <- subset(du,is.finite(location))
du <- subset(du,location>0)
du <- subset(du,index_cases==1)

# just keep the predicted origins for the number of individuals in emergent subgraphs
du <- merge(du,loc,by='location')
du[,chains := 'emergent']

cat('\nCombine subgraph origins...')
do <- merge(do,du,by=c('iteration','subgraph','subtype','index_cases','location','loc_label','chains'),all=T)
do <- merge(do, pars.basic$ds,by.x='subtype',by.y='subtypes')

#	summarise subgraphs by origins
cat('\nSummarise origins...')
ans <- do[, list(N=length(subgraph)),by=c('iteration','loc_label','subtype')]
total <- do[, list(total=length(subgraph)),by=c('iteration','subtype')]
ans <- merge(ans,total,by=c('iteration','subtype'))
ans[,origin:= N/total]
ans <- ans[, list(q= quantile(origin, prob=ps),
									q_label=p_labs), 
					 by=c('loc_label','subtype')]		

ans <- dcast.data.table(ans, loc_label+subtype~q_label, value.var='q')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,pars.basic$ds,by.x='subtype',by.y='subtypes')

make_gqs_plots_st(ans,
									ans$loc_label,
									'Proportion of chains originating \nfrom geographic region \n(posterior median and 95% credibility intervals)\n',
									outfile.base,
									'origins_sbt_gqs')

saveRDS(ans,file=paste0(outfile.base,'-','origins_sbt_gqs.RDS'))

# sum over Ams/non-Ams origin to get subgraphs of external origins
do[, loc_Ams:= as.character(factor(grepl('Ams',loc_label), levels=c(TRUE,FALSE), labels=c('Ams','External')))]

tot <- do[, list(total=length(subgraph)),by=c('iteration','subtype')]
do <- do[, list(N=length(subgraph)),by=c('iteration','subtype','loc_Ams')]
do <- merge(do, tot, by=c('iteration','subtype'),all.x=T)
do[, p:=N/total]

##### add origins to estimate % local infections
#ex <- merge(dt,subset(ex,loc_Ams=='External'),by=c('iteration','subtype'))
ex <- merge(dt,subset(do,loc_Ams=='External', select=c('subtype','iteration','p')),by=c('iteration','subtype'))
ex[, ext_imp:= p * importations]
ex[, inf_Ams:= 1 - ext_imp]

saveRDS(ex,file=paste0(outfile.base,'-','external_importations_sbt_samples.RDS'))

ans <- ex[, list(ext_imp= quantile(ext_imp, prob=ps,na.rm=T),inf_Ams = quantile(inf_Ams, prob=ps,na.rm=T),
								q_label=p_labs),by=c('TRM','subtype')]		

ans <- dcast.data.table(ans, TRM+subtype~q_label, value.var='ext_imp')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,pars.basic$ds,by.x='subtype',by.y='subtypes')

make_gqs_plots_st(ans,
									ans$TRM,
									paste0('External importations between ',args$start_d, ' and ', args$end_d,' \n(posterior median and 95% credibility intervals)\n'),
									outfile.base,
									'ext_importations_sbt_gqs')

ans <- ex[, list(ext_imp= quantile(ext_imp, prob=ps,na.rm=T),inf_Ams = quantile(inf_Ams, prob=ps,na.rm=T),
								 q_label=p_labs),by=c('TRM','subtype')]		

ans <- dcast.data.table(ans, TRM+subtype~q_label, value.var='inf_Ams')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]
ans <- merge(ans,pars.basic$ds,by.x='subtype',by.y='subtypes')

make_gqs_plots_st(ans,
									ans$TRM,
									paste0('Infections acquired in Amsterdam between ',args$start_d, ' and ', args$end_d,' \n(posterior median and 95% credibility intervals)\n'),
									outfile.base,
									'infections_acquired_in_Amsterdam_sbt_gqs')

cat('\nSave outputs...')
saveRDS(dt,file=paste0(outfile.base,'-','estimate_importations_sbt_gqs.RDS'))


cat('\nCalculate local infections by migrant groups...')

geo <- data.table(read.csv('/rds/general/project/ratmann_roadmap_data_analysis/live/misc/NEWGEO.csv'))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']

load(file='/rds/general/project/ratmann_roadmap_data_analysis/live/analysis_200917/misc/200917_sequence_labels.rda')
dseq <- merge(dseq,geo,by.x='ORIGIN',by.y='Alpha_2_code',all.x=T)
setnames(dseq,'ORIGIN','BIRTH_COUNTRY_ISO')
dsubgraphtaxa <- merge(dsubgraphtaxa,unique(subset(dseq,select=c('PATIENT','BIRTH_COUNTRY_ISO','WRLD'))),by.x='ID',by.y='PATIENT',all.x=T)

## msm
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
dsubgraphtaxa[TRANSM=='HSX' & LOC_BIRTH %in% c("Africa"), mwmb:="G4"]
# caribbean and south america
dsubgraphtaxa[TRANSM=='HSX' & LOC_BIRTH %in% c("LaAmCar","FormerCurrDutchColonies"), mwmb:="G5"]


#dsubgraphtaxa <- subset(dsubgraphtaxa,TRANSM==args_dir$trsm & SELECT!='Ams' & REP=='000')
#if(args$infdate==1){
#	dsubgraphtaxa <- subset(dsubgraphtaxa,INF_D>=start_d & INF_D<end_d)
#}else{
#	dsubgraphtaxa <- subset(dsubgraphtaxa,HIV1_POS_D>=start_d & HIV1_POS_D<end_d)
#}

# for each strata, what is the proportion of each subtype
st_mwmb <- dsubgraphtaxa[, list(p=length(ID)),by=c('mwmb','ST')]
st_mwmb <- subset(st_mwmb,!is.na(mwmb))
st_mwmb <- st_mwmb[, list(ST=ST,p_st=p/sum(p)),by=c('mwmb')]

## load birth places from GQs
bp <- as.data.table( reshape2::melt( fit.bplace ) )
setnames(bp, 1:6, c('iteration','subtype','bplace','index_cases','N','chaintype'))
bp <- subset(bp,is.finite(N))
bp <- bp[, list(N=sum(N)),by=c('iteration','bplace','subtype')]
bp <- bp[, list(subtype=subtype,p_st=N/sum(N)),by=c('iteration','bplace')]
# add subtype labs
bp <- merge(bp, pars.basic$ds,by.x='subtype',by.y='subtypes')
# add birthplace labs (groups are ordered alphabetically)
labs <- data.table(mwmb=sort(unique(st_mwmb$mwmb)),bplace=seq(1:length(unique(st_mwmb$mwmb))))
bp <- merge(bp, labs,by='bplace',all=T)

ex <- merge(ex,bp,by=c('iteration','subtypes_name','subtype'),all=T,allow.cartesian=T)
#ex <- merge(ex,st_mwmb,by.x=c('subtypes_name'),by.y=('ST'),all=T,allow.cartesian=T)
#ex <- merge(ex,st,by=c('iteration','subtype'),all=T)
ex[is.nan(inf_Ams),inf_Ams:=0]

mwmb <- ex[, list(inf_Ams=sum(inf_Ams*p_st)),by=c('iteration','mwmb')]

ans <- mwmb[, list(inf_Ams = quantile(inf_Ams, prob=ps,na.rm=T),
								 q_label=p_labs),by=c('mwmb')]		

ans <- dcast.data.table(ans, mwmb~q_label, value.var='inf_Ams')
ans[, L:= paste0( round(M*100, d=1), '% [',  round(CL*100, d=1),'-', round(CU*100, d=1),'%]')]

g <- ggplot(ans, aes(x=mwmb)) +
	geom_errorbar(aes(ymin=CL, ymax=CU), width=0.3) +
	geom_point(aes(y=M)) +			
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
	coord_cartesian(ylim=c(0,1)) +
	theme_bw() +
	theme(axis.text.x=element_text(size=12,angle=90, vjust = 0.5, hjust=0),
				legend.position="bottom") +
	labs(y="% infections acquired in Amsterdam (median and 95% credible interval)", x='Migrant group')
ggsave(paste0(outfile.base,'-','locally_acquired_mwmb','.png'), w=6, h=6)

saveRDS(st_mwmb,file=paste0(outfile.base,'-','subtypes_prop_mwmb.RDS'))
saveRDS(mwmb,file=paste0(outfile.base,'-','local_infections_samples_mwmb.RDS'))
saveRDS(ans,file=paste0(outfile.base,'-','local_infections_mwmb.RDS'))

cat('\nCalculate local infections by risk groups...')

ex <- readRDS(file=paste0(outfile.base,'-','external_importations_sbt_samples.RDS'))
#st <- dsubgraphtaxa[, list(p=length(ID)),by=c('TRANSM','ST')]
#st <- st[, list(ST=ST,p_st=p/sum(p)),by=c('TRANSM')]

#ex <- merge(ex,st,by.x=c('subtypes_name'),by.y=('ST'),all=T,allow.cartesian=T)
ex <- merge(ex,st,by=c('iteration','subtype'),all=T)
ex[is.nan(inf_Ams),inf_Ams:=0]

local <- ex[, list(inf_Ams=sum(inf_Ams*prop)),by=c('iteration','TRM')]

saveRDS(st,file=paste0(outfile.base,'-','subtypes_prop_TRANSM.RDS'))
saveRDS(local,file=paste0(outfile.base,'-','local_infections_samples_TRANSM.RDS'))


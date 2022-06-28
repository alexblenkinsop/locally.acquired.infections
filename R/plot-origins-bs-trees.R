
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
library(ggsci)


outdir <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model'
outfile.base <- outdir

# plot origins

dsubgraphtaxa <- readRDS('/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/subgraphs_withmetadata_allreps.RDS')

dat0 <- subset(dsubgraphtaxa,REP=='000')
setnames(dat0,'ORIGINHOST','ORIGIN0')
dat1 <- unique(subset(dsubgraphtaxa,select=c(ST,ST_CLADE,REP,SELECT,NAME,
																						FULL_NAME,ORIGINHOST)))
set(dat0,NULL,c('REP','TAXA','ID','SEQ_YEAR','FULL_NAME'),NULL)
dat1 <- merge(dat1,dat0,by=c("ST","ST_CLADE","SELECT","NAME"),all.x=T)
dat1 <- subset(dat1,!is.na(ORIGIN0))

outfile <- file.path(args$outdir,'subgraphs_withmetadata_allreps_merged.RDS')
saveRDS(dsubgraphtaxa, file=outfile)

dat <- dsubgraphtaxa[, list(ORIGIN=ORIGIN0,pct=sum(ORIGINHOST==ORIGIN)),
										 by=c("ST","ST_CLADE","SELECT","NAME",
										 		 "TAXA","ID","SEQ_YEAR","FULL_NAME")]


dat <- dsubgraphtaxa[, list(REP=REP,FULL_NAME=FULL_NAME,ORIGIN=ORIGINHOST[REP=='000']),
										 by=c("ST","ST_CLADE","SELECT","NAME",
										 		 "TAXA","ID","SEQ_YEAR")]
outfile <- file.path(outdir,'subgraphs_withmetadata_allreps_origin0.RDS')
saveRDS(dat, file=outfile)

# can reduce dataset to not include all taxa - only need origin of subgraph ID


## counts by origin
ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

dat0 <- dsubgraphtaxa[REP=='000', list(N_0=length(FULL_NAME)),
										 by=c("ST","SELECT","REP","ORIGINHOST")]

dat <- dsubgraphtaxa[, list(N=length(FULL_NAME)),
										 by=c("ST","SELECT","REP","ORIGINHOST")]

dat <- dsubgraphtaxa[, list(N=length(FULL_NAME)),
										 by=c("ST","SELECT","REP","ORIGINHOST")]
tmp <- dat[, list(N = quantile(N, prob=ps,na.rm=T),
									q_label=p_labs),by=c("ST","SELECT","ORIGINHOST")]		

tmp <- dcast(tmp,ST+SELECT+ORIGINHOST~q_label,value.var='N')

dat <-merge(dat0,tmp,by=c("ST","SELECT","ORIGINHOST"),all=T)


g <- ggplot(subset(dat,SELECT!='Ams' & SELECT!='withmetadata'), aes(x=ORIGINHOST, y=N_0,fill=ORIGINHOST)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_grid(ST ~ SELECT,scales="free_y") +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
	#									 limits = c(0, 1),labels=scales::percent) +
	#coord_flip() +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_text( size=26,vjust = 0.5, hjust=0.5,angle=45),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Number of subgraphs with origin \n(median and 95% credible interval)", x='') #+
	#scale_fill_manual(values=col_pal) +
	#scale_colour_manual(values=col_pal)
#ggsci::scale_color_lancet() +
#ggsci::scale_fill_lancet()
ggsave(paste0(outfile.base,'/','bar_origins_BStrees','.png'), g, w=22, h=16)



##### for each subgraph in central tree, find subgraph in each bootstrap tree with most common tips to match

# test on one subgraph

dsubgraphtaxa <- subset(dsubgraphtaxa,ST!='subgraphs_withmetadata.rda' & SELECT!='Ams')

dt <- list()
dt_x <- list()
subtypes <- unique(dsubgraphtaxa$ST)
subtypes <- subtypes[!grepl('B',subtypes)]
for(s in subtypes){

	# find the unique subgraphs in the central tree and distinguish their names/origins for later
	#tmp <- subset(dsubgraphtaxa,ST==s & REP=='000' & NAME=='Ams-SPLIT1')
	tmp <- subset(dsubgraphtaxa,ST==s & REP=='000')
	#tmp[, FULL_NAME_0:=FULL_NAME]
	setnames(tmp,c('FULL_NAME','ORIGINHOST'),c('FULL_NAME_0','ORIGINHOST_0'))
	tmp[ST=='B',FULL_NAME_0:=paste0(ST,'_c',ST_CLADE,'_',REP,'_',SELECT,'_',NAME)]
	
	# find the subgraphs and taxa from the bootstrap trees
	tmp2 <- subset(dsubgraphtaxa,ST==s & REP!='000')
	#tmp2 <- subset(dsubgraphtaxa, REP!='000')
	
	# for each subgraph in the central tree, find every subgraph ID in the bootstrap tree which could be a match
	#tmp2 <- merge(tmp2,unique(subset(tmp,select=c('ST','FULL_NAME_0'))),by='ST',all.x=T)
	tmp2 <- merge(tmp2,unique(subset(tmp,select=c('SELECT','ST','FULL_NAME_0'))),by=c('SELECT','ST'),all.x=T,allow.cartesian = T)
	#tmp2 <- merge(tmp2,unique(subset(tmp,select=c('SELECT','ST','FULL_NAME_0'))),by=c('SELECT','ST'),all.x=T)
	#tmp3 <- tmp2[, list(FULL_NAME_0=FULL_NAME_0, Ntips=sum(TAXA %in% tmp$TAXA[tmp$FULL_NAME_0==FULL_NAME_0])),by=c('SELECT','ST','REP','FULL_NAME')]
	# for each subgraph from the central tree, count the number of tips in the bootstrap trees which are in the central tree
	tmp3 <- tmp2[, list(Ntips=sum(TAXA %in% tmp$TAXA[tmp$FULL_NAME_0==FULL_NAME_0])),by=c('SELECT','ST','REP','FULL_NAME','FULL_NAME_0')]

	# for each BS replicate, find the subgraph ID of the subgraph which has the most common tips with the central tree
	tmp4 <- tmp3[, list(Ntips=max(Ntips),FULL_NAME=FULL_NAME[Ntips==max(Ntips)]),by=c('SELECT','ST','REP','FULL_NAME_0')]

	# add origins from central and bootstrap trees
	tmp4 <- merge(tmp4,subset(tmp,select=c('SELECT','ST','FULL_NAME_0','ORIGINHOST_0')),by=c('SELECT','ST','FULL_NAME_0'),all.x=T,allow.cartesian=T)
	tmp4 <- merge(tmp4,unique(subset(dsubgraphtaxa,select=c('SELECT','ST','FULL_NAME','ORIGINHOST'))),by=c('SELECT','ST','FULL_NAME'),all.x=T)
	
	# calculate % of closest BS trees that share the same origin as the central tree
	dat <- tmp4[, list(pct=sum(ORIGINHOST==ORIGINHOST_0,na.rm=T)/length(ORIGINHOST[!is.na(ORIGINHOST)])),
											 by=c("SELECT","ST","FULL_NAME_0",'ORIGINHOST_0')]

	dt[[s]] <- dat
	
	g <- ggplot(subset(dat,!is.na(ORIGINHOST_0)), aes(x=FULL_NAME_0, y=pct,colour=ORIGINHOST_0,size=5)) +
		geom_point(stat='identity', position = "dodge") +
		#geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
		facet_grid(. ~ SELECT,scales="free",space="free") +
		#coord_cartesian(ylim=c(0,1)) +
		#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
		scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
											 limits = c(0, 1),labels=scales::percent) +
		#coord_flip() +
		theme_bw(base_size = 20) +
		theme(axis.text.x=element_blank(),
					axis.text.y=element_text(size=24),
					axis.title=element_text(size=26),
					legend.position="none",
					strip.text.x=element_text(size=24),
					strip.background=element_blank(),
					strip.placement = "outside") +
		labs(y="% subgraphs from boostrap trees \nwith same origin as central tree", x='subgraph') #+
	#ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_subtype_',s,'.png'), g, w=22, h=16)
	
	# relabel Ams/nonAms
	tmp4[grepl('Ams',ORIGINHOST_0) & !is.na(ORIGINHOST_0), ORIGIN0:='Ams']
	tmp4[!grepl('Ams',ORIGINHOST_0) & !is.na(ORIGINHOST), ORIGIN0:='Ext']
	tmp4[grepl('Ams',ORIGINHOST) & !is.na(ORIGINHOST_0), ORIGIN:='Ams']
	tmp4[!grepl('Ams',ORIGINHOST) & !is.na(ORIGINHOST), ORIGIN:='Ext']
	
	dt_x[[s]] <- tmp4[, list(pct=sum(ORIGIN==ORIGIN0,na.rm=T)/length(ORIGIN[!is.na(ORIGIN)])),
									by=c("SELECT","ST","FULL_NAME_0",'ORIGIN0')]

}


s <- 'B'
clades <- c(1,2,3,4)
dt_b <- list()
dtb_x <- list()
for(c in clades){
	
	#tmp <- subset(dsubgraphtaxa,ST==s & REP=='000' & NAME=='Ams-SPLIT1')
	tmp <- subset(dsubgraphtaxa,ST==s & ST_CLADE==c & REP=='000')
	#tmp[, FULL_NAME_0:=FULL_NAME]
	setnames(tmp,c('FULL_NAME','ORIGINHOST'),c('FULL_NAME_0','ORIGINHOST_0'))
	tmp[ST=='B',FULL_NAME_0:=paste0(ST,'_c',ST_CLADE,'_',REP,'_',SELECT,'_',NAME)]
	
	tmp2 <- subset(dsubgraphtaxa,ST==s & REP!='000')
	tmp2[ST=='B',FULL_NAME:=paste0(ST,'_c',ST_CLADE,'_',REP,'_',SELECT,'_',NAME)]
	reps <- tmp2
	#tmp2 <- subset(dsubgraphtaxa, REP!='000')
	
	#tmp2 <- merge(tmp2,unique(subset(tmp,select=c('ST','FULL_NAME_0'))),by='ST',all.x=T)
	tmp2 <- merge(tmp2,unique(subset(tmp,select=c('SELECT','ST','FULL_NAME_0'))),by=c('SELECT','ST'),all.x=T,allow.cartesian = T)
	#tmp2 <- merge(tmp2,unique(subset(tmp,select=c('SELECT','ST','FULL_NAME_0'))),by=c('SELECT','ST'),all.x=T)
	#tmp3 <- tmp2[, list(FULL_NAME_0=FULL_NAME_0, Ntips=sum(TAXA %in% tmp$TAXA[tmp$FULL_NAME_0==FULL_NAME_0])),by=c('SELECT','ST','REP','FULL_NAME')]
	tmp3 <- tmp2[, list(Ntips=sum(TAXA %in% tmp$TAXA[tmp$FULL_NAME_0==FULL_NAME_0])),by=c('SELECT','ST','REP','FULL_NAME','FULL_NAME_0')]
	
	tmp4 <- tmp3[, list(Ntips=max(Ntips),FULL_NAME=FULL_NAME[Ntips==max(Ntips)]),by=c('SELECT','ST','REP','FULL_NAME_0')]
	
	# add origins from central and bootstrap trees
	tmp4 <- merge(tmp4,subset(tmp,select=c('SELECT','ST','FULL_NAME_0','ORIGINHOST_0')),by=c('SELECT','ST','FULL_NAME_0'),all.x=T,allow.cartesian=T)
	tmp4 <- merge(tmp4,unique(subset(reps,select=c('SELECT','ST','FULL_NAME','ORIGINHOST'))),by=c('SELECT','ST','FULL_NAME'),all.x=T)
	
	dat <- tmp4[, list(pct=sum(ORIGINHOST==ORIGINHOST_0,na.rm=T)/length(ORIGINHOST[!is.na(ORIGINHOST)])),
							by=c("SELECT","ST","FULL_NAME_0",'ORIGINHOST_0')]
	
	dt_b[[c]] <- dat
	
	g <- ggplot(subset(dt_b[[c]],!is.na(ORIGINHOST_0)), aes(x=FULL_NAME_0, y=pct,colour=ORIGINHOST_0,size=5)) +
		geom_point(stat='identity', position = "dodge") +
		#geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
		facet_grid(. ~ SELECT,scales="free",space="free") +
		#coord_cartesian(ylim=c(0,1)) +
		#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
		scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
											 limits = c(0, 1),labels=scales::percent) +
		#coord_flip() +
		theme_bw(base_size = 20) +
		theme(axis.text.x=element_blank(),
					axis.text.y=element_text(size=24),
					axis.title=element_text(size=26),
					legend.position="none",
					strip.text.x=element_text(size=24),
					strip.background=element_blank(),
					strip.placement = "outside") +
		labs(y="% subgraphs from boostrap trees \nwith same origin as central tree", x='subgraph') #+
	#ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_subtype_',s,'_c',c,'.png'), g, w=22, h=16)
	
	# relabel Ams/nonAms
	tmp4[grepl('Ams',ORIGINHOST_0) & !is.na(ORIGINHOST_0), ORIGIN0:='Ams']
	tmp4[!grepl('Ams',ORIGINHOST_0) & !is.na(ORIGINHOST), ORIGIN0:='Ext']
	tmp4[grepl('Ams',ORIGINHOST) & !is.na(ORIGINHOST_0), ORIGIN:='Ams']
	tmp4[!grepl('Ams',ORIGINHOST) & !is.na(ORIGINHOST), ORIGIN:='Ext']
	
	dtb_x[[s]] <- tmp4[, list(pct=sum(ORIGIN==ORIGIN0,na.rm=T)/length(ORIGIN[!is.na(ORIGIN)])),
									by=c("SELECT","ST","FULL_NAME_0",'ORIGIN0')]
	
}

# plot all subtypes
dt_ball <- do.call('rbind',dt_b)

dat <- do.call('rbind',dt)
dat <- rbind(dat,dt_ball)

saveRDS(dat,file=paste0(outfile.base,'/','origins_bootstrap_trees_agreement.RDS'))

g <- ggplot(subset(dat,!is.na(ORIGINHOST_0)), aes(x=FULL_NAME_0, y=pct,colour=ORIGINHOST_0,size=5)) +
	geom_point(stat='identity', position = "dodge") +
	#geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#facet_grid(ST ~ SELECT,scales="free",space="free") +
	facet_wrap(ST ~ SELECT,scales="free",ncol=2) +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	#coord_flip() +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_blank(),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="none",
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% subgraphs from boostrap trees \nwith same origin as central tree", x='subgraph') +
#scale_fill_manual(values=col_pal) +
#scale_colour_manual(values=col_pal)
ggsci::scale_color_npg()
#ggsci::scale_fill_lancet()
#ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes','.png'), g, w=22, h=16)
ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes','.png'), g, w=22, h=24)

## MSM
g <- ggplot(subset(dat,!is.na(ORIGINHOST_0) & SELECT=='AmsMSM'), aes(x=FULL_NAME_0, y=pct,colour=ORIGINHOST_0)) +
	geom_point(size=5, position = "dodge") +
	#geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#facet_grid(ST ~ .,scales="free",space="free") +
	facet_wrap(ST ~ .,scales="free",ncol=2) +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	#coord_flip() +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_blank(),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="bottom",
				legend.text=element_text(size=24),
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% subgraphs from boostrap trees \nwith same origin as central tree", x='subgraph',colour='origin',size='') +
	#scale_fill_manual(values=col_pal) +
	#scale_colour_manual(values=col_pal)
	ggsci::scale_color_npg()
#ggsci::scale_fill_lancet()
#ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes','.png'), g, w=22, h=16)
ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes_MSM','.png'), g, w=22, h=24)



## HSX
g <- ggplot(subset(dat,!is.na(ORIGINHOST_0) & SELECT=='AmsHSX'), aes(x=FULL_NAME_0, y=pct,colour=ORIGINHOST_0)) +
	geom_point(size=5, position = "dodge") +
	#geom_errorbar(aes(x=ORIGINHOST,ymin=CL, ymax=CU,fill=ORIGINHOST),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	#facet_grid(ST ~ .,scales="free",space="free") +
	facet_wrap(ST ~ .,scales="free",ncol=2) +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales::percent) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	#coord_flip() +
	theme_bw(base_size = 20) +
	theme(axis.text.x=element_blank(),
				axis.text.y=element_text(size=24),
				axis.title=element_text(size=26),
				legend.position="bottom",
				legend.text=element_text(size=24),
				strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% subgraphs from boostrap trees \nwith same origin as central tree", x='subgraph',colour='origin',size='') +
	#scale_fill_manual(values=col_pal) +
	#scale_colour_manual(values=col_pal)
	ggsci::scale_color_npg()
#ggsci::scale_fill_lancet()
#ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes','.png'), g, w=22, h=16)
ggsave(paste0(outfile.base,'/','pct_subgraph_origins_BStrees_allsubtypes_HSX','.png'), g, w=22, h=24)


ds <- subset(dat,!is.na(ORIGINHOST_0))
ds <- ds[, list(pct80=length(pct[pct>0.75])/length(pct))]


hist(subset(dat,SELECT=='AmsHSX' & ST=='01AE')$pct)

st_N <- dat[,list(N=length(FULL_NAME_0)),by=c('ST')]
st_N <- st_N[order(-N),]

dat$ST <- factor(dat$ST,levels=st_N$ST)
brk <- function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)
dat$TRSM=factor(dat$SELECT,levels=c('AmsMSM','AmsHSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
g <- ggplot(subset(dat,!is.na(ORIGINHOST_0))) + geom_histogram(aes(x=pct,fill=ST)) +
	facet_grid(ST ~ TRSM,scales="free") +
	theme_bw(base_size = 34) +
	theme(#axis.text.x=element_blank(),
				#axis.text.y=element_text(size=24),
				#axis.title=element_text(size=26),
				legend.position="none",
				#legend.text=element_text(size=24),
				#strip.text.x=element_text(size=24),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="Number subgraphs", x='% of bootstrap trees with same origin as central tree',fill='') +
	ggsci::scale_fill_npg()
ggsave(paste0(outfile.base,'/','hist_pct_subgraph_origins_BStrees_allsubtypes','.png'), g, w=22, h=24)


## histogram by Ams/ext ----

# plot all subtypes
dt_bxall <- do.call('rbind',dtb_x)

#dat <- do.call('rbind',dt_x)
dat <- do.call('rbind',dt)
dat <- rbind(dat,dt_bxall)

saveRDS(dat,file=paste0(outfile.base,'/','origins_bootstrap_trees_agreement_Ams_ext.RDS'))

st_N <- dat[,list(N=length(FULL_NAME_0)),by=c('ST')]
st_N <- st_N[order(-N),]

dat$ST <- factor(dat$ST,levels=st_N$ST)
brk <- function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)
dat$TRSM=factor(dat$SELECT,levels=c('AmsMSM','AmsHSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
g <- ggplot(subset(dat,!is.na(ORIGIN0))) + geom_histogram(aes(x=pct,fill=ST)) +
	facet_grid(ST ~ TRSM,scales="free") +
	theme_bw(base_size = 34) +
	theme(#axis.text.x=element_blank(),
		#axis.text.y=element_text(size=24),
		#axis.title=element_text(size=26),
		legend.position="none",
		#legend.text=element_text(size=24),
		#strip.text.x=element_text(size=24),
		strip.background=element_blank(),
		strip.placement = "outside") +
	labs(y="Number subgraphs", x='% of bootstrap trees with same origin (Ams/external) \nas central tree',fill='') +
	ggsci::scale_fill_npg()
ggsave(paste0(outfile.base,'/','hist_pct_subgraph_origins_Ams_ext_BStrees_allsubtypes','.png'), g, w=22, h=24)


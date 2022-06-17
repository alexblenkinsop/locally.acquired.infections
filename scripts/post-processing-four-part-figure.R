require(ggpubr)
require(gridExtra)
require(ggsci)
col_pal = c(pal_npg("nrc")(9)[c(1:5)],pal_npg("nrc")(9)[c(8,6,7,9)])
col_pal_c = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

outfile.base <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/branching_process_210810m_cmdstan-undiagnosed_untilmay_2014-2018_HSX/branching_process_210810m_cmdstan-undiagnosed_untilmay_2014-2018_HSX'

cat(" \n -------------------------------- plot A -------------------------------- \n")
# from post-processing-make-figure-inf-diag-seq-CIs.R

dd <- readRDS(paste0(outfile.base,'-','N_undiagnosed_samples_MSM_HSX.RDS'))

dz <- dd[, list(N_diag=quantile(diag,prob=c(0.025,0.5,0.975)),
									N_inf=quantile(N_inf,prob=c(0.025,0.5,0.975)),
									#undiag=quantile(undiag,prob=c(0.025,0.5,0.975)),
									qlabel=c('L','M','U')),
					 by=c('trsm','mg','variable')] # summarise quantiles for Jan of each year
dz <- subset(dz,qlabel=='M',select=c('trsm','mg','variable','N_diag'))
dz$variable <- factor(dz$variable,levels=c('N_diagnoses_NL_MSM','N_diagnoses_G1_MSM','N_diagnoses_G2_MSM',
																						 'N_diagnoses_G3_MSM','N_diagnoses_Oth_MSM',
																						 'N_diagnoses_NL_HSX', 'N_diagnoses_G4_HSX',
																						 'N_diagnoses_G5_HSX','N_diagnoses_Oth_HSX'))
dz[, mg_fb:=1]
dz[trsm=='MSM' & mg==4, mg_fb:=0]
dz[trsm=='HSX' & mg==3, mg_fb:=0]
dz$mg_fb <- as.character(dz$mg_fb)
dz$trsm <- factor(dz$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
setnames(dz,'variable','lab')

pct <- subset(pct,variable=='Diagnosed' & bplace!='All')
pct[, pct:=round(value/sum(value)*100,0)]
pct[, lab:=paste0(value,'\n(',pct,'%)')]

dat4 <- merge(dz,subset(dat3,variable=='Diagnosed'& bplace!='All'),by=c('trsm','mg_fb'),all=T)

#load(file=paste0(outfile.base,'-','plot_A_new.Rdata'))

plot <- ggplot(data=subset(dat4,variable=='Diagnosed' & bplace!='All')) +
	geom_bar(aes(x=bplace,y=N_diag,fill=lab),stat='identity', position = "stack") +
	geom_errorbar(aes(x=bplace,ymin=L, ymax=U), width=0.5, colour="black",size=1.5)	+
	#geom_text(data=pct,aes(label=lab,x=bplace,y=value),  vjust=-0.2,size = 10) +
	#scale_y_sqrt(breaks=seq(0,600,50)) +
	facet_grid(.~trsm,scales="free_y") +
	scale_y_continuous(expand = c(0,0))  +
	coord_cartesian(ylim=c(0,400)) +
	labs(x='',y="Number of new \n diagnosed infections\n(median and 95% credible interval)",fill="") +
	theme_bw(base_size=36) +
	theme(legend.position="none",
				strip.background=element_blank(),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_pal) 
ggsave(file=paste0(outfile.base,'-inf-bplace_ci_2014-2018.png'),plot,w=12, h=5)
ggsave(file=paste0(outfile.base,'-inf-bplace_ci_2014-2018.pdf'),plot,w=12, h=5)

plot_A <- ggarrange(plot,ncol=1,labels=c('A'),font.label=list(size=40))
#ggsave(paste0(outfile.base,'-','A_inf_2014-2018','.png'), plot_A, w=22, h=10)

save(dat4,pct,file=paste0(outfile.base,'-','plot_A_new.Rdata'))

cat(" \n -------------------------------- plot B -------------------------------- \n")
# from undiagnosed_bytrmgroup.R
#load(file=paste0(outfile.base,'-','plot_B_av.Rdata')) # do I need this?

outdir <- '/rds/general/project/ratmann_roadmap_data_analysis/live/branching_process_model/undiagnosed_211102-undiag_untilmay2019_weights'
job_tag <- 'undiag_untilmay2019_weights'

do <- readRDS(file.path(outdir,'metadata_with_inftime_mgroup.RDS'))

## bar plot 
msm_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_MSM",'.RDS')))
hsx_s <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_HSX",'.RDS')))
du <- rbind(msm_s,hsx_s)
#du <- readRDS(dd,file=paste0(outfile.base,'-','p_undiagnosed_samples_MSM_HSX.RDS'))
du <- data.table(du)

msm_all <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",'MSM','.RDS')))
hsx_all <- readRDS(file=file.path(outdir,paste0('p_undiagnosed_average_',job_tag,"_all_",'HSX','.RDS')))
da <- rbind(msm_all,hsx_all)
da <- data.table(da)
da[, mlab:='All']
da[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
da[trsm=='MSM', trsm := 'Amsterdam MSM']


du[trsm=='HSX', trsm := 'Amsterdam heterosexuals']
du[trsm=='MSM', trsm := 'Amsterdam MSM']
du <- merge(du, unique(subset(do,select=c('trsm','mgid','mlab'))),by.x=c('trsm','mg'),by.y=c('trsm','mgid'),all=T)
du <- du[!is.na(du$trsm),]

du <- merge(du,da,by=c('trsm','mlab','p0.025','p0.5','p0.975'),all=T)

du$trsm <- factor(du$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexuals'))

du$mlab <- factor(du$mlab,levels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. Europe',
																	 'S. America &\n Caribbean','Sub-Saharan\nAfrica','Other'),
									labels=c('All','NL',' W.Europe,\nN.America,Oceania','E. & C. Europe',
													 'S. America &\n Caribbean','Sub-Saharan Africa','Other'))
du[, mwmb_lab_trsm:=paste0(trsm,' ',mlab)]
order <- unique(subset(du,select=c('trsm','mlab','mwmb_lab_trsm')))
order <- order[order(trsm,mlab),]
du$mwmb_lab_trsm <- factor(du$mwmb_lab_trsm,levels=unique(order$mwmb_lab_trsm))

plot_msmhsx <- ggplot(data=subset(du,!is.na(mwmb_lab_trsm)),aes(x=mlab,y=p0.5,fill=mwmb_lab_trsm)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mwmb_lab_trsm),position=position_dodge(width=0.9), width=0.5, colour="black",size=1)	+
	facet_wrap(. ~ trsm,scales="free_y") +
	#scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	coord_flip() +
	labs(x='',y="% infections undiagnosed by May 2019\n(median and 95% credible interval)",fill="Ethnicity") +
	theme_bw(base_size=36) +
	theme(legend.position="none",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_pal_c) 
ggsave(file.path(outdir,'undiagnosedav_mwmb_plotB_nolab.png'),plot_msmhsx,w=34, h=16)
ggsave(file.path(outdir,'undiagnosedav_mwmb_plotB_nolab.pdf'),plot_msmhsx,w=34, h=16)

plot_B <- ggarrange(plot_msmhsx,ncol=1,labels=c('B'),font.label=list(size=40))
#plot_B <- ggarrange(plot_msmhsx,ncol=1,labels=c('A'),font.label=list(size=40))

ggsave(file.path(outdir,'undiagnosedav_mwmb_plotB.png'),plot_B,w=34, h=16)

save(du,file=paste0(outfile.base,'-','plot_B_av.Rdata'))

# get leg
col_msm = pal_npg("nrc")(9)[c(1:5)]
col_hsx = pal_npg("nrc")(9)[c(8,6,7,9)]

plot_msm <- ggplot(data=subset(du,trsm=='Amsterdam MSM' & !is.na(mwmb_lab_trsm) & mlab!='All' ),aes(x=mlab,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_wrap(. ~ trsm,scales="free_y") +
	#scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	coord_flip() +
	labs(x='',y="% infections undiagnosed \n by May 2019",fill="MSM, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_msm) + guides(fill = guide_legend(title.position = "top"))
#ggsave(file.path(outdir,'undiagnosedav_mwmb_MSM_plotB.png'),plot_msm,w=34, h=16)

plot_hsx <- ggplot(data=subset(du,trsm=='Amsterdam heterosexuals' & !is.na(mwmb_lab_trsm) & mlab!='All'),aes(x=mlab,y=p0.5,fill=mlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mlab,ymin=p0.025, ymax=p0.975,fill=mlab),position=position_dodge(width=0.9), width=0.5, colour="black")	+
	facet_wrap(. ~ trsm,scales="free_y") +
	#scale_y_continuous(expand = c(0,0),labels = scales::percent)  +
	coord_cartesian(ylim=c(0,1)) +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	coord_flip() +
	labs(x='',y="% infections undiagnosed \n by May 2019",fill="Heterosexuals, area of birth") +
	theme_bw(base_size=36) +
	theme(legend.position="bottom",
				legend.text=element_text(size=36),
				strip.background=element_blank(),
				#strip.text = element_blank(),
				plot.title = element_text(hjust = 0.5),
				axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				strip.text.x=element_text(size=38)
	) +
	scale_fill_manual(values=col_hsx) + guides(fill = guide_legend(title.position = "top"))
#ggsave(file.path(outdir,'undiagnosedav_mwmb_HSX_plotB.png'),plot_hsx,w=34, h=16)

leg_msm <- get_legend(plot_msm, position = NULL)
leg_hsx <- get_legend(plot_hsx, position = NULL)
leg <- ggarrange(leg_msm,leg_hsx,ncol=2,align="hv")
#leg <- annotate_figure(leg,
#								top = text_grob("Area of birth", color = "black", face = "bold", size = 42))

ggsave(file.path(outdir,'undiagnosedav_mwmb_legend_plotB.png'),leg,w=34, h=4)
ggsave(file.path(outdir,'undiagnosedav_mwmb_legend_plotB.pdf'),leg,w=34, h=4)

cat(" \n -------------------------------- plot C -------------------------------- \n")
# from post-processing-estimates-by-migrant-groups.R

#load(file=paste0(outfile.base,'-','plot_C.Rdata'))

col_pal_c = c("grey50",pal_npg("nrc")(9)[c(1:5)],"grey50",pal_npg("nrc")(9)[c(8,6,7,9)])

#tab <- dat[, list(qs= quantile(inf_Ams, prob=ps), qlab=p_labs), by=c('trsm','mwmb_lab','xlab','time')]
#ds <- dcast(tab,trsm+mwmb_lab+time+xlab~qlab,value.var='qs')

ds <- readRDS(file=paste0(outfile.base,'-','locally_acquired_bplace_MSM-HSX_2014-2018_cases.RDS'))

ds$trsm <- factor(ds$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexual'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
g_mwmb <- ggplot(subset(ds), aes(x=mwmb_lab, y=M,fill=xlab)) +
	geom_bar(stat='identity', position = "dodge") +
	geom_errorbar(aes(x=mwmb_lab,ymin=CL, ymax=CU,fill=xlab),position=position_dodge(width=0.9), width=0.5, colour="black",size=1)	+
	facet_wrap(. ~ trsm,scales="free_y") +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),
										 limits = c(0, 1),labels=scales::percent) +
	coord_flip() +
	theme_bw(base_size = 36) +
	theme(axis.text.x=element_text( size=36,vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				legend.position="none",
				strip.text.x=element_text(size=38),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% infections acquired in Amsterdam \n(median and 95% credible interval)", x='') +
	scale_fill_manual(values=col_pal_c) +
	scale_colour_manual(values=col_pal_c)
plot_C <- ggarrange(g_mwmb,ncol=1,labels=c('C'),font.label=list(size=40))
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_plotC','.png'), g_mwmb, w=22, h=10)
ggsave(paste0(outfile.base,'-','locally_acquired_birthplace_plotC','.pdf'), g_mwmb, w=22, h=10)

save(ds,file=paste0(outfile.base,'-','plot_C.Rdata'))

cat(" \n -------------------------------- plot D -------------------------------- \n")
#load(file=paste0(outfile.base,'-','plot_D.Rdata'))

dat_s <- readRDS(,file=paste0(outfile.base,'-','p_local_infections_migrantgroup_samples.RDS'))
bp_s <- readRDS(file=paste0(outfile.base,'-','N_infections_migrantgroup_samples.RDS'))

col_pal = c(pal_npg("nrc")(9)[c(1:5)],pal_npg("nrc")(9)[c(8,6,7,9)])

dat <- merge(dat_s,bp_s,by=c('trsm','mwmb','iteration'),all=T)
dat[, N_local:=floor(inf_Ams*N)]

dat$mwmb_lab <- factor(dat$mwmb,levels=c('All','NL','G1','G2','G3','G4','G5','Other'),labels=c('All','NL','W.Europe,\nN.America,\nOceania','E. & C. \nEurope','S. America \n& Caribbean','Sub-Saharan \nAfrica','S. America \n& Caribbean','Other'))
dat[, fb:='Foreign-born']
dat[mwmb=='NL', fb:='NL']
dat_c <- subset(dat,mwmb_lab!='All')

dats <- dat_c[, list(N_local=sum(N_local)),by=c('trsm','fb','iteration')]
dats <- dats[, list(trsm=trsm,fb=fb,pct_gp=N_local/sum(N_local)),by=c('iteration')]

tab_f <- dats[, list(pct_gp= quantile(pct_gp, prob=ps),
									qlab=p_labs), by=c('trsm','fb')]

dat <- dat_c[, list(trsm=trsm,mwmb=mwmb,mwmb_lab=mwmb_lab,fb=fb,pct_gp=N_local/sum(N_local)),by=c('iteration')]

tab <- dat[, list(pct_gp= quantile(pct_gp, prob=ps),
									qlab=p_labs), by=c('trsm','mwmb_lab','fb')]

dat <- dcast(tab_f,trsm+fb~qlab,value.var=c('pct_gp'))
dat_mg <- dcast(tab,trsm+mwmb_lab+fb~qlab,value.var=c('pct_gp'))
setnames(dat_mg,c('M','CL','CU'),c('fb_M','fb_CL','fb_CU'))
dat <- merge(dat,dat_mg,by=c('trsm','fb'),all=T)

dat$trsm <- factor(dat$trsm,levels=c('MSM','HSX'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
dat$fb <- factor(dat$fb,levels=c('NL','Foreign-born'),labels=c('Dutch-born','Foreign-born'))
dat[, mwmb_lab_trsm:=paste0(trsm,' ',mwmb_lab)]

order <- unique(subset(dat,select=c('trsm','fb')))
order <- unique(subset(dat,select=c('trsm','fb','mwmb_lab','mwmb_lab_trsm')))
order <- order[order(trsm,fb,mwmb_lab),]
dat$mwmb_lab_trsm <- factor(dat$mwmb_lab_trsm,levels=unique(order$mwmb_lab_trsm))

#dat$trsm <- factor(dat$trsm,levels=c('Amsterdam MSM','Amsterdam heterosexual'),labels=c('Amsterdam MSM','Amsterdam heterosexuals'))
	g_mwmb_N <- ggplot(subset(dat), aes(x=fb, y=fb_M,fill=mwmb_lab_trsm)) +
	geom_bar(stat='identity', position = "stack") +
	geom_errorbar(aes(x=fb,ymin=CL, ymax=CU), width=0.5, colour="black",size=1.5)	+
	facet_wrap(. ~ trsm,scales="free_x") +
	scale_y_continuous(expand=expansion(mult = c(0, .05)), breaks=seq(0,0.5,0.1), minor_breaks=seq(0,0.5,0.05),
											 limits = c(0, 0.5),labels=scales::percent_format(accuracy = 1L)) +
	theme_bw(base_size = 36) +
	theme(axis.text.x=element_text( size=36, vjust = 0.5, hjust=0.5),
				axis.text.y=element_text(size=32),
				axis.title=element_text(size=38),
				legend.position="none",
				strip.text.x=element_text(size=38),
				strip.background=element_blank(),
				strip.placement = "outside") +
	labs(y="% of infections acquired in \n Amsterdam by ethnicity \n(median and 95% credible interval)", x='',fill='Ethnicity') +
	scale_fill_manual(values=col_pal,labels=c('NL (MSM)','W.Europe, N.America,\nOceania (MSM)','E. & C. Europe (MSM)','S. America & Caribbean (MSM)',
																						'Other (MSM)','NL (heterosexual)','S. America & \nCaribbean (heterosexual)','Sub-Saharan \nAfrica (heterosexual)','Other (heterosexual)')) +
	scale_colour_manual(values=col_pal)
#plot_D <- ggarrange(g_mwmb_N,ncol=1,labels=c('D'),font.label=list(size=40))
plot_D <- ggarrange(g_mwmb_N,ncol=1,labels=c('B'),font.label=list(size=40))
ggsave(paste0(outfile.base,'-','N_local_inf_plot_D','.png'), plot_D, w=22, h=10)
ggsave(paste0(outfile.base,'-','N_local_inf_plot_D','.pdf'), plot_D, w=22, h=10)

save(dat,file=paste0(outfile.base,'-','plot_D.Rdata'))

### combine panels
plot_all <- grid.arrange(
	plot_A,
	plot_B,
	plot_C,
	plot_D,
	widths = c(0.6, 0.4),
	layout_matrix = rbind(c(1,NA),
												c(2,2),
												c(3,3),
												c(4,NA)))
plot_all <- ggarrange(plot_all,leg,ncol=1,heights=c(1,0.1))
ggsave(paste0(outfile.base,'-','four-part-plot_final','.png'), plot_all, w=32, h=42)
ggsave(paste0(outfile.base,'-','four-part-plot_final','.pdf'), plot_all, w=32, h=42)

ggsave(paste0(outfile.base,'-','four-part-plot_final_hires2','.jpeg'), plot_all,dpi=300, w=32, h=42)

ggsave(paste0(outfile.base,'-','four-part-plot_final_hires3','.png'), plot_all,dpi=600, w=32, h=42)


ggsave(paste0(outfile.base,'-','four-part-plot_final_hires300','.png'), plot_all,dpi=300, w=32, h=42)
ggsave(paste0(outfile.base,'-','four-part-plot_final_hires200','.png'), plot_all,dpi=200, w=32, h=42)


### combine panels
plot_all <- grid.arrange(
	plot_B,
	plot_D,
	widths = c(0.7, 0.3),
	layout_matrix = rbind(c(1,1),
												c(2,NA)))
plot_all <- ggarrange(plot_all,leg,ncol=1,heights=c(1,0.1))
ggsave(paste0(outfile.base,'-','two-part-plot_final_hires300','.png'), plot_all,dpi=300, w=34, h=34)

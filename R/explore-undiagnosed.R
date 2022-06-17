


dat <- data.table(TRSM=c(rep('MSM',9),rep('HSX',9)),
									year=c(rep(seq(2010,2018,1),2)),
									pr_undiag=rep(0,18))

for(trsm in c('MSM','HSX')){
	for(i in 2010:2018){
		#print(trsm)
		#print(i)
		infile.undiag <- paste0('~/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/undiagnosed_AMS_',trsm,'.csv')
		pr <- est_prop_undiagnosed(infile.undiag,i,i+1)
		print(pr)
		dat[TRSM==trsm & year==i,pr_undiag:=pr]
	}
}

average <- dat[, list(av=mean(pr_undiag)),by=c('TRSM')]

av <- dat %>%                                        # Specify data frame
	group_by(trsm) %>%                         # Specify group indicator
	summarise_at(vars(Sepal.Length),              # Specify column
							 list(name = mean)) 
dat <- merge(dat,average,by=c('TRSM'),all.x=T)

g <- ggplot(data=dat) + geom_line(aes(x=year,y=pr_undiag)) +
	facet_grid(.~TRSM) +
	geom_hline(data= average, aes(yintercept=av),colour="black",size=0.5,linetype=2) +
	theme_bw(base_size=25) +
	labs(x='year',y='proportion undiagnosed') +
	theme(strip.background=element_blank())
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/undiagnosed.png',g,w=15, h=10)


infile.indinfo <- file.path(args$indir,'Data','data_191223_Amsterdam','SHM_1902_ROADMAP_191223_tblBAS.csv')
dbas <- read.csv(infile.indinfo,header=T)
dbas <- subset(dbas,select=c(PATIENT,RECART_D,HIV1_POS_D,HIV1_POS_D_A,MODE))
dbas[,'HIV1_POS_D'] <- as.Date(dbas[,'HIV1_POS_D'],format="%d/%m/%Y")
dbas <- data.table(dbas)
dbas[, TRSM:='Other']
dbas[MODE==1, TRSM:='MSM']
dbas[MODE==6, TRSM:='HSX']

dat <- read.csv(infile.cd4,header=T)
dat <- subset(dat,CD4_U==1)
dat$CD4_D <- as.Date(dat$CD4_D,format=c("%Y-%m-%d"))

#diagd <- subset(dsubgraphtaxa,select=c('ID','HIV1_POS_D'))
dat <- merge(dat,dbas,by='PATIENT',all.x=T)
dat <- data.table(dat)
dat[CD4_D>=HIV1_POS_D, A:=1]
dat <- subset(dat,A==1)

tmp <- dat %>%
	group_by(PATIENT) %>%
	filter(CD4_D == min(CD4_D))
setnames(tmp,c('CD4_D','CD4_V'),c('CD4_first_D','CD4_first_V'))

tmp <- data.table(tmp)

pr <- tmp[, list(pr_200=length(PATIENT[CD4_first_V<200]),pr_350=length(PATIENT[CD4_first_V<350]),pr_500=length(PATIENT[CD4_first_V<500])),by=c('TRSM')]
pr <- subset(pr,TRSM!='Other')

pr <- melt(pr)
pr[, cd4:=gsub('pr_','',variable)]
pr[cd4==200, time:=7.93]
pr[cd4==350, time:=4.19]
pr[cd4==500, time:=1.19]
pr[, N:=nrow(tmp[tmp$TRSM=='MSM'])]
pr[TRSM=='HSX', N:=nrow(tmp[tmp$TRSM=='HSX'])]
pr[, p:= value/N]

cd4 <- ggplot(data=pr,aes(x=time,y=p)) + geom_point() +
	geom_smooth(method=lm)	+
	labs(x='time to cd4 count below threshold',y='pr(diagnoses w/cd4 count below threshold)') +
	theme_bw(base_size=25)
	ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/cd4_time_prop.png',cd4,w=15, h=10)
	
pr[, nvalue:=N-value]	
pr$TRSM <- as.factor(pr$TRSM)
lm <- glm(as.matrix(cbind(value,nvalue)) ~ TRSM + time ,data=pr,family = binomial(link = "logit"))
#lm <- glm(p ~ TRSM + time ,data=pr,family = binomial(link = "logit"))
years <- seq(1,10,1)
glm.pr <- predict(lm, newdata=data.frame(time=rep(years,2),TRSM=c(rep('MSM',length(years)),rep('HSX',length(years)))), se.fit=TRUE)
glm.pr <- data.table(time=rep(years,2), TRSM=c(rep('MSM',length(years)),rep('HSX',length(years))), log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

#glm.pr <- predict(lm, newdata=pr, se.fit=TRUE)
#glm.pr <- data.table(time=pr$time, TRSM=pr$TRSM, log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

require(boot)
glm.pr[, pr_low:=inv.logit(log_p)]
glm.pr[, pr_low_cl:=inv.logit(log_p-2*log_p_sd)]
glm.pr[, pr_low_cu:=inv.logit(log_p+2*log_p_sd)]

#dat <- merge(pr,glm.pr,by=c('time'),all=T)

g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=pr_low_cl, ymax=pr_low_cu ),alpha=0.6) + 
	geom_point(data=glm.pr, aes(x=time, y=pr_low),col="red",alpha=0.6) + 
	geom_point(data=pr, aes(x=time, y=p), colour='blue') +
	#geom_smooth(data=pr, aes(x=time, y= p)) +
	labs(x='time to cd4 count below threshold',y="pr(diagnoses w/cd4 count below threshold)") +
	facet_wrap(.~TRSM) +
	theme_bw(base_size=25) +
	theme(strip.background = element_blank())
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/cd4_time_prop_glm.png',g2,w=15, h=10)

g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=log_p-2*log_p_sd, ymax=log_p+2*log_p_sd),alpha=0.6) + 
	geom_point(aes(x=time, y=log_p),col="red",alpha=0.6) + 
	geom_point(data=pr, aes(x=time, y=logit(p)), colour='blue') +
	labs(x='time to cd4 count below threshold',y="pr(diagnoses w/cd4 count below threshold)") +
	facet_wrap(.~TRSM) +
	theme_bw() +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/cd4_time_prop_glm_logit.png',g2,w=15, h=10)


###############
lm <- glm(as.matrix(cbind(value,N)) ~ time ,data=subset(pr,TRSM=='HSX'),family = binomial)
years <- seq(1,10,1)
glm.pr <- predict(lm, newdata=data.frame(time=years,length(years)), se.fit=TRUE)
glm.pr <- data.table(time=years, length(years), log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

glm.pr[, pr_low:=inv.logit(log_p)]
glm.pr[, pr_low_cl:=inv.logit(log_p-2*log_p_sd)]
glm.pr[, pr_low_cu:=inv.logit(log_p+2*log_p_sd)]

hsx <- copy(glm.pr)
hsx[, TRSM:='HSX']

lm <- glm(as.matrix(cbind(value,N)) ~ time ,data=subset(pr,TRSM=='MSM'),family = binomial)
years <- seq(1,10,1)
glm.pr <- predict(lm, newdata=data.frame(time=years,length(years)), se.fit=TRUE)
glm.pr <- data.table(time=years, length(years), log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

glm.pr[, pr_low:=inv.logit(log_p)]
glm.pr[, pr_low_cl:=inv.logit(log_p-2*log_p_sd)]
glm.pr[, pr_low_cu:=inv.logit(log_p+2*log_p_sd)]

msm <- copy(glm.pr)
msm[, TRSM:='MSM']

glm.pr <- rbind(hsx,msm)

g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=pr_low_cl, ymax=pr_low_cu ),alpha=0.6) + 
	geom_point(aes(x=time, y=pr_low),col="red",alpha=0.6) + 
	geom_point(data=pr, aes(x=time, y=p), colour='blue') +
	labs(x='time to cd4 count below threshold',y="pr(diagnoses w/cd4 count below threshold)") +
	facet_wrap(.~TRSM) +
	theme_bw() +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/cd4_time_prop_glm.png',g2,w=15, h=10)

####################
lm <- glm(as.matrix(cbind(value,N)) ~ TRSM + time ,data=pr,family = binomial)
years <- seq(1,10,1)
glm.pr <- predict(lm, newdata=data.frame(time=years,TRSM=rep(c('MSM','HSX'),length(years))), se.fit=TRUE)
glm.pr <- data.table(time=years, TRSM=rep(c('MSM','HSX'),length(years)), log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

glm.pr <- predict(lm, newdata=pr, se.fit=TRUE)
glm.pr <- data.table(time=pr$time, TRSM=pr$TRSM, log_p=glm.pr$fit, log_p_sd=glm.pr$se.fit)

require(boot)
glm.pr[, pr_low:=inv.logit(log_p)]
glm.pr[, pr_low_cl:=inv.logit(log_p-2*log_p_sd)]
glm.pr[, pr_low_cu:=inv.logit(log_p+2*log_p_sd)]

g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=pr_low_cl, ymax=pr_low_cu ),alpha=0.6) + 
	geom_point(aes(x=time, y=pr_low),col="red",alpha=0.6) + 
	geom_point(data=pr, aes(x=time, y=p), colour='blue') +
	labs(x='time to cd4 count below threshold',y="pr(diagnoses w/cd4 count below threshold)") +
	facet_wrap(.~TRSM) +
	theme_bw() +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file='/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Data/Undiagnosed/cd4_time_prop_glm.png',g2,w=15, h=10)


lm(p~time,data=pr)

pr <- tmp[, list(pr_200=length(PATIENT[CD4_first_V<200])),by='A']
dsubgraphsize <- dsubgraphtaxa[, list(icasesart=length(ID[INF_D<args$start_d & RECART_D>=args$start_d]),icases=length(ID[INF_D<args$start_d & supp==0]),jcases=length(ID[INF_D>=args$start_d & INF_D<args$end_d])), by=c('ST','ST_CLADE','REP','SELECT','NAME')]


dsubgraphtaxa <- merge(dsubgraphtaxa,subset(tmp,select=c('PATIENT','CD4_first_D','CD4_first_V')),by.x=c('ID'),by.y=c('PATIENT'),all.x=T)




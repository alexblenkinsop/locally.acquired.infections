

dat <- subset(dsubgraphtaxa,REP=='000' & SELECT!='Ams')

dat[, posd:=cut(HIV1_POS_D, breaks=seq(1980,2019.5,0.5), right = FALSE)]
#dat[, t1:= as.numeric( sub("\\((.+),.*", "\\1", posd))]
#dat[, t2:= as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", posd))]

dat[, t1:= as.numeric( sub("\\[(.+),.*", "\\1", posd))]
dat[, t2:= as.numeric( sub("[^,]*,([^]]*)\\)", "\\1", posd))]

lp <- dat[, list(diag=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2]),
								 late=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2 & CD4_first_V<350])),by=c('SELECT','posd','t1','t2')]
lp[, notlate:=diag-late]
lp[, plate:=late/diag]
lp[, loglatep:=logit(plate)]

lp <- subset(lp,t1>=2010)
lp <- subset(lp,t1>=2000)

g <- ggplot(lp, aes(x=posd, y= plate)) + 
	geom_point() +
	geom_smooth(aes(x=posd, y= plate)) +
	#stat_smooth(method = "lm") +
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	facet_wrap(.~SELECT) +
	labs(x='time',y='% late presenters') +
	theme_bw(base_size=24) +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_6months.png",g,w=20,h=15)

##### by subtype
lst <- dat[, list(diag=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2]),
									late=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2 & CD4_first_V<350])),by=c('SELECT','ST','posd','t1','t2')]
lst[, notlate:=diag-late]
lst[, plate:=late/diag]
lst[, loglatep:=logit(plate)]

lst <- subset(lst,t1>=2010)
lst <- subset(lst,t1>=2000)
lst[,time_discrete:= (t1+t2)/2] 

g <- ggplot(lst, aes(x=time_discrete, y= plate)) + 
	geom_point() +
	geom_smooth(aes(x=time_discrete, y= plate)) +
	#stat_smooth(method = "lm") +
	scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	facet_grid(ST~SELECT) +
	labs(x='time',y='% late presenters') +
	theme_bw(base_size=24) +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_bysubtype_trmgroup2.png",g,w=25,h=20)

# logit glm
glm_st <- glm(as.matrix(cbind(late,notlate)) ~ SELECT + time_discrete + ST ,data=lst,family = binomial(link = "logit"))
years <- seq(2010,2020,0.01)
dat.pr <- data.frame(time_discrete=years,SELECT=rep(c('AmsMSM','AmsHSX'),length(years)),ST=rep(rep(unique(lst$ST),2),length(years)))
#dat.pr <- dat.pr %>% expand(time_discrete,SELECT,ST)
#glm.pr <- predict(glm_st, newdata=data.frame(time_discrete=years,SELECT=rep(rep(c('AmsMSM','AmsHSX'),length(years)),length((unique(lst$ST)))),ST=rep(rep(unique(lst$ST),length(years)),2)), se.fit=TRUE)
glm.pr <- predict(glm_st, newdata=dat.pr, se.fit=TRUE)
glm.pr <- data.table(dat.pr, log_latep=glm.pr$fit, log_latep_sd=glm.pr$se.fit)
glm.pr[, l_late_cl:=log_latep-2*log_latep_sd]
glm.pr[, l_late_cu:=log_latep+2*log_latep_sd]
glm.pr[, pr_latep:=inv.logit(log_latep)]
glm.pr[, late_cl:=inv.logit(l_late_cl)]
glm.pr[, late_cu:=inv.logit(l_late_cu)]
lst_ss <- subset(lst,time_discrete>=2010)
st.order <- lst_ss[, list(diag=sum(diag)),by=c('ST')]
st.order <- st.order[order(-st.order$diag),]
st.labs <- lst_ss[, list(diag=sum(diag)),by=c('ST','SELECT')]
glm.pr$ST <- factor(glm.pr$ST,levels=c(st.order$ST))
lst$ST <- factor(lst$ST,levels=c(st.order$ST))
st.labs$ST <- factor(st.labs$ST,levels=c(st.order$ST))
g2 <- ggplot(subset(glm.pr,time_discrete>=2010)) + 
	geom_ribbon(aes(x=time_discrete, ymin=late_cl, ymax=late_cu ),alpha=0.6) + 
	geom_line(aes(x=time_discrete, y=pr_latep)) +
	geom_point(data=subset(lst,time_discrete>=2010), aes(x=time_discrete, y=plate), colour='blue') +
	labs(x='year', y='pr(late presentation)') +
	geom_text(data=st.labs,aes(2019, 0.75,label=paste("Diagnoses= ", diag))) +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	#coord_cartesian(ylim=c(0,1)) +
	scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_grid(ST~SELECT) +
	theme_bw(base_size=24)
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_ST_glm_pred.png",g2,w=20,h=15)

#################
lm <- glm(plate ~ SELECT + t1 ,data=lp,family = binomial(link = "logit"))
years <- seq(2010,2020,0.01)
glm.pr <- predict(lm, newdata=data.frame(t1=years,SELECT=rep(c('AmsMSM','AmsHSX'),length(years))), se.fit=TRUE)
glm.pr <- data.table(time=years, SELECT=rep(c('AmsMSM','AmsHSX'),length(years)), log_latep=glm.pr$fit, log_latep_sd=glm.pr$se.fit)
g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=log_latep-2*log_latep_sd, ymax=log_latep+2*log_latep_sd ),alpha=0.6) + 
	geom_line(aes(x=time, y=log_latep)) +
	geom_point(data=lp, aes(x=t1, y=loglatep), colour='blue') +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	#coord_cartesian(ylim=c(0,1)) +
	scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_wrap(.~SELECT) +
	theme_bw()
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_pred.png",g2,w=20,h=15)

lm <- glm(as.matrix(cbind(late,notlate)) ~ SELECT + posd ,data=lp,family = binomial)
years <- sort(unique(lp$posd))
glm.pr <- predict(lm, newdata=data.frame(posd=years,SELECT=rep(c('AmsMSM','AmsHSX'),length(years))), se.fit=TRUE)
glm.pr <- data.table(time=years, SELECT=rep(c('AmsMSM','AmsHSX'),length(years)), log_latep=glm.pr$fit, log_latep_sd=glm.pr$se.fit)

glm.pr[, pr_latep:=inv.logit(log_latep)]

g2 <- ggplot(glm.pr) + 
	geom_ribbon(aes(x=time, ymin=log_latep-2*log_latep_sd, ymax=log_latep+2*log_latep_sd ),alpha=0.6) + 
	geom_line(aes(x=time, y=log_latep)) +
	geom_point(aes(x=time, y=log_latep),col="red",alpha=0.6) + 
	geom_point(data=lp, aes(x=posd, y=loglatep), colour='blue') +
	labs(y="logit pr(late presentation)") +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_wrap(.~SELECT) +
	theme_bw() +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_pred_6month2.png",g2,w=20,h=15)

g2 <- ggplot(glm.pr) + 
	geom_point(aes(x=time, y=pr_latep),col="red",alpha=0.6) + 
	geom_point(data=lp, aes(x=posd, y=plate), colour='blue') +
	labs(y="pr(late presentation)") +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	#coord_cartesian(ylim=c(0,1)) +
	#scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_wrap(.~SELECT) +
	theme_bw(base_size=24) +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_pred_6months.png",g2,w=20,h=15)

###############
library(mgcv, quietly = TRUE)
require(splines)
st_b <- subset(lp,ST=='B')
spline_fit_st <- gam(
	formula = as.matrix(cbind(late,notlate)) ~
		# cr is cubic spline with k knots
		s(t2, bs = "cr", k = 2, by = as.factor(SELECT)),
	family = binomial, data = lp)
pdf("/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_spline_2k.pdf",w=10,h=7) 
plot(spline_fit_st, scale = 0, page = 1, rug = FALSE)
dev.off() 
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_spline_ST.png",gsp,w=20,h=15)

glm_fit_big$coefficients
spline_fit$coefficients["(Intercept)"]

##########################################################
#### GLM WITH SPLINE ON TIME
lp[,time_discrete:= (t1+t2)/2] 
num_knots <- 3
knots = seq(min(lp$t1), max(lp$t2), length.out = num_knots+2)[-c(1, num_knots+2)]
spline_degree <- 2

hsx <- subset(lp,SELECT=='AmsHSX')
hsx <- hsx[order(hsx$t1),]
spline_fit_hsx <- glm(
	formula = as.matrix(cbind(late,notlate)) ~
	bs(time_discrete, df= num_knots + spline_degree - 1, degree= spline_degree, intercept = F),
	family = binomial, data = hsx)
spline_fit_hsx

msm <- subset(lp,SELECT=='AmsMSM')
msm <- msm[order(msm$t1),]
spline_fit_msm <- glm(
	formula = as.matrix(cbind(late,notlate)) ~
		bs(time_discrete, df= num_knots + spline_degree - 1, degree= spline_degree, intercept = FALSE),
	family = binomial, data = msm)
spline_fit_msm

summary(spline_fit_hsx)
model.matrix(spline_fit_hsx)
est <- data.table(model.matrix(spline_fit_hsx) %*% coef(spline_fit_hsx))
# est <- spline_fit_hsx$fitted.values # is equivalent to the above
est[, time_discrete:=sort(hsx$time_discrete)]
est[, pt:= inv.logit(V1)]
est[, plate:=hsx$plate]
hsx <- merge(hsx,subset(est,select=c('time_discrete','pt')),by=c('time_discrete'),all.x=T)

summary(spline_fit_msm)
model.matrix(spline_fit_msm)
est <- data.table(model.matrix(spline_fit_msm) %*% coef(spline_fit_msm))
est[, time_discrete:=sort(msm$time_discrete)]
est[, pt:= inv.logit(V1)]
est[, plate:=msm$plate]
msm <- merge(msm,subset(est,select=c('time_discrete','pt')),by=c('time_discrete'),all.x=T)

dat <- rbind(hsx,msm)

X <- model.matrix(spline_fit_msm)
sigma <- summary(spline_fit_msm)$sigma
var.Yhat <- (diag(X %*% solve(t(X) %*% X) %*% t(X)) + 1) * sigma^2

g2 <- ggplot(dat) + 
	#geom_point(aes(x=time_discrete, y=pt),col="red",alpha=0.6) + 
	geom_point(aes(x=time_discrete, y=plate), colour='blue') +
	geom_line(aes(x=time_discrete, y=pt),col="red",alpha=0.6) + 
	labs(y="pr(late presentation)",x="year") +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_wrap(.~SELECT) +
	theme_bw(base_size=24) +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_sp_hsxmsm.png",g2,w=15,h=15)

###########################################################
#### GLM WITH SPLINE ON TIME bY SUBTYPE
lst[,time_discrete:= (t1+t2)/2] 
num_knots <- 3
knots = seq(min(lst$t1), max(lst$t2), length.out = num_knots+2)[-c(1, num_knots+2)]
spline_degree <- 2

#hsx <- subset(lst,SELECT=='AmsHSX')
#msm <- subset(lst,SELECT=='AmsMSM')
hsx <- list()
msm <- list()

for(st in unique(lst$ST)){
	hsx[[st]] <- subset(lst,SELECT=='AmsHSX' & ST==st)
	hsx[[st]] <- hsx[[st]][order(hsx[[st]]$t1),]
	spline_fit_hsx <- glm(
		formula = as.matrix(cbind(late,notlate)) ~
			bs(time_discrete, df= num_knots + spline_degree - 1, degree= spline_degree, intercept = F),
		family = binomial, data = hsx[[st]])

	msm[[st]] <- subset(lst,SELECT=='AmsMSM' & ST==st)
	msm[[st]] <- msm[[st]][order(msm[[st]]$t1),]
	spline_fit_msm <- glm(
		formula = as.matrix(cbind(late,notlate)) ~
			bs(time_discrete, df= num_knots + spline_degree - 1, degree= spline_degree, intercept = FALSE),
		family = binomial, data = msm[[st]])

	#summary(spline_fit_hsx)
	#model.matrix(spline_fit_hsx)
	est <- data.table(model.matrix(spline_fit_hsx) %*% coef(spline_fit_hsx))
	# est <- spline_fit_hsx$fitted.values # is equivalent to the above
	est[, time_discrete:=sort(hsx[[st]]$time_discrete)]
	est[, pt:= inv.logit(V1)]
	est[, plate:=hsx[[st]]$plate]
	hsx[[st]] <- merge(hsx[[st]],subset(est,select=c('time_discrete','pt')),by=c('time_discrete'),all.x=T)
	
	#summary(spline_fit_msm)
	#model.matrix(spline_fit_msm)
	est <- data.table(model.matrix(spline_fit_msm) %*% coef(spline_fit_msm))
	est[, time_discrete:=sort(msm[[st]]$time_discrete)]
	est[, pt:= inv.logit(V1)]
	est[, plate:=msm[[st]]$plate]
	msm[[st]] <- merge(msm[[st]],subset(est,select=c('time_discrete','pt')),by=c('time_discrete'),all.x=T)
	
}
hsx <- do.call('rbind',hsx)
msm <- do.call('rbind',msm)
dat <- rbind(hsx,msm)

#X <- model.matrix(spline_fit_msm)
#sigma <- summary(spline_fit_msm)$sigma
#var.Yhat <- (diag(X %*% solve(t(X) %*% X) %*% t(X)) + 1) * sigma^2

g2 <- ggplot(dat) + 
	#geom_point(aes(x=time_discrete, y=pt),col="red",alpha=0.6) + 
	geom_point(aes(x=time_discrete, y=plate), colour='blue') +
	geom_line(aes(x=time_discrete, y=pt),col="red",alpha=0.6) + 
	labs(y="pr(late presentation)",x="year") +
	#scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05),labels=scales:::percent) +
	coord_cartesian(ylim=c(0,1)) +
	#scale_x_continuous(breaks=seq(2010,2020,5)) +
	facet_grid(ST~SELECT) +
	theme_bw(base_size=24) +
	theme(strip.background = element_blank(),
				axis.text.x=element_text(size=18,angle=45, vjust = 0.5, hjust=0.5))
ggsave(file="/Users/alexb/Box Sync/Roadmap/RQ1 Estimating introductions/Manuscript/latepresenters_sp_subtype_hsxmsm.png",g2,w=25,h=20)

###########################################################

matplot(hsx$time_discrete, model.matrix(spline_fit_hsx)[,-1], type='l', 
				+        col=rainbow(5), lwd=2)
legend("top", legend = dimnames(model.matrix(spline_fit_hsx))[[2]][-1], 
				 +        col=rainbow(5), lty=1:5, bty="n", lwd=2)

glm.pr <- predict(spline_fit, type="response")
plot(spline_fit)
plot(hsx$time_discrete, hsx$log_latep)
#lines(x, exp(mu), lwd = 2, col = "grey")
lines(x, glm.pr, col = "red", lwd = 2)

years <- sort(unique(lp$posd))
years <- knots
glm.pr <- predict(spline_fit, newdata=data.frame(posd=years,length(years)), se.fit=TRUE)
glm.pr <- data.table(time=years, log_latep=glm.pr$fit, log_latep_sd=glm.pr$se.fit)

glm.pr[, pr_latep:=inv.logit(log_latep)]

B <- t(bs(X, knots=seq(-5,5,1), degree=3, intercept = F)) 
B <- t(bs(lp$time_discrete, df= num_knots + spline_degree - 1, degree= spline_degree, intercept = F))

intervals = find_intervals(knots, spline_degree)
bs(sort(unique(lp$time_discrete)), knots, spline_degree)

###############
lp <- dat[, list(t1=seq(2010,2019,0.5),t2=seq(2010.5,2019.5,0.5)),by=c('SELECT')]
lp <- dat[, {t1=seq(2010,2019,0.5)
						t2=seq(2010.5,2019.5,0.5)
						list(diag=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2]),
																											 late=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2 & CD4_first_V<350]))},
					by=c('SELECT')]


dt[,{tmp1=mean(mpg); tmp2=mean(abs(mpg-tmp1)); tmp3=round(tmp2, 2); list(tmp2=tmp2, tmp3=tmp3)}, by=cyl]
lp[,{tmp1=mean(mpg); tmp2=mean(abs(mpg-tmp1)); tmp3=round(tmp2, 2); list(tmp2=tmp2, tmp3=tmp3)}, by=cyl]

								 
								 diag=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2]),
								 late=length(ID[HIV1_POS_D>=t1 & HIV1_POS_D<t2 & CD4_first_V<350])),by=c('SELECT')]
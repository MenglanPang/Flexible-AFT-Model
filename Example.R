#Example for Flexible AFT model using septic shock dataset
library(survival)
library(splines)
sepsis<-read.csv("sepsis.csv")
source("FlexAFT.R")
Vars<-c("immunodep","knauss","urinaire","cirrhose","infection","age","sofa")


#NL<-c(rep(0,5),rep(1,2))
#nknot.NL<-c(rep(NA,5),rep(1,2))
#degree.NL<-c(rep(NA,5),rep(2,2))
#TD<-c(0,0,0,rep(1,4))
#nknot.TD<-c(NA,NA,NA,rep(1,4))
#degree.TD<-c(NA,NA,NA,rep(2,4))


##Input:
#Data: Dataset (no missing data;categorial variables have to be coded as dummy variables)
#Var: names of the variable to be included in the spline-based AFT model
#NL: Indicate for each variable if the NL effect is modeled (0/1). 
#    Can be 1 only for continuous variables
#TD: Indicate for each variable if the TD effect is modeled (0/1)
#nknot.NL: number of interior knots of the splines for NL effects
#nknot.TD: number of interior knots of the splines for TD effects
#degree.NL: degree of the splines for NL effects
#degree.TD: degree of the splines for TD effects
#nknot.bh: number of the interior knots of the splines for modeling the baseline hazard
#degree.bh: degree of the splines for modeling the baseline hazard
#Time.Obs: name of the observed time variable
#Delta: name of the event indicator
#knot_time: "eventtime" or "alltime"; whether the knots are placed based on the quantiles of the uncensored event time or the overal observed (event or censoring) time 
#tol: convergece criter evaluating two consecutive log-likelihood;default is 1^-5
#ndivision: the number of intervals defining the granularity of the numerical computation of integrals; 
##default is "maxtime" such that the number of interval=100*floor(Observed event/censoring time) 

sink("septicshock_flexible_AFT.txt")
sepsis_fit<-FlexAFT(Data=sepsis,Var=Vars,
                    NL=c(rep(0,5),rep(1,2)),TD=c(0,0,0,rep(1,4)),
                    nknot.NL=c(rep(NA,5),rep(1,2)),nknot.TD=c(NA,NA,NA,rep(1,4)),
                    degree.NL=c(rep(NA,5),rep(2,2)),degree.TD=c(NA,NA,NA,rep(2,4)),
                    nknot.bh=2,degree.bh=3,
                    Time.Obs="time.obs",Delta="event",knot_time="eventtime",tol=1e-5,ndivision=600)
sink()
save(sepsis_fit,file="sepsis_flexaft.rda")

cov.val<-rep(NA,7)
for (i in 1:5){
  cov.val[i]<-0
}
cov.val[6]<-min(sepsis$age)
cov.val[7]<-min(sepsis$sofa)


par(mfrow=c(2,3),mgp=c(2,1,0),mar=c(3.5,3.5,2,1),
    oma=c(0.25,0.25,0.25,0.25))
Hazard.est<-HazardEst(fit=sepsis_fit,time=seq(0.3,60,0.1),cov=cov.val,Data=sepsis)
plot(0,0,type="n",xlim=c(0,60),ylim=c(0,0.006),xlab="Time (day)",ylab="Hazard",cex.main=0.75,main="Hazard")
lines(Hazard.est$time,Hazard.est$hazard)
Data<-sepsis
Time.Obs<-"time.obs"
Delta<-"event"

Time.obs<-Data[Data[,Delta]==1,Time.Obs]
Time.obs.freq<-table(Time.obs)
Time.obs.uniq<-as.numeric(names(Time.obs.freq))
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
mtext('(a)', side=1, line=2, at=1,cex=0.75)


plot(0,0,type="n",xlim=c(0,60),ylim=c(0,6),ylab=expression(e^beta(t)),xlab="Time(Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

TD_effect<-TDest_bin(fit=sepsis_fit,Data=sepsis,var.name="cirrhose",time=seq(0,60,0.1))

lines(TD_effect$time,exp(TD_effect$est.TD),lty=1,lwd=1)

TD_effect<-TDest_bin(fit=sepsis_fit,Data=sepsis,var.name="infection",time=seq(0,60,0.1))

lines(TD_effect$time,exp(TD_effect$est.TD),lty=2)

legend("topright",c("Cirrhosis","Infection type(nosocomial)"),
       lty=1:2,
       #fill=c(rgb(0.192,0.192,0.192,0.1),rgb(0.128,0.128,0.128,0.3)),
       cex=0.65)
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
#abline(h=0,lty=2,col="grey")
mtext('(b)', side=1, line=2, at=1,cex=0.75)

###age
plot(0,0,type="n",xlim=c(20,80),ylim=c(-0.8,0.6),xlab="Age",ylab=expression(g(X)),cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="NL effect")
NL_est<-NLest(fit=sepsis_fit,Data=sepsis,var.name="age",
                 ref.val=66,Q.high=0.9,Q.low=0.1,var.val=seq(20,80,0.1))


lines(NL_est$X,NL_est$est.NL)

X.freq<-table(sepsis[,"age"])
X.uniq<-as.numeric(names(X.freq))
for (j in 1:length(X.uniq)){
  rug(X.uniq[j],ticksize=X.freq[j]/50)
  
}
mtext('(c)', side=1, line=2, at=20,cex=0.75)

plot(0,0,type="n",xlim=c(0,60),ylim=c(0.45,5),ylab=expression(beta(t)),xlab="Time (Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

TD_est<-TDest_con(fit=sepsis_fit,Data=sepsis,var.name="age",
                 Q.high=0.9,Q.low=0.1,time=seq(0,60,0.1))

lines(TD_est$time,TD_est$est.TD)
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
legend("bottomright","Age",lty=1,cex=0.65)
mtext('(d)', side=1, line=2, at=1,cex=0.75)

plot(0,0,type="n",xlim=c(3,20),ylim=c(-0.55,0.75),xlab="SOFA",ylab=expression(g(X)),cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="NL effect")
NL_est<-NLest(fit=sepsis_fit,Data=sepsis,var.name="sofa",
                 ref.val=11,Q.high=0.9,Q.low=0.1,var.val=seq(3,20,0.1))

lines(NL_est$X,NL_est$est.NL)

X.freq<-table(sepsis[,"sofa"])
X.uniq<-as.numeric(names(X.freq))
for (j in 1:length(X.uniq)){
  rug(X.uniq[j],ticksize=X.freq[j]/1000)
  
}
mtext('(e)', side=1, line=2, at=4,cex=0.75)

plot(0,0,type="n",xlim=c(0,60),ylim=c(0,11),ylab=expression(beta(t)),xlab="Time (Day)",cex.axis=0.8,cex.lab=0.8,cex.main=0.75,main="TD effect")

TD_est<-TDest_con(fit=sepsis_fit,Data=sepsis,var.name="sofa",
                 Q.high=0.9,Q.low=0.1,time=seq(0,60,0.1))

lines(TD_est$time,TD_est$est.TD)
for (j in 1:length(Time.obs.uniq)){
  rug(Time.obs.uniq[j],ticksize=0.005*Time.obs.freq[j])
}
legend("bottomright","SOFA",lty=1,cex=0.65)
mtext('(f)', side=1, line=2, at=1,cex=0.75)


###survival
par(mfrow=c(1,1))
Surv.est<-SurvEst(fit=sepsis_fit,time=seq(0.3,60,0.1),cov=cov.val,Data=sepsis)
plot(0,0,type="n",xlim=c(0,60),ylim=c(0,1),xlab="Time (day)",ylab="Survival",cex.main=0.75,main="Survival")
lines(Surv.est$time,Surv.est$survival)

##most common values
#immunodep=0; knauss=0; urinaire=0;infection=0;
#cirrhose=0;age=65;sofa=10
X.ref<-rep(NA,7)
for (i in 1:5){
  X.ref[i]<-0
}
X.ref[6]<-65
X.ref[7]<-11

X.index<-X.ref
X.index[5]<-1
pcttime<-c(0.95,0.9,0.85,0.8,0.75,0.7)


infection.tr<-TDTRest(fit=sepsis_fit,Data=sepsis,time=seq(0.01,90,0.01),cov0=X.ref,cov1=X.index,pct=pcttime)

X.index<-X.ref
X.index[4]<-1
cirrhosis.tr<-TDTRest(fit=sepsis_fit,Data=sepsis,time=seq(0.01,90,0.01),cov0=X.ref,cov1=X.index,pct=pcttime)


X.ref<-c(0,0,0,0,0,65,11)
X.index<-c(0,0,0,0,0,65,12)
sofa.tr<-TDTRest(fit=sepsis_fit,Data=sepsis,time=seq(0.01,90,0.01),cov0=X.ref,cov1=X.index,pct=pcttime)


X.index<-X.ref
X.index[6]<-70
age.tr<-TDTRest(fit=sepsis_fit,Data=sepsis,time=seq(0.01,90,0.01),cov0=X.ref,cov1=X.index,pct=pcttime)


###linear predictor
predlp<-pred.lp(fit=sepsis_fit,Data=sepsis,time=10,newdata=sepsis)

#Loglikelihood
loglik<-AFTLoglike(fit=sepsis_fit,Data=sepsis,newdata=sepsis)

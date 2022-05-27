################################################################################
#Date: Oct 26 2021
#Project: CARE Part 3 Tuberculosis project
#Author: Kuan Peng
################################################################################

#Load required package 
pacman::p_load(data.table,tidyverse,tableone,writexl,readxl,ggplot2,WeightIt,mlogit,survey,lubridate,AER,broom,MatchIt)

#Main analysis IP primary----

TB <- readRDS("TB_all.RDS")
#correctify the datastructure
fac <-subset(colnames(TB),str_detect(colnames(TB),"dx|rx|vaccinated|sex|TB.bet|TB.bet.w"))
TB[,(fac):=lapply(.SD,as.factor),.SDcol=fac]
str(TB)#check variable class

TB[Vaccine.Brand.1st=="BioNTech/Fosun"&Vaccine.Brand.2nd=="BioNTech/Fosun"&vaccinated==1,vaccine:="BNT162b2"]
TB[Vaccine.Brand.1st=="Sinovac"&Vaccine.Brand.2nd=="Sinovac"&vaccinated==1,vaccine:="CoronaVac"]
TB[vaccinated==0,vaccine:="None"]

#Adjust previous COVID-19 history
COVID <- readRDS("D:/Batch 5/R/LAB_ALL_COVID.RDS")# Identify COVID-19 case using 213,215,216,217
setDT(COVID)
covid <- COVID[str_detect(T_NUM,"^213|^215|^216|^217")&result=="detected"]
setorder(covid,date,by=patient_pssn)
covid <- covid[,head(.SD,1),by=patient_pssn]
covid <- left_join(covid,TB[,.(patient_pssn,index.date2)])
cohis <- covid[date<=index.date2]$patient_pssn
TB[,covid_his:=fifelse(patient_pssn%in%cohis,1,0)]
coaft <- covid[date>index.date2]
setnames(coaft,"date","date.cov")
TB <- left_join(TB,coaft[,.(patient_pssn,date.cov)])
# Add risk window 14 days 

TB[TB==1,onset_time:=as_date(date.x)-as_date(index.date2)]
TB[onset_time>=30,date.z:=date.x]
TB[,TB.z:=as.numeric(onset_time>=30)]
TB[is.na(TB.z),TB.z:=0]

#For those unvaccinated before 2021/9/30 but vaccinated later, we will need to censor them 
TB[vaccinated==0&!is.na(Vaccine.Brand.1st),date.re:=Date.of.vaccination.1st]
#For those vaccinated two dose before 2021/9/30, we will censor them if they receive third dose
# TB[vaccinated==1&!is.na(Date.of.vaccination.3rd),date.3rd:=Date.of.vaccination.3rd]
TB[,censordate:=as_date("2022-1-31")]
TB[,censor1:=pmin(as_date(censordate),as_date(date.x),as_date(death_date_ymd),as_date(date.re),as_date(date.cov),na.rm = T)]
TB[,censor2:=pmin(as_date(censordate),as_date(date.y),as_date(death_date_ymd),as_date(date.re),as_date(date.cov),na.rm = T)]
TB[,censor3:=pmin(as_date(censordate),as_date(death_date_ymd),as_date(date.re),as_date(date.cov),na.rm = T)]
TB[,censor4:=pmin(as_date(censordate),as_date(date.ap),as_date(death_date_ymd),as_date(date.re),as_date(date.cov),na.rm = T)]
TB[,censor5:=pmin(as_date(censordate),as_date(date.z),as_date(death_date_ymd),as_date(date.re),as_date(date.cov),na.rm = T)]



TB[,time2nd_tb:=as.integer(difftime(censor1,index.date2,units = "days")+1)]
TB[,time2nd_tb.w:=as.integer(difftime(censor2,index.date2,units = "days")+1)]
TB[,time2nd_de:=as.integer(difftime(censor3,index.date2,units = "days")+1)]
TB[,time2nd_ap:=as.integer(difftime(censor4,index.date2,units = "days")+1)]
TB[,time2nd_risktb:=as.integer(difftime(censor5,index.date2,units = "days")+1)]



TB$vaccine <- factor(TB$vaccine, levels = c("None", "BNT162b2", "CoronaVac"))
TB[death_date_ymd>index.date2,death:=1]
TB[is.na(death),death:=0]
TB$covid_his <- as.factor(TB$covid_his)
# all
tb.w <- weightit(vaccine~Age+sex+
                   dx.stroke_embo+ dx.asthma+ dx.infection_resp+
                   dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
                   ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,
                 data = TB,
                 method = "ps",
                 estimand = "ATT",
  focal = "None",
  use.mlogit = FALSE)
w.test <- trim(tb.w,at = .99,lower = TRUE )
summary(w.test)



tb.tableone <- copy(TB)[,.(vaccine,Age,sex,dx.mi,dx.pvd,dx.cbd,dx.copd,dx.dementia,dx.paralysis,
                           dx.liver_mild,dx.liver_modsev,dx.ulcers,dx.ra,dx.cancer,dx.dm_com0,dx.dm_com1,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,
                           dx.infection_viral,dx.sle,dx.spa,dx.psa,dx.ms,dx.ibd,
                           rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,ip,ae,sopc,gopc,cci,covid_his,TB.bet)]
s1 <- print(CreateTableOne(data=tb.tableone,strata = "vaccine",test = FALSE),
        showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)

weighted.trim <- svydesign(ids=~1,data=tb.tableone,weights=~w.test$weights)
weighted2.trim <- print(svyCreateTableOne(data=weighted.trim,strata = "vaccine",test = FALSE),
                   showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)

bw <- data.table(s1,keep.rownames = T)
aw <- data.table(weighted2.trim,keep.rownames=T)
write_xlsx(list("beforeweighting"=bw,"afterweighting"=aw),"table1.xlsx")


# >=60
tb.w.l4 <- weightit(vaccine~Age+sex+ 
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,
                 data = TB[Age>=60],
  method = "ps",
  estimand = "ATT",
  focal = "None",
  use.mlogit = FALSE)
w.test.l4 <- trim(tb.w.l4,at = .99,lower = TRUE )

tb.tableone.l4 <- copy(TB[Age>=60])[,.(vaccine,Age,sex,dx.mi,dx.pvd,dx.cbd,dx.copd,dx.dementia,dx.paralysis,
  dx.liver_mild,dx.liver_modsev,dx.ulcers,dx.ra,dx.cancer,dx.dm_com0,dx.dm_com1,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,
  dx.infection_viral,dx.sle,dx.spa,dx.psa,dx.ms,dx.ibd,
  rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,ip,ae,sopc,gopc,cci,covid_his,TB.bet)]

weighted.l4.trim <- svydesign(ids=~1,data=tb.tableone.l4,weights=~w.test.l4$weights)
weighted2.l4.trim <- print(svyCreateTableOne(data=weighted.l4.trim,strata = "vaccine",test = FALSE),
                      showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)




# tb.w.l4.w <- weightit(vaccine~Age+sex+ dx.chf+
#     dx.stroke_embo+ dx.asthma+ dx.infection_resp+
#     dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
#     ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+covid_his+TB.bet.w,
#   data = TB[Age>=60],
#   method = "ps",
#   estimand = "ATT",
#   stabilize = FALSE,
#   focal = "None",
#   by = NULL,
#   s.weights = NULL,
#   moments = NULL,
#   int = FALSE,
#   subclass = NULL,
#   missing = NULL,
#   verbose = FALSE,
#   include.obj = FALSE,
#   use.mlogit=FALSE)
# w.test.l4.w <- trim(tb.w.l4.w,at = .99,lower = TRUE )

#<60
tb.w.l3 <- weightit(vaccine~Age+sex+ 
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,
                    data = TB[Age<60],
  method = "ps",
  estimand = "ATT",
  focal = "None",
  use.mlogit=FALSE)
w.test.l3 <- trim(tb.w.l3,at = .99,lower = TRUE )
summary(w.test.l3)

tb.tableone.l3 <- copy(TB[Age<60])[,.(vaccine,Age,sex,dx.mi,dx.pvd,dx.cbd,dx.copd,dx.dementia,dx.paralysis,
  dx.liver_mild,dx.liver_modsev,dx.ulcers,dx.ra,dx.cancer,dx.dm_com0,dx.dm_com1,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,
  dx.infection_viral,dx.sle,dx.spa,dx.psa,dx.ms,dx.ibd,
  rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,ip,ae,sopc,gopc,cci,covid_his,TB.bet)]

weighted.l3.trim <- svydesign(ids=~1,data=tb.tableone.l3,weights=~w.test.l3$weights)
weighted2.l3.trim <- print(svyCreateTableOne(data=weighted.l3.trim,strata = "vaccine",test = FALSE),
                           showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)


#Male

tb.w.l2 <- weightit(vaccine~Age+
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,
  data = TB[sex=="M"],
  method = "ps",
  estimand = "ATT",
  focal = "None",
  use.mlogit=FALSE)
w.test.l2 <- trim(tb.w.l2,at = .99,lower = TRUE )
summary(w.test.l2)

tb.tableone.l2 <- copy(TB[sex=="M"])[,.(vaccine,Age,dx.mi,dx.pvd,dx.cbd,dx.copd,dx.dementia,dx.paralysis,
  dx.liver_mild,dx.liver_modsev,dx.ulcers,dx.ra,dx.cancer,dx.dm_com0,dx.dm_com1,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,
  dx.infection_viral,dx.sle,dx.spa,dx.psa,dx.ms,dx.ibd,
  rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,ip,ae,sopc,gopc,cci,covid_his,TB.bet)]

weighted.l2.trim <- svydesign(ids=~1,data=tb.tableone.l2,weights=~w.test.l2$weights)
weighted2.l2.trim <- print(svyCreateTableOne(data=weighted.l2.trim,strata = "vaccine",test = FALSE),
  showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)


#Female

tb.w.l1 <- weightit(vaccine~Age+ 
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,
  data = TB[sex=="F"],
  method = "ps",
  estimand = "ATT",
  focal = "None",
  use.mlogit=FALSE)
w.test.l1 <- trim(tb.w.l1,at = .99,lower = TRUE )
summary(w.test.l1)

tb.tableone.l1 <- copy(TB[sex=="F"])[,.(vaccine,Age,dx.mi,dx.pvd,dx.cbd,dx.copd,dx.dementia,dx.paralysis,
  dx.liver_mild,dx.liver_modsev,dx.ulcers,dx.ra,dx.cancer,dx.dm_com0,dx.dm_com1,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,
  dx.infection_viral,dx.sle,dx.spa,dx.psa,dx.ms,dx.ibd,
  rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,ip,ae,sopc,gopc,cci,covid_his,TB.bet)]

weighted.l1.trim <- svydesign(ids=~1,data=tb.tableone.l1,weights=~w.test.l1$weights)
weighted2.l1.trim <- print(svyCreateTableOne(data=weighted.l1.trim,strata = "vaccine",test = FALSE),
  showAllLevels = F, quote = T, noSpaces = T,smd=TRUE)




#additionally adjusted variable with SMD > 0.100

m1 <- coxph(Surv(time2nd_tb,TB)~vaccine+Age,data=TB,weights =w.test$weights)
m2 <- coxph(Surv(time2nd_ap,ap)~vaccine+Age,data=TB,weights =w.test$weights)
# m3 <- coxph(Surv(time2nd_de,death)~vaccine+Age,data=TB,weights =w.test$weights)
m4 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine+Age,data=TB,weights =w.test$weights)
m31 <- coxph(Surv(time2nd_risktb,TB.z)~vaccine+Age,data=TB,weights =w.test$weights)

m5 <- coxph(Surv(time2nd_tb,TB)~vaccine+ Age,data=TB[Age>=60],weights =w.test.l4$weights)
# m6 <- coxph(Surv(time2nd_ap,ap)~vaccine+Age,data=TB[Age>=60],weights =w.test.l4$weights)
# m7 <- coxph(Surv(time2nd_de,death)~vaccine+Age,data=TB[Age>=60],weights =w.test.l4$weights)
# m8 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine+Age,data=TB[Age>=60],weights =w.test.l4$weights)

m9 <- coxph(Surv(time2nd_tb,TB)~vaccine,data=TB[Age<60],weights =w.test.l3$weights)
# m10 <- coxph(Surv(time2nd_ap,ap)~vaccine,data=TB[Age<60],weights =w.test.l3$weights)
# m11 <- coxph(Surv(time2nd_de,death)~vaccine,data=TB[Age<60],weights =w.test.l3$weights)
# m12 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine,data=TB[Age<60],weights =w.test.l3$weights)

m13 <- coxph(Surv(time2nd_tb,TB)~vaccine,data=TB[sex=="M"],weights =w.test.l2$weights)
# m14 <- coxph(Surv(time2nd_ap,ap)~vaccine,data=TB[sex=="M"],weights =w.test.l2$weights)
# m15 <- coxph(Surv(time2nd_de,death)~vaccine,data=TB[sex=="M"],weights =w.test.l2$weights)
# m16 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine,data=TB[sex=="M"],weights =w.test.l2$weights)

m17 <- coxph(Surv(time2nd_tb,TB)~vaccine+Age,data=TB[sex=="F"],weights =w.test.l1$weights)
# m18 <- coxph(Surv(time2nd_ap,ap)~vaccine+Age,data=TB[sex=="F"],weights =w.test.l1$weights)
# m19 <- coxph(Surv(time2nd_de,death)~vaccine+Age,data=TB[sex=="F"],weights =w.test.l1$weights)
# m20 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine+Age,data=TB[sex=="F"],weights =w.test.l1$weights)





check1 <- cox.zph(m1,terms = F)
cox.zph(m2)
cox.zph(m4)
check2 <- cox.zph(m5,terms = F)
check3 <- cox.zph(m9,terms = F)
check4 <- cox.zph(m13,terms = F)
check5 <- cox.zph(m17,terms = F)
cox.zph(m21)
cox.zph(m22)
cox.zph(m23)
cox.zph(m31)


plot(check1,hr=T)
plot(check2,hr=T)
plot(check3,hr=T)
plot(check4,hr=T)

jpeg("All.jpg",width = 2400,height = 1600,res=300)
par(mfrow=c(1,2))
plot(check1,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineBNT162b2",se = F)
plot(check1,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineCoronaVac",se = F)
dev.off()

jpeg("Old.jpg",width = 2400,height = 1600,res=300)
par(mfrow=c(1,2))
plot(check2,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineBNT162b2",se = F)
plot(check2,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineCoronaVac",se = F)
dev.off()

jpeg("Adult.jpg",width = 2400,height = 1600,res=300)
par(mfrow=c(1,2))
plot(check3,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineBNT162b2",se = F)
plot(check3,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineCoronaVac",se = F)
dev.off()


jpeg("Male.jpg",width = 2400,height = 1600,res=300)
par(mfrow=c(1,2))
plot(check4,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineBNT162b2",se = F)
plot(check4,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineCoronaVac",se = F)
dev.off()

jpeg("Female.jpg",width = 2400,height = 1600,res=300)
par(mfrow=c(1,2))
plot(check5,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineBNT162b2",se = F)
plot(check5,hr=T,col = "lightseagreen",lty = c(1,2),lwd = c(2,1),nsmo = 100,xlab="Time after complete 2nd dose (days)",var = "vaccineCoronaVac",se = F)
dev.off()


# TB2 <- TB[Age<60]
# TB2$weightsl3 <- w.test.l3$weights
# TB2 <- TB2[!vaccine=="BNT162b2"]
# TB2$vaccine <- factor(TB2$vaccine, levels = c("None", "CoronaVac"))
# fit <- survreg(Surv(time2nd_tb,TB)~vaccine,data=TB2,dist="weibull",weights =weightsl3)

# tb2 <- survSplit(Surv(time2nd_tb,TB)~vaccine+weightsl3,data=TB2,cut = c(90,180),episode = "tgroup",id = "id")
# tb3 <- survSplit(Surv(time2nd_tb,TB)~vaccine+weightsl3,data=TB2,cut = c(50,100,150),episode = "tgroup",id = "id")
# tb2 <- data.table(tb2)
# tb2 <- tb2[!vaccine=="BNT162b2"]
# tb2$vaccine <- factor(tb2$vaccine, levels = c("None", "CoronaVac"))
# tb3 <- data.table(tb3)
# tb3 <- tb3[!vaccine=="BNT162b2"]
# tb3$vaccine <- factor(tb3$vaccine, levels = c("None", "CoronaVac"))
# 
# m21 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine*strata(tgroup),data=tb2,weights =weightsl3)
# m21.1 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb2[time2nd_tb<91],weights =weightsl3)
# m21.2 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb2[time2nd_tb<181],weights =weightsl3)
# cox.zph(m21.1)


# m22 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine*strata(tgroup),data=tb3,weights =weightsl3)
m23 <- coxph(Surv(time2nd_tb,TB)~vaccine+ strata(Age),data=TB[Age>=60],weights =w.test.l4$weights)
# m24 <- coxph(Surv(time2nd_tb.w,TB.w)~vaccine+strata(Age),data=TB[Age>=60],weights =w.test.l4$weights)

# tb3 <- data.table(tb3)
# m22.1 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb3[time2nd_tb<51],weights =weightsl3)
# m22.2 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb3[time2nd_tb<101&time2nd_tb>=51],weights =weightsl3)
# m22.3 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb3[time2nd_tb<151&time2nd_tb>=101],weights =weightsl3)
# m22.4 <- coxph(Surv(tstart,time2nd_tb,TB)~vaccine,data=tb3[time2nd_tb>=151],weights =weightsl3)
# tidy(m22.1) %>% mutate(exp.estimate=exp(estimate),.after=estimate) %>% mutate(upper.CI=exp(estimate+1.96*std.error),.after=std.error) %>% mutate(lower.CI=exp(estimate-1.96*std.error),.after=std.error)
# tidy(m22.2) %>% mutate(exp.estimate=exp(estimate),.after=estimate) %>% mutate(upper.CI=exp(estimate+1.96*std.error),.after=std.error) %>% mutate(lower.CI=exp(estimate-1.96*std.error),.after=std.error)
# tidy(m22.3) %>% mutate(exp.estimate=exp(estimate),.after=estimate) %>% mutate(upper.CI=exp(estimate+1.96*std.error),.after=std.error) %>% mutate(lower.CI=exp(estimate-1.96*std.error),.after=std.error)
# tidy(m22.4) %>% mutate(exp.estimate=exp(estimate),.after=estimate) %>% mutate(upper.CI=exp(estimate+1.96*std.error),.after=std.error) %>% mutate(lower.CI=exp(estimate-1.96*std.error),.after=std.error)


# check the competiting risk of death 
etime <- with(TB[!duplicated(patient_pssn)],ifelse(TB==0,time2nd_de,time2nd_tb))
event <- with(TB[!duplicated(patient_pssn)],ifelse(TB==0,2*death,1))
event <- factor(event,0:2,labels=c("censor","tb","death"))
# cfit1 <- coxph(Surv(etime,event)~vaccine+Age,data=TB[!duplicated(patient_pssn)],id=patient_pssn) #too slow
tbdat <- finegray(Surv(etime,event)~vaccine+Age+sex+
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=TB[!duplicated(patient_pssn)],etype = "tb")
tbdat <- data.table(tbdat)
fgfit <- coxph(Surv(fgstart,fgstop,fgstatus)~vaccine+Age+ 
    dx.stroke_embo+ dx.asthma+ dx.infection_resp+
    dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
    ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=tbdat,weight=fgwt)
summary(fgfit)

# fgfit2 <- coxph(Surv(fgstart,fgstop,fgstatus)~vaccine+Age+ 
#     dx.stroke_embo+ dx.asthma+ dx.infection_resp+
#     dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
#     ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=tbdat[Age>=60],weight=fgwt)
# summary(fgfit2)
# 
# fgfit3 <- coxph(Surv(fgstart,fgstop,fgstatus)~vaccine+Age+ 
#     dx.stroke_embo+ dx.asthma+ dx.infection_resp+
#     dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
#     ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=tbdat[Age<60],weight=fgwt)
# summary(fgfit3)
# 
# fgfit4 <- coxph(Surv(fgstart,fgstop,fgstatus)~vaccine+Age+ 
#     dx.stroke_embo+ dx.asthma+ dx.infection_resp+
#     dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
#     ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=tbdat[sex=="M"],weight=fgwt)
# summary(fgfit4)
# 
# fgfit5 <- coxph(Surv(fgstart,fgstop,fgstatus)~vaccine+Age+ 
#     dx.stroke_embo+ dx.asthma+ dx.infection_resp+
#     dx.infection_viral+dx.sle +dx.spa+dx.psa+dx.ibd+dx.ms+
#     ip+ae+sopc+gopc+cci+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+TB.bet+covid_his,data=tbdat[sex=="F"],weight=fgwt)
# summary(fgfit5)


# get_fk_result <- function(x,y,z){
#   temp <- merge(z[get(x)==1,.(case=.N),by=c("vaccine")],
#                 z[,.(total=.N),by=c("vaccine")])
#   temp <- merge(temp,
#                 z[get(x)==1,.(median_FU=paste0(median(get(y))," (",paste(quantile(get(y),0.25),quantile(get(y),0.75),sep = ", "),")")),vaccine])
#   
#   temp <- merge(temp,
#                 z[,.(fu_py=round(sum(get(y))/365.25,1)),by=vaccine])
#   
#   
#   frac <- left_join(z[get(x)==1,.N,by=c("vaccine")],
#                     z[,round(sum(get(y))/365.25,1),by=vaccine])
#   
#   ci_dt <- data.table()
#   for (i in 1:3) {
#     n <- frac$N[i]
#     d <- frac$V1[i]
#     fit <- glm(n ~ offset(log(d)), family=poisson)
#     vaccine <- frac$vaccine[i]
#     est_95 <- paste0(round(n/d*10000,2),"(",round(exp(confint(fit))[1]*10000,2),", ",round(exp(confint(fit))[2]*10000,2),")")
#     if(i ==1){
#       ci_dt <- data.table(vaccine,est_95)
#     }else{
#       ci_dt <- rbind(ci_dt,data.table(vaccine,est_95))
#     }
#   }
#   temp <- merge(temp,ci_dt)
#   return(temp)
# }

get_fk_result <- function(x,y,z){
  temp <- merge(z[get(x)==1,.(case=.N),by=c("vaccine")],
    z[,.(total=.N),by=c("vaccine")])
  temp <- merge(temp,
    z[,.(median_FU=paste0(median(get(y))," (",paste(quantile(get(y),0.25),quantile(get(y),0.75),sep = ", "),")")),vaccine])
  
  temp <- merge(temp,
    z[,.(fu_py=round(sum(get(y))/365.25,1)),by=vaccine])
  
  
  frac <- left_join(z[get(x)==1,.N,by=c("vaccine")],
    z[,round(sum(get(y))/365.25,1),by=vaccine])
  
  ci_dt <- data.table()
  for (i in 1:3) {
    n <- frac$N[i]
    d <- frac$V1[i]
    fit <- glm(n ~ offset(log(d)), family=poisson)
    vaccine <- frac$vaccine[i]
    est_95 <- paste0(round(n/d*10000,2),"(",round(exp(confint(fit))[1]*10000,2),", ",round(exp(confint(fit))[2]*10000,2),")")
    if(i ==1){
      ci_dt <- data.table(vaccine,est_95)
    }else{
      ci_dt <- rbind(ci_dt,data.table(vaccine,est_95))
    }
  }
  temp <- merge(temp,ci_dt)
  return(temp)
}


g1 <- rbind(get_fk_result("TB","time2nd_tb",TB),get_fk_result("TB.z","time2nd_risktb",TB),get_fk_result("ap","time2nd_ap",TB),get_fk_result("death","time2nd_de",TB),get_fk_result("TB.w","time2nd_tb.w",TB))
g2 <- rbind(get_fk_result("TB","time2nd_tb",TB[Age>=60]),get_fk_result("ap","time2nd_ap",TB[Age>=60]),get_fk_result("death","time2nd_de",TB[Age>=60]),get_fk_result("TB.w","time2nd_tb.w",TB[Age>=60]))
g3 <- rbind(get_fk_result("TB","time2nd_tb",TB[Age<60]),get_fk_result("ap","time2nd_ap",TB[Age<60]),get_fk_result("death","time2nd_de",TB[Age<60]),get_fk_result("TB.w","time2nd_tb.w",TB[Age<60]))
g4 <- rbind(get_fk_result("TB","time2nd_tb",TB[sex=="M"]),get_fk_result("ap","time2nd_ap",TB[sex=="M"]),get_fk_result("death","time2nd_de",TB[sex=="M"]),get_fk_result("TB.w","time2nd_tb.w",TB[sex=="M"]))
g5 <- rbind(get_fk_result("TB","time2nd_tb",TB[sex=="F"]),get_fk_result("ap","time2nd_ap",TB[sex=="F"]),get_fk_result("death","time2nd_de",TB[sex=="F"]),get_fk_result("TB.w","time2nd_tb.w",TB[sex=="F"]))
# g6 <- rbind(get_fk_result("TB","time2nd_tb",tb3[tgroup=="1"]),get_fk_result("TB","time2nd_tb",tb3[tgroup=="2"]),get_fk_result("TB","time2nd_tb",tb3[tgroup=="3"]),get_fk_result("TB","time2nd_tb",tb3[tgroup=="4"]))

# g <- list(g1,g2,g3,g4,g5,g6)
g <- list(g1,g2,g3,g4,g5)

write_xlsx(g,"5.table2.xlsx")

r1 <- merge( tidy(m1) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m1)),keep.rownames=T),by.x = "term",by.y = "rn")
r5 <- merge( tidy(m5) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m5)),keep.rownames=T),by.x = "term",by.y = "rn")
r9 <- merge( tidy(m9) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m9)),keep.rownames=T),by.x = "term",by.y = "rn")
r13 <- merge( tidy(m13) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m13)),keep.rownames=T),by.x = "term",by.y = "rn")
r17 <- merge( tidy(m17) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m17)),keep.rownames=T),by.x = "term",by.y = "rn")
r2 <- merge( tidy(m2) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m2)),keep.rownames=T),by.x = "term",by.y = "rn")
r4 <- merge( tidy(m4) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m4)),keep.rownames=T),by.x = "term",by.y = "rn")
r23 <- merge( tidy(m23) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m23)),keep.rownames=T),by.x = "term",by.y = "rn")
r25 <- merge( tidy(fgfit) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(fgfit)),keep.rownames=T),by.x = "term",by.y = "rn")
r27 <- merge( tidy(m22.1) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m22.1)),keep.rownames=T),by.x = "term",by.y = "rn")
r28 <- merge( tidy(m22.2) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m22.2)),keep.rownames=T),by.x = "term",by.y = "rn")
r29 <- merge( tidy(m22.3) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m22.3)),keep.rownames=T),by.x = "term",by.y = "rn")
r30 <- merge( tidy(m22.4) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m22.4)),keep.rownames=T),by.x = "term",by.y = "rn")
r31 <- merge( tidy(m31) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(m31)),keep.rownames=T),by.x = "term",by.y = "rn")

# r <- list("all"=r1,">=60"=r5,"<60"=r9,"male" = r13,"female"= r17, "ap"=r2, "tb.sec"=r4, "tb.strata"= r23,"competing risk"= r25,"group1"= r27,"group2"= r28,"group3"= r29,"group4"= r30,"window30"= r31)
r <- list("all"=r1,">=60"=r5,"<60"=r9,"male" = r13,"female"= r17, "ap"=r2, "tb.sec"=r4, "tb.strata"= r23,"competing risk"= r25,"window30"= r31)

write_xlsx(r,"5.all.trimed.xlsx")



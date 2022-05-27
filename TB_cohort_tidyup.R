library(data.table)
library(readxl)
library(tidyverse)
library(tableone)
library(writexl)
library(lubridate)
memory.limit(1000000000000)
cohort.matched <- readRDS("cohort_matched_kuan_TB.RDS")
cohort.death <- cohort.matched[death_date_ymd<index.date2]#18475 #remove death before index date
cohort.matched <- cohort.matched[!patient_pssn%in%cohort.death$patient_pssn]#4068145
cohort.matched[,Age:=Age-1]
cohort.matched <- cohort.matched[Age>=18]#3954539
# DX
load("D:/Batch 5/R/DX.RData")

#remove TB and metastatic cancer diagnosis before index date 1 
TB1 <- unique(dx_clean[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^01[0-8]")]$patient_pssn)
dx_latest_index <- merge(dx_latest, cohort.matched[, .(patient_pssn, index.date)], by="patient_pssn")
TB2 <- unique(dx_latest_index[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^01[0-8]")&as_date(date)<as_date(index.date)]$patient_pssn)
TB <- union(TB1,TB2)  #6974

metastatic1 <-unique(dx_clean[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^19[6-9]")]$patient_pssn) 
metastatic2 <- unique(dx_latest_index[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^19[6-9]")&as_date(date)<as_date(index.date)]$patient_pssn)
metastatic <- union(metastatic1,metastatic2)#14734 

cohort.matched <- cohort.matched[!patient_pssn%in%TB]#3947565
cohort.matched <- cohort.matched[!patient_pssn%in%metastatic]#3932982

#remove TB related drugs before index date 1

drug_clean <- readRDS("D:/Batch 5/R/RX_clean.RDS")
drug_latest <- readRDS("D:/Batch 5/R/RX_latest.RDS")
drug_clean_combined <- rbind(drug_clean$`2018`,drug_clean$`2019`,drug_clean$`2020`)


drug_latest_cohort <- drug_latest[patient_pssn%in%cohort.matched$patient_pssn]
drug_clean_cohort <- drug_clean_combined[patient_pssn%in%cohort.matched$patient_pssn]
drug_clean_cohort$presc_start_date_ym <-  as_date(paste0(drug_clean_cohort$presc_start_date_ym,"-15"))
drug_clean_cohort[,presc_end:=presc_start_date_ym+presc_duration_day]
drug_latest_cohort[,presc_end:=as.Date(presc_start_date_ymd)+presc_duration_day]
drug_latest_cohort <-left_join(drug_latest_cohort,cohort.matched[,.(patient_pssn,index.date)])


tbrx_before <- drug_latest_cohort[str_detect(bnfno_p,"^5.1.9")&as_date(disp_date_ymd)<index.date]
tbrx_before <- tbrx_before[str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05",ignore_case = T))]
tbrx1 <- unique(tbrx_before$patient_pssn)
tbrx_before2 <- drug_clean_cohort[str_detect(bnfno_p,"^5.1.9")]
tbrx_before2 <- tbrx_before2[str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05",ignore_case = T))]
tbrx2 <- unique(tbrx_before2$patient_pssn)
tbrx <- union(tbrx1,tbrx2)#

cohort.matched <- cohort.matched[!patient_pssn%in%tbrx]#3927864

#exclude subjects with hospital admission start within 30 days before index date and discharge date after index date
# load("D:/Batch 5/R/A10X.RData")
# iptb <- admin_latest[patient_pssn%in%cohort.matched$patient_pssn]
# iptb <- left_join(iptb,cohort.matched[,.(patient_pssn, index.date)])
# iptb$date <- as_date(iptb$date)
# iptb$dischg_date_ymd <- as_date(iptb$dischg_date_ymd)
# iptb <-iptb[date>=(index.date-30)&date<index.date&dischg_date_ymd>index.date]
# cohort.matched <- cohort.matched[!patient_pssn%in%iptb$patient_pssn]#3942604  removed 8647

dx_latest_index2 <- merge(dx_latest, cohort.matched[, .(patient_pssn, index.date2)], by="patient_pssn")
dx.case <- dx_latest_index2[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^01[0-8]")][str_detect(ranking,"P")]
dx.case <- left_join(dx.case,cohort.matched[, .(patient_pssn, death_date_ymd,index.date)])
dx.case <- dx.case[date>=index.date]
setorder(dx.case,date,by=patient_pssn)
dx.case <- dx.case[,head(.SD,1),by=patient_pssn]#keep the earliest records


drug_latest_tb <- drug_latest[patient_pssn%in%dx.case$patient_pssn][str_detect(bnfno_p,"^5.1.9")][str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05|^pyra0[346]|^etha0[239]|^etha11|^stre01",ignore_case = T))]
drug_latest_tb <- left_join(drug_latest_tb,dx.case[,.(patient_pssn,date)])
drug_latest_tb$date <- as_date(drug_latest_tb$date)
drug_latest_tb$disp_date_ymd <- as_date(drug_latest_tb$disp_date_ymd)
drug_latest_tb <- drug_latest_tb[( disp_date_ymd<date+14)&disp_date_ymd>=date]
drug_latest_tb <- drug_latest_tb[str_detect(item_cd,"^rifi0[12]|^rifa05"),item_cd:="RIFI"]
drug_latest_tb[,item_cd:=substr(item_cd,1,4)]
tb_px_n <-distinct(drug_latest_tb,patient_pssn,item_cd,.keep_all = T)


dx.case[,PYRA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"PYRA")]$patient_pssn,1,0)]
dx.case[,ETHA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"ETHA")]$patient_pssn,1,0)]
dx.case[,RIFA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"RIFA")]$patient_pssn,1,0)]
dx.case[,ISON:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"ISON")]$patient_pssn,1,0)]
dx.case[,STRE:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"STRE")]$patient_pssn,1,0)]
dx.case[,RIFI:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"RIFI")]$patient_pssn,1,0)]
dx.case.v <- rbind(dx.case[PYRA==1&ISON==1&RIFA==1][STRE==1|ETHA==1],dx.case[PYRA==1&RIFI==1][STRE==1|ETHA==1])#532 out of 913
# rbind(dx.case[ISON==1&RIFA==1],dx.case[RIFI==1])

dx.case.1.v <- dx.case.v[is.na(death_date_ymd)][date>=index.date2&date<=as.Date("2022-1-31")][str_detect(Source,"1.IP")&ranking=="P"]#442
dx.case.2.v <- dx.case.v[!is.na(death_date_ymd)][date>=index.date2&date<=death_date_ymd][str_detect(Source,"1.IP")&ranking=="P"]#32
dx.case.new.v <- rbind(dx.case.1.v,dx.case.2.v)
dx.case.3.v <- dx.case.v[date<index.date2&date>=index.date][str_detect(Source,"1.IP")&ranking=="P"]
newcase.v <- dx.case.new.v$patient_pssn
cohort.matched[,TB:=ifelse(patient_pssn%in%newcase.v,1,0)]#474

# ip_admin_latest <- readRDS("D:/Batch 5/R/ip/ip_admin_latest.RDS")
# ip_admin_latest <- data.table(ip_admin_latest)
# k <- left_join(dx.case,ip_admin_latest[,.(patient_pssn,adm_date_ymd,dischg_date_ymd)])
# kk <- k[date==adm_date_ymd]
# kk <- kk[!(patient_pssn%in%kk[duplicated(patient_pssn)]$patient_pssn&adm_date_ymd==dischg_date_ymd)]
# kk[, LOS:=as.integer(as_date(dischg_date_ymd)-as_date(adm_date_ymd))]
# summary(kk[!patient_pssn%in%newcase.v,LOS])
# hist(kk[!patient_pssn%in%newcase.v,LOS],xlim = c(0,100),nclass = 500)
# kk[!patient_pssn%in%newcase.v&LOS<=14]
# summary(kk[,LOS])
# hist(kk[,LOS],xlim = c(0,100),nclass = 500)
# kk[LOS<=14]

dx.case.1 <- dx.case[is.na(death_date_ymd)][date>=index.date2&date<=as.Date("2022-1-31")][str_detect(Source,"1.IP")&ranking=="P"]#726
dx.case.2 <- dx.case[!is.na(death_date_ymd)][date>=index.date2&date<=death_date_ymd][str_detect(Source,"1.IP")&ranking=="P"]#90
dx.case.new <- rbind(dx.case.1,dx.case.2)
dx.case.3 <- dx.case[date<index.date2&date>=index.date][str_detect(Source,"1.IP")&ranking=="P"]
newcase <- dx.case.new$patient_pssn#816
cohort.matched[,TB.w:=ifelse(patient_pssn%in%newcase,1,0)]

cohort.matched <- left_join(cohort.matched,dx.case.new.v[,.(patient_pssn,date)],by="patient_pssn")
cohort.matched <- left_join(cohort.matched,dx.case.new[,.(patient_pssn,date)],by="patient_pssn")
# use incident MI as negative control
#remove case with MI history
# mihis <- dx_clean[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^410")]$patient_pssn
# 
# fx.case <- dx_latest_index2[patient_pssn%in%cohort.matched$patient_pssn][str_detect(code,"^410")]
# fx.case <- left_join(fx.case,cohort.matched[, .(patient_pssn, death_date_ymd,index.date)])
# fx.case <- fx.case[date>=index.date]
# setorder(fx.case,date,by=patient_pssn)
# fx.case <- fx.case[,head(.SD,1),by=patient_pssn]#keep the earliest records
# fx.case <-fx.case[!patient_pssn%in%mihis]
# 
# fx.case.1 <- fx.case[is.na(death_date_ymd)][date>=index.date2&date<=as.Date("2022-1-31")][str_detect(Source,"1.IP")&ranking=="P"]
# fx.case.2 <- fx.case[!is.na(death_date_ymd)][date>=index.date2&date<=death_date_ymd][str_detect(Source,"1.IP")&ranking=="P"]
# fx.case.new <- rbind(fx.case.1,fx.case.2)
# fx.case.3 <- fx.case[date<index.date2&date>=index.date][str_detect(Source,"1.IP")&ranking=="P"]
# fx.between <- fx.case.3[,head(.SD,1),by=patient_pssn]
# newcase.f <- fx.case.new$patient_pssn
# cohort.matched[,mi:=ifelse(patient_pssn%in%newcase.f,1,0)]
# cohort.matched <- left_join(cohort.matched,fx.case.new[,.(patient_pssn,date.mi=date)],by="patient_pssn")
# 

cohort.matched[,TB.bet:=ifelse(patient_pssn%in%dx.case.3.v$patient_pssn,1,0)]
cohort.matched <- left_join(cohort.matched,dx.case.3.v[,.(patient_pssn,date.bet=date)])
cohort.matched[,TB.bet.w:=ifelse(patient_pssn%in%dx.case.3$patient_pssn,1,0)]
cohort.matched <- left_join(cohort.matched,dx.case.3[,.(patient_pssn,date.bet.w=date)])



tb_ex <- cohort.matched
drug_clean_tb <- drug_clean_cohort
drug_latest_tb <- drug_latest_cohort
remove(drug_clean_cohort)
remove(drug_latest_cohort)
#1 exclusion----


drug_latest_tb [,index.date:=NULL]
colnames(drug_clean_tb) <- colnames(drug_latest_tb)
drug_latest_tb$presc_start_date_ymd <- as_date(drug_latest_tb$presc_start_date_ymd)
drug_latest_tb$disp_date_ymd <- as_date(drug_latest_tb$disp_date_ymd)
drug_clean_tb$disp_date_ymd <- as_date(paste0(drug_clean_tb$disp_date_ymd,"-15"))
drug_tb <- rbind(drug_clean_tb,drug_latest_tb)
drug_tb <- left_join(drug_tb,tb_ex[,.(patient_pssn,index.date)])
drug_tb$disp_date_ymd <- as_date(drug_tb$disp_date_ymd)
drug_tb[,presc_end:=disp_date_ymd +presc_duration_day]

setDT(drug_tb)

dx_clean_tb <- dx_clean[patient_pssn%in%tb_ex$patient_pssn]
dx_latest_tb <- dx_latest[patient_pssn%in%tb_ex$patient_pssn]

dx_clean_tb$date <- as_date(paste0(dx_clean_tb$date,"-15"))
dx_latest_tb$date <- as_date(dx_latest_tb$date)
dx_tb <- rbind(dx_clean_tb,dx_latest_tb)

# 2 generate variables ----

adj <- setDT(read_xlsx("codes_TB.xlsx",sheet = 2))
adj_drug <- adj[str_detect(Name,"rx")]
adj_dx <- adj[str_detect(Name,"dx")]

#dx variable

for (i in 1:length(adj_dx$Regex)) {
  dx.tb <- dx_tb[str_detect(code,regex(adj_dx$Regex[i],ignore_case = T))]
  dx.tb <- merge(dx.tb,tb_ex[,.(patient_pssn,index.date)])
  dx.tb <- dx.tb[as_date(date)<as_date(index.date)]$patient_pssn
  tb_ex[,adj_dx$Name[i]:=as.numeric(tb_ex$patient_pssn%in%dx.tb)]
  print(i)
}
# rx variable

for (i in 1:length(adj_drug$Regex)) {
  rx.drug <- drug_tb[str_detect(bnfno_p,regex(adj_drug$Regex[i],ignore_case = T))]
  rx.drug <- merge(rx.drug,tb_ex[,.(patient_pssn,index.date)])
  rx.drug <- rx.drug[disp_date_ymd<as_date(index.date)&presc_end>=as_date(index.date-180)]$patient_pssn
  tb_ex[,adj_drug$Name[i]:=as.numeric(tb_ex$patient_pssn%in%rx.drug)]
  print(i)
}

tb_ex[,cci:=dx.mi+dx.chf+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+(dx.dm_com0&!dx.dm_com1)+dx.dm_com1*2+dx.crf*2+(dx.liver_mild&!dx.liver_modsev)+dx.liver_modsev*3+dx.ulcers+dx.ra+dx.aids*6+dx.cancer*2+dx.cancer_mets*6]

adj_rx <- setDT(read_xlsx("CARE-TB-Crosscheck.xlsx",sheet = 5))
remove(rx.drug)
for (i in 1:length(adj_rx$Regex)) {
  rx.drug <- drug_tb[str_detect(item_cd,regex(adj_rx$Regex[i],ignore_case = T))]
  rx.drug <- merge(rx.drug,tb_ex[,.(patient_pssn,index.date)])
  rx.drug <- rx.drug[disp_date_ymd <as_date(index.date)&presc_end>=as_date(index.date-180)]$patient_pssn
  tb_ex[,adj_rx$Name[i]:=as.numeric(tb_ex$patient_pssn%in%rx.drug)]
  print(i)
}

# hospital services utilization within one year of diagnosis
load("D:/Batch 5/R/A10X.RData")
admin_clean_tb <- admin_clean[patient_pssn%in%tb_ex$patient_pssn]
admin_clean_tb[,date:=as_date(paste0(date,"-15"))]
admin_latest_tb <- admin_latest[patient_pssn%in%tb_ex$patient_pssn]
admin_latest_tb$date <- as_date(admin_latest_tb$date)
admin_combine <- rbind(admin_clean_tb[,dischg_date_ym:=NULL],admin_latest_tb[,dischg_date_ymd:=NULL])
admin_combine<- left_join(admin_combine,tb_ex[,.(patient_pssn,index.date)])
admin_tb <- admin_combine[date<index.date&date>as_date(index.date-365)]

admin_tb_ip <-admin_tb[,.N,by=c("patient_pssn","Source")][str_detect(Source,"1.IP"),.(ip=sum(N)),by=patient_pssn]
admin_tb_ae <-admin_tb[,.N,by=c("patient_pssn","Source")][str_detect(Source,"2.AE"),.(ae=sum(N)),by=patient_pssn]
admin_tb_gopc <-admin_tb[,.N,by=c("patient_pssn","Source")][str_detect(Source,"4.GO"),.(gopc=sum(N)),by=patient_pssn]
admin_tb_sopc <-admin_tb[,.N,by=c("patient_pssn","Source")][str_detect(Source,"3.SO"),.(sopc=sum(N)),by=patient_pssn]
tb_ex <- left_join(tb_ex,admin_tb_ae)
tb_ex <- left_join(tb_ex,admin_tb_ip)
tb_ex <- left_join(tb_ex,admin_tb_gopc)
tb_ex <- left_join(tb_ex,admin_tb_sopc)
# NA value means no ipae/op visits
tb_ex[is.na(sopc),sopc:=0]
tb_ex[is.na(ip),ip:=0]
tb_ex[is.na(gopc),gopc:=0]
tb_ex[is.na(ae),ae:=0]

# #use appendicitis as negative control
aphis <- dx_clean[patient_pssn%in%tb_ex$patient_pssn][str_detect(code,"^54[0-3]|^47[.]")]$patient_pssn
dx_latest_index2 <- merge(dx_latest, tb_ex[, .(patient_pssn, index.date,index.date2,death_date_ymd)], by="patient_pssn")

ap.case <- dx_latest_index2[patient_pssn%in%tb_ex$patient_pssn][str_detect(code,"^54[0-3]|^47[.]")]
ap.case <- ap.case[date>=index.date]
setorder(ap.case,date,by=patient_pssn)
ap.case <- ap.case[,head(.SD,1),by=patient_pssn]#keep the earliest records
ap.case <-ap.case[!patient_pssn%in%aphis]

ap.case.1 <- ap.case[is.na(death_date_ymd)][date>=index.date2&date<=as.Date("2022-1-31")][str_detect(Source,"1.IP")&ranking=="P"]
ap.case.2 <- ap.case[!is.na(death_date_ymd)][date>=index.date2&date<=death_date_ymd][str_detect(Source,"1.IP")&ranking=="P"]
ap.case.new <- rbind(ap.case.1,ap.case.2)


tb_ex[,ap:=ifelse(patient_pssn%in%ap.case.new$patient_pssn,1,0)]
tb_ex <- left_join(tb_ex,ap.case.new[,.(patient_pssn,date.ap=date)],by="patient_pssn")

saveRDS(tb_ex,"TB_all.RDS")#3927864

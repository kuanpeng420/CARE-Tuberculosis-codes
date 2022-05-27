#Load required package 
  pacman::p_load(data.table,tidyverse,writexl,readxl,lubridate)

# 1. Generating cohort ----

  memory.limit(1000000000000)
  #read combined DH data
  cohort <- readRDS("D:/Batch 5/R/4.cohort_full.RDS")
  #if dod is "" then NA
  cohort[death_date_ymd=="", death_date_ymd := NA]
  setDT(cohort)
  cohort[,Age:=Age-1]
  # 1. Identify cohort (all active patients in HA during study period) 
  colnames(cohort) <- make.names(colnames(cohort))
  cohort1 <- cohort %>% filter(!is.na(patient_pssn)) #4275346
  cohort_two_dif_vaccine <- cohort1 %>% 
    filter(!is.na(Vaccine.Brand.1st),!is.na(Vaccine.Brand.2nd),Vaccine.Brand.1st!=Vaccine.Brand.2nd) %>% 
    distinct(patient_pssn) 
  
  exclude_pssn<-cohort_two_dif_vaccine #98
  cohort2 <- cohort1 %>% 
    filter(!patient_pssn%in%exclude_pssn$patient_pssn) #4275248

  load("D:/Batch 5/R/DX.RData")
  dx_latest[,date:=as_date(date)]
  #exclude diagnosis of TB from any setting before 2.23 and prescription of any drugs
  TB1 <- unique(dx_clean[patient_pssn%in%cohort2$patient_pssn][str_detect(code,"^01[0-8]")]$patient_pssn)
  TB2 <- unique(dx_latest[patient_pssn%in%cohort2$patient_pssn][str_detect(code,"^01[0-8]")&as_date(date)<as_date("2021-2-23")]$patient_pssn)
  oldtb <- union(TB1,TB2)  #6878
  
  drug_clean <- readRDS("D:/Batch 5/R/RX_clean.RDS")
  drug_latest <- readRDS("D:/Batch 5/R/RX_latest.RDS")
  drug_clean_combined <- rbind(drug_clean$`2018`,drug_clean$`2019`,drug_clean$`2020`)
  drug_latest[,disp_date_ymd:=as_date(disp_date_ymd)]
  drug_latest[,presc_end_date_ymd:=disp_date_ymd+presc_duration_day]
  
  tbrx_before <- drug_latest[str_detect(bnfno_p,"^5.1.9")&disp_date_ymd<as.Date("2021-2-23")]
  tbrx_before <- tbrx_before[str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05",ignore_case = T))]
  tbrx1 <- unique(tbrx_before$patient_pssn)
  tbrx_before2 <- drug_clean_combined[str_detect(bnfno_p,"^5.1.9")]
  tbrx_before2 <- tbrx_before2[str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05",ignore_case = T))]
  tbrx2 <- unique(tbrx_before2$patient_pssn)
  tbrx <- union(tbrx1,tbrx2)#10529
  
  cohort3 <- cohort2[!(patient_pssn%in%oldtb)]#4268370
  cohort3 <- cohort3[!(patient_pssn%in%tbrx)]#4263220
  
  #exclude diagnosis of metastatic cancer  from any setting before 2.23
  metastatic1 <-unique(dx_clean[patient_pssn%in%cohort3$patient_pssn][str_detect(code,"^19[6-9]")]$patient_pssn) 
  metastatic2 <- unique(dx_latest[patient_pssn%in%cohort3$patient_pssn][str_detect(code,"^19[6-9]")&as_date(date)<as_date("2021-2-23")]$patient_pssn)
  metastatic <- union(metastatic1,metastatic2)#14872
  cohort3 <- cohort3[!(patient_pssn%in%metastatic)]#4248348
  
  #Identify new cases based on IP diagnosis primary
  tbnew <- dx_latest[patient_pssn%in%cohort3$patient_pssn][str_detect(code,"^01[0-8]")&as_date(date)>=as_date("2021-2-23")&Source=="1.IP"&ranking=="P"]
  tbnew <-tbnew[order(date),.SD,by=patient_pssn][,head(.SD,1),by=patient_pssn]
  #verify by rx records
  drug_latest_tb <- drug_latest[patient_pssn%in%tbnew$patient_pssn][str_detect(bnfno_p,"^5.1.9")][str_detect(item_cd,regex("^ison0[1-9]|^ISON1[012]|^rifa0[12389]|^RIFA1[16]|^rifi0[12]|^rifa05|^pyra0[346]|^etha0[239]|^etha11|^stre01",ignore_case = T))][disp_date_ymd>as.Date("2021-2-23")]
  drug_latest_tb <- left_join(drug_latest_tb,tbnew[,.(patient_pssn,date)])
  drug_latest_tb <- drug_latest_tb[(disp_date_ymd<date+14)&disp_date_ymd>=date]
  drug_latest_tb <- drug_latest_tb[str_detect(item_cd,"^rifi0[12]|^rifa05"),item_cd:="RIFI"]
  drug_latest_tb[,item_cd:=substr(item_cd,1,4)]
  tb_px_n <-distinct(drug_latest_tb,patient_pssn,item_cd,.keep_all = T)

  cohort3[,PYRA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"PYRA")]$patient_pssn,1,0)]
  cohort3[,ETHA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"ETHA")]$patient_pssn,1,0)]
  cohort3[,RIFA:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"RIFA")]$patient_pssn,1,0)]
  cohort3[,ISON:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"ISON")]$patient_pssn,1,0)]
  cohort3[,STRE:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"STRE")]$patient_pssn,1,0)]
  cohort3[,RIFI:=fifelse(patient_pssn%in%tb_px_n[str_detect(item_cd,"RIFI")]$patient_pssn,1,0)]
  
  cohort4 <- cohort3[patient_pssn%in%tbnew$patient_pssn]
  selected <- rbind(cohort4[PYRA==1&ISON==1&RIFA==1][STRE==1|ETHA==1],cohort4[PYRA==1&RIFI==1][STRE==1|ETHA==1])
  cases <- tbnew[patient_pssn%in%selected$patient_pssn]#908 out of 1503 have antibiotics within 14 days after diagnosis 1521-1503=18 metastic cancer

  Allip <- dx_latest[as_date(date)>=as_date("2021-2-23")&Source=="1.IP"]  #3348137
  setorder(Allip,"patient_pssn","date")
  ctrls <- Allip[patient_pssn%in%cohort3$patient_pssn]#3271676/3071193 
#remove records of cases and records after first tuberculosis diagnosis
  ctrls <- ctrls[!patient_pssn%in%cases$patient_pssn]#3064137
  ctrls_cases <- Allip[patient_pssn%in%cases$patient_pssn]
  ctrls_cases[str_detect(code,"^01[0-8]"),indexdate:=min(date),by=patient_pssn]
  ctrls_cases[,indexdate:=max(indexdate,na.rm = TRUE),by=patient_pssn]
  ctrls_bef_cases <- ctrls_cases[date<indexdate][,!c("indexdate")]
  ctrls <- rbind(ctrls,ctrls_bef_cases)#3065399
  
  cohort1[!is.na(Vaccine.Brand.1st)&!is.na(Vaccine.Brand.2nd),vaccinated:=1]
  cohort1[is.na(vaccinated),vaccinated:=0]
  cases <- left_join(cases,cohort1[,.(patient_pssn,Age,sex,Vaccine.Brand.1st,Date.of.vaccination.1st,Vaccine.Brand.2nd,Date.of.vaccination.2nd,death_date_ymd,death_diag_cd,vaccinated,Vaccine.Brand.3rd,Date.of.vaccination.3rd )])
  ctrls <- left_join(ctrls,cohort1[,.(patient_pssn,Age,sex,Vaccine.Brand.1st,Date.of.vaccination.1st,Vaccine.Brand.2nd,Date.of.vaccination.2nd,death_date_ymd,death_diag_cd,vaccinated,Vaccine.Brand.3rd,Date.of.vaccination.3rd )])
#remove missing value in age and sex
  ctrls <- ctrls[(sex=="M" | sex=="F") & !is.na(Age)]#3065399
  cases <- cases[(sex=="M" | sex=="F") & !is.na(Age)]#908
  cases[,date:=as_date(date)]
## Correct age by reducing one year

  cases <- cases[Age>=18]#890
  ctrls <- ctrls[Age>=18]#3016632

  # remove admission records between first and second dose, and after third dose
  ctrls[, vaccinated := as.numeric((!is.na(Date.of.vaccination.1st)) & Date.of.vaccination.1st <= date)]
  ctrls[, `2nd_dose` := as.numeric(vaccinated & (!is.na(Date.of.vaccination.2nd)) & date >= Date.of.vaccination.2nd)]
  ctrls <- ctrls[!(vaccinated==1&`2nd_dose`==0)]#2925549
  # ctrls <- ctrls[date<as_date(Date.of.vaccination.3rd)|is.na(Date.of.vaccination.3rd)]#2906266
  cases[, vaccinated := as.numeric((!is.na(Date.of.vaccination.1st)) & Date.of.vaccination.1st <= date)]
  cases[, `2nd_dose` := as.numeric(vaccinated & (!is.na(Date.of.vaccination.2nd)) & date >= Date.of.vaccination.2nd)]
  cases <- cases[!(vaccinated==1&`2nd_dose`==0)]#863
  # cases <- cases[date<as_date(Date.of.vaccination.3rd)|is.na(Date.of.vaccination.3rd)]#858
  #generate cohort for matching
  cohort_tb <- rbind(cases[, tb:=1], ctrls[, tb:=0], fill=T)
  setorder(cohort_tb,"patient_pssn","date")
  # 2. Matching ----

 cc.match <- function(cohort, case.var, age.var, age.gap, date.var, date.gap, match.var, K, seed) {
  setorderv(cohort, c(match.var, age.var, "patient_pssn"))
  cases <- cohort[get(case.var)==1]
  ctrls <- cohort[get(case.var)==0]
  
  matched <- copy(cases[, c("patient_pssn", age.var, date.var, match.var), with=F])
  .Random.seed <- seed
  
  for(i in 1:nrow(cases)) {
    case <- cases[i]
    ctrls_ <- ctrls[abs(get(age.var)-case[[age.var]]) <= age.gap]
    ctrls_ <- ctrls_[abs(get(date.var)-case[[date.var]]) <= date.gap]
    pool <- Reduce(intersect, lapply(match.var, function(v) ctrls_[get(v)==case[[v]], unique(patient_pssn)]))
    if(length(pool)==0) {cat("No match for", case$patient_pssn, "\n"); ctrl <- rep(NA, K)}
    else ctrl <- c(sample(as.character(pool), min(length(pool),K), replace=F),rep(NA,K-min(length(pool),K)))
    matched[i, paste0("M",1:K) := as.list(ctrl)]
    matched[i, pool.size := length(pool)]
  }
  
  return(matched)
 }

# Generate matches 1)No age difference is allowed 2) allow one day deviation in admission date
set.seed(475)
seed <- .Random.seed
matched <- cc.match(copy(cohort_tb), "tb", "Age", 0, "date", 1, c("sex"), 10, seed)
# matched[pool.size<10]#8 cases fails to match 10
  setnames(matched,"patient_pssn","Case")
  matched <-  melt(matched,id.vars = c("date","Case","Age","sex"), measure.vars = c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10"),
          variable.name = "Index", value.name = "patient_pssn")
  matched <- matched[!Case==patient_pssn]# remove 1 self-match case and NA
  uniqueN(matched$patient_pssn)##8195 unique reference key
  
  matched.check <- matched[patient_pssn%in%matched$Case]# 6 cases in both control and cases group
  matched[,patient_pssn:=as.integer(patient_pssn)]
  
  matched.ctrls <- left_join(matched,cohort1,by=c("Age","sex","patient_pssn"))
  setnames(matched.ctrls,"Case","match")
  matched.cases <- cases[,match:=patient_pssn]
  matched.cases <- matched.cases[,.(patient_pssn,date,Age,sex,match,Vaccine.Brand.1st,Date.of.vaccination.1st,Vaccine.Brand.2nd,Date.of.vaccination.2nd,death_date_ymd,death_diag_cd,vaccinated)]
  matched.ctrls <-matched.ctrls[,.(patient_pssn,date,Age,sex,match,Vaccine.Brand.1st,Date.of.vaccination.1st,Vaccine.Brand.2nd,Date.of.vaccination.2nd,death_date_ymd,death_diag_cd,vaccinated)]
  
  saveRDS(matched,"Matched.RDS")
  saveRDS(matched.cases,"cases.RDS")
  saveRDS(matched.ctrls,"ctrls.RDS")


  
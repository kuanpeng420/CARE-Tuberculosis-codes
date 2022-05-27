#Load required package 
  pacman::p_load(data.table,tidyverse,tableone,writexl,readxl,comorbidity,broom,survival,ggplot2,powerSurvEpi,elrm,lubridate)

#1. Generate covariates for cohorts----
#Load data  

  case <- readRDS("cases.RDS")
  control <- readRDS("ctrls.RDS")
  cc <- rbind(case,control)
  setnames(cc,"date","admission_date")
  load("D:/Batch 5/R/DX.RData")
  dx_latest[,date:=as_date(date)]
  dx_clean[,date:=as_date(paste0(date,"-15"))]
  dx <- rbind(dx_latest,dx_clean)
#Identify cohort disease history
  
  adj <- setDT(read_xlsx("codes_TB.xlsx",sheet = 2))
  adj_drug <- adj[str_detect(Name,"rx")]
  adj_dx <- adj[str_detect(Name,"dx")]
  for (i in 1:length(adj_dx$Regex)) {
    dx.tb <- dx[str_detect(code,regex(adj_dx$Regex[i],ignore_case = T))]
    dx.tb <- left_join(dx.tb,cc[,.(patient_pssn,admission_date)],by="patient_pssn")
    dx.tb <- dx.tb[as_date(date)<as_date(admission_date)]
    dx.tb <- unique(full_join(cc[,.(patient_pssn,admission_date)],dx.tb,by=c("patient_pssn","admission_date")),by=c("patient_pssn","admission_date"))
    dx.tb <- dx.tb[,adj_dx$Name[i]:=as.numeric(!is.na(code))]
    cc <- left_join(cc,dx.tb[,c("date","order","code","code_ext","Source","ranking"):=NULL],by=c("patient_pssn","admission_date"))
    print(i)
  }
  drug_clean <- readRDS("D:/Batch 5/R/RX_clean.RDS")
  drug_latest <- readRDS("D:/Batch 5/R/RX_latest.RDS")
  drug_clean_combined <- rbind(drug_clean$`2018`,drug_clean$`2019`,drug_clean$`2020`)

  drug_latest[,disp_date_ymd:=as_date(disp_date_ymd)]
  drug_latest[,presc_end_date_ymd:=disp_date_ymd+presc_duration_day]
  drug_latest_tb <- drug_latest[patient_pssn%in%cc$patient_pssn]
  drug_clean_combined[,disp_date_ymd:=as_date(paste0(disp_date_ym,"-15"))]
  drug_clean_combined[,presc_end_date_ymd:=disp_date_ymd+presc_duration_day]
  drug_clean_combined[,presc_start_date_ym:=NULL]
  drug_clean_tb <- drug_clean_combined[patient_pssn%in%cc$patient_pssn]
  drug_tb <- rbind(drug_clean_tb,drug_latest_tb,fill=TRUE)
  for (i in 1:length(adj_drug$Regex)) {
    rx.drug <- drug_tb[str_detect(bnfno_p,regex(adj_drug$Regex[i],ignore_case = T))]
    rx.drug <- left_join(rx.drug,cc[,.(patient_pssn,admission_date)],by="patient_pssn")
    rx.drug <- rx.drug[disp_date_ymd<=as_date(admission_date)&presc_end_date_ymd>=as_date(admission_date-180)]
    rx.drug <- unique(full_join(cc[,.(patient_pssn,admission_date)],rx.drug,by=c("patient_pssn","admission_date")),by=c("patient_pssn","admission_date"))
    rx.drug <- rx.drug[,adj_drug$Name[i]:=as.numeric(!is.na(item_cd))]
    cc <- left_join(cc,rx.drug[,c(1,2,16)],by=c("patient_pssn","admission_date"))
    print(i)
  }
  cc[,cci:=dx.mi+dx.chf+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+(dx.dm_com0&!dx.dm_com1)+dx.dm_com1*2+dx.crf*2+(dx.liver_mild&!dx.liver_modsev)+dx.liver_modsev*3+dx.ulcers+dx.ra+dx.aids*6+dx.cancer*2+dx.cancer_mets*6]
  
  adj_rx <- setDT(read_xlsx("CARE-TB-Crosscheck.xlsx",sheet = 5))
  remove(rx.drug)
  for (i in 1:length(adj_rx$Regex)) {
    rx.drug <- drug_tb[str_detect(item_cd,regex(adj_rx$Regex[i],ignore_case = T))]
    rx.drug <- left_join(rx.drug,cc[,.(patient_pssn,admission_date)],by="patient_pssn")
    rx.drug <- rx.drug[disp_date_ymd<=as_date(admission_date)&presc_end_date_ymd>=as_date(admission_date-180)]
    rx.drug <- unique(full_join(cc[,.(patient_pssn,admission_date)],rx.drug,by=c("patient_pssn","admission_date")),by=c("patient_pssn","admission_date"))
    rx.drug <- rx.drug[,adj_rx$Name[i]:=as.numeric(!is.na(item_cd))]
    cc <- left_join(cc,rx.drug[,c(1,2,16)],by=c("patient_pssn","admission_date"))
    print(i)
  }
  
  load("D:/Batch 5/R/A10X.RData")
  admin_clean_tb <- admin_clean[patient_pssn%in%cc$patient_pssn]
  admin_clean_tb[,date:=as_date(paste0(date,"-15"))]
  admin_latest_tb <- admin_latest[patient_pssn%in%cc$patient_pssn]
  admin_latest_tb$date <- as_date(admin_latest_tb$date)
  admin_combine <- rbind(admin_clean_tb[,dischg_date_ym:=NULL],admin_latest_tb[,dischg_date_ymd:=NULL])
  admin_combine<- left_join(admin_combine,cc[,.(patient_pssn,admission_date,match)])
  admin_tb <- admin_combine[date<admission_date&date>as_date(admission_date-365)]

  admin_tb_ip <-admin_tb[,.N,by=c("patient_pssn","Source","match")][str_detect(Source,"1.IP"),.(ip=sum(N)),by=.(patient_pssn,match)]
  admin_tb_ae <-admin_tb[,.N,by=c("patient_pssn","Source","match")][str_detect(Source,"2.AE"),.(ae=sum(N)),by=.(patient_pssn,match)]
  admin_tb_gopc <-admin_tb[,.N,by=c("patient_pssn","Source","match")][str_detect(Source,"4.GO"),.(gopc=sum(N)),by=.(patient_pssn,match)]
  admin_tb_sopc <-admin_tb[,.N,by=c("patient_pssn","Source","match")][str_detect(Source,"3.SO"),.(sopc=sum(N)),by=.(patient_pssn,match)]
  cc <- left_join(cc,admin_tb_ae)
  cc <- left_join(cc,admin_tb_ip)
  cc <- left_join(cc,admin_tb_gopc)
  cc <- left_join(cc,admin_tb_sopc)
  # NA value means no ipae/op visits
  cc[is.na(sopc),sopc:=0]
  cc[is.na(ip),ip:=0]
  cc[is.na(gopc),gopc:=0]
  cc[is.na(ae),ae:=0]

  cc[is.na(Date.of.vaccination.1st) & !is.na(Date.of.vaccination.2nd)] # None
  cc[,`:=`(Date.of.vaccination.1st=as_date(Date.of.vaccination.1st),Date.of.vaccination.2nd=as_date(Date.of.vaccination.2nd))]
  cc[, vaccinated := as.numeric((!is.na(Date.of.vaccination.1st)) & Date.of.vaccination.1st < admission_date)]
  
  cc[, vaccinated.sinovac := as.numeric(vaccinated & Vaccine.Brand.1st=="Sinovac")]
  cc[, vaccinated.biontech := as.numeric(vaccinated & Vaccine.Brand.1st=="BioNTech/Fosun")]
  cc[, vaccinated.brand := ifelse(vaccinated, Vaccine.Brand.1st, "None")]
  cc[, vaccinated.brand := factor(vaccinated.brand, c("None","BioNTech/Fosun","Sinovac"))]
  
# 2nd dose and onset
  cc[, `2nd_dose` := as.numeric(vaccinated & (!is.na(Date.of.vaccination.2nd)) & admission_date >= Date.of.vaccination.2nd)]
  cc[vaccinated==1, onset_time2 := as.numeric(admission_date-ifelse(`2nd_dose`,Date.of.vaccination.2nd,Date.of.vaccination.1st))]
  cc[match==patient_pssn,tb:=as.integer(1)]
  cc[!match==patient_pssn,tb:=as.integer(0)]
  cc$`2nd_dose` <- as.factor(cc$`2nd_dose`)
  cc$vaccinated <- as.factor(cc$vaccinated)
  cc$vaccinated.sinovac <- as.factor(cc$vaccinated.sinovac)
  cc$vaccinated.biontech <- as.factor(cc$vaccinated.biontech)
  
  factorvar <- c(subset(colnames(cc),str_detect(colnames(cc),"rx|dx|vaccinated")))
  cc[, (factorvar):= lapply(.SD,as.factor),.SDcol=factorvar]#convert to factor for clogit

#Categorize cohort by agegroup
  agegroup <- function(x){
         if (x>=60) {
                 return("Old")
         }  else{
                 return("Adult")
         }
  }
  cc$agegroup <- sapply(cc$Age, agegroup)
  
  #adjust covide history before admission date
  COVID <- readRDS("D:/CARE Part 3/R/LAB_ALL_COVID.RDS")# Identify COVID-19 case using 213,215,216,217
  setDT(COVID)
  covid <- COVID[str_detect(T_NUM,"^213|^215|^216|^217")&result=="detected"]
  setorder(covid,date,by=patient_pssn)
  covid <- covid[,head(.SD,1),by=patient_pssn]
  setnames(covid,"date","covid_date")
  cc <- left_join(cc,covid[,.(patient_pssn, covid_date)])
  cc[!is.na(covid_date)&covid_date<admission_date,covid:=1]
  cc[is.na(covid),covid:=0]
  cc$covid <- as.factor(cc$covid)
  
  saveRDS(cc, "casectrl.tb.RDS")


 
#3. Generate table one----
  
# Table one
  data <- readRDS("casectrl.tb.RDS")
  data30 <- data[!match%in%data[(tb==1&onset_time2<30)]$match]
  
  table.cohorts <- data[,-c("patient_pssn","match","admission_date",
                                    "onset_time2","2nd_dose","Date.of.vaccination.1st","Date.of.vaccination.2nd","death_date_ymd","death_diag_cd")]
  table.cohorts30 <- data30[,-c("patient_pssn","match","admission_date",
    "onset_time2","2nd_dose","Date.of.vaccination.1st","Date.of.vaccination.2nd","death_date_ymd","death_diag_cd")]
  
  s1 <- print(CreateTableOne(data=table.cohorts[tb==1][,.(Age,sex,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,dx.infection_viral,dx.sle,dx.spa,dx.psa, dx.ms, dx.ibd,covid,cci,ip,ae,
    sopc,gopc, dx.mi, dx.pvd, dx.cbd, dx.copd, dx.dementia, dx.paralysis, dx.crf, dx.liver_mild, dx.liver_modsev,dx.ulcers,dx.ra, dx.cancer,
    dx.cancer_mets, dx.dm_com0,dx.dm_com1,rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,vaccinated.brand)],strata=c("vaccinated.brand")),
          showAllLevels = F, quote = T, noSpaces = T,test = F,smd = T)
  s2 <- print(CreateTableOne(data=table.cohorts[tb==0][,.(Age,sex,dx.chf,dx.stroke_embo,dx.asthma,dx.infection_resp,dx.infection_viral,dx.sle,dx.spa,dx.psa, dx.ms, dx.ibd,covid,cci,ip,ae,
    sopc,gopc, dx.mi, dx.pvd, dx.cbd, dx.copd, dx.dementia, dx.paralysis, dx.crf, dx.liver_mild, dx.liver_modsev,dx.ulcers,dx.ra, dx.cancer,
    dx.cancer_mets, dx.dm_com0,dx.dm_com1,rx.av,rx.ab,rx.gc,rx.immunsupp,rx.ac,rx.am,rx.bb,rx.ccb,rx.statin,rx.acei,rx.arb,rx.dig,vaccinated.brand)],strata=c("vaccinated.brand")),
              showAllLevels = F, quote = T, noSpaces = T,test = F,smd = T)
  tableone <- list(s1,s2)
  write.csv(tableone,"Tableone.new.csv")
#Table one with two factor stratification
  
  t1 <- print(CreateTableOne(data=table.cohorts,vars = c("vaccinated.brand"),strata=c("tb")),
    showAllLevels = F, quote = T, noSpaces = T)
  t2 <- print(CreateTableOne(data=table.cohorts,vars = c("vaccinated.brand"),strata=c("tb","sex")),
    showAllLevels = F, quote = T, noSpaces = T)
  
  t3 <- print(CreateTableOne(data=table.cohorts,vars = c("vaccinated.brand"),strata=c("tb","agegroup")),
          showAllLevels = F, quote = T, noSpaces = T)
  
  tabletwo <- list(t1,t2,t3)
  write.csv(tabletwo,"Tabletwo.new.csv")
  
  t1s <- print(CreateTableOne(data=table.cohorts30,vars = c("vaccinated.brand"),strata=c("tb")),
    showAllLevels = F, quote = T, noSpaces = T)
  t2s <- print(CreateTableOne(data=table.cohorts30,vars = c("vaccinated.brand"),strata=c("tb","sex")),
    showAllLevels = F, quote = T, noSpaces = T)
  
  t3s <- print(CreateTableOne(data=table.cohorts30,vars = c("vaccinated.brand"),strata=c("tb","agegroup")),
    showAllLevels = F, quote = T, noSpaces = T)
  
  tablethree <- list(t1s,t2s,t3s)
  write.csv(tablethree,"Tablethree.new.csv")  

#4.1 Conditional logistic regression----

  data <- readRDS("casectrl.tb.RDS")# Still there are 3 people vaccinated but 2nd_dose is 0, this is because 1 day deviation in admission date.
  data30 <- data[!match%in%data[(tb==1&onset_time2<30)]$match]
  
  
  clm1 <- clogit(tb ~ vaccinated.sinovac + vaccinated.biontech + strata(match),data)
  clm2 <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid,data)
  #stratified by sex
  clm3 <- clogit(tb ~ vaccinated.sinovac + vaccinated.biontech + strata(match), data[sex=="M"])
  clm4 <- clogit(tb ~ vaccinated.sinovac + vaccinated.biontech + strata(match), data[sex=="F"])
  clm5 <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data[sex=="M"])#remove ms
  clm6 <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data[sex=="F"])#remove ms
  #stratified by agegroup
  clm9 <- clogit(tb ~ vaccinated.sinovac+vaccinated.biontech + strata(match), data[Age>=18&Age<=59])
  clm10 <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data[Age>=18&Age<=59])
  
  clm11 <- clogit(tb ~ vaccinated.sinovac+vaccinated.biontech + strata(match), data[Age>=60])
  clm12 <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data[Age>=60])
  #data30
  clm1s <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid,data30)
  #stratified by sex
  clm2s <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data30[sex=="M"])#remove ms
  clm3s <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data30[sex=="F"])#remove ms
  #stratified by agegroup
  clm4s <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data30[Age>=18&Age<=59])
  clm5s <- clogit(tb ~ vaccinated.sinovac +
      vaccinated.biontech+ strata(match) +dx.stroke_embo+dx.asthma+dx.infection_resp+dx.infection_viral+dx.ms+
      dx.sle+dx.spa+dx.psa+dx.ibd+rx.av+rx.ab+rx.gc+rx.immunsupp+rx.ac+rx.am+rx.bb+rx.ccb+rx.statin+rx.acei+rx.arb+rx.dig+ip+ae+sopc+gopc+cci+covid, data30[Age>=60])
  
  r1 <- merge( tidy(clm1) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm1)),keep.rownames=T),by.x = "term",by.y = "rn")
  r2 <- merge( tidy(clm2) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm2)),keep.rownames=T),by.x = "term",by.y = "rn")
  r3 <- merge( tidy(clm3) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm3)),keep.rownames=T),by.x = "term",by.y = "rn")
  r4 <- merge( tidy(clm4) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm4)),keep.rownames=T),by.x = "term",by.y = "rn")
  r5 <- merge( tidy(clm5) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm5)),keep.rownames=T),by.x = "term",by.y = "rn")
  r6 <- merge( tidy(clm6) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm6)),keep.rownames=T),by.x = "term",by.y = "rn")
  # r7 <- merge( tidy(clm7) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm7)),keep.rownames=T),by.x = "term",by.y = "rn")
  # r8 <- merge( tidy(clm8) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm8)),keep.rownames=T),by.x = "term",by.y = "rn")
  r9 <- merge( tidy(clm9) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm9)),keep.rownames=T),by.x = "term",by.y = "rn")
  r10 <- merge( tidy(clm10) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm10)),keep.rownames=T),by.x = "term",by.y = "rn")
  r11 <- merge( tidy(clm11) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm11)),keep.rownames=T),by.x = "term",by.y = "rn")
  r12 <- merge( tidy(clm12) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm12)),keep.rownames=T),by.x = "term",by.y = "rn")
  r13 <- merge( tidy(clm1s) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm1s)),keep.rownames=T),by.x = "term",by.y = "rn")
  r14 <- merge( tidy(clm2s) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm2s)),keep.rownames=T),by.x = "term",by.y = "rn")
  r15 <- merge( tidy(clm3s) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm3s)),keep.rownames=T),by.x = "term",by.y = "rn")
  r16 <- merge( tidy(clm4s) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm4s)),keep.rownames=T),by.x = "term",by.y = "rn")
  r17 <- merge( tidy(clm5s) %>% mutate(exp.estimate=exp(estimate),.after=estimate), data.table(exp(confint.default(clm5s)),keep.rownames=T),by.x = "term",by.y = "rn")
  
  # r <- list("Main.crude"=r1,"Main.adj"=r2,"Male.crude"=r3,"Female.crude"=r4,"Male.adjust"=r5,
  #           "Female.adjust"=r6,"Adult.crude"=r9,
  #           "Adult.adj"=r10,"old.crude"=r11,"old.adj"=r12,"bef.crude"=r13,"bef.adj"=r14,"aft.crude"=r15,"aft.adj"=r16)
  r <- list("Main.crude"=r1,"Main.adj"=r2,"Male.crude"=r3,"Female.crude"=r4,"Male.adjust"=r5,
    "Female.adjust"=r6,"Adult.crude"=r9,
    "Adult.adj"=r10,"old.crude"=r11,"old.adj"=r12,"s.all"=r13,"s.male"=r14,"s.female"=r15,"s.adult"=r16,"s.old"=r17)
  write_xlsx(r,"casectrl.tb.xlsx")
  
  

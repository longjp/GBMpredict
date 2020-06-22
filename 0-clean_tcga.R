## clean tcga GBM data
## create variables with death time and indicator for censoring
rm(list=ls())
##getwd()
##setwd("~/Desktop/fjord/2019-12-ferguson-lts/analysis")
library(openxlsx)

dat <- read.csv("../data/clinical.tsv",
                sep="\t",
                na.strings=c("NA","--"),
                stringsAsFactors=FALSE)

## patients are all entered twice, 1 row for different treatment
## types
nrow(dat)
length(unique(dat$case_id))
length(unique(paste0(dat$treatment_type,dat$case_id)))
table(dat$treatment_type)
## only interested in survival times, so just choose radiation therapy
##dat <- dat[dat$treatment_type=="Pharmaceutical Therapy, NOS",]
dat <- dat[dat$treatment_type=="Radiation Therapy, NOS",]
nrow(dat)

## # individuals with no days to death and no days to last follow up
## remove patients?
sum(is.na(dat$days_to_death) & is.na(dat$days_to_last_follow_up))
dat <- dat[!(is.na(dat$days_to_death) & is.na(dat$days_to_last_follow_up)),]

## remove patients with unknown vital status
table(dat$vital_status)
dat <- dat[dat$vital_status!="Not Reported",]
table(dat$vital_status)

## death is always at least as far out as follow-up, GOOD
mean(dat$days_to_death >= dat$days_to_last_follow_up,na.rm=TRUE)

## ti is observed time, either death or censoring
dat$ti <- dat$days_to_death
dat$ti[is.na(dat$ti)] <- dat$days_to_last_follow_up[is.na(dat$ti)]
dat$ti_obs <- rep(1,nrow(dat))
dat$ti_obs[dat$vital_status=="Alive"] <- 0
table(dat$ti_obs,dat$vital_status)


### read in KPS data
### see email from Alfaro-Munoz,Kristin D <KDAlfaro@mdanderson.org>
### on February 17
##kps <- read.xlsx("1-s2.0-S009286741501692X-mmc2.xlsx",startRow=2)
kps <- read.xlsx("../data/TCGA pten tp53.xlsx",startRow=2)
## contains many non GBM data, eg when kps$Case is TCGA-WY
## https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
## so we just grab the cases that have matching submitter_id in dat
length(unique(kps$Case))
nrow(kps)
mean(dat$submitter_id %in% kps$Case)
rownames(kps) <- kps$Case
kps_sub <- kps[dat$submitter_id,]
tcga <- cbind(dat,kps_sub)
identical(tcga$submitter_id,tcga$Case)

#### sanity checks
## age looks good
plot(tcga$age_at_diagnosis,tcga$`Age.(years.at.diagnosis)`)
## gender good
table(tcga$gender,tcga$Gender)
## survival usually matches
plot(tcga$`Survival.(months)`*365/12,tcga$ti)
## survival usually matches
table(tcga$`Vital.status.(1=dead)`,tcga$vital_status,useNA="always")
## when vital status does not match, the patient is living longer in
## tcga portal data, suggesting that discrepancy is caused by portal being 
## more up to date than spreadsheet
ix <- tcga$`Vital.status.(1=dead)`==0 & tcga$vital_status=="Dead" & !is.na(tcga$`Vital.status.(1=dead)`==0)
cols <- rep("#00000030",length(ix))
cols[ix] <- "red"
plot(tcga$`Survival.(months)`*365/12,tcga$ti,col=cols)

## convert survival time from days to years
tcga$ti <- tcga$ti / 365.25



### read in TML data
### see email from Roel on March 9
tml <- read.xlsx("../data/Supplementary_Table_3a.xlsx")
tml <- tml[!is.na(tml$Cancer.Type),]
tml <- tml[tml$Cancer.Type=="gbm",]
dim(tml) ## missing a lot of TML measures 
head(tml)
colnames(tml)[2:5] <- c("Mutations","BP","TML..per.Mb.","Cancer.Type")


tml$TCGA.ID <- substr(tml$TCGA.ID,1,12)
tml <- tml[tml$TCGA.ID %in%  tcga$submitter_id,]
colnames(tml)[1] <- "submitter_id"
tcga <- merge(tcga,tml,by="submitter_id",all=TRUE)

### read in completed mutation data
dat <- read.csv("../data/sample_matrix.txt",sep="\t")
dat_names <- substr(gsub("lgggbm_tcga_pub:","",dat[,1],fixed=TRUE),1,12)
length(dat_names)
length(unique(dat_names))
mean(dat_names %in% tcga$submitter_id)
mean(tcga$submitter_id %in% dat_names)

dat[,1] <- dat_names
colnames(dat)[1] <- "submitter_id"
tcga <- merge(tcga,dat,by="submitter_id")

## compare previous genomic features (lots of missingness) with new ones
table(tcga$IDH1,tcga$IDH.status,useNA="always")
table(tcga$`PTEN.Status.(1=Altered)`,tcga$PTEN,useNA="always")
table(tcga$`TP53.Status.(1=Altered)`,tcga$TP53,useNA="always")
head(dat)
table(tcga$BRAF.V600E.status,useNA="always")

## from current model: age, kps, tml, t1/t2, idh1, mgmt 
ix <- c("ti",
        "ti_obs",
        "age_at_diagnosis",
        "IDH1",
        "Karnofsky.Performance.Score",
        "MGMT.promoter.status",
        "gender",
        "BRAF.V600E.status",
        "PTEN",
        "TP53",
        "TML..per.Mb.")
tcga_sub <- tcga[,ix]
rownames(tcga_sub) <- tcga$submitter_id
head(tcga_sub)
colnames(tcga_sub)[3:ncol(tcga_sub)] <- c("Age","IDH1","KPS",
                                          "MGMT.Meth","Sex","BRAF","PTEN",
                                          "TP53","TML..per.Mb.")

#### convert all variable categories to match MDA cohort
## Age
summary(tcga_sub$Age)
tcga_sub$Age <- tcga_sub$Age/365.25
## IDH1
table(tcga_sub$IDH1,useNA="always")
## MGMT.Meth
table(tcga_sub$MGMT.Meth,useNA="always")
tcga_sub$MGMT.Meth[is.na(tcga_sub$MGMT.Meth)] <- "Unknown"
recode <- c("Methylated"="Y","Unmethylated"="N","Unknown"="Unknown")
tcga_sub$MGMT.Meth <- recode[tcga_sub$MGMT.Meth]
tcga_sub$MGMT.Meth <- as.factor(tcga_sub$MGMT.Meth)
table(tcga_sub$MGMT.Meth)
## Sex
table(tcga_sub$Sex,useNA="always")
recode <- c("female"="F","male"="M")
tcga_sub$Sex <- recode[tcga_sub$Sex]
tcga_sub$Sex <- as.factor(tcga_sub$Sex)
table(tcga_sub$Sex)
## BRAF
table(tcga_sub$BRAF,useNA="always")
tcga_sub$BRAF[is.na(tcga_sub$BRAF)] <- "unknown"
recode <- c("Mutant"=1,"WT"=0,"unknown"="unknown")
tcga_sub$BRAF <- recode[tcga_sub$BRAF]
table(tcga_sub$BRAF)
## PTEN and TP53
table(tcga_sub$PTEN,useNA="always")
table(tcga_sub$TP53,useNA="always")

head(tcga_sub)
lapply(tcga_sub,class)

## convert genomic to factors
genomic_feats <- c("IDH1","BRAF","PTEN","TP53")
for(ii in genomic_feats){
  tcga_sub[,ii] <- as.factor(tcga_sub[,ii])
}

## rename objects
tcga_full <- tcga
tcga <- tcga_sub


## save cleaned data
save(tcga,tcga_full,file="0-clean_tcga.RData")


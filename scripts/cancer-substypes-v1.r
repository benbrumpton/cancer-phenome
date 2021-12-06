# Next statement
# source('~/projects/cancer_phenome/scripts/cancer-phenome-v3.r')

#skipper <- TRUE
#if (skipper) {

rm(list=ls())

## load required packages
library("foreign")
library("data.table")
library("dplyr")

## load descriptions
d1 <- fread("/home/benb/projects/cancer_phenome/data/ICDO3-codes-topography.txt",header=F,sep=";")
setnames(d1,c("V1","V2"),c("ICDO3","Description"))

## add any cancer as a separate entry
d2 <- data.table(ICDO3 = "CX", Description = "Any type of cancer")
topo.des <- rbind(d1,d2)

## load key
key.master <- fread("/mnt/work/bridge/master-key/allin-master-key-20160702.csv",header=T,colClasses="character")
key.105880.105118 <- data.table(read.spss("/mnt/work/bridge/bridges-from-hunt/allin-endo-thyroidea-2015_14040/PID@105880-PID@105118.sav", to.data.frame = TRUE))

## Other key file for HUNT-genes
# Add Genotype ID and Cov
cov = data.table(read.table(gzfile("/mnt/work/genotypes/DATASET_20161002/SAMPLE_QC/Masterkey_DATASET.20161002.txt.gz"), header=T))
cov[, PID.105118 := gid.current] # Rename key
cov$gid.current <- as.character(cov$gid.current) # Match type
key.master.qc <- merge(key.master, cov, by="SentrixID") #69424 passed QC with SentrixID

# Check duplicates
dup <- key.master.qc %>% group_by(gid.current.x) %>% mutate(n = n()) %>% filter(n > 1)

# Remove duplicate
key.master.qc <- key.master.qc[!duplicated(key.master.qc), ]

## load data
d <- fread("/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/kilde/cancer/hveem_thyroidea_utlevert_1.csv",header=T,colClasses="character")
d[,DS := as.integer(DS)]
dat <- d[DS>=3]  ## malignant cancers only
dat[,TOPOGRAFI_ICDO3_letter := paste0("C",TOPOGRAFI_ICDO3)]
dat[,TOPOGRAFI_ICDO3_letter := paste0(substr(TOPOGRAFI_ICDO3_letter,1,3),".",substr(TOPOGRAFI_ICDO3_letter,4,4))]
dat[, ANY_CANCER := "CX"]

## roll up levels based on ICD-O3 codes
dat[, TOPOGRAFI_ICDO3_letter_upper := substr(TOPOGRAFI_ICDO3_letter,1,3)]
dat.m <- melt(dat, id.vars = c("PID_105880","DIAGNOSEAAR","ALDER"),
                measure.vars =  c("TOPOGRAFI_ICDO3_letter","TOPOGRAFI_ICDO3_letter_upper","ANY_CANCER"))
dat.m.u <- unique(dat.m,by=c("PID_105880","DIAGNOSEAAR","ALDER","value"))
sumByID <- data.table(as.data.frame(dat.m.u[, list(nr.appearances=.N),by = c("PID_105880","value")]))
sumByID[nr.appearances>=1,case := 1]
sumByID.dcast <- dcast(sumByID[,c("PID_105880","value","case"),with=F], PID_105880 ~ value, value.var = "case")

## merge back to descriptions
dat2 <- merge(key.105880.105118,sumByID.dcast,by.x="PID.105880",by.y="PID_105880",all.x=FALSE,all.y=TRUE)

#        } else {print ("next section")
dat2$PID.105118 <- as.character(dat2$PID.105118)
key.master.qc$PID.105118 <- as.character(key.master.qc$PID.105118)
dat3 <- merge(key.master.qc[,c("IID","FID","Sex","BirthYear","PID.105118", "batch",paste0("PC",1:20)),with=F],dat2,by="PID.105118",all.x=TRUE)
dat3[, PID.105118 := NULL]
dat3[, PID.105880 := NULL]

#         } # End of if statment

start_col <- 26

## set all persone without the disease to 'controls'
for (i in names(dat3)[start_col:NCOL(dat3)])      ## -- HI BEN -- YOU HAVE TO CHECK AT WHICH COLUMN YOUR PHENOTYPES START -- OR OTHERWISE MODIFY THE SCRIPT SO YOU ONLY TAKE THE "C COLUMS"
    dat3[is.na(get(i)), (i) := 0]

## exclude controls if they have cancer in same anatomical location
#for (i in names(dat3)[25:NCOL(dat3)])
#		dat3[get(i)==0 & get(substr(i,1,3))==1, (i) := NA] 

## exclude controls if they have any other cancer
for (i in names(dat3)[start_col:NCOL(dat3)])
		dat3[get(i)==0 & CX==1, (i) := NA] 

# add PATID and MATID
dat3$PATID <- 0
dat3$MATID <- 0

# order columns IID, FID, PATID, MATID
dat3 <- dat3 %>%
  select("FID", "IID", "PATID", "MATID", everything())

# add sub.group
dat3$sub.group <- NA
dat3$sub.group[dat3$C34==1 | dat3$C00==1| dat3$C15==1 | dat3$C16==1 | dat3$C25==1 | dat3$C67==1 ] <- 1 # Missing C90 but they are in ICD10_gr?
dat3$sub.group[dat3$CX==0] <- 0

# subset dataframe
dat4 <- dat3[,c("FID", "IID", "PATID", "MATID", "Sex", "BirthYear", "batch", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "sub.group")]

## output phenotypes
write.table(dat4, file="/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-subtypes-v1.txt", quote=F, sep="\t", row.names=F)

#######Not using this#######

## summarize
dat3m = melt(dat3, id.vars = c("IID"),
                measure.vars = names(dat3[,start_col:NCOL(dat3)]))
dat3m[value==1,case := 1]
dat3m[value==0,control := 1]
dat3m[is.na(value),Exclude := 1]

## sum cases, controls, and exclusions
dat3m.caseSum <- dat3m[,list(Ncase=sum(case,na.rm = TRUE)),by=variable]
dat3m.controlSum <- dat3m[,list(Ncontrol=sum(control,na.rm = TRUE)),by=variable]
dat3m.caseExcludeSum <- dat3m[,list(CaseControlExclude=sum(Exclude,na.rm = TRUE)),by=variable]
m1 <- merge(dat3m.caseSum,dat3m.controlSum,by="variable")
m2 <- merge(m1,dat3m.caseExcludeSum,by="variable")
m2[,Nsum := Ncase+Ncontrol+CaseControlExclude]
m3 <- merge(m2,topo.des,by.x="variable",by.y="ICDO3")
setorderv(m3,"Ncase",-1L)

## output summary
#write.table(m3,file="/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-phenome-by-ICDO3-summary-v1.txt",quote=F,row.names=F,sep="\t")

print('Done')

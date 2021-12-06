rm(list=ls())

## load required packages
require(data.table)
#library("foreign")

## load descriptions
d1 <- fread("/Users/jbn/Documents/RESEARCH/HUNT/HUNT_Analysis2/Cancer_phenome/ICDO3-codes-topography.txt",header=F,sep=";")
setnames(d1,c("V1","V2"),c("ICDO3","Description"))

## add any cancer as a separate entry
d2 <- data.table(ICDO3 = "CX", Description = "Any type of cancer")
topo.des <- rbind(d1,d2)

## load key
key.master <- fread("/Users/jbn/Documents/RESEARCH/HUNT/HUNT_Analysis2/AAA/key23.txt",header=T,colClasses="character")
key.105880.105118 <- data.table(read.spss("/Users/jbn/HUNT_Data/Thyroide_data/PID@105880-PID@105118.sav", to.data.frame = TRUE))

## load data
d <- fread("/Users/jbn/HUNT_Data/Thyroide_data/Thyroid_cancer/hveem_thyroidea_utlevert_1.csv",header=T,colClasses="character")
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
dat3 <- merge(key.master[,c("IID","FID","Sex","BirthYear","PID.105118",paste0("PC",1:20)),with=F],dat2,by="PID.105118",all.x=TRUE)
dat3[, PID.105118 := NULL]
dat3[, PID.105880 := NULL]

## set all persone without the disease to 'controls'
for (i in names(dat3)[25:NCOL(dat3)])      ## -- HI BEN -- YOU HAVE TO CHECK AT WHICH COLUMN YOUR PHENOTYPES START -- OR OTHERWISE MODIFY THE SCRIPT SO YOU ONLY TAKE THE "C COLUMS"
    dat3[is.na(get(i)), (i) := 0]

## exclude controls if they have cancer in same anatomical location
#for (i in names(dat3)[25:NCOL(dat3)])
#		dat3[get(i)==0 & get(substr(i,1,3))==1, (i) := NA] 

## exclude controls if they have any other cancer
for (i in names(dat3)[25:NCOL(dat3)])
		dat3[get(i)==0 & CX==1, (i) := NA] 

## output phenotypes
write.table(...)
	
## summarize
dat3m = melt(dat3, id.vars = c("IID"),
                measure.vars = names(dat3[,25:NCOL(dat3)]))
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
write.table(m3,file="/Users/jbn/Documents/RESEARCH/HUNT/HUNT_Analysis2/Cancer-phenome-by-ICDO3-summary-v1.tbl",quote=F,row.names=F,sep="\t")

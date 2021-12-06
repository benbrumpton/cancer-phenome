# Create key tables
#rm(list=ls())
#source('~/projects/cancer_phenome/scripts/create_key_file.R')

pheno <- read.table("/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-phenome-v1.txt", header=T)

## load required packages
require(data.table)

## load descriptions
d1 <- fread("/home/benb/projects/cancer_phenome/data/ICDO3-codes-topography.txt",header=F,sep=";")
setnames(d1,c("V1","V2"),c("ICDO3","Description"))

## add any cancer as a separate entry
d2 <- data.table(ICDO3 = "CX", Description = "Any type of cancer")
topo.des <- rbind(d1,d2)

# add the other column names from pheno
ICDO3 <- names(pheno[,1:24])
Description <- names(pheno[,1:24])
df <- data.frame(ICDO3, Description)
topo.des.header <- rbind(topo.des, df)

pheno2 <- pheno

names(pheno2) <- topo.des.header$Description[match(names(pheno2), topo.des.header$ICDO3)]

# check for NA's
idx <- is.na(names(pheno2))
names(pheno2[idx])
names(pheno[idx])

# add missing column names
ICDO3 <- names(pheno[idx])
Description <- names(pheno[idx])
df <- data.frame(ICDO3, Description)

topo.des.header.miss <- rbind(topo.des.header, df)

# add new column names again
pheno3 <- pheno
names(pheno3) <- topo.des.header.miss$Description[match(names(pheno3), topo.des.header.miss$ICDO3)]
table(is.na(names(pheno3)))

pheno <- pheno3

pheno_names <- names(pheno[25:344])

n_pheno <- length(pheno_names)

key <- data.frame(matrix(ncol=11,nrow=n_pheno))
names(key)=c("name", "version", "type", "baseQcovar", "baseBcovar", "addQcovar", "addBcovar", "description", "N", "N_cases", "N_controls")

key$name <- pheno_names
key$version <- "09.11.18"

for (i in 1:n_pheno){

#n_values <- length(which(!is.na(unique(pheno$snoring_H3))))
n_values <- length(which(!is.na(unique(pheno[, pheno_names[i]]))))

if (n_values==2) {

type <- "binary"
cases <- table(pheno[, pheno_names[i]])[2]
controls <- table(pheno[, pheno_names[i]])[1]
individuals <- addmargins(table(pheno[, pheno_names[i]]))[3]

} else {

type <- "quantitative"
cases <- "NA"
controls <- "NA"
individuals <- length(which(!is.na(pheno[, pheno_names[i]])))
} # End of else statment

key[i,3] <- type
key[i,11] <- controls
key[i, 10] <- cases
key[i, 9] <- individuals

} #End of loop

key$baseQcovar <- "BirthYear|PC1|PC2|PC3|PC4"
key$baseBcovar <- "batch"
key$addQcovar <- ""
key$addBcovar <- "Sex"
key$description <- "ICDO3-codes"

# remove cases with zero or NA
key <- key[key$N_cases != 'NA',]

# order the key by case numbers
key$N_cases <- as.numeric(key$N_cases)
key <- key[order(key$N_cases, decreasing=T),]

write.table(key, "/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-phenome-key-file.txt", col.names=T, row.names=F,quote=F, sep = "\t")

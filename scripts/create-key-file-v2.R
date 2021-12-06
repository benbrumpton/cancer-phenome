# Create key tables
#rm(list=ls())
#source('~/projects/cancer_phenome/scripts/create-key-file-v2.R')

pheno <- read.table("/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-phenome-v1.txt", header=T)

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

# lets just start with the top 3
a <- head(key, 1)

#b <- key[3:10,]
b <- key[11:15,]
#key <- rbind(a,b)

#Key for pipeline 1
b <- key[16:20,]
#key <- rbind(a,b)

#key for hunt-genes
b <- key[21:35,]
#key <- rbind(a,b)

# How many between 94 and 50
tt <- key[key$N_cases>=50,]

# key for pipeline1
b <- key[36:51,]

# key for hunt-genes
b <- key[52:65,]

write.table(b, "/mnt/work/phenotypes/phenotypes-from-data-owners/allin-endo-thyroidea-2015_14040/constructs/custom/cancer-phenome-key-file.txt", col.names=T, row.names=F,quote=F, sep = "\t")
#write.table(b, "/mnt/cargo/benb_out/cancer-phenome-key-file.txt", col.names=T, row.names=F,quote=F, sep = "\t")

#write.table(key, "~/projects/cancer_phenome/data/cancer-phenome-key-file.txt", col.names=T, row.names=F,quote=F, sep = "\t")

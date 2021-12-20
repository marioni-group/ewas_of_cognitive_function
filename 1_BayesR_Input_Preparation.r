#### Combined Model ~ 10k ####
## DNA methylation file 
require(EpiSmokEr)
require(lumi)
setwd("")
## Wave 1
meth_w1 = readRDS("norm_mvals_5087.rds")

# Mean impute missings
w1_nas = apply(meth_w1,1,function(x){which(is.na(x))})
nas = which(unlist(lapply(w1_nas, length))>0)
for(i in names(nas)){
	imp = mean(meth_w1[i,], na.rm=T)
	ind = which(is.na(meth_w1[i,]))
	meth_w1[i,ind] = rep(imp, length(ind))
}
which(is.na(meth_w1))
# numeric(0)
which(is.infinite(meth_w1))
# numeric(0)
# Calculate EpiSmokEr score with NA-imputed data
w1_epismk = epismoker(m2beta(meth_w1), method="SSc")

# Filter probes (crosshyb/snp/xchr)
list_w1 =  readRDS("Cog_w1_cpgs.rds") 
list_w1 = list_w1$cog_w1
meth_w1 = meth_w1[which(rownames(meth_w1) %in% list_w1),]


## Wave 3
meth_w3 = readRDS("wave3_mvals.rds")
meth_w3[which(is.nan(meth_w3))] <- NA
meth_w3[which(is.infinite(meth_w3))] <- NA

# Mean impute missings
w3_nas = apply(meth_w3,1,function(x){which(is.na(x))})
nas = which(unlist(lapply(w3_nas, length))>0)
for(i in names(nas)){
	imp = mean(meth_w3[i,], na.rm=T)
	ind = which(is.na(meth_w3[i,]))
	meth_w3[i,ind] = rep(imp, length(ind))
}
which(is.na(meth_w3))
# numeric(0)
which(is.infinite(meth_w3))
# numeric(0)

# Calculate EpiSmokEr score with NA-imputed data
w3_epismk = epismoker(m2beta(meth_w3), method="SSc")

list_w3 = read.table("w3_probes_to_keep.txt")
meth_w3 = meth_w3[which(row.names(meth_w3) %in% list_w3$V1),] 


## Find overlap between CpG sites 
cpgs = intersect(rownames(meth_w1), rownames(meth_w3))


## Make sure order of CpG sites match in both files 
meth_w1 = meth_w1[cpgs,]
meth_w3 = meth_w3[cpgs,]


## Phenotype file 

## Read in w3 phenotypes - already prepared 

pheno_w3 = readRDS("w3_phenos.rds") 

# Remove cog scores and calculate from scratch
pheno_w3 = pheno_w3[,-(20:42)]

epi_score_w3 <- w3_epismk
names(epi_score_w3)[1] <- "Basename"
names(epi_score_w3)[2] <- "EpiScore"
pheno_w3 <- merge(pheno_w3, epi_score_w3, by.x = "Sample_Sentrix_ID",  by.y = "Basename")

## read in w1 - need all 5087 
 	## need cognitive variables + age, sex, bmi, pack_years and ever_smoke 

## Age, Sex 
demo = readRDS("stradl-samples-5087.rds") 

## Combine with BMI 
bmi = read.csv("body.csv")
names(bmi)[1] <- "Sample_Name"
pheno_w1 <- merge(demo, bmi, by = "Sample_Name")

## Combine with Smoking info 
pck <- read.csv("../Bex_EWAS/pack_years.csv")
smk <- read.csv("../Bex_EWAS/ever_smoke.csv")

smoke <- merge(pck,smk,by="Sample_Name")

pheno_w1 <- merge(pheno_w1, smoke, by = "Sample_Name")


## Combine cognitive information 
source("pheno_prep_nov2020.r")
cog1$Sample_Name = cog1$ID

pheno_w1 <- merge(pheno_w1, cog1, by = "Sample_Name")
pheno_w3 <- merge(pheno_w3, cog1, by = "Sample_Name")

# Get Age months for pheno_w1
load("Age_Sex_5101_GS_02Feb2018.RData")
pheno_w1$Basename <- paste(pheno_w1$Sentrix_ID, pheno_w1$Sentrix_Position, sep = "_")

pheno_w1$age = age[match(pheno_w1$Basename, age$ID),"Age"]




## Get Basename and subset to the necessary columns - w1
epi_score_w1 <- w1_epismk
names(epi_score_w1)[1] <- "Basename"
names(epi_score_w1)[2] <- "EpiScore"
pheno_w1 <- merge(pheno_w1, epi_score_w1, by = "Basename") 


pheno_w1 <- pheno_w1[,c("Basename", "Sample_Name", "digit_symbol", "LM", "vocabulary", "verbal_total", "g", "gf", "age", "sex", "bmi", "EpiScore", "plate_processing_batch")]

## Rename batch so its the same as pheno_w3 
names(pheno_w1)[ncol(pheno_w1)] <- "Batch"

## Get Basename and subset to the necessary columns - w3

pheno_w3$Basename <- paste(pheno_w3$Sentrix_ID, pheno_w3$Sentrix_Position, sep = "_")
pheno_w3 <- pheno_w3[,c("Basename", "Sample_Name", "digit_symbol", "LM", "vocabulary", "verbal_total", "g", "gf", "age", "sex", "bmi", "EpiScore", "Batch")] 


## Combine the two waves 

pheno_w1$Set <- "W1"
pheno_w3$Set <- "W3"
phenotypes = rbind(pheno_w1, pheno_w3)

# Batches
phenotypes$Batch = paste0(phenotypes$Set, "_", phenotypes$Batch)

# Remove self-report AD individuals
ad = read.csv("disease.csv")
ad_y = ad[which(ad$alzheimers_Y==1), "ID"]
phenotypes = phenotypes[-which(phenotypes$Sample_Name %in% ad_y), ]

# complete cases
phenotypes = phenotypes[complete.cases(phenotypes), ] # 9172

# Subset to those also in GWAS data
samps_9162 = read.csv("Cognitive_BayesR_10k_phenotypes.csv")
phenotypes = phenotypes[which(phenotypes$Sample_Name %in% samps_9162$Sample_Name), ]

## Save phenotypes file 
write.csv(phenotypes, "Cognitive_BayesR_10k_phenotypes.csv",row.names = F)

# dir.create("Cog_Inputs")
## Prepare Phenotypes for 10k
phenotypes = read.csv("Cognitive_BayesR_10k_phenotypes.csv")
rownames(phenotypes) = phenotypes$Basename


# Calculate resduals, regressing out age, sex, BMI, epismoker
phenotypes$verbal_total <- resid(lm(phenotypes$verbal_total ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
verbal_total = phenotypes$verbal_total
verbal_total = scale(verbal_total) 
write.table(x = t(as.matrix(as.numeric(verbal_total))), "Cog_Inputs/Verbal_Total_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$vocabulary <- resid(lm(phenotypes$vocabulary ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
vocabulary = phenotypes$vocabulary
vocabulary = scale(vocabulary)
write.table(x = t(as.matrix(as.numeric(vocabulary))), "Cog_Inputs/Vocabulary_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$g <- resid(lm(phenotypes$g ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
g = phenotypes$g 
g = scale(g)
write.table(x = t(as.matrix(as.numeric(g))), "Cog_Inputs/G_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$gf <- resid(lm(phenotypes$gf ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
gf = phenotypes$gf 
gf = scale(gf)
write.table(x = t(as.matrix(as.numeric(gf))), "Cog_Inputs/GF_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 

phenotypes$LM <- resid(lm(phenotypes$LM ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
lm = phenotypes$LM 
lm = scale(lm)
write.table(x = t(as.matrix(as.numeric(lm))), "LM_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$digit_symbol <- resid(lm(phenotypes$digit_symbol ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$bmi + phenotypes$EpiScore, na.action = na.exclude))
digit_symbol = phenotypes$digit_symbol 
digit_symbol = scale(digit_symbol) 
write.table(x = t(as.matrix(as.numeric(digit_symbol))), "Cog_Inputs/Digit_Symbol_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 





## Save out CpG IDs 

cpgs = as.data.frame(cpgs)
names(cpgs) <- "cog_10k"
saveRDS(cpgs, "Cog_10k_cpgs.rds")

## Combine two waves 

total = cbind(meth_w1, meth_w3) 

## Write out till this far 


## Make sure order of IDs matches up with phenotypic files 
id = phenotypes$Basename 
meth = total[,id]

meth = t(meth)

# Regress age, sex, batch, epigenetic smoking from DNAm data
for(i in 1:ncol(meth)){ 
meth[,i] <- resid(lm(meth[,i] ~ phenotypes$age + as.factor(phenotypes$sex) + as.factor(phenotypes$Batch) + as.factor(phenotypes$Set) + phenotypes$EpiScore, na.action = na.exclude))
print(i)
} 

## Scale CpGs 
tmeth1 = t(meth) 
tmeth1 = as.data.frame(tmeth1) 
tmeth1 = apply(tmeth1, 1, scale)


write.table(x = t(tmeth1), "Cog_Inputs/DNAm_9162.csv", sep = ",", row.names = F, col.names = F, quote = F) 

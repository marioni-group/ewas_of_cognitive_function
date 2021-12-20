
# Write table of samples for GWAS data
phenos = read.csv("Cog_GWAS_BayesR/Cognitive_BayesR_10k_phenotypes.csv")
x = read.table("GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015.fam")
x = x[which(x$V2 %in% phenos$Sample_Name), ]
fam = x[,1:2]
write.table(fam[,1:2], file="ewas_samps.txt", sep='\t', quote=F, row.names=F)

# In terminal 
# plink19 --recodeA --bfile /GWAS_Source/GS_GWAS/GS20K_QC\'d/GS20K_QCpsychprotocol_SPH_04112015 --keep ewas_samps.txt --maf 0.01 --out GS_9162_gwas

# Back in R
setwd("Cog_GWAS_BayesR")
require(data.table)
snps = fread("GS_9162_gwas.raw", header=T, stringsAsFactors=F)
snps = as.data.frame(snps)

snp_df = as.data.frame(names(snps)[7:ncol(snps)])
names(snp_df) <- "cog_10k"
saveRDS(snp_df, "Cog_GWAS_BayesR/Cog_10k_snps.rds")



rownames(snps) = snps$IID
snps = snps[,7:ncol(snps)]

snps_scaled = apply(snps, 2, scale)

# mean impute
for(i in 1:ncol(snps_scaled)){
	na = which(is.na(snps_scaled[,i]))
	snps_scaled[na,i] = mean(snps_scaled[,i], na.rm=T)
	}

phenos = read.csv("../Cog_EWAS_BayesR/Cognitive_BayesR_10k_phenotypes.csv")
phenos = phenos[which(phenos$Sample_Name %in% rownames(snps)), ]

# Line up phenos and snps

snps_scaled = snps_scaled[match(phenos$Sample_Name, rownames(snps)),]

write.table(t(snps_scaled), "Cog_GWAS_BayesR/Cog_Inputs/GWAS_9162.csv", sep = ",", row.names = F, col.names = F, quote = F) 
		
phenotypes=phenos
pcs = read.table("GS20K_ALL_MAF5_PCA.eigenvec")
pcs = pcs[which(pcs$V2 %in% phenos$Sample_Name),]
pcs = pcs[match(phenos$Sample_Name, pcs$V2),]
stopifnot(isTRUE(all.equal(pcs$V2, phenos$Sample_Name)))

		phenotypes = cbind(phenotypes, pcs[,3:22])
write.csv(phenotypes, "Cog_GWAS_BayesR/Cognitive_BayesR_10k_phenotypes.csv",row.names = F)

phenotypes$verbal_total <- resid(lm(phenotypes$verbal_total ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
verbal_total = phenotypes$verbal_total
verbal_total = scale(verbal_total) 
write.table(x = t(as.matrix(as.numeric(verbal_total))), "Cog_GWAS_BayesR/Cog_Inputs/Verbal_Total.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$vocabulary <- resid(lm(phenotypes$vocabulary ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
vocabulary = phenotypes$vocabulary
vocabulary = scale(vocabulary)
write.table(x = t(as.matrix(as.numeric(vocabulary))), "Cog_GWAS_BayesR/Cog_Inputs/Vocabulary_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$g <- resid(lm(phenotypes$g ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
g = phenotypes$g 
g = scale(g)
write.table(x = t(as.matrix(as.numeric(g))), "Cog_GWAS_BayesR/Cog_Inputs/G_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$gf <- resid(lm(phenotypes$gf ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
gf = phenotypes$gf 
gf = scale(gf)
write.table(x = t(as.matrix(as.numeric(gf))), "Cog_GWAS_BayesR/Cog_Inputs/GF_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 

phenotypes$LM <- resid(lm(phenotypes$LM ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
lm = phenotypes$LM 
lm = scale(lm)
write.table(x = t(as.matrix(as.numeric(lm))), "Cog_GWAS_BayesR/Cog_Inputs/LM_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 


phenotypes$digit_symbol <- resid(lm(phenotypes$digit_symbol ~ phenotypes$age + factor(phenotypes$sex) + phenotypes$V3 + phenotypes$V4 + phenotypes$V5 + phenotypes$V6 + phenotypes$V7 + phenotypes$V8 + phenotypes$V9 + phenotypes$V10 + phenotypes$V11 + phenotypes$V12 + phenotypes$V13 + phenotypes$V14 + phenotypes$V15 + phenotypes$V16 + phenotypes$V17 + phenotypes$V18 + phenotypes$V19 + phenotypes$V20 + phenotypes$V21 + phenotypes$V22, na.action = na.exclude))
digit_symbol = phenotypes$digit_symbol 
digit_symbol = scale(digit_symbol) 
write.table(x = t(as.matrix(as.numeric(digit_symbol))), "Cog_GWAS_BayesR/Cog_Inputs/Digit_Symbol_10k.csvphen", sep = ",", row.names = F, col.names = F, quote = F) 

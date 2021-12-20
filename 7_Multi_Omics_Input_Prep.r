# In R
setwd("Cog_Combined_BayesR")
# Make Group file
cpgs = readRDS("../Cog_EWAS_BayesR/Cog_10k_cpgs.rds")
snps = readRDS("../Cog_GWAS_BayesR/Cog_10k_snps.rds")


cpgs$group = 0
snps$group = 1

features = rbind(cpgs, snps)

dir.create("Cog_Inputs")
write.table(features, file="Cog_Inputs/group_file.txt", sep=' ', row.names=F, col.names=F)

require(data.table)
ewas = fread("../Cog_EWAS_BayesR/Cog_Inputs/DNAm_9162.csv")
gwas = fread("../Cog_GWAS_BayesR/Cog_Inputs/GWAS_9162.csv")

ewas = as.data.frame(ewas)
gwas = as.data.frame(gwas)

data = rbind(ewas, gwas)


write.table(x = data, "Cog_Inputs/DNAm_GWAS_9162.csv", sep = ",", row.names = F, col.names = F, quote = F) 



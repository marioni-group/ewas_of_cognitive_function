# Calculate variance explained (GWAS)
## Open R 
	setwd("Cog_GWAS_BayesR/") 
	library(data.table) 
# dir.create("Cog_Summstats")
	loop = list.files("Cog_Comp/", pattern = "_processed.csv") 

## Calculate Mean Variance explained by all probes and credible intervals
names = readRDS("Cog_10k_snps.rds")  
names = names$cog_10k
out2 = list()
for(i in loop){  
  
  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  
  sigma <- read.csv(paste("Cog_Sigma/", i, sep=""))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T)
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  A <-gsub(".csv*","",i)
  output$Biomarker <- A 
  out2[[i]] = output
  write.csv(output, paste0("Cog_Summary/", A, "_output_10k.csv"), row.names = F) 
  
  
# Calculate variance explained (Multi-Omics)
## Open R 
	setwd("Cog_Combined_BayesR/") 

	loop = list.files("Cog_Comp/", pattern = "_processed.csv") 

  ## Calculate Mean Variance explained by all probes and credible intervals
out2 = list()
for(i in loop){  
  
  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  
  # For each iteration, sigma1/sigma1+sigma2+sigmaE
  sigma <- read.csv(paste("Cog_Sigma/", i, sep=""))  
  output$DNAm_Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:3]),na.rm = T) 
  output$DNAm_Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:3]),na.rm = T, prob =0.025)
  output$DNAm_Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:3]),na.rm = T, prob =0.975)

  output$GWAS_Mean_Variance_Explained = mean(sigma[,3]/rowSums(sigma[,1:3]),na.rm = T) 
  output$GWAS_Variance_LCI = quantile(sigma[,3]/rowSums(sigma[,1:3]),na.rm = T, prob =0.025)
  output$GWAS_Variance_HCI = quantile(sigma[,3]/rowSums(sigma[,1:3]),na.rm = T, prob =0.975)

  output$Mean_Variance_Explained = mean(rowSums(sigma[,2:3])/rowSums(sigma[,1:3]),na.rm = T) 
  output$Variance_LCI = quantile(rowSums(sigma[,2:3])/rowSums(sigma[,1:3]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(rowSums(sigma[,2:3])/rowSums(sigma[,1:3]),na.rm = T, prob =0.975)
  output$Data = "Combined"

  A <-gsub(".csv*","",i)
  output$Biomarker <- A 
  out2[[i]] = output
}

write.table(do.call("rbind", out2), file="Cog_Summary/Var_Explained_EWAS_GWAS_Combined.csv", quote=F, row.names=F)


  loop_dnam = list.files("../Cog_EWAS_BayesR/Cog_Comp/", pattern = "_processed.csv") 

## Step 1 - Calculate Mean Variance explained by all probes and credible intervals
out_dnam = list()
for(i in loop_dnam){  
  
  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  
  # For each iteration, sigma1/sigma1+sigma2+sigmaE
  sigma <- read.csv(paste("../Cog_EWAS_BayesR/Cog_Sigma/", i, sep=""))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T) 
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  output$Data = "DNAm"


  A <-gsub(".csv*","",i)
  output$Biomarker <- A 
  out_dnam[[i]] = output
}


loop_gwas = list.files("../Cog_GWAS_BayesR/Cog_Comp/", pattern = "_processed.csv") 
## Step 1 - Calculate Mean Variance explained by all probes and credible intervals
out_gwas = list()
for(i in loop_gwas){  
  
  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  
  # For each iteration, sigma1/sigma1+sigma2+sigmaE
  sigma <- read.csv(paste("../Cog_GWAS_BayesR/Cog_Sigma/", i, sep=""))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T) 
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  output$Data = "GWAS"

  A <-gsub(".csv*","",i)
  output$Biomarker <- A 
  out_gwas[[i]] = output
}

combine_var_exp = do.call("rbind", out2)
dnam_var_exp = do.call("rbind", out_dnam)
gwas_var_exp = do.call("rbind", out_gwas)

 write.csv(do.call("rbind", out2,), file="Cog_Combined_BayesR/Cog_Summary/Var_Explained_Combined_output_10k.csv", row.names=F, quote=F)


plotdata = rbind(combine_var_exp[,c("Biomarker", "Mean_Variance_Explained", "Variance_LCI", "Variance_HCI", "Data")], 
                 dnam_var_exp[,c("Biomarker", "Mean_Variance_Explained", "Variance_LCI", "Variance_HCI", "Data")],
                 gwas_var_exp[,c("Biomarker", "Mean_Variance_Explained", "Variance_LCI", "Variance_HCI", "Data")])
plotdata$Biomarker[grep("igit", plotdata$Biomarker)] = "Digit Symbol"
plotdata$Biomarker[grep("LM|lm", plotdata$Biomarker)] = "Logical Memory"
plotdata$Biomarker[grep("g_", plotdata$Biomarker)] = "g"
plotdata$Biomarker[grep("gf_", plotdata$Biomarker)] = "gf"
plotdata$Biomarker[grep("ocab", plotdata$Biomarker)] = "Vocabulary"
plotdata$Biomarker[grep("erbal", plotdata$Biomarker)] = "Verbal Total"

plotdata$Data = factor(plotdata$Data, levels=c("DNAm", "GWAS", "Combined"))
names(plotdata) = gsub("Biomarker", "Cognitive_Phenotype", names(plotdata))

plotdata$Cognitive_Phenotype = factor(plotdata$Cognitive_Phenotype, 
  levels = c("Logical Memory", "Digit Symbol", "Verbal Total", "Vocabulary", "gf", "g"))

saveRDS(plotdata, file="Cog_Combined_BayesR/plotdata.rds")

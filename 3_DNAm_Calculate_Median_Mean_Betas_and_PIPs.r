# Calculate variance explained (R)
setwd("Cog_EWAS_BayesR") 
library(data.table) 

loop = list.files("Cog_Comp/", pattern = "_processed.csv") 

## Step 1 - Calculate Mean Variance explained by all probes and credible intervals
names = readRDS("Cog_10k_cpgs.rds")  
names = names$cog_10k
out2 = list()
for(i in loop){  
  
  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
# Mean Variance Explained  
  sigma <- read.csv(paste("Cog_Sigma/", i, sep=""))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T)
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  A <-gsub(".csv*","",i)
  output$Biomarker <- A 
  out2[[i]] = output
  write.csv(output, paste0("Cog_Summary/", A, "_output_10k.csv"), row.names = F) 
  }
  
## Step 2 - Calculate the proportion of variance that is attributable to small, medium and large effects - 1,2,3
for(file in loop){    
  A <-gsub(".csv*","",file)
  betas <- fread(paste("Cog_Beta/", file, sep=""))
  betas = as.data.frame(betas)
  comp <- fread(paste("Cog_Comp/", file, sep=""))  
  comp = as.data.frame(comp) 
  names(comp) <- names
  names(betas) <- names
  list = apply(comp,1,function(x)which(!x%in% c(1,2,3)))  
  
  
  x <- as.matrix(0, ncol = 1, nrow = 1000) 
  x <-as.data.frame(x) 
  for(i in 1:length(list)){ 
    x[[i]]<-length(list[[i]]) == ncol(comp)
  } 
  
  if(length(which(x[1,] %in% "TRUE")) > 0){ 
    comp = comp[-which(x %in% "TRUE"),] 
  } else { 
    comp = comp} 
  
  
  t <- vector() 
  list = apply(comp,1,function(x)which(x %in% 3)) 
  for(i in 1:length(list)){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m1_true <- matrix(0, ncol = 1, nrow = n)
  m1_true <- as.data.frame(m1_true) 
  m1_true$ind <- ind_true
  x<-vector()
  for(j in m1_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m1= as.data.frame(x) 
  m1$x[is.na(m1$x)] <- 0 
  names(m1) <- "Variance_Small_Effects" 
  
  list = apply(comp,1,function(x)which(x %in% 2)) 
  for(i in 1:length(list)){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m2_true <- matrix(0, ncol = 1, nrow = n)
  m2_true <- as.data.frame(m2_true) 
  m2_true$ind <- ind_true
  x<-vector()
  for(j in m2_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m2= as.data.frame(x) 
  m2$x[is.na(m2$x)] <- 0 
  names(m2) <- "Variance_Medium_Effects"
  
  list = apply(comp,1,function(x)which(x %in% 1)) 
  for(i in 1:length(list)){ 
    t[[i]] <- length(list[[i]]) > 0 
  } 
  ind_true = which(t %in% "TRUE")
  ind_false = which(t %in% "FALSE")
  list_true = list[ind_true]
  list_false = list[ind_false] 
  n = length(list_true) 
  m3_true <- matrix(0, ncol = 1, nrow = n)
  m3_true <- as.data.frame(m3_true) 
  m3_true$ind <- ind_true
  x<-vector()
  for(j in m3_true$ind){ 
    x[j] <- sum((betas[j,list[[j]]])^2) 
  } 
  m3= as.data.frame(x) 
  m3$x[is.na(m3$x)] <- 0 
  names(m3) <- "Variance_Large_Effects"
  
  m1$num <- row.names(m1) 
  m2$num <- row.names(m2) 
  m3$num <- row.names(m3) 
  all = merge(m1, m2, by = "num", all = T) 
  var = merge(all, m3, by = "num", all = T) 
  var[is.na(var)] <- 0 
  var$num <- NULL
  var$Total_Variance <- var[,1] + var[,2] + var[,3]
  var$Proportion_Small_Effects <- var[,1]/var[,4]
  var$Proportion_Medium_Effects <- var[,2]/var[,4]
  var$Proportion_Large_Effects <- var[,3]/var[,4]
  out2[[file]]$Proportion_Small_Effects <- mean(var$Proportion_Small_Effects, na.rm=T) 
  out2[[file]]$Proportion_Medium_Effects <- mean(var$Proportion_Medium_Effects, na.rm=T) 
  out2[[file]]$Proportion_Large_Effects <- mean(var$Proportion_Large_Effects, na.rm=T) 
  
  out2[[file]]$Biomarker <- A
  print(A)
  write.csv(output, file = paste("Cog_Summary/", A, "_output_10k.csv", sep = ""), row.names = F) 
  
}
 write.csv(do.call("rbind", out2), paste0("Cog_Summary/", "Var_Explained_Combined_output_10k.csv"), row.names = F) 

out2 = read.csv("Cog_Summary/Var_Explained_Combined_output_10k.csv")


## Calculate Median Betas and PIPs
names = readRDS("Cog_10k_cpgs.rds")  
names = names$cog_10k

cpgs = names

for(i in loop){ 
  tmp = fread(paste0("Cog_Beta/", i))
  tmp = as.data.frame(tmp)
  print(i)
  A <-gsub(".csv*","",i)
  median = apply(tmp, 2, median)
  median_LLCI = apply(tmp, 2, function(x) quantile(x, probs =  0.025))
  median_LCI = apply(tmp, 2, function(x) quantile(x, probs =  0.05))
  median_HCI = apply(tmp, 2, function(x) quantile(x, probs =  0.95))
  median_HHCI = apply(tmp, 2, function(x) quantile(x, probs =  0.975))
  names(median) <- cpgs
  names(median_LLCI) <- cpgs
  names(median_LCI) <- cpgs
  names(median_HCI) <- cpgs
  names(median_HHCI) <- cpgs
  median = as.matrix(median)
  median_LCI = as.matrix(median_LCI)
  median_LLCI = as.matrix(median_LLCI)
  median_HCI = as.matrix(median_HCI)
  median_HHCI = as.matrix(median_HHCI)
  median_LCI <- cbind(median_LLCI, median_LCI)
  median_HCI <- cbind(median_HCI, median_HHCI)
  betas <- cbind(median,cbind(median_LCI,median_HCI))
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(6,1,2,3,4,5)]
  names(betas)[2:6] <- c("Median_Beta", "Beta_2.5", "Beta_5", "Beta_95", "Beta_97.5")
  comp = fread(paste0("./Cog_Comp/", i))
  comp = as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(betas, pip[,2])
  print(A)
  write.csv(betas1, paste0("Cog_EWAS_BayesR/Cog_Summstats/", A), row.names = F)
} 


# Calculate Mean Betas (PRS Prediction and summstat upload)
for(i in loop){ 
  tmp = fread(paste0("Cog_Beta/", i))
  tmp = as.data.frame(tmp)
  print(i)
  A <-gsub(".csv*","",i)
  means = apply(tmp, 2, mean)
  ses = apply(tmp, 2, function(x){sd(x)/sqrt(length(x))})
  names(means) <- cpgs
  names(ses) <- cpgs
  means = as.matrix(means)
  ses = as.matrix(ses)  
  betas <- cbind(means,ses)
  betas <- as.data.frame(betas)
  betas$CpG <- row.names(betas)
  betas <- betas[,c(3,1,2)]
  names(betas)[2:3] <- c("Mean_Beta", "SE")
  comp = fread(paste0("./Cog_Comp/", i))
  comp = as.data.frame(comp)
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- cpgs
  betas1 <- cbind(betas, pip[,2])
  print(A)
  write.table(betas1, paste0("Cog_EWAS_BayesR/Cog_Summstats/", A, "_Mean_Beta_PIP.txt"), row.names = F)
} 

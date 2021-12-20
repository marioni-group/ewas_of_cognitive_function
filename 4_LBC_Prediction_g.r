########################################################
### GS Age demographics for results in the main text ###
########################################################

d <- read.csv("Cognitive_BayesR_10k_phenotypes.csv")
summary(d$age)
sd(d$age)

###########################################################################################################
### Script to generate DNAm predictor of g in LBC1936 at wave 1 (age 70) and LBC1921 at wave 1 (age 79) ###
###########################################################################################################

### read in weights ###
g_wts <- read.table("Cog_Summstats/g_processed_Mean_Beta_PIP.txt", header=T, stringsAsFactors=F)
meanBetas <- g_wts[,1:2]

### load LBC target file ###
d = readRDS("targets_3489_bloodonly.rds")
d36_w1 <- d[d$cohort=="LBC36" & d$WAVE==1 & d$set==1,]
d21_w1 <- d[d$cohort=="LBC21" & d$WAVE==1 & d$set==1,]
 
### load in methylation data ###
dat <- readRDS("LBC_betas_3489_bloodonly.rds")
meth = t(dat)
meth1 = as.data.frame(meth)
meth1$id = as.character(rownames(meth1))

### subset DNAm object to LBC36 wave 1 data ###
tmp36 <- which(rownames(meth1) %in% d36_w1$Basename)
meth36 <- meth1[tmp36,]

tmp21 <- which(rownames(meth1) %in% d21_w1$Basename)
meth21 <- meth1[tmp21,]

### subset to relevant CpG sites ###
a36 = which(names(meth36) %in% meanBetas$CpG)
meth36a <- meth36[,a36]

a21 = which(names(meth21) %in% meanBetas$CpG)
meth21a <- meth21[,a21]

### scale methylation data ###
meth36b <- scale(meth36a)
meth21b <- scale(meth21a)

### replace missing values with 0 ###
meth36b[is.na(meth36b)] <- 0
meth21b[is.na(meth21b)] <- 0

### line up weights and CpGs ###
b36 = which(meanBetas$CpG %in% colnames(meth36b)) 
mean_betas <- meanBetas[b36,]
meth36c <- meth36b[,match(colnames(meth36b), mean_betas$CpG)]
meth36d <- meth36c[,mean_betas$CpG]

b21 = which(meanBetas$CpG %in% colnames(meth21b)) 
mean_betas <- meanBetas[b21,]
meth21c <- meth21b[,match(colnames(meth21b), mean_betas$CpG)]
meth21d <- meth21c[,mean_betas$CpG]

### create predictor ###
g_pred36 <- meth36d %*% mean_betas$Mean_Beta
pred36 <- as.data.frame(g_pred36)
names(pred36) <- c("g_pred")
pred36$Basename <- rownames(pred36)
d36 <- merge(d36_w1, pred36, by="Basename")

g_pred21 <- meth21d %*% mean_betas$Mean_Beta
pred21 <- as.data.frame(g_pred21)
names(pred21) <- c("g_pred")
pred21$Basename <- rownames(pred21)
d21 <- merge(d21_w1, pred21, by="Basename")


#######################################
### LBC1936 phenotypes and analyses ###
#######################################

### load in phenotypes ###
library("foreign")
ph36 <- read.spss("../Riccardo/Cognitive_EWAS_GS_LBC/LBC1936_BloodBasedEWASofCognitiveAbilities_RM_06MAY2021.sav", to.data.frame=T, use.value.labels=F)

### take max grip strength, recode NA SIMD (99999) value ###
ph36$grip <- pmax(ph36$griplh_w1, ph36$griprh_w1)
ph36$depind_w1[ph36$depind_w1==99999] <- NA

### remove outliers (based on visual inspection of histograms) and impossible values ###
ph36$symsear_w1[ph36$symsear_w1<0] <- NA

### merge phenos with DNAm predictor ###
d36a = merge(d36, ph36[,c(1,4:17,20:23)], by.x="ID_raw", by.y="lbc36no")

### Merge in PRS derived from Davies et al. 2018 Nat Comms Z-scores - retain predictor with P<1 ###
prs36 = read.table("Cog_EWAS_BayesR/PRS/cog_prs.all.score", header=T)
tmp = merge(d36a, prs36[,c("IID","X1.000000")], by.x="ID_raw", by.y="IID", all.x=T)
tmp1 <- tmp[,c(4,6,14,19:38)]

### remove anyone with missing cognitive data ###
a <- rowSums(tmp1[,5:11])
w1_36 <- tmp1[which(!is.na(a)), ]

### adjust lung function for age, sex, and height ###
w1_36$lung <- resid(lm(fev_w1 ~ age + sex + height_w1, na.action=na.exclude, data=w1_36))

### PCA of cognitive test data at wave 1 of LBC36 ###
library(psych)
w1_36$g <- scale(principal(w1_36[,5:11], rotate="none", nfactors=1)$scores)

### PCA output for variance explained by g ###
principal(w1_36[,5:11], rotate="none", nfactors=7)


#######################################
### Main LBC1936 Prediction Results ###
#######################################

### correlation and variance explained ###
r <- cor(w1_36$g, w1_36$g_pred, use="pairwise.complete.obs")
r

### incremental DNAm R2 ###
null <- summary(lm(g ~ age + sex, data=w1_36))$r.squared
full <- summary(lm(g ~ age + sex + g_pred, data=w1_36))$r.squared
round(100*(full - null), 3)
# [1] 3.426

### incremental PRS R2 ###
full2 <- summary(lm(g ~ age + sex + X1.000000, data=w1_36))$r.squared
round(100*(full2 - null), 3)
# [1] 7.266


### incremental DNAm + PRS R2 ###
full3 <- summary(lm(g ~ age + sex + g_pred + X1.000000, data=w1_36))$r.squared
round(100*(full3 - null), 3)
# [1] 10.672

### get pvals ###
summary(lm(g ~ age + sex + g_pred, data=w1_36))
summary(lm(g ~ age + sex + g_pred + X1.000000, data=w1_36))

### difference in g between top and bottom decile of DNAm_g ###
w1_36$g_pred_dec <- as.numeric(cut(w1_36$g_pred, quantile(w1_36$g_pred, prob = seq(0, 1, length = 11)), include.lowest=T ))
summary(lm(g ~ age + sex + factor(g_pred_dec), data=w1_36))


### Fully adjusted incremental DNAm R2 ###
nullf <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung, data=w1_36))$r.squared
fullf <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + g_pred, data=w1_36))$r.squared
round(100*(fullf - nullf), 3)

### Fully adjusted incremental PRS R2 ###
full2f <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + X1.000000, data=w1_36))$r.squared
round(100*(full2f - nullf), 3)

### Fully adjusted incremental DNAm + PRS R2 ###
full3f <- summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + g_pred + X1.000000, data=w1_36))$r.squared
round(100*(full3f - nullf), 3)

### Get pvals for EpiScore and PGRS ###
summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + g_pred, data=w1_36))
summary(lm(g ~ age + sex + log(bmi_w1) + yrsedu_w1 + factor(smokcat_w1) + log(alcunitwk_w1+1) + hibp_w1 + hadsd_w1 + depind_w1 + log(sixmwk_w1) + grip + lung + X1.000000, data=w1_36))

#######################################
### LBC1921 phenotypes and analyses ###
#######################################

### load in phenotypes ###
pheno21 <- read.spss("../Riccardo/Cognitive_EWAS_GS_LBC/LBC1921_BloodBasedEWASofCognitiveAbilities_RM_06MAY2021.sav", to.data.frame=T, use.value.labels=F)

### recode unsure history of hibp as missing (n=4) ###
pheno21$hyphist[pheno21$hyphist==3] <- NA

### merge phenos with DNAm predictor ###
d21a = merge(d21, pheno21[,c(1,4:10,12,13,15:20)], by.x="ID_raw", by.y="studyno")

### Merge in PRS derived from Davies et al. 2018 Nat Comms Z-scores - retain predictor with P<1 ###
prs21 = read.table("Cog_EWAS_BayesR/PRS/cog_prs_lbc21.all.score", header=T)
lbc_bld = read.table("../Riccardo/LBC1921_FullSampleSex_RM_15JUL2013.dat", header=T, sep="\t")
prs21a <- merge(prs21, lbc_bld[,c(1,3)], by.x="IID", by.y="bloodnum")

tmp = merge(d21a, prs21a[,c("studyno","X1.000000")], by.x="ID_raw", by.y="studyno", all.x=T)

tmp1 <- tmp[,c(4,6,14,19:35)]

### remove anyone with missing cognitive data ###
a <- rowSums(tmp1[,5:8])
w1_21 <- tmp1[which(!is.na(a)), ]

### adjust lung function for age, sex, and height ###
w1_21$lung <- resid(lm(fev ~ age + sex + height, na.action=na.exclude, data=w1_21))

### PCA of cognitive test data at wave 1 of LBC36 ###
library(psych)
w1_21$g <- scale(principal(w1_21[,5:8], rotate="none", nfactors=1)$scores)

### PCA output for variance explained by g ###
principal(w1_21[,5:8], rotate="none", nfactors=4)

#######################################
### Main LBC1921 Prediction Results ###
#######################################

### correlation and variance explained ###
r <- cor(w1_21$g, w1_21$g_pred, use="pairwise.complete.obs")
r

### incremental DNAm R2 ###
null <- summary(lm(g ~ age + sex, data=w1_21))$r.squared
full <- summary(lm(g ~ age + sex + g_pred, data=w1_21))$r.squared
round(100*(full - null), 3)
# [1] 4.478

### incremental PRS R2 ###
full2 <- summary(lm(g ~ age + sex + X1.000000, data=w1_21))$r.squared
round(100*(full2 - null), 3)
# [1] 6.928

### incremental DNAm + PRS R2 ###
full3 <- summary(lm(g ~ age + sex + g_pred + X1.000000, data=w1_21))$r.squared
round(100*(full3 - null), 3)
# [1] 10.546

### get pvals ###
summary(lm(g ~ age + sex + g_pred, data=w1_21))
summary(lm(g ~ age + sex + X1.000000, data=w1_21))

### difference in g between top and bottom quintile of DNAm_g ###
w1_21$g_pred_dec <- as.numeric(cut(w1_21$g_pred, quantile(w1_21$g_pred, prob = seq(0, 1, length = 6)), include.lowest=T ))
summary(lm(g ~ age + sex + factor(g_pred_dec), data=w1_21))



### Fully adjusted incremental DNAm R2 ###
nullf <- summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung, data=w1_21))$r.squared
fullf <- summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung + g_pred, data=w1_21))$r.squared
round(100*(fullf - nullf), 3)

### Fully adjusted incremental PRS R2 ###
full2f <- summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung + X1.000000, data=w1_21))$r.squared
round(100*(full2f - nullf), 3)

### Fully adjusted incremental DNAm + PRS R2 ###
full3f <- summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung + g_pred + X1.000000, data=w1_21))$r.squared
round(100*(full3f - nullf), 3)


### get pvals ###
summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung + g_pred, data=w1_21))
summary(lm(g ~ age + sex + log(bmi) + yrseduc + factor(smoker) + log(alcpw+1) + factor(hyphist) + HADSD + soclcode + log(sixmtime) + gripstr + lung + X1.000000, data=w1_21))






#############################
### Supplementary Table 8 ###
#############################

### Regression of covariates on DNAm g score ###
summary(lm(scale(g_pred) ~ scale(age) + sex + scale(yrsedu_w1) + scale(log(bmi_w1)) + factor(smokcat_w1) + scale(log(alcunitwk_w1+1)) + hibp_w1 + scale(hadsd_w1) + scale(depind_w1) + scale(log(sixmwk_w1)) + scale(grip) + scale(lung), data=w1_36))$coefficients[,c(1,2,4)]
nobs(lm(scale(g_pred) ~ scale(age) + sex + yrsedu_w1 + scale(bmi_w1) + factor(smokcat_w1) + scale(alcunitwk_w1) + hibp_w1 + scale(hadsd_w1) + scale(depind_w1) + scale(sixmwk_w1) + scale(grip) + scale(lung), data=w1_36))

summary(lm(scale(g_pred) ~ scale(age) + sex + scale(yrsedu_w1) + scale(log(bmi_w1)) + factor(smokcat_w1) + scale(log(alcunitwk_w1+1)) + hibp_w1 + scale(hadsd_w1) + scale(depind_w1) + scale(log(sixmwk_w1)) + scale(grip) + scale(lung), data=w1_36))$r.squared
summary(lm(scale(g_pred) ~ scale(age) + sex, data=w1_36))$r.squared

w1_21$smoker1 <- relevel(as.factor(ifelse(w1_21$smoker==1, "current", ifelse(w1_21$smoker==2, "never", ifelse(w1_21$smoker==3, "ex", "NA")))), ref="never")
summary(lm(scale(g_pred) ~ scale(age) + sex + scale(yrseduc) + scale(log(bmi)) + factor(smoker1) + scale(log(alcpw+1)) + factor(hyphist) + scale(HADSD) + scale(soclcode) + scale(log(sixmtime)) + scale(gripstr) + scale(lung), data=w1_21))$coefficients[,c(1,2,4)]
nobs(lm(g ~ age + sex + bmi + yrseduc + factor(smoker1) + alcpw + factor(hyphist) + HADSD + soclcode + sixmtime + gripstr + lung + g_pred + X1.000000, data=w1_21))

summary(lm(scale(g_pred) ~ scale(age) + sex + scale(yrseduc) + scale(log(bmi)) + factor(smoker) + scale(log(alcpw+1)) + factor(hyphist) + scale(HADSD) + scale(soclcode) + scale(log(sixmtime)) + scale(gripstr) + scale(lung), data=w1_21))$r.squared
summary(lm(scale(g_pred) ~ scale(age) + sex, data=w1_21))$r.squared



#############################
### Supplementary Table X ###
#############################

### Correlation between DNAm g score and measured blood cell proportions ###
targets = readRDS("targets_3489_bloodonly.rds")
targets = targets[which(targets$WAVE==1),]

w1_36 = merge(targets[,c("ID_raw", "neut", "lymph", "mono", "eosin", "baso")], w1_36, by.x="ID_raw", by.y="ID")

blood_cor = cor(scale(w1_36[,c("neut", "lymph", "mono", "eosin", "baso", "g_pred")]), use="pairwise.complete.obs")
blood_cor[upper.tri(blood_cor)] = NA
write.table(blood_cor, file="Supplementary_Table_BloodCount_Correlations.xls", quote=F, sep='\t', col.names=NA)

blood_cor36 = data.frame(Cell= c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils"), cor=rep(NA, 5), p = rep(NA, 5))
rownames(blood_cor36) = c("neut", "lymph", "mono", "eosin", "baso")
for(i in rownames(blood_cor36)){
blood_cor36[i, "cor"] = signif(cor.test(w1_36[,i], w1_36[,"g_pred"], use="pairwise.complete.obs")$estimate,2)
blood_cor36[i, "p"] = signif(cor.test(w1_36[,i], w1_36[,"g_pred"], use="pairwise.complete.obs")$p.value,2)
}
# write.table(blood_cor36, file="Supplementary_Table_BloodCount_Correlations_v2.xls", quote=F, sep='\t', col.names=NA)
blood_cor36
             # Cell     cor    p
# neut  Neutrophils -0.0160 0.65
# lymph Lymphocytes -0.0390 0.26
# mono    Monocytes  0.0027 0.94
# eosin Eosinophils  0.0170 0.61
# baso    Basophils -0.0680 0.05



w1_21 = merge(targets[,c("ID_raw", "neut", "lymph", "mono", "eosin", "baso")], w1_21, by.x="ID_raw", by.y="ID")
for(i in c("neut", "lymph", "mono", "eosin", "baso")){
w1_21[which(w1_21[,i] == 9999),i] = NA # Remove unmeasured props
}
blood_cor = cor(scale(w1_21[,c("neut", "lymph", "mono", "eosin", "baso", "g_pred")]), use="pairwise.complete.obs")
blood_cor[upper.tri(blood_cor)] = NA
write.table(blood_cor, file="Supplementary_Table_BloodCount_Correlations.xls", quote=F, sep='\t', col.names=NA)

blood_cor21 = data.frame(Cell= c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils"), cor=rep(NA, 5), p = rep(NA, 5))
rownames(blood_cor21) = c("neut", "lymph", "mono", "eosin", "baso")
for(i in rownames(blood_cor21)){
blood_cor21[i, "cor"] = signif(cor.test(w1_21[,i], w1_21[,"g_pred"], use="pairwise.complete.obs")$estimate,2)
blood_cor21[i, "p"] = signif(cor.test(w1_21[,i], w1_21[,"g_pred"], use="pairwise.complete.obs")$p.value,2)
}
blood_cor21
             # Cell     cor    p
# neut  Neutrophils  0.0720 0.15
# lymph Lymphocytes -0.0029 0.95
# mono    Monocytes -0.0630 0.20
# eosin Eosinophils -0.0620 0.21
# baso    Basophils  0.0440 0.37


# write.table(rbind(blood_cor36 blood_cor21), file="Supplementary_Table_BloodCount_Correlations_v2.xls", quote=F, sep='\t', col.names=NA)


###############################################
### Supplementary Table 8: Cohort summaries ###
###############################################
lbc21_demos1 = matrix(nrow=9, ncol=3)
rownames(lbc21_demos1) = c("age", "bmi", "alcpw" ,"soclcode" ,"sixmtime", "yrseduc", "HADSD", "gripstr", "fev")
colnames(lbc21_demos1) = c("mean", "sd", "n")
for(i in rownames(lbc21_demos1)){
 lbc21_demos1[i,"mean"] = signif(mean(w1_21[,i], na.rm=T),4)
 lbc21_demos1[i,"sd"] = signif(sd(w1_21[,i], na.rm=T),4)
 lbc21_demos1[i,"n"] = length(which(!is.na(w1_21[,i])))
 }
 lbc21_demos1
 
            # mean     sd   n
# age      79.130 0.5751 427
# bmi      26.190 4.0470 424
# alcpw     5.468 8.9910 427
# soclcode  2.216 0.8840 426
# sixmtime  4.712 1.9240 423
# yrseduc  11.000 2.5170 425
# HADSD     3.522 2.2920 425
# gripstr  26.200 8.9080 424
# fev       1.884 0.6226 424

 
 props = list()
 for(i in c("hyphist", "smoker", "sex")){
 props[[i]] = rbind(prop.table(table(w1_21[,i])), table(w1_21[,i]))
 }
props
# $hyphist
               # 0           1
# [1,]   0.5947867   0.4052133
# [2,] 251.0000000 171.0000000

# $smoker
               # 1           2           3
# [1,]  0.06807512   0.4507042   0.4812207
# [2,] 29.00000000 192.0000000 205.0000000

# $sex
               # F           M
# [1,]   0.6042155   0.3957845
# [2,] 258.0000000 169.0000000




lbc36_demos1 = matrix(nrow=9, ncol=3)
rownames(lbc36_demos1) = c("age", "bmi_w1", "alcunitwk_w1" ,"depind_w1" ,"sixmwk_w1", "yrsedu_w1", "hadsd_w1", "grip", "fev_w1")
colnames(lbc36_demos1) = c("mean", "sd", "n")
for(i in rownames(lbc36_demos1)){
 lbc36_demos1[i,"mean"] = signif(mean(w1_36[,i], na.rm=T),4)
 lbc36_demos1[i,"sd"] = signif(sd(w1_36[,i], na.rm=T),4)
 lbc36_demos1[i,"n"] = length(which(!is.na(w1_36[,i])))
 }
 
 props = list()
 for(i in c("hibp_w1", "smokcat_w1", "sex")){
 props[[i]] = rbind(prop.table(table(w1_36[,i])), table(w1_36[,i]))
 }

lbc36_demos1
                 # mean        sd   n
# age            69.590    0.8301 844
# bmi_w1         27.780    4.3400 843
# alcunitwk_w1   10.190   13.7100 844
# depind_w1    4652.000 1882.0000 838
# sixmwk_w1       3.844    1.2070 840
# yrsedu_w1      10.740    1.1110 844
# hadsd_w1        2.728    2.1810 842
# grip           29.550   10.1500 840
# fev_w1          2.364    0.6718 842

props
# $hibp_w1
               # 0           1
# [1,]   0.5876777   0.4123223
# [2,] 496.0000000 348.0000000

# $smokcat_w1
               # 0           1          2
# [1,]   0.4668246   0.4206161  0.1125592
# [2,] 394.0000000 355.0000000 95.0000000

# $sex
               # F           M
# [1,]   0.4952607   0.5047393
# [2,] 418.0000000 426.0000000


### Table S10: Model Outputs ##
w1_36$ever_smoke = ifelse(w1_36$smokcat_w1==0, 0, 1)
phenos_36 = c("bmi_w1", "alcunitwk_w1", "hadsd_w1", "depind_w1", "sixmwk_w1", "yrsedu_w1", "grip", "lung", "ever_smoke", "hibp_w1")
s10_36_epi = s10_36_g = data.frame(Trait=phenos_36, Effect = NA, SE = NA, P= NA)
for(i in phenos_36){
# if(i == "ever_smoke" | i == "hibp_w1"){
# s10_36_epi[which(s10_36_epi$Trait==i), c("Effect" ,"SE", "P")] = 
# summary(glm(get(i) ~ scale(g_pred) + (age) + sex, data=w1_36, family="binomial"))$coef[2,c(1,2,4)]
# } else { 
s10_36_epi[which(s10_36_epi$Trait==i), c("Effect" ,"SE", "P")] = 
summary(lm(scale(get(i)) ~ scale(g_pred) + (age) + sex, data=w1_36))$coef[2,c(1,2,4)]
}
# } # Using standard lm for ever smoke/htn


for(i in phenos_36){
# if(i == "ever_smoke" | i == "hibp_w1"){
# s10_36_g[which(s10_36_epi$Trait==i), c("Effect" ,"SE", "P")] = 
# summary(glm(get(i) ~ scale(g) + (age) + sex, data=w1_36, family="binomial"))$coef[2,c(1,2,4)]} 
# else {
s10_36_g[which(s10_36_epi$Trait==i), c("Effect" ,"SE", "P")] = 
summary(lm(scale(get(i)) ~ scale(g) + (age) + sex, data=w1_36))$coef[2,c(1,2,4)]
}
# } # Using standard lm for ever smoke/htn


w1_21$ever_smoke = ifelse(w1_21$smoker==3, 0, 1)
phenos_21 = c("bmi", "alcpw", "HADSD", "soclcode", "sixmtime", "yrseduc", "gripstr", "lung", "ever_smoke", "hyphist")
s10_21_epi = s10_21_g = data.frame(Trait=phenos_21, Effect = NA, SE = NA, P= NA)
for(i in phenos_21){
# if(i == "ever_smoke" | i == "hyphist"){
# s10_21_epi[which(s10_21_epi$Trait==i), c("Effect" ,"SE", "P")] = 
# summary(glm(get(i) ~ scale(g_pred) + (age) + sex, data=w1_21, family="binomial"))$coef[2,c(1,2,4)]
# } else{
s10_21_epi[which(s10_21_epi$Trait==i), c("Effect" ,"SE", "P")] = 
summary(lm(scale(get(i)) ~ scale(g_pred) + (age) + sex, data=w1_21))$coef[2,c(1,2,4)]
}
# } # Using standard lm for ever smoke/htn for plotting purposes 

for(i in phenos_21){
# if(i == "ever_smoke" | i == "hyphist"){
# s10_21_g[which(s10_21_epi$Trait==i), c("Effect" ,"SE", "P")] = 
# summary(glm(get(i) ~ scale(g) + age + sex, data=w1_21, family="binomial"))$coef[2,c(1,2,4)]
# } else {
s10_21_g[which(s10_21_epi$Trait==i), c("Effect" ,"SE", "P")] = 
summary(lm(scale(get(i)) ~ scale(g) + (age) + sex, data=w1_21))$coef[2,c(1,2,4)]
}
# } # Using standard lm for ever smoke/htn for plotting purposes  


saveRDS(list(g_pred=s10_36_epi, measured=s10_36_g), file="Figure_1B_plotdata.rds")


s10_21_epi
        # Trait       Effect         SE            P
# 1         bmi -0.069551912 0.04884927 0.1552449339
# 2       alcpw  0.111105028 0.04667710 0.0177411663
# 3       HADSD -0.058205400 0.04900844 0.2356367578
# 4    soclcode -0.144236134 0.04796327 0.0027941103
# 5    sixmtime -0.115692476 0.04815796 0.0167241696
# 6     yrseduc  0.097762359 0.04802758 0.0424216503
# 7     gripstr  0.108702546 0.03195502 0.0007337646
# 8        lung  0.102176826 0.04872908 0.0366057146
# 9  ever_smoke -0.112049819 0.04821127 0.0205918726
# 10    hyphist -0.003447182 0.04882353 0.9437459576

s10_21_g
        # Trait      Effect         SE            P
# 1         bmi -0.07348844 0.04954561 1.387585e-01
# 2       alcpw  0.08382065 0.04675228 7.370844e-02
# 3       HADSD -0.08682352 0.04886691 7.633394e-02
# 4    soclcode -0.43365200 0.04359534 4.394887e-21
# 5    sixmtime -0.16680624 0.04847676 6.378775e-04
# 6     yrseduc  0.47086876 0.04238993 2.576224e-25
# 7     gripstr  0.14459004 0.03209439 8.608497e-06
# 8        lung  0.19241150 0.04879689 9.422033e-05
# 9  ever_smoke -0.01613459 0.04852715 7.396875e-01
# 10    hyphist  0.02698143 0.04917759 5.835365e-01

s10_36_epi
          # Trait      Effect         SE            P
# 1        bmi_w1 -0.09478764 0.03450534 6.142426e-03
# 2  alcunitwk_w1  0.08896867 0.03314304 7.409462e-03
# 3      hadsd_w1 -0.12192227 0.03455909 4.415901e-04
# 4     depind_w1  0.13300955 0.03387158 9.316508e-05
# 5     sixmwk_w1 -0.09035396 0.03392959 7.894030e-03
# 6     yrsedu_w1  0.12846452 0.03426174 1.893285e-04
# 7          grip  0.07169225 0.02129540 7.958581e-04
# 8          lung  0.06663276 0.03473194 5.538883e-02
# 9    ever_smoke -0.03257296 0.03443427 3.444481e-01
# 10      hibp_w1 -0.03233753 0.03464007 3.508148e-01

s10_36_g
          # Trait      Effect         SE            P
# 1        bmi_w1 -0.15317978 0.03527350 1.579948e-05
# 2  alcunitwk_w1  0.06291547 0.03420490 6.621351e-02
# 3      hadsd_w1 -0.19817179 0.03515368 2.360477e-08
# 4     depind_w1  0.34160281 0.03309180 1.360961e-23
# 5     sixmwk_w1 -0.24598546 0.03417659 1.365886e-12
# 6     yrsedu_w1  0.45491377 0.03192359 1.992752e-41
# 7          grip  0.14196279 0.02145202 6.513466e-11
# 8          lung  0.17596748 0.03534844 7.796731e-07
# 9    ever_smoke -0.09089817 0.03533721 1.027306e-02
# 10      hibp_w1 -0.09824538 0.03552629 5.809516e-03

# Saturated model(Table S11)
summary(lm(scale(g) ~ scale(age) + 
                      sex + 
					  scale(g_pred) +
					  scale(yrsedu_w1) + 
					  scale(log(bmi_w1)) + 
					  factor(ever_smoke) + 
					  scale(log(alcunitwk_w1+1)) + 
					  factor(hibp_w1) + 
					  scale(hadsd_w1) + 
					  scale(depind_w1) + 
					  scale(log(sixmwk_w1)) + 
					  scale(grip) + 
					  scale(lung), data=w1_36))$coefficients

                               # Estimate Std. Error    t value     Pr(>|t|)
# (Intercept)                   0.16607076 0.06864364  2.4193174 1.576871e-02
# scale(age)                   -0.17194102 0.02965858 -5.7973446 9.646459e-09
# sexM                         -0.26466185 0.09559022 -2.7687126 5.756251e-03
# scale(g_pred)                 0.09252673 0.02955899  3.1302400 1.809213e-03
# scale(yrsedu_w1)              0.32282799 0.03099240 10.4163589 6.238455e-24
# scale(log(bmi_w1))           -0.03333052 0.03057144 -1.0902502 2.759268e-01
# factor(ever_smoke)1           0.01164714 0.06064919  0.1920412 8.477580e-01
# scale(log(alcunitwk_w1 + 1))  0.01170504 0.03110484  0.3763093 7.067854e-01
# factor(hibp_w1)1             -0.08615049 0.05948416 -1.4482928 1.479217e-01
# scale(hadsd_w1)              -0.07692980 0.03025147 -2.5430101 1.117493e-02
# scale(depind_w1)              0.17671204 0.03192641  5.5349799 4.203865e-08
# scale(log(sixmwk_w1))        -0.06693662 0.03344467 -2.0014138 4.568055e-02
# scale(grip)                   0.14226201 0.04979344  2.8570434 4.385472e-03
# scale(lung)                   0.03171145 0.03113443  1.0185330 3.087285e-01


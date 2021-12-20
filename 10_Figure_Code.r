##############################
### Supplementary Figure 1 ###
##############################

d = read.csv("Cog_EWAS_BayesR/Cognitive_BayesR_10k_phenotypes.csv")
mean(d$EpiScore)
sd(d$EpiScore)
length(!is.na(d$EpiScore))

smk <- read.csv("../Marioni_Group/GS/GS_dataset/ever_smoke.csv")
smk1 <- merge(d, smk, by="Sample_Name")
smk1$smk_cat <- smk1$ever_smoke
smk1$smk_cat[smk1$smk_cat==5] <- NA
smk1$smk_cat[smk1$smk_cat==1] <- "Current \nN=1,530"
smk1$smk_cat[smk1$smk_cat==2] <- "Former (<12 months)\nN=218"
smk1$smk_cat[smk1$smk_cat==3] <- "Former (>12 months)\nN=2,513"
smk1$smk_cat[smk1$smk_cat==4] <- "Never\nN=4,650"
smk1$smk_cat[is.na(smk1$smk_cat)] <- "Unknown\nN=251"
table(smk1$smk_cat, smk1$ever_smoke)
  
library(ggplot2)

ggplot(smk1, aes(x=smk_cat, y=EpiScore)) + 
	geom_violin(trim=FALSE, fill="gray")+
	labs(title="EpiSmoker vs Self-report Smoking",x="Self-reported Smoking", y = "EpiSmoker Score")+
	geom_boxplot(width=0.1)+
	theme_classic() +
	theme(axis.title.x = element_text(size = 16),
		axis.title.y = element_text(size = 16), 
		axis.text=element_text(size=14),
		plot.title=element_text(size = 16))
  

##############################
### Supplementary Figure 2 ###
##############################

group_pips = readRDS("Cog_EWAS_BayesR/group_pips.rds")
pips = readRDS("pips.rds")

# Left Panel
sig_probes = unique(unlist(lapply(pips, function(x){x$CpG[which(x$PIP>0.8)]})))
sigs = list()
for(i in names(pips)){
  tmp = pips[[i]][which(pips[[i]]$CpG %in% sig_probes),] 
  sigs[[i]] = pips[[i]][which(pips[[i]]$CpG %in% tmp$CpG),]
}
sigs = do.call("rbind", sigs)
plotdat = sigs
plotdat$Trait = gsub("_processed*.*", "", rownames(plotdat))
plotdat$Beta = plotdat$Median_Beta
plotdat$LCI = plotdat$Beta_2.5
plotdat$HCI = plotdat$Beta_97.5
require(RColorBrewer)
colours <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # Colourblind-friendly
colvalues=colours[1:6]

plotdat$lci <- plotdat$Beta_5
plotdat$uci <- plotdat$Beta_95

plotdat$CpG1 <- factor(plotdat$CpG, levels=c("cg05325763", "cg13108341", "cg09984392"))
plotdat$Trait1 <- factor(plotdat$Trait, levels=rev(c("Digit", "LM", "Verbal", "Vocabulary","gf","g")))
levels(plotdat$Trait1)[levels(plotdat$Trait1)=="Digit"] <- "Digit Symbol"
levels(plotdat$Trait1)[levels(plotdat$Trait1)=="LM"] <- "Logical Memory"
levels(plotdat$Trait1)[levels(plotdat$Trait1)=="Verbal"] <- "Verbal Fluency"

s2a = ggplot(plotdat, aes(x=CpG1, y=Beta, ymin=lci, ymax=uci, col=Trait1, fill=Trait1)) + 
         geom_point(size=4, position=position_dodge(width = 0.9)) + 
         geom_linerange(size=2, position=position_dodge(width = 0.9), show.legend=FALSE) +
         geom_linerange(aes(ymin=LCI, ymax=HCI), size=1, position=position_dodge(width = 0.9), show.legend=FALSE) +
         scale_fill_manual(values=colvalues, breaks=rev(levels(plotdat$Trait1))) +
         scale_color_manual(values=colvalues, breaks=rev(levels(plotdat$Trait1))) +
         coord_flip() + 
         ylim(-0.13,0.11) + 
         ggtitle("Median Betas and 95% Credible Intervals") + 
         theme_light() + 
		 ylab("Effect Size (Beta)") +
		 xlab("CpG") +
		 theme(axis.title.x = element_text(size = 16),
		       axis.title.y = element_text(size = 16), 
               plot.title=element_text(size = 16),
			   axis.text=element_text(size=16),
               legend.text=element_text(size=16),
			   legend.position=c(0.8,0.2),
			   legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
			   legend.title=element_blank())

# Right Panel
sigs = list()
for(i in names(pips)){
  tmp = pips[[i]][which(pips[[i]]$PIP > 0.8),] # Lead CpGs with group pips>0.8
  sigs[[i]] = pips[[i]][which(pips[[i]]$CpG %in% tmp$CpG),]
}
sigs = do.call("rbind", sigs)
plotdat = sigs
plotdat$Trait = gsub("_processed*.*", "", rownames(plotdat))
plotdat$Beta = plotdat$Median_Beta
plotdat$LCI = plotdat$Beta_2.5
plotdat$HCI = plotdat$Beta_97.5
require(RColorBrewer)
colours <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # Colourblind-friendly
colvalues=colours[c(1,3:6)]
require(ggplot2)

plotdat$lci <- plotdat$Beta_5
plotdat$uci <- plotdat$Beta_95

plotdat$CpG1 <- factor(plotdat$CpG, levels=c("cg05325763", "cg13108341", "cg09984392"))
plotdat$Trait1 <- factor(plotdat$Trait, levels=rev(c("Digit", "Verbal", "Vocabulary","gf","g")))
levels(plotdat$Trait1)[levels(plotdat$Trait1)=="Digit"] <- "Digit Symbol"
levels(plotdat$Trait1)[levels(plotdat$Trait1)=="Verbal"] <- "Verbal Fluency"

s2b = ggplot(plotdat, aes(x=CpG1, y=Beta, ymin=lci, ymax=uci, col=Trait1, fill=Trait1)) + 
         geom_point(size=4, position=position_dodge(width = 0.9)) + 
         geom_linerange(size=2, position=position_dodge(width = 0.9), show.legend=FALSE) +
         geom_linerange(aes(ymin=LCI, ymax=HCI), size=1, position=position_dodge(width = 0.9), show.legend=FALSE) +
         scale_fill_manual(values=colvalues, breaks=rev(levels(plotdat$Trait1))) +
         scale_color_manual(values=colvalues, breaks=rev(levels(plotdat$Trait1))) +
         coord_flip() + 
         ylim(-0.13,0.11) + 
         ggtitle("Median Betas and 95% Credible Intervals") + 
         theme_light() + 
		 ylab("Effect Size (Beta)") +
		 xlab("CpG") +
		 theme(axis.title.x = element_text(size = 16),
		       axis.title.y = element_text(size = 16), 
               plot.title=element_text(size = 16),
			   axis.text=element_text(size=16),
               legend.text=element_text(size=16),
			   legend.position=c(0.8,0.2),
			   legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
			   legend.title=element_blank())

require(gridExtra)
png("Supplementary_Figure2_ColourBlindFriendly.png", width=1200, height=600)
grid.arrange(s2a, s2b, nrow=1)
dev.off()

##############################
### Supplementary Figure 3 ###
##############################
w1_36 <- readRDS("../Riccardo/Cog_EWAS_BayesR_GS_PGRS_to_LBC/w1_36.rds")
p36 <- ggplot(w1_36, aes(x=g_pred, y=g)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() +
  labs(title="LBC1936",x="Epigenetic g Score", y = "Measured g Score") +
  theme(axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text=element_text(size=14),
  plot.title=element_text(size = 16))

w1_21 <- readRDS("../Riccardo/Cog_EWAS_BayesR_GS_PGRS_to_LBC/w1_21.rds")
p21 <- ggplot(w1_21, aes(x=g_pred, y=g)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme_classic() +
  labs(title="LBC1921",x="Epigenetic g Score", y = "Measured g Score") +
  theme(axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16), 
  axis.text=element_text(size=14),
  plot.title=element_text(size = 16))
 
require(gridExtra)
#tiff("Figure_2_April2021.tiff", width=800, height=800)
grid.arrange(p36, p21, nrow=2)
#dev.off() 





 
  
###  Figure 1:
# 1A:
r2_36 <- c(3.4, 7.3, 10.7) # Obtained from LBC_Prediction_g.R (Incremental_DNAm_R2)
ph <- c("DNAm","SNP","DNAm + SNP")
fig_dat36 <- data.frame(cbind(ph, r2_36))
fig_dat36$r2_36 <- as.numeric(fig_dat36$r2_36)
fig_dat36$ph1 <- factor(fig_dat36$ph, levels=c("DNAm","SNP","DNAm + SNP"))
colours = brewer.pal(12, "Paired")

figure_1a <- ggplot(fig_dat36, aes(x=ph1, y=r2_36, fill=ph)) + 
	geom_bar(stat="identity") +
	scale_fill_manual(values=colours[c(3,7,9)])+
	ylab("Incremental R2 (%)") +
	xlab("Omics data type") +
	ylim(0,15) + 
#	ggtitle("EpiScore and polygenic prediction of g in LBC1936") + 
	labs(title = bquote(bold('A:') ~ "EpiScore and polygenic prediction of g in LBC1936")) +
	theme_light() +
	theme(axis.title.x = element_text(size = 16),
	axis.title.y = element_text(size = 16), 
        axis.text=element_text(size=14),
	plot.title=element_text(size = 16),
	legend.position = "none") 

# 1B:
fig1b = readRDS("Cog_EWAS_BayesR/Figure_1B_plotdata.rds")
fig1b = do.call("rbind", fig1b)
fig1b$Type =  gsub("\\..*", "", rownames(fig1b))
fig1b$Type = gsub("g_pred", "EpiScore", fig1b$Type)
fig1b$Type = gsub("measured", "Measured g", fig1b$Type)

fig1b$Trait = gsub("yrsedu_w1", "Education", fig1b$Trait)
fig1b$Trait = gsub("bmi_w1", "BMI", fig1b$Trait)
fig1b$Trait = gsub("sixmwk_w1", "6m Walk", fig1b$Trait)
fig1b$Trait = gsub("lung", "Lung function", fig1b$Trait)
fig1b$Trait = gsub("hibp_w1", "Hypertension", fig1b$Trait)
fig1b$Trait = gsub("hadsd_w1", "Depression", fig1b$Trait)
fig1b$Trait = gsub("grip", "Grip strength", fig1b$Trait)
fig1b$Trait = gsub("ever_smoke", "Ever smoked", fig1b$Trait)
fig1b$Trait = gsub("depind_w1", "Deprivation", fig1b$Trait)
fig1b$Trait = gsub("bmi_w1", "BMI", fig1b$Trait)
fig1b$Trait = gsub("alcunitwk_w1", "Alcohol", fig1b$Trait)
fig1b$Trait = factor(fig1b$Trait, levels=rev(c("Lung function", "Hypertension", "Grip strength", "Ever smoked", "Education", "Deprivation", "Depression", "BMI", "Alcohol", "6m Walk")))

fig1b[which(fig1b$Trait %in% c("Deprivation", "6m Walk")), "Effect"] = -fig1b[which(fig1b$Trait %in% c("Deprivation", "6m Walk")), "Effect"] # reverse code these variables for interpretability
fig1b$LCI = fig1b$Effect - 1.96*fig1b$SE
fig1b$UCI = fig1b$Effect + 1.96*fig1b$SE

figure_1b = ggplot(fig1b, aes(x=Effect, y=Trait, group=Type)) + 
geom_point(aes(colour=Type), position = position_dodge(width=0.75)) + 
geom_errorbar(aes(xmin=LCI, xmax=UCI, colour=Type), position = position_dodge(width=0.75)) + 
xlab("Effect Size (Beta)") + 
ylab("") +
		theme_light() + 

theme(axis.title.x = element_text(size = 16),
     	axis.title.y = element_text(size = 16), 
        axis.text=element_text(size=14),
	    plot.title=element_text(size = 16),
        legend.text=element_text(size=16),
	    legend.position=c(0.85,0.3),
	    legend.background = element_rect(size=0.5, linetype="solid", colour="black"),
	    legend.title=element_blank()) + 
labs(title = bquote(bold('B:') ~ "EpiScore and measured g associations with ageing traits in LBC1936"))
	  



#######################################
### Plot Inflam/Neuro (1C/1D)       ###
#######################################
library(ggrepel)

# Read in 70 proteins 
proteins_inf <- read.csv("../Danni/LBC_proteins/protein_data/Normalised_Inflammatory_Proteins_W1.csv")

# Remove plate ID (not needed as feature) and merge

d36 <- merge(w1_36, proteins_inf[,-72], by.x="ID", by.y="lbc36no", all.x=T)

EpiScore_beta <- NA
EpiScore_p <- NA
g_beta <- NA
g_p <- NA

for(i in 26:95){
EpiScore_beta[i-25] <- summary(lm(scale(d36[,i]) ~ scale(age) + sex + scale(g_pred), data=d36))$coefficients[4,1]
EpiScore_p[i-25] <- summary(lm(scale(d36[,i]) ~ scale(age) + sex + scale(g_pred), data=d36))$coefficients[4,4]
g_beta[i-25] <- summary(lm(scale(d36[,i]) ~ scale(age) + sex + scale(g), data=d36))$coefficients[4,1]
g_p[i-25] <- summary(lm(scale(d36[,i]) ~ scale(age) + sex + scale(g), data=d36))$coefficients[4,4]
}

vars <- names(d36)[26:95]

out <- data.frame(cbind(EpiScore_beta, EpiScore_p, g_beta, g_p))
out$protein <- vars

### protein infl plot object ###
out$Name = NA
out$Name[which(out$g_p<(0.05/70) & out$EpiScore_p<(0.05/70))] = out$protein[which(out$g_p<(0.05/70) & out$EpiScore_p<(0.05/70))]
figure_1c <- ggplot(out, aes(EpiScore_beta, g_beta), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_light() +
      xlab("EpiScore g - Protein Effect Size") +
      ylab("Mesaured g - Protein Effect Size") +
	  labs(title = bquote(bold('C:') ~ "Inflammation Proteins with EpiScore and Measured g in LBC1936")) +
      theme(axis.title.x = element_text(size = 16),
	axis.title.y = element_text(size = 16), 
        axis.text=element_text(size=14),
	plot.title=element_text(size = 16),
	legend.position = "none") +
  geom_abline(slope=1, intercept=0, linetype = 2) +
  # geom_smooth(method="lm", colour="darkgrey") + 
  geom_point(data=out[out$EpiScore_p<(0.05/70),], shape=21, fill="red", size=2)+
  geom_point(data=out[out$g_p<(0.05/70),], colour='blue') +
  geom_point(data=out[out$g_p<(0.05/70) & out$EpiScore_p<(0.05/70),], colour='green') + 
  geom_point(data=out[out$g_p>(0.05/70) & out$EpiScore_p>(0.05/70),], colour='white') + 
  geom_label_repel(aes(label = Name),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')



### protein neuro plot object ###

d36_w2a <- readRDS("../Riccardo/Cog_EWAS_BayesR_GS_PGRS_to_LBC/d36_w2a.rds")

# Read in proteins # 751 people
proteins_neu <- read.csv("../Danni/LBC_proteins/protein_data/W2_Normalised_Neurology_Proteins.csv")
names(proteins_neu)[1] <- "lbc36no"

d36_w2b <- merge(d36_w2a, proteins_neu, by.x="ID", by.y="lbc36no", all.x=T)

EpiScore_beta <- NA
EpiScore_p <- NA

for(i in 20:111){
EpiScore_beta[i-19] <- summary(lm(scale(d36_w2b[,i]) ~ scale(age) + sex + scale(g_pred), data=d36_w2b))$coefficients[4,1]
EpiScore_p[i-19] <- summary(lm(scale(d36_w2b[,i]) ~ scale(age) + sex + scale(g_pred), data=d36_w2b))$coefficients[4,4]
}


vars <- names(d36_w2b)[20:111]

out <- data.frame(cbind(EpiScore_beta, EpiScore_p))
out$Protein <- vars

sarah <- read.csv("../Riccardo/Cog_EWAS_BayesR_GS_PGRS_to_LBC/Harris_et_al_protein_cog_assocs.csv")
sarah$Protein <- gsub("-", ".", sarah$Protein)
sarah$Protein <- gsub(" ", ".", sarah$Protein)
out_n <- merge(out, sarah, by="Protein")
out_n$Name = NA
out_n$Name[which(out_n$P_LBC1936<(0.05/90) & out_n$EpiScore_p<(0.05/90))] = out_n$Protein[which(out_n$P_LBC1936<(0.05/90) & out_n$EpiScore_p<(0.05/90))]

figure_1d <- ggplot(out_n, aes(EpiScore_beta, Beta_LBC1936), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_light() +
      xlab("EpiScore g - Protein Effect Size") +
      ylab("Mesaured g - Protein Effect Size") +
	  labs(title = bquote(bold('D:') ~ "Neurology Proteins with EpiScore and Measured g in LBC1936")) +
      theme(axis.title.x = element_text(size = 16),
	axis.title.y = element_text(size = 16), 
        axis.text=element_text(size=14),
	plot.title=element_text(size = 16),
	legend.position = "none") +
  geom_abline(slope=1, intercept=0, linetype = 2) +
  # geom_smooth(method="lm", colour="darkgrey") + 
  geom_point(data=out_n[out_n$EpiScore_p<(0.05/90),], shape=21, fill="red", size=2)+
  geom_point(data=out_n[out_n$P_LBC1936<(0.05/90),], colour='blue') +
  geom_point(data=out_n[out_n$P_LBC1936<(0.05/90) & out_n$EpiScore_p<(0.05/90),], colour='green') + 
  geom_point(data=out_n[out_n$P_LBC1936>(0.05/90) & out_n$EpiScore_p>(0.05/90),], colour='white') +
  geom_point(data=out_n[out_n$g_p>(0.05/90) & out$EpiScore_p>(0.05/90),], colour='black') +
  geom_label_repel(aes(label = Name),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')




fig1 = grid.arrange(figure_1a, figure_1b, figure_1c, figure_1d, nrow=2)

png("Cog_EWAS_BayesR/Figure_1_PostReview.png", width=1200, height=600)
grid.arrange(figure_1a, figure_1b, figure_1c, figure_1d, nrow=2)
dev.off()


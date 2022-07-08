setwd('enter working directory')

# Packages required
library(readxl)
library(tidyverse)
library(limma)
library(superheat)
library(reshape2)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(ggrepel)
library(sva)
library('variancePartition')
library('edgeR')
library('BiocParallel')
library(lme4)
library(lmerTest)
library(GGally)
library(mice)
library(miceadds)
library(pheatmap)
library(cluster)
library(purrr)
library(NbClust)
library(tscR)
library(multiMiR)
library(clusterProfiler)
library("org.Hs.eg.db")
library(enrichplot)
library(TissueEnrich)
library("miRBaseConverter")

#Read in miR Data
dat1= read_excel("Raw data_final.xlsx", 
                sheet = 2, 
                na = c("","NA", "?", "#VALUE!", "Undetermined"),
                trim_ws = TRUE)
# Read in Phenotype Table
pheno=read_excel("Raw data_final.xlsx",
                 sheet = 1,
                 na = c("","NA", "?"),
                 trim_ws = TRUE)

# RELEVEL Factors
pheno$Timepoint <- factor(pheno$Timepoint,
                          levels = c("T0","T1", "T2"))
pheno$ID <- as.factor(pheno$ID)

#Conduct Imputation (requires time to impute)
pred=quickpred(M_wide, mincor=0.4)
imputed <- mice(data = M_wide, pred=pred, m = 3, print= F, seed = 3) 
impute_long=mice::complete(imputed,"long")
write.csv(impute_long, "imputed_long.csv")

# Read in imputed dataframe from above (imputed once as it takes hours and stored it as CSV)
dat1= read_excel("Raw data_final.xlsx", 
                 sheet = 4, 
                 na = c("","NA", "?", "#VALUE!", "Undetermined"),
                 trim_ws = TRUE)

# Select Imputation 3 for all further analysis analysis 
M_imp=  dat1 %>%
  filter(IMP_NUM == "3") %>%
  dplyr:: select(2:476)

#Convert to wide 
M_long= pivot_longer(data =M_imp,
                     cols = -c("ID"),
                     names_to = c(".value","Timepoint"), 
                     names_pattern = "(.+)_(T0|T1|T2)$",
                     values_drop_na = 'FALSE')


# Data processing for Heatmap
dat1_T0 = M_long %>%
  filter(Timepoint == "T0")%>%
  summarise(across(where(is.numeric), ~psych::geometric.mean(., na.rm=TRUE)))
write.csv(dat1_T0, "Timepoint 0.csv")

dat1_T1 = M_long %>%
  filter(Timepoint == "T1")%>%
  summarise(across(where(is.numeric), ~psych::geometric.mean(., na.rm=TRUE)))
write.csv(dat1_T1, "Timepoint 1.csv")

dat1_T2 = M_long %>%
  filter(Timepoint == "T2")%>%
summarise(across(where(is.numeric), ~psych::geometric.mean(., na.rm=TRUE)))
write.csv(dat1_T2, "Timepoint 2.csv")

# This data has been averaged for each microRNA for each timepoint
# Sheet 1 = mean
# Sheet 2 = Geomean
Zdat= read_excel("Phase.xlsx", 
                 sheet = 2, 
                 na = c("","NA", "?", "#VALUE!", "Undetermined"),
                 trim_ws = TRUE)

#Relabel Time
levels(Zdat$Time) <- list(EF = "0",
                          O = "1",
                          ML = "2")

data.heatmap <- column_to_rownames(Zdat, 'Time') #pheatmap needs  row and column names
rownames(data.heatmap)[1]="EF"
rownames(data.heatmap)[2]="O"
rownames(data.heatmap)[3]="ML"
data.heatmap <- as.matrix(data.heatmap) # make a matrix
df3=t(scale(data.heatmap))
row.names= c("OP", "FP", "LP")
zscore.heatmap = pheatmap(mat= df3,
                         legend= TRUE,
                         cluster_cols=T,
                         angle_col = "315",
                         show_rownames = F,
                         filename = "Heatmap.tiff")
dev.off()

# Create Matrix for analysis
M_long2= M_long %>%
  unite("ID", ID:Timepoint, remove= T)
M_impute=M_long2 %>% remove_rownames() %>%
  column_to_rownames(var = "ID") #Make miR names become row names

# MIXED MODELLING with DREAM PACKAGE
# The rowname of phenoTable must be the same as column name of Matrix
M_impute= t(M_impute)
P2=pheno %>% remove_rownames() %>%
  column_to_rownames(var = "SampleName") #Make miR names become row names 

# Q1- Do c-miRNAs change throughout the menstrual cycle?
form<- ~  Time_NUM + Age + (1|ID)
fitmm = dream(M_impute, formula= form, data = P2)
fitmm = eBayes(fitmm)

results=topTable(fitmm, 
                 coef='Time_NUM',
                 adjust.method = "BH",
                 number=nrow(M_impute))

#violin plot of contribution of each variable to total variance
vp = fitExtractVarPartModel(M_impute, form, P2)
fig=plotVarPart(vp)
vp <- sortCols(vp)

# CONTRASTS TO COMPARE COEFFICIENTS
# Does expression of miRNAs change at different phases of the cycle?
form <- ~ 0 + Timepoint + (1|ID)
L = makeContrastsDream(form, P2, 
                       contrasts = c("TimepointT1 - TimepointT2",
                                     "TimepointT0 - TimepointT1",
                                     "TimepointT0 - TimepointT2"))
# fit dream model with contrasts
fit = dream(M_impute, form, P2, L)
fit = eBayes(fit)
results=topTable(fit,
                 coef = "TimepointT1 - TimepointT2",
                 number = nrow(M_impute),
                 adjust.method = "BH")

# Do any of the c-miRs cluster similarily?
M_long = pivot_longer(data = Zdat,
                      cols= -c("Time"))
M_wide = pivot_wider(data = M_long, 
                     names_from = Time,
                     values_from = value)
M_wide2=M_wide %>% remove_rownames() %>%
  column_to_rownames(var = "name") #Make miR names become row names 
scaleM = scale(M_wide2)

#How many clusters?
# function to compute total within-cluster sum of squares
fviz_nbclust(scaleM, kmeans, method = "wss", k.max = 24) + 
  theme_minimal() + ggtitle("the Elbow Method")

gap_stat <- clusGap(scaleM, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

### CLUSTERING ANALYSIS
time <- c(1,2,3)
#grouping the trajectories based on similarities in their slopes regardless of the distance between them
sDist2=slopeDist(scaleM, time)
sClust2= getClusters(sDist2, k=3)
plotCluster(scaleM, sClust2, "all")

endc=data.frame(miRNAs = rownames(scaleM), cluster = sClust2$clustering)
View(endc)

Cmerge= merge(endc, M_wide2, by= "row.names")
Cmerge_long= pivot_longer(data = Cmerge,
                          cols= -c("miRNAs", "cluster", "Row.names" ))
Cmerge_long$name <- factor(Cmerge_long$name,
                          levels = c("1","2", "3")) 

# Create Cluster Plot
cluster=ggplot(Cmerge_long, aes(x=name, y=value, color=cluster, group= miRNAs)) +
  geom_point(size=2)+
  geom_line(size=0.8)+
  facet_grid(~cluster)+
  scale_y_continuous(trans='log10')+
  scale_color_gradient(low = "blue", high = "red")+
  labs(y=paste("cf-miRNA expression (AU)"))+
  theme(text = element_text(size=16), 
        legend.position = "none",
        axis.title.x = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text= element_blank(), 
        axis.text = element_text(size = 14),
        axis.line = element_line(colour = "black", size=1))
ggplot2::ggsave(file = "c-miRNA_cluster.tiff", plot = cluster, 
                width = 10, height = 10)

Cmerge_long= read.csv2(file = "Cluster Analysis.csv", sep = ",")

Cmerge_long2= Cmerge_long %>%
  mutate(miRNAs = str_replace_all(miRNAs, "_", "-"))

# CLUSTER ANALYSIS
C1=  Cmerge_long2 %>%
  filter(cluster == "1") 

C2=  Cmerge_long2 %>%
  filter(cluster == "2") 

C3=  Cmerge_long2 %>%
  filter(cluster == "3") 

#Are clusters changing across the cycle?
Cmerge_long$name <- as.numeric(Cmerge_long$name)
Model1 = lmer(value ~ name + (1|miRNAs), data = C2)
summary(Model1)

#Imported saved results of DREAM analysis
dat1= read_excel("Results_DREAM_Imputed.xlsx", 
                 sheet = 3, 
                 na = c("","NA", "?", "#VALUE!", "Undetermined"),
                 trim_ws = TRUE)

# multiMiR:Prediction tool
# Replace the underscores (_) in data with dashes (-)
TargetmiR = (dat1$`MIR-contrasts`)
View(TargetmiR)
C1_GENE= get_multimir(org = "hsa",
                      mirna = TargetmiR,
                      table = "validated",
                      summary = T)

validated=C1_GENE@data
head(C1_GENE@data)
# Which interactions that are supported by Luciferase assay or western blotting
Luciferase=C1_GENE@data[grep("Luciferase", C1_GENE@data[, "experiment"]), ]
WB=C1_GENE@data[grep("Western", C1_GENE@data[, "experiment"]), ]
Val= rbind(Luciferase, WB)

### TISSUE ENRICHMENT
genes_ensembl= c(Val$target_ensembl)
genes_2=unique(genes_ensembl) 
#Background List
genes_BG= c(Val$target_ensembl)
genes_BG2=unique(bg)

gs<-GeneSet(geneIds=genes_2,
            organism="Homo Sapiens",
            geneIdType=ENSEMBLIdentifier())

output<-teEnrichment(inputGenes = gs,
                     tissueSpecificGeneType = 1)

seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)

tENRICH= ggplot(enrichmentOutput,aes(x=reorder(Tissue,-fold.change),fold.change,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  geom_hline(yintercept = 1.301)+
  theme_bw()+
  theme(legend.position="none")+
  theme(axis.text = element_text(size = 10))        +
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

ggpubr::ggarrange(tENRICH)+
  tiff('EnrichTissue_miRNAs contrasts.tiff', width = 10, height = 7.5, units = 'in', res=600)
dev.off() 

# Pathway Enrichment
#Create Background Genes
Val= read.csv2("miRNAs BACKGROUND unique.csv", sep = ",")
gene_BG <- Val$BG_Genes
GO <- enrichGO(gene         = genes_2,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont           = "all",
               pAdjustMethod = "BH",
               universe = gene_BG, #Background Genes
               pvalueCutoff= 0.05,
               qvalueCutoff  = 0.1)
edox2 <- pairwise_termsim(GO)
treeplot(edox2, 
         #label_format = 3, 
         nWords = 5, #Number of words for Cluster Tags 
         #offset= rel(10), #Moves the dendogram clustering
         #offset_tiplab =rel(1), # Moves GO pathways
         #label_format_cladelab=2,
         extend=0.5,
         hexpand =0.5, #Moves position of cluster tags
         group_color = c("#FF0000", "#000099", "#336600", "#990099", "#990033"),
         fontsize=5)+ 
  ggplot2::theme(legend.position=c(0.9,#across  
                                   0.5), #Up/down 
                   text = element_text(size = 15))
                   
ggpubr::ggarrange(TREE)+
  tiff('EnrichTREE.tiff', width = 20, height = 10, units = 'in', res=600)
dev.off()

#Q2: Are miRNAs altered by differing levels of the hormones?
glimpse(pheno)
collinearity=ggpairs(pheno, columns =  5:10) # Can I run hormones in same model? No collinearity between hormones
ggplot2::ggsave(file = "c-miRNA_collinearity.tiff", plot = collinearity, 
                width = 10, height = 10)
dev.off()

# Check Distribution of hormones
hist(pheno$FSH)
pheno = pheno%>%
  mutate(E2log=log(E2)) %>%
  mutate(Plog= log(Progesterone))%>%
  mutate(LHlog= log(LH))%>%
  mutate(FSHlog= log(FSH))

# The rowname of phenoTable must be the same as column name of Matrix
P2=pheno %>% remove_rownames() %>%
  column_to_rownames(var = "SampleName") #Make miR names become row names 

M_impute= t(M_impute)
form<- ~  FSHlog + Time_NUM + Age + (1|ID)
fitmm = dream(M_impute, formula= form, data = P2)
fitmm = eBayes(fitmm)
summary(decideTests(fitmm))
results=topTable(fitmm, 
                 coef='FSHlog',
                 adjust.method = "BH",
                 number=nrow(M_impute))

# CONTRASTS 
# Does expression of miRNAs change at different phases of the cycle?
form <- ~ 0 + E2log * Timepoint + (1|ID)
L = makeContrastsDream(form, P2, 
                       contrasts = c("TimepointT1 - TimepointT2",
                                     "TimepointT0 - TimepointT1",
                                     "TimepointT0 - TimepointT2"))
# fit dream model with contrasts
fit = dream(M_impute, form, P2, L)
fit = eBayes(fit)
results=topTable(fit,
                 coef = "E2log:TimepointT1 - E2log:TimepointT2",
                 number = nrow(M_impute),
                 adjust.method = "BH")

# GRAPH RESULTS
miR <- rownames(results)
miR <- which(rownames(M_impute) == "hsa_miR_19a_3p")
which(rownames(M_impute) == "hsa_miR_92a_3p") 
pheno_with_miR <- cbind(pheno,
                        M_impute = as.numeric(M_impute[1,]))
# Box plot of timepoints
ggplot(pheno_with_miR, aes(x=Timepoint, y=M_impute, color=ID, group=ID)) +
  geom_point()+
  geom_line()+
# geom_boxplot(aes(group=Timepoint))+
  labs(y=paste("Expression level", miR ))+
  theme(text = element_text(size=12), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text= element_blank(), 
        axis.line = element_line(colour = "black", size=1))
dev.off()
# correlation plot
dev.off()
ggplot(pheno_with_miR, aes(x=LHlog, y=M_impute )) +
  geom_point()+
  geom_smooth(method = "rlm", se=T, 
              formula = y ~ x)+
   labs(y=paste("Expression level", miR ))+
  theme(text = element_text(size=12), 
        panel.background = element_blank(), 
        strip.background = element_blank(),
        strip.text= element_blank(), 
        axis.line = element_line(colour = "black", size=1))

#Table One
colnames(pheno)
pheno_subset= dplyr::select(pheno, 
                            c("ID", "Timepoint","Age", "E2", "Progesterone","LH","FSH","LH.FSH_Ratio"))
colnames(pheno_subset)
aov(E2 ~ Timepoint, data = pheno_subset)
Model1 = lmer(E2 ~ 0+ Timepoint + (1|ID), data = pheno_subset)
summary(Model1)$coefficients
em<- emmeans(Model1,  ~Timepoint) #Compare levels of Timepoint within Group
summary(em)

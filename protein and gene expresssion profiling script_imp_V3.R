######mRNA transcriptOMICS profiling###################
#########--DIVYA AGRAWAL--####################
#######################################################
setwd("E:/data")  ##modified on 2/9/2025
getwd()
getProjectSummary('TCGA-BRCA')
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("ggpubr")
BiocManager::install("TCGAbiolinks")
BiocManager::install("TCGAbiolinksGUI.data")
BiocManager::install("glmnet")
BiocManager::install("gProfileR")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("magrittr")
BiocManager::install("ggplotify")
BiocManager::install("GenomicDataCommons")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("tidyverse")
BiocManager::install("SummarizedExperiment")
BiocManager::install("dplyr")
BiocManager::install("data.table")
BiocManager::install("ggplot2")
BiocManager::install("RColorBrewer")
BiocManager::install("genefilter")
BiocManager::install("clusterProfiler")
BiocManager::install("BiocStyle")
BiocManager::install("hpar")
BiocManager::install("HPAanalyze")
BiocManager::install("FactoMineR")
BiocManager::install("factoextra")
devtools::install_github("trannhatanh89/HPAanalyze")
devtools::install_github("anhtr/HPAanalyze")
install.packages("httr")
BiocManager::install("caret")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("tidyselect")
BiocManager::install("plotly")
BiocManager::install("UniprotR")
BiocManager::install("openxlsx")
BiocManager::install("TCGAbiolinks")
BiocManager::install("ggVennDiagram")
BiocManager::install("dplyr")
BiocManager::install("tidyverse")
BiocManager::install("limma")
BiocManager::install("preprocessCore")
BiocManager::install("venn")
BiocManager::install("ggVennDiagram")
BiocManager::install("pathview")
BiocManager::install("corrplot")
BiocManager::install("rlang")
BiocManager::install("EnhancedVolcano")
BiocManager::install("ggrepel")


library(corrplot)
library(venn)
library("EnhancedVolcano")
library(rlang)
library(ggVennDiagram)
library(pathview)
library(preprocessCore)
library(TCGAbiolinks)
library(ggVennDiagram)
library(dplyr)
library(tidyverse)
library(limma)
library(HPAanalyze)
library(BiocStyle)
library(gProfileR)
library(ggpubr)
library(glmnet)
library(ggplotify)
library(hpar)
library(plotly)
library(org.Hs.eg.db)
library(openxlsx)
library(TCGAbiolinks)
library(TCGAbiolinksGUI.data)
library(GenomicDataCommons)
library(GenomeInfoDbData)
library(magrittr)
library(UniprotR)
library(limma)
library(edgeR)
library(ggVennDiagram)
library(dplyr)
library(data.table)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(glmnet)
library(factoextra)
library(FactoMineR)
library(caret)
library(RColorBrewer)
library(gProfileR)
library(genefilter)
library(clusterProfiler)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(tidyselect)

q_mRNA<- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts"
)
m_mRNA <- q_mRNA[[1]][[1]]
m_mRNA$sample_type   ##1231
table(m_mRNA$sample_type) #Metastatic 7 #Primary Tumor 1111 #Solid tissue normal 113

## query for Primary Tumor patient###
q_PT_T<- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts",
                    sample.type = c("Primary Tumor"),
                  access = 'open'
)
meta_PT_T <- q_PT_T[[1]][[1]]
meta_PT_T[duplicated(meta_PT_T$cases.submitter_id), ]  ##16 replicates

## query for Solid Tissue Normal patient###
q_STN_T<- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  sample.type = c("Solid Tissue Normal"),
                  access = 'open'
)
meta_STN_T <- q_STN_T[[1]][[1]]
meta_STN_T[duplicated(meta_STN_T$cases.submitter_id), ]  ##0 replicates

##intersect both the above metadata#
i_mRNA<-intersect(meta_PT_T$cases.submitter_id, meta_STN_T$cases.submitter_id) #113+113
i_mRNA <-as.data.frame(i_mRNA)
colnames(i_mRNA)<-"cases.submitter_id"  ##name the column

##venn diagram representation##
seta<- meta_PT_T$cases.submitter_id
setb<-meta_STN_T$cases.submitter_id
x = list(seta, setb)
####metastatic tumor patients
q_MT<- GDCquery(project = "TCGA-BRCA",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                experimental.strategy = "RNA-Seq",
                workflow.type = "STAR - Counts",
                sample.type = c("Metastatic"),
                access ='open'
)
meta_MT <- q_MT[[1]][[1]]
meta_MT$sample_type   ##1231
table(meta_MT$sample_type) #Metastatic 7 
setc<- meta_MT$cases.submitter_id
ggVennDiagram(x = list(seta, setb, setc),label_alpha = 1, category.names = c("PT", "STN", "MT"), label_size = 5, set_color = c("red","lightgreen","lightblue"))  ##Not working

##ggVennDiagram(x = list(seta, setb),label_alpha = 1, category.names = c("PT", "STN"), set_color = c("black"))  ##Not working
venn = Venn(x)
plot_venn(process_data(venn, shape_id = "201"), set_color = c("red","lightgreen"), label_alpha = 5, label_size = 7, edge_size = 2, edge_lty = "solid")

##merge the dataframes
meta_PT_T<-meta_PT_T[,c(3,29)]
meta_STN_T<-meta_STN_T[,c(3,29)]
mer<- merge(x = meta_STN_T, y= meta_PT_T, by="cases.submitter_id", all = TRUE)
colnames(mer)[2]<-'STN'
colnames(mer)[3]<-'PT'
mer<-na.omit(mer)    #119+119

###overlapped patients IDs##BARCODES##
barcodes<-unique(c(mer[,2], mer[,3]))  #113 STN #119 PT
barcodes <- barcodes[c((114:232), (1:113))]   ##232 total samples

##########Actuall filtered#### pateint samples####
query_mRNA<- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts",
                    sample.type = c("Primary Tumor","Solid Tissue Normal"),
                    barcode = barcodes,
                    access = 'open')
                  
m_bar <- query_mRNA[[1]][[1]]

GDCdownload(query_mRNA, method = "api", files.per.chunk = 6)    ###232 files
data_mRNA <- GDCprepare(query_mRNA,
                      save = T,
                      save.filename = "./brca_pt_trans.rda")
load("./brca_pt_trans.rda")

######################################################
##########edge R Differential expression analysis########
#######################################################
PT_matrix <- assay(data_mRNA)
dgeBRCA <- DGEList(PT_matrix) ###don't follow the next code##always

##How many patients are having the zero counts##
df <-data.frame(PT_matrix)
colSums(df==0)
mean(colSums(df==0))    ##27477.53

###gene annotation
dgeBRCA <-DGEList(counts = PT_matrix, genes = mcols(data_mRNA)$gene_name)

###filtered to protein_coding genes in the summarized experiment object
hgnc_id<-rowData(data_mRNA)$hgnc_id
dgeBRCA$genes$hgnc_id <-hgnc_id
gene_type<- rowData(data_mRNA)$gene_type
dgeBRCA$genes$gene_type <- gene_type
dgeBRCA <- dgeBRCA[dgeBRCA$genes$gene_type == "protein_coding", ]

dim(dgeBRCA)    ##19962  232
group <- factor(data_mRNA$sample_type_id)
dgeBRCA$samples$group<-group
dgeBRCA$samples$group
rownames(dgeBRCA) <- lapply(rownames(dgeBRCA),  sub, pattern = "\\.\\d+$", replacement = "")

###to count the number of zeros in the single rows in summarized object##
table(rowSums(dgeBRCA$counts==0)==232)    #FALSE 19446 #TRUE 516

####density plot of the non-filtered protein_coding genes#####
plot(density(cpm((dgeBRCA),log=TRUE)) , main="Raw data", xlab="log CPM", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, font.axis = 2)

##one filtering method####USE THIS ##
cpmBRCA <-cpm(dgeBRCA)
thrBRCA <- cpmBRCA>0.5
table(rowSums(thrBRCA))
keepBRCA <- rowSums(thrBRCA) >= 10
summary(keepBRCA)
dgeBRCA <- dgeBRCA[keepBRCA,]
plot(density(cpm((dgeBRCA),log=TRUE)), main="Filtered data", xlab="log CPM", cex.main = 1.5, cex.lab =1.5, cex.axis = 1.5, font.axis = 2)

######################################
###before normalization ##MA plot ### boxplot####
nsamples <- ncol(dgeBRCA)
col <- colorRampPalette(c("blue","red"))(232)
lcpm <-cpm(dgeBRCA, log = TRUE)
boxplot(lcpm, las=2, col=col, main="", cex.axis = 0.3)
title(main="Raw", ylab="Log-cpm", cex.main = 1)

###MA plot before normalisation 
average <- rowMeans(lcpm)
log2FC <- rowMeans(lcpm[, c(1:119)]) - rowMeans(lcpm[, c(120:232)])
correlation<-cor(average, log2FC, method = "pearson")
print(correlation)    ###0.07792364

# create MA plot
plot(average, log2FC, pch = 16, col = "blue", main = "Raw",
     xlab = "Average Expression(A)", ylab = "Log2 Fold Change(M)")
abline(h = 0, lty = 2, col="yellow")

##normalization method###TMM method
dgeBRCA <-calcNormFactors(dgeBRCA, method = "TMM")   ##TMM method###
dgeBRCA$samples$norm.factors
lcpm2<-cpm(dgeBRCA, log = TRUE)
boxplot(lcpm2, las=2, col=col, main="", cex.axis = 0.3)
title(main="Normalized", ylab="Log-cpm", cex.main = 1)

###MA plot after normalization###
average <- rowMeans(lcpm2)
log2FC <- rowMeans(lcpm2[, c(1:119)]) - rowMeans(lcpm2[, c(120:232)])
correlation<-cor(average, log2FC, method = "pearson")
print(correlation)   ###0.07746663

# create MA plot
plot(average, log2FC, pch = 16, col = "black", main = "Normalized",
     xlab = "Average Expression(A)", ylab = "Log2 Fold Change(M)", cex = 0.01, cex.main = 1.5, cex.lab =1.5, cex.axis = 1.5, font.axis = 2)
abline(h = 0, lty = 10, col="blue", lwd = 2)

#####plotsmear##
plotSmear(dgeBRCA, de.tags=rownames(dgeBRCA), main="Normalized", cex=0.5)
abline(h=0, col="black", lty=2, lwd=2) ###plotsmear is based on the limma model

##MD plot
plotMD(cpm(dgeBRCA, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

###MDS plot
points <- c(0,1)
colors <- rep(c("black", "red"), 2)
plotMDS(lcpm2, col=colors[group], pch=points[group])
legend("bottomleft", legend=levels(group), pch=points, col=colors, ncol=1)

##estimating the dispersion###
##############################
group <- factor(data_mRNA$sample_type_id)
design <-model.matrix(~0+group)
dgeBRCA <-estimateDisp(dgeBRCA, design = design, robust = TRUE)
dim(dgeBRCA)     ##16404   ##232
plotBCV(dgeBRCA) ##Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
fit <- glmQLFit(dgeBRCA,design, robust = TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)
Gexp<- makeContrasts(group01-group11, levels=design)

##glmQLFtest### ###use this object for volcano plot only##
res <- glmQLFTest(fit, contrast=Gexp)
all_genes <- rownames(res$coefficients)
top_tags <- topTags(res, nrow(res$table), adjust.method = "BH")

is.de <- decideTests(res)
summary(is.de)  #down 4711 #3718 #7975  total 16404  ###down 4791 #not sig 3573 ##up 8040
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright", cex = 0.1)

###downloaded as a result
top_tags<-as.data.frame(top_tags)
top_tags$expression<-ifelse(top_tags$logFC >0.1 & top_tags$FDR <=0.05, "Higher abundance",
                            ifelse(top_tags$logFC < -0.1 & top_tags$FDR <=0.05, "Lower abundance", "Non significant"))
table(sum(top_tags$expression == "Higher abundance"))  ##7963 #4708 Down #3733 non sig now 8023 up
top_tags <- top_tags[order(top_tags$FDR),]

# Create the plot
volcano_plot <- ggplot(top_tags, aes(x=logFC, y=-log10(FDR), color=expression)) + 
  geom_point(size=1) +
  scale_color_manual(values=c("green", "red", "blue")) +
  theme_classic() +
  xlim(c(-14, 14)) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color="grey50", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", alpha=0.5) +
  labs(x="Log2 Fold Change", y="-Log10(FDR)", color="Expression") +
  theme(legend.position="bottom")

# Add the top 20 genes
volcano_plot + geom_label_repel(data=top_tags[1:20,], aes(x=logFC, y=-log10(FDR), label=top_tags$genes[1:20]), 
                         size=3, color = "black", box.padding = 0.5, force =1, max.overlaps = Inf) +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))


############################################################
###################PROTEOMICS PROFILING#####################
############################################################

query<- GDCquery(project = "TCGA-BRCA",
                 data.category = "Proteome Profiling",
                 data.type = "Protein Expression Quantification",
                 access = 'open')
# Get metadata matrix
meta <- query[[1]][[1]]

meta$sample_type ##919
table(meta$sample_type)  #metastatic 5 # Primary Tumor 881 #Solid Tissue Normal 33

##Intersecting the common donor ids between primary and normal sample patients##
#extracting the primary tumor samples patients IDs
q_PT<- GDCquery(project = "TCGA-BRCA",
                data.category = "Proteome Profiling",
                data.type = "Protein Expression Quantification",
                sample.type = c("Primary Tumor"),
                access = 'open')

##get matadata matrix
m_PT<-q_PT[[1]][[1]]

table(m_PT$sample_type)  #Primary Tumor 881

#extracting the solid tissue normal samples patients IDs
q_STN<- GDCquery(project = "TCGA-BRCA",
                 data.category = "Proteome Profiling",
                 data.type = "Protein Expression Quantification",
                 sample.type = c("Solid Tissue Normal"),
                 access = 'open')

##get the metadata matrix
m_STN <-q_STN[[1]][[1]]

table(m_STN$sample_type)  #Solid Tissue normal 33
#extracting the metastatic samples patients IDs
q_MTp<- GDCquery(project = "TCGA-BRCA",
                 data.category = "Proteome Profiling",
                 data.type = "Protein Expression Quantification",
                 sample.type = c("Metastatic"),
                 access = 'open')

##get the metadata matrix
m_MTp <-q_MTp[[1]][[1]]

table(m_MTp$sample_type)  #Metastatic 5

##intersecting/common overlapped donor IDs (the above both sample_types)
i<-intersect(m_PT$cases.submitter_id, m_STN$cases.submitter_id)  ##33 from both
i <-as.data.frame(i)
colnames(i)<-"cases.submitter_id"  ##name the column

#now intersect the above i vector to the m_PT cases
intersect(i$cases.submitter_id, m_PT$cases.submitter_id)  ##confirmed intersected
##venn diagram representation##
set1<- m_PT$cases.submitter_id
set2<-m_STN$cases.submitter_id
set3<- m_MTp$cases.submitter_id

ggVennDiagram(x = list(set1, set2, set3),label_alpha = 1, label_size= 5,set_color = c("red","lightgreen","lightblue"))
###
x<-list(set1, set2)
venn_p = Venn(x)
plot_venn(process_data(venn_p, shape_id = "201"), set_color = c("red", "lightgreen"), label_alpha = 1, edge_size = 2, edge_lty = "solid", label_size = 7)

##merge the dataframes
m_STN<-m_STN[,-c(1:22)]
m_PT<-m_PT[,-c(1:22)]
m<- merge(x = m_STN, y= m_PT, by="cases.submitter_id", all = TRUE)
m<-na.omit(m)
colnames(m)[2]<-'STN'
colnames(m)[3]<-'PT'

###overlapped patients IDs##BARCODES##
my_barcodes<-unique(c(m[,2], m[,3]))
my_barcodes <- my_barcodes[c((34:66), (1:33))]
query.prot<-GDCquery(project = "TCGA-BRCA",
                     data.category = "Proteome Profiling",
                     sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                     data.type = "Protein Expression Quantification",
                     barcode = my_barcodes,
                     access = 'open')
                     
# Get metadata matrix
metadata1 <- query.prot[[1]][[1]]

# Download data using api
GDCdownload(query.prot, method = "api", files.per.chunk = 6)

data <- GDCprepare(query.prot,
                   save = T,
                   save.filename = "./brca_prot.rda",
                   summarizedExperiment = TRUE)
load("./brca_prot.rda")
class(data)

###
d<-data.frame(data)
d<-d[,-c(1:4)]
d <- column_to_rownames(d, var = "peptide_target")
sum(is.na(d))           #1797
d<-na.omit(d)
dim(d)                  #457 rows #66 col

##MA PLOT###############################raw data##
# calculate average expression and fold change
raw_d <-d[,c(1:66)]
raw_average <- rowMeans(raw_d)
raw_log2FC <- rowMeans(raw_d[, c(1:33)]) - rowMeans(raw_d[, c(34:66)])
correlation<-cor(raw_average, raw_log2FC, method = "pearson")
print(correlation)   ## -0.7765531
# create MA plot
plot(raw_average, raw_log2FC, pch = 16, col = "blue", main = "Raw data",
     xlab = "Average Expression(A)", ylab = "Log2 Fold Change(M)", cex= 0.5,
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, font.lab = 1, font.main = 2, font.axis = 2)
abline(h = 0, lty = 2, lwd = 2)
par(font.lab = 2)


##special visualizaiton of MA plot##
f <- data.frame(raw_average, raw_log2FC)
ggplot(f, aes(x = raw_average, y = raw_log2FC)) + 
  geom_point(color = "#0072B2", size = 1) +
  scale_x_continuous(limits = c(0, 0.9)) + 
  coord_cartesian(ylim = c(-4, 4)) +
  labs(title = "MA Plot", x = "Average Expression(A)", y = "Log2 Fold Change(M)") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

####quantile normalization##

# Define a set of candidate reference genes
ref_genes <- rownames(d)
ref_genes

# Load the necessary libraries
d<-as.matrix(d) 
head(d)

# Perform normalization using the normalize.quantiles function
norm_data <- normalize.quantiles(d)
head(norm_data)

# Combine candidate reference genes with the original data
norm_data <- cbind(norm_data, d[,-which(colnames(d) %in% ref_genes)])

# Check the normalization
boxplot(norm_data)

#MA plot ###
# calculate average expression and fold change
m <-norm_data[,c(1:66)]
average <- rowMeans(m)
log2FC <- rowMeans(m[, c(1:33)]) - rowMeans(m[, c(34:66)])
# create MA plot
plot(average, log2FC, pch = 16, col = "blue", main = "Normalized",
     xlab = "Average Expression(A)", ylab = "Log2 Fold Change(M)", cex = 0.5,      cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, font.lab = 1, font.main = 2, font.axis = 2)

abline(h = 0, lty = 2, lwd =2)
correlation<-cor(average, log2FC, method = "pearson")
print(correlation)   ##-0.5926424

##visually more appealing##
# create data frame
df <- data.frame(average, log2FC)

# create MA plot with ggplot2
ggplot(df, aes(x = average, y = log2FC)) + 
  geom_point(color = "#0072B2", size = 1) +
  scale_x_continuous(limits = c(0, 0.9)) +
  coord_cartesian(ylim = c(-4, 4)) +
  labs(title = "MA Plot", x = "Average Expression(A)", y = "Log2 Fold Change(M)") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


##################LIMMA MODEL#############
# make the group for limma analysis in the end##
colnames(d)[1:33] <- "PT"
colnames(d)[34:66] <- "STN"
group <-colnames(d)
group

# Construct a design matrix with two groups of samples
norm_data<-as.data.frame(norm_data)
design <- model.matrix(~ 0 + group, data = norm_data)

##voom.y.d<-voom(norm_data, design, plot = T)##cant be performed as negative counts are not allowed
# Fit a linear model to the data
fit <- lmFit(norm_data, design)
coef.fit<-fit$coefficients
head(coef(fit))

# Define the contrast matrix to compare Group1 and Group2
contrast_matrix <- makeContrasts(groupPT - groupSTN, levels = colnames(coef(fit)))
contrast_matrix
design

# Perform hypothesis testing using empirical Bayes moderation
results <- contrasts.fit(fit, contrast_matrix)
results<-eBayes(results)

#perform which proteins are significantly differential expressed proteins
res_final <- topTable(results, n = Inf, adjust = "BH", sort.by = "p")
head(res_final, 5)


###compared to coefficient##
coef_DEP<-coef.fit[rownames(coef.fit) %in% rownames(res_final)[1:5],]
coef_DEP
length(which(res_final$adj.P.Val <= 0.05))   ##267

library(dplyr)
DEPs<-res_final %>% arrange(logFC) %>% filter(adj.P.Val <0.05)
head(DEPs)
res_final$Proteins<-rownames(res_final)
res_final<-res_final[,c("Proteins",names(res_final)[1:6])]

# View the differentially expressed proteins
print(results)

###downloded as a result
res_final<-as.data.frame(res_final)

res_final$expression<-ifelse(res_final$logFC >0.1 & res_final$adj.P.Val <=0.05, "Higher abundance",
                             ifelse(res_final$logFC < -0.1 & res_final$adj.P.Val <=0.05, "Lower abundance", "Non significant"))
res_final <- res_final[order(res_final$adj.P.Val),]
head(res_final)

table(sum(res_final$expression == "Non significant"))  #190 #156 Down #111 up

# Create the volcano plot
volcano_plot <- ggplot(res_final, aes(x=logFC, y=-log10(adj.P.Val), color=expression)) + 
  geom_point(size=1) +
  scale_color_manual(values=c("green", "red", "blue")) +
  theme_classic() +
  xlim(c(-10, 10)) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color="grey50", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50", alpha=0.5) +
  labs(x="Log2Fold Change", y="-Log10(FDR)", color="Expression") +
  theme(legend.position="bottom")

# Add the top 20 genes
volcano_plot + geom_label_repel(data=res_final[1:20,], aes(x=logFC, y=-log10(adj.P.Val), label=res_final$Proteins[1:20]), 
                         size=3, color = "black", box.padding = 0.5, force = 1, max.overlaps = Inf) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 14))
###DEP list###ready#
res_final<-res_final[,-1]

#order the rownames#
res_final <- res_final[order(rownames(res_final)), ]

##row to column shifting##n(var = "peptide_targets")

res_final.new <- res_final %>% rownames_to_column
##MAPPING AND ADDITION OF THE GENE SYMBOL LIST###OVERLAPPS
map.d<-res_final.new$rowname    ##457 

#mapping
exp.list <- read.table("Expanded_Ab_List_Updated.txt", sep="\t", header=T)
pep.exp <- exp.list[,3]     ##499

mrna.list <-dgeBRCA$genes$genes   ##16404

pep.dat <- toupper(map.d)   ##converts all characters in string to uppercase map.d=457
pep.exp <- toupper(pep.exp)   ##pep.exp=499

length(intersect(pep.dat, pep.exp)) # 353 

length(pep.dat)                  # 457

x <- setdiff(pep.dat, pep.exp)   #104 ##these are not present in exp.list##also some are present but not the format that can be mapped to pep.dat

new.pep.exp <- gsub("-", "", pep.exp)  ##here removing all hyphens from the string of pep.exp
new.pep.dat <- gsub("-", "", pep.dat)

pep.exp <- new.pep.exp 
pep.dat <- new.pep.dat 

length(intersect(pep.dat, pep.exp)) # 397

pep.dat[grep("_P", pep.dat)] <- NA  ##grep searches for the _P for the pattern and returns the element as NA
pep.exp[grep("_P", pep.exp)] <- NA  ##grep searches for the _P for the pattern and returns the element as NA

length(intersect(pep.dat, pep.exp)) # 329

x <- setdiff(pep.dat, pep.exp)
length(x)                           # 37

rr1<- read_xlsx("rr1.xlsx") 
write.table(rr1, "yyy.txt", sep = "\t", row.names = FALSE)
rr <- read.table("yyy.txt", head=T, sep="\t")
new.pep.dat <- pep.dat
for(i in 1:nrow(rr) ){
  new.pep.dat[new.pep.dat %in% rr[i,1]] <- rr[i,2]
}

length(intersect(new.pep.dat, pep.exp))  #337

length(setdiff(new.pep.dat, pep.exp))     #0
setdiff(new.pep.dat, pep.exp)

new.pep.dat.sym <- new.pep.dat
for( i in 1:length(new.pep.dat)){
  if( !is.na(new.pep.dat[i])){
    new.pep.dat.sym[i] <- exp.list$Gene.Name[pep.exp %in% new.pep.dat[i]]
  }
  else{
    new.pep.dat.sym[i] <- NA
  }
}

length(new.pep.dat.sym) # 457    ##122 NA ##320 gene_symbols with mRNA are overlapped##
length(pep.dat)         # 457

#for venn diagram of mapping removing NA and duplicates if exist in proteins.
fra <- na.omit(new.pep.dat.sym)
sum(duplicated(fra))  # > 0 means multiple peptides mapped to same gene symbol
fra[duplicated(fra)]  ##2 #"CASP8"  "NOTCH1"
fr <- unique(fra) ##333
#to investigates any duplicates if exist in mRNA.
dups <- top_tags$genes[duplicated(top_tags$genes)]
print(dups)  ##POLR2J3

##venn diagram 
#set.11<-new.pep.dat.sym   
set.11<-fr
set.12<-top_tags$genes   ##16404
set_target<- list(set.11, set.12)
#ggVennDiagram(x = list(set.11, set.12),label_alpha = 0.8, category.names = c("Prot_gene.symbol", "mrna_gene.symbol") )
venn_p = Venn(set_target)
plot_venn(process_data(venn_p, shape_id = "201"), set_color = c("red", "cyan"), label_alpha = 1, edge_size = 2, edge_lty = "solid", label_size =7)
####

new.pep.dat.sym<-as.data.frame(new.pep.dat.sym)
colnames(new.pep.dat.sym)[1]<-'gene_symbol'
res_final.new <- cbind(res_final.new[, 1], new.pep.dat.sym$gene_symbol, res_final.new[, 2:ncol(res_final.new)])
colnames(res_final.new)[1]<-'peptide_target'
colnames(res_final.new)[2]<-'genes'

##Now merge with mRNA final result list#
#335
res_final.new1<-na.omit(res_final.new)  ##122 NA

# get the index of the column to exclude 
exclude_idx <- grep("genes", colnames(res_final.new1))

# add "_prot" to the remaining column names
new_colnames <- colnames(res_final.new1)[-exclude_idx]
new_colnames <- paste(new_colnames, "_prot", sep = "")

# replace the old column names with the new ones
colnames(res_final.new1)[-exclude_idx] <- new_colnames

###Total-overlapped genes in both expression analysis  ##322 17
res_final.new1<-merge( x = res_final.new1, y = top_tags, by = "genes")

###RESULT WITH 3 TABLES 
##CELL SURFACE TARGET TABLE##
##uni-prot##sub-cellular location##phospho proteins and cell surface targets#
final.list<-res_final.new1
##final.list<-read.csv("final.list.csv")
final.list$col3 <- ifelse(final.list$expression_prot == final.list$expression, final.list$expression_prot, paste(final.list$expression_prot, "vs", final.list$expression))
colnames(final.list)[18]<-'prot|mrna'
require(biomaRt)
ensembl <-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
description<- getBM(
  mart=ensembl,
  attributes=c("uniprotswissprot", "hgnc_id"),
  filter=c("hgnc_id"),
  values= unlist(final.list$hgnc_id),
  uniqueRows=FALSE)
description<-unique(description)
description$uniprotswissprot <- ifelse(description$uniprotswissprot == "", NA, description$uniprotswissprot)
description<-na.omit(description)
final.list.new<-merge(final.list, description, by= "hgnc_id")


##whole list##map/subcellular location##GO terms ids
library("UniprotR")
final.list.new1<-GetSubcellular_location(ProteinAccList = final.list.new$uniprotswissprot, directorypath = NULL)
go<-GetProteinGOInfo(ProteinAccList = final.list.new$uniprotswissprot, directorypath = NULL)
target.list<- cbind(final.list.new, final.list.new1, go)
target.list <- target.list %>% rownames_to_column(var = "uniprotIDs")
target.list <- target.list[!duplicated(target.list$hgnc_id),]  ##320 total##
table(target.list$`prot|mrna`)   #up 59

#removing the not required columns for downloading the result##
target.list <- subset(target.list, select = -c(uniprotIDs, uniprotswissprot, AveExpr_prot, t_prot, B_prot, P.Value_prot, gene_type, logCPM, F, PValue, Gene.Ontology.IDs))

##membrane proteins##55###ER system is also included##
ta <- target.list[!is.na(target.list$Topological.domain) & !is.na(target.list$Transmembrane), ]
ta<-ta[, -12]  ##remove the intramembrane column ##NA present so garbage##
table(ta$`prot|mrna`)

##up##upregulated proteins
up_ta <- subset(ta, `prot|mrna` == "up") ##11
Down_ta<-subset(ta, `prot|mrna` == "Down")
up_nonsig_ta <-subset(ta, `prot|mrna`== "up vs Non sig")
nonsig_up_ta <-subset(ta, `prot|mrna` == "Non sig vs up")


##downloading of the files_final mapping file###Target file##FOR SAVING FILES
write.table(target.list, file = "target.list.txt", row.names = F, sep = "\t",quote = T)
target.list<-read.table("target.list.txt")
target.list<-read_xlsx("F:/target.list.xlsx")

write.table(ta, file = "ta.txt", row.names = F, sep = "\t",quote = T)
ta<-read.table("ta.txt")
ta<-read_xlsx("F:/ta.xlsx")

write.table(up_ta, file = "up_ta.txt", row.names = F, sep = "\t",quote = T)
up_ta<-read.table("up_ta.txt")

write.table(up_nonsig_ta, file = "up_nonsig_ta.txt", row.names = F, sep = "\t",quote = T)
up_nonsig_ta<-read.table("up_nonsig_ta.txt")
 
write.table(nonsig_up_ta, file = "nonsig_up_ta.txt", row.names = F, sep = "\t",quote = T)
nonsig_up_ta<-read.table("nonsig_up_ta.txt")

##tissue/patho tissue expression results/--hpa###
hpa_data <- allHparData()
h<-hpaDownload(downloadList = 'histology')
h<- data("hpa_histology_data")
hpaListParam(h)
##Open tsv file
h<- read.table(file = "C:/Users/Divya Agrawal/Downloads/normal_ihc_tissues.tsv/normal_ihc_tissues.tsv", sep = '\t', header = TRUE)
#head(normtissue <- hpaNormalTissue())
library("ExperimentHub")
eh <- ExperimentHub()
query(eh, "hpar")
subcell <- hpaSubcellularLoc()
vt<-hpaSubset(h, targetGene = up_ta$V2,targetTissue = "breast", targetCancer = "breat cancer")

##These plots are not required##HPA
##
hpaVis(h,
       targetGene = c("DDR1"),
       targetTissue = c("Breast"),
       targetCancer = c("breast cancer"))    ###ploted this.
##to plot it##
hpaVisSubcell(                  #hpa
  data = NULL,
  targetGene = c("DDR1"),
  reliability = c("enhanced", "supported", "approved", "uncertain"),
  color = c("#ffffb2"),
  customTheme = FALSE
)

##other tissue expression ##hpa  ###I will use GEPIA database for extracting the expression
hpaVisTissue(
  data = h,
  targetGene = c("DDR1"),
  targetTissue = "breast",
  targetCellType = NULL,
  color = c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"),
  customTheme = FALSE
)

##other patho tissue expression##hpa
hpaVisPatho(
  data = h,
  targetGene = c("DDR1"),
  targetCancer = NULL,
  facetBy = "cancer",
  color = c("#FCFDBF", "#FE9F6D", "#DE4968", "#8C2981"),
  customTheme = FALSE
)

##biomart ### uniprot##
##uniprot##
BiocManager::install("UniprotR")
GetSubcellular_location(ProteinAccList = "Q08345",directorypath = NULL)

go<-GetProteinGOInfo(ProteinAccList = "Q08345",directorypath = NULL)

Plot.GOSubCellular(go, Top = 10, directorypath = NULL)
Plot.GOMolecular(go, Top = 20, directorypath = NULL)
PlotGOAll(go, Top = 50, directorypath = NULL)

GetExpression(ProteinAccList = "Q08345",directorypath = NULL)
Accs <- "Q08345"
GetPTM_Processing(ProteinAccList = "Q08345", directorypath = NULL)
GetSeqLength(Accs, directorypath = NULL)
GetpdbStructure(Accs, directorypath = NULL)
GetSubcellular_location(Accs, directorypath = NULL)
GetStructureInfo(Accs, directorypath = NULL)
##commercial antibody search##
###hpaxmlAntibody object code to extract the information about the proteins
##used for antibody formation##
GCH1xml <- hpaXmlGet('ENSG00000131979')
hpaXmlAntibody(GCH1xml)
##########################################################################
# Connect to the Ensembl database in BioMart
ensembl <- useMart("ensembl")

# Choose the dataset and attributes (e.g., human dataset and gene attributes)
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
attributes <- c("entrezgene_id","hgnc_symbol","external_gene_name")
listAttributes()
# Retrieve Entrez IDs for gene names using BioMart
genes_of_interest <- final.list.new$genes  # Replace with your gene names

# Get Entrez IDs for your genes of interest
gene_info <- getBM(attributes, filters = "hgnc_symbol", values = genes_of_interest, mart = ensembl_dataset)
# Sort the dataframe in decreasing order by column A
gene_info <- gene_info %>% arrange(gene_info$hgnc_symbol)
final.list.new <- final.list.new %>% arrange(final.list.new$genes)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = final.list.new$logFC

## feature 2: named vector
names(geneList) = as.character(gene_info$entrezgene_id)

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

data(geneList, package="DOSE")
gene <- names(geneList)

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.01)
gg<- enrichGO(gene         = gene,
             OrgDb=org.Hs.eg.db,
             ont = "BP",
              pvalueCutoff = 0.05)
dotplot(gg)
head(kk)
dotplot(kk)
barplot(kk)
browseKEGG(kk, 'hsa04110')
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


###HPA###
  
library(HPAanalyze)
h<-hpaDownload(downloadList = 'all')
hpaListParam(h)
vt<-hpaSubset(h, targetGene = "DDR1",targetTissue = "breast", targetCancer = "breat cancer")
hpaVis(visType = c("Patho", "Tissue", "Subcell"),
       targetGene = c("DDR1"),
       targetCancer = c("breast cancer"), customTheme = FALSE)

##cell surface target##
hpaVisTissue(
  data = h,
  targetGene = "DDR1",
  targetTissue = "breast",
  targetCellType = NULL,
  color = c("#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"),
  customTheme = FALSE
)
hpaVisPatho(
  data = h,
  targetGene = "DDR1",
  targetCancer = NULL,
  facetBy = "cancer",
  color = c("#FCFDBF", "#FE9F6D", "#DE4968", "#8C2981"),
  customTheme = FALSE
)

hpaVisSubcell(
  data = NULL,
  targetGene = "DDR1",
  reliability = c("enhanced", "supported", "approved", "uncertain"),
  color = c("#ffffb2", "#e31a1c"),
  customTheme = FALSE
)   ###hpaxmlAntibody object code to extract the information about the proteins
##used for antibody formation##
GCH1xml <- hpaXmlGet('ENSG00000131979')
hpaXmlAntibody(GCH1xml)


##biomart ### uniprot##
##uniprot##
BiocManager::install("UniprotR")
GetSubcellular_location(ProteinAccList = "Q08345",directorypath = NULL)

go<-GetProteinGOInfo(ProteinAccList = "Q08345",directorypath = NULL)

Plot.GOSubCellular(go, Top = 10, directorypath = NULL)
Plot.GOMolecular(go, Top = 20, directorypath = NULL)
PlotGOAll(go, Top = 20, directorypath = NULL)

GetExpression(ProteinAccList = "Q08345",directorypath = NULL)
Accs <- "Q08345"
GetPTM_Processing(ProteinAccList = "Q08345", directorypath = NULL)
GetSeqLength(Accs, directorypath = NULL)
GetpdbStructure(Accs, directorypath = NULL)
GetSubcellular_location(Accs, directorypath = NULL)
GetStructureInfo(Accs, directorypath = NULL)


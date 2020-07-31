# Load libraries
library(edgeR)
library(dplyr)
library(tibble)
library(openxlsx)

#### Data Pre-processing ####
# Import data
seqdata <- read.delim("./Data/exon_counts_new_clean.txt", stringsAsFactors = FALSE)

# Format the data 
countdata <- seqdata[,-(1:2)]
rownames(countdata) <- seqdata[,1]
colnames(countdata) <- c("Red._","Red.P2","Red.P4","RedGreen._","RedGreen.P2","RedGreen.P4")
countdata <- countdata[,c(6,4,5,3,1,2)]
head(countdata)
dim(countdata)

# experimental design
DataGroups <- c("RACM","RACM","RACM","PC","PC","PC")

# create DGE object of edgeR
dge <- DGEList(counts=countdata,group=factor(DataGroups))

# Filter by the filterByExpr function in the edgeR package.
keep.exprs <- filterByExpr(dge, group=DataGroups)
dge_filt <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge_filt$counts)
# sanity check for library sizes
dge_filt$samples$lib.size==colSums(dge_filt$counts)

# Extract Gene IDs
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
ensembl_list <- rownames(dge_filt)
gene_list <- getBM(filters= "ensembl_gene_id", 
                   attributes= c("ensembl_gene_id", "external_gene_name"),
                   values=ensembl_list,mart=mart)

# normalization using TMM method
dge_filt <- calcNormFactors(dge_filt, method="TMM")

#### Dispersion estimation ####
# Create the contrast matrix
design.mat <- model.matrix(~ 0 + dge_filt$samples$group)
colnames(design.mat) <- levels(dge_filt$samples$group)
# First, get common dispersion estimates - the overall BCV of the dataset, averaged over all genes:
dge_filt <- estimateGLMCommonDisp(dge_filt,design.mat)
# Then estimate gene-wise dispersion estimates
dge_filt <- estimateGLMTrendedDisp(dge_filt,design.mat, method="power")
dge_filt <- estimateGLMTagwiseDisp(dge_filt,design.mat)
# Plot the estimated dispersions:
plotBCV(dge_filt)

#### Get Normalized Counts, FPKM, TPM ####
# Output TMM normalized counts with edgeR (via code from https://www.biostars.org/p/317701/#317704)
dge_filt <- estimateCommonDisp(dge_filt)
dge_filt <- estimateTagwiseDisp(dge_filt)
norm_counts <- t(t(dge_filt$pseudo.counts)*(dge_filt$samples$norm.factors))
# get the gene names
norm_counts <- as.data.frame(norm_counts)
norm_counts <- add_column(norm_counts, gene_list[match(rownames(norm_counts),gene_list$ensembl_gene_id),], .before = 1)
#write.csv(norm_counts, file="./Data/normalized_counts_edger.csv", row.names = FALSE)

CPM <- cpm(dge_filt)

## Use normalized DGEList object to get FPKM ###
# get gene lengths
gene_lengths <- seqdata %>% filter(Geneid %in% rownames(dge_filt)) %>% pull(Length)
logRPKM <- rpkm(dge_filt, gene.length=gene_lengths, log=TRUE)
logRPKM_sorted <- logRPKM[order(row.names(logRPKM)), ]
# get the gene names
logRPKM_sorted <- as.data.frame(logRPKM_sorted)
logRPKM_sorted <- add_column(logRPKM_sorted, gene_list[match(rownames(logRPKM_sorted),gene_list$ensembl_gene_id),], .before = 1)

#write.csv(logRPKM_sorted, file="logFPKM.csv", row.names = FALSE)

## Calculate TPM from normalized DGEList object ##
TPM_table <- as.data.frame(dge_filt$counts)
gene_lengths <- seqdata %>% filter(Geneid %in% rownames(dge_filt)) %>% pull(Length)
TPM_table <- TPM_table/gene_lengths
# Absolute values
TPM <- TPM_table %>% apply(2,function(x){x/(sum(x)/1000000)})
colSums(TPM) # check if all equal to 10^6
TPM <- TPM[order(row.names(TPM)), ]
TPM <- as.data.frame(TPM)
TPM <- add_column(TPM, gene_list[match(rownames(TPM),gene_list$ensembl_gene_id),], .before = 1)
# log2TPM
logTPM <- log2(TPM[,3:8])
logTPM <- add_column(logTPM,TPM[,1:2],.before = 1)

#### Differential expression analysis ####

design.mat

contrasts <- makeContrasts(PC - RACM, levels=design.mat)

contrasts

# fit genewise glms
fit <- glmFit(dge_filt, design.mat)
# conduct likelihood ratio tests and show the top DE genes
lrt <- glmLRT(fit, contrast=contrasts)
edgeR_results <- topTags(lrt, n=Inf, adjust.method="BH", sort.by="PValue")
results <- edgeR_results$table
# add gene names
results <- add_column(results, gene_list[match(rownames(results),gene_list$ensembl_gene_id),], .before = 1)

# Show the total number of DE genes identified at an FDR of 0.05 
is.de <- decideTestsDGE(lrt)
summary(is.de)
#1*PC -1*RACM
#Down           1127
#NotSig        15915
#Up              973

#### Export results ####
df_list <- list("edger_results" = results, "TPM" = TPM, "norm_counts" = norm_counts)
write.xlsx(df_list, file = "./Data/results_by_Pvalue.xlsx")


#### Heatmap of genes ####
library(pheatmap)
logCPM <- cpm(dge_filt, prior.count=2, log=TRUE)
G_list_match <- gene_list[match(rownames(logCPM),gene_list$ensembl_gene_id),]
rownames(logCPM) <- G_list_match$external_gene_name
colnames(logCPM) <- c("RACM1", "RACM2", "RACM3", "PC1", "PC2", "PC3")
# select the logCPM values for the 20 top genes:
#logCPM <- logCPM[match(results$external_gene_name[1:20],rownames(logCPM)),]
# scale each row (each gene) to have mean zero and standard deviation one:
logCPM_norm <- t(scale(t(logCPM)))

genes_plot <- c("Isl1","Hcn4","Tbx3","Shox2","Bmp4","Tbx18",     #SAN>RA
                "Tbx5","Mef2c","Gata4",                          #SAN~RA
                "Nkx2-5","Nppa","Nppb","Bmp10","Gja1","Gja5")    #SAN<RA

heat <- pheatmap(
  logCPM_norm,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = FALSE
)
png("./Heatmap_VV.png", width = 30, height = 20, units = "cm", res = 300)
grid.draw(heat)
dev.off()

#### Data Exploration - MDS plot ####
# Set up colour schemes for groups
col.cell <- c("red","gold")[dge_filt$samples$group]
png("MDSplot.png")
plotMDS(dge_filt, col=col.cell)
legend("topleft",fill=c("red","gold"),legend=levels(dge_filt$samples$group))
title("Cell Populations")
dev.off()
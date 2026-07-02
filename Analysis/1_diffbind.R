library(DiffBind)
library(liteRnaSeqFanc)
library(DESeq2)

setwd("/scratch/mmoore/Epitherapy3D/analysis/HiC/paper_final/cutrun")
samples = read.csv("epitherapy_3D_diffbind.csv")
samples = samples[1:12,]
samples$Peaks = gsub("0.000","0.0000",samples$Peaks)

samples_dba = dba(sampleSheet=samples) %>% dba.count(summits = TRUE)
peaks = samples_dba$peaks[[1]]
summits = peaks[,c(1,9)]
summits[,3] = summits$Summits+1
write.table(summits,"1_diffbind/peaks.summits", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

samples_dba = dba(sampleSheet=samples) %>% dba.count(summits = FALSE)
dbo <- dba(samples_dba) %>% dba.contrast(minMembers = 2) %>% dba.analyze(method = DBA_DESEQ2, bRetrieveAnalysis=FALSE)
ds2 <- dba.analyze(dbo, bRetrieveAnalysis = TRUE)
pn <- dba.peakset(dbo, bRetrieve = T, bRemoveM = T, bRemoveRandom = T) %>% utilsFanc::gr.get.loci()
ct <- ds2@assays@data$counts
rownames(ct) <- pn
saveRDS(ct, file = "1_diffbind/GSC_deseq2_mat.rds")

count_df <- readRDS(file = "1_diffbind/GSC_deseq2_mat.rds")
# Create a sample information data frame by separating the column names of the counts data frame
sample_info <- strsplit(colnames(count_df), split = "_")
sample_df <- do.call(rbind, sample_info)
colnames(sample_df) <- c("Cell", "Treatment", "Rep")
sample_df <- as.data.frame(sample_df)
rownames(sample_df) <- paste0(sample_df$Cell, "_", sample_df$Treatment, "_", sample_df$Rep)
colnames(count_df) <- rownames(sample_df)
sample_df$group <- paste0("GSC_",sample_df$Treatment)

dds <- DESeqDataSetFromMatrix(countData = count_df,
                              colData = sample_df,
                              design = ~group )
dds1 <- DESeq(dds)
saveRDS(dds1, file = "1_diffbind/GSC_dds_obj.rds")

#sizeFactors(dds1)
#B36_DMSO_BRep1 B36_DMSO_BRep2   B36_DP_BRep1   B36_DP_BRep2 B49_DMSO_BRep1 B49_DMSO_BRep2   B49_DP_BRep1   B49_DP_BRep2 B66_DMSO_BRep1 
#2.2464908      2.6753272      1.8682791      1.6629605      0.9540908      1.1810070      0.7529448      0.7253656      0.5066440 
#B66_DMSO_BRep2   B66_DP_BRep1   B66_DP_BRep2 
#0.5134536      0.7009572      0.6437500 

# Comp1: (treatment) = DP v. DMSO
res1 <- results(dds1, name = "group_GSC_DP_vs_GSC_DMSO")
saveRDS(res1, file = "1_diffbind/GSC_DPvDMSO.rds")

res1 <- readRDS(file = "1_diffbind/GSC_DPvDMSO.rds") %>% as.data.frame()
colnames(res1) <- paste("GSC_DPvDMSO", colnames(res1), sep = "_")
res1 <- res1 %>% rownames_to_column(var = "region_name")

# Join this data frame with the counts information from the DDS object
dds <- readRDS(file = "1_diffbind/GSC_dds_obj.rds")
ct <- counts(dds, normalized = TRUE) %>% as.data.frame %>% rownames_to_column(var = "region_name")
# Join the counts with the joint data frame
joint_ct <- left_join(ct, res1, by = "region_name")
# Write the joint counts file with DE information to a .csv
write.csv(joint_ct, file = "1_diffbind/GSC_joint_counts_DEG.csv")
joint_ct_bed = joint_ct %>% separate(col = region_name, into = c("chr","coord"), sep=":", remove=TRUE)
joint_ct_bed = joint_ct_bed %>% separate(col = coord,into = c("start", "end"),sep = "-",remove = TRUE)
write.table(joint_ct_bed, file = "1_diffbind/GSC_joint_counts_DEG.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

DP_peaks = joint_ct[joint_ct$GSC_DPvDMSO_padj<0.05 & joint_ct$GSC_DPvDMSO_log2FoldChange > 0.5,]
DP_peaks = DP_peaks %>%
     separate(col = region_name,into = c("chr", "coord"),sep = ":",remove = TRUE)
DP_peaks = DP_peaks %>%
     separate(col = coord,into = c("start", "end"),sep = "-",remove = TRUE)
write.table(DP_peaks,"1_diffbind/GSC_CTCF_upSig_peaks.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

DMSO_peaks = joint_ct[joint_ct$GSC_DPvDMSO_padj<0.05 & joint_ct$GSC_DPvDMSO_log2FoldChange < -0.5,]
DMSO_peaks = DMSO_peaks %>%
  separate(col = region_name,into = c("chr", "coord"),sep = ":",remove = TRUE)
DMSO_peaks = DMSO_peaks %>%
  separate(col = coord,into = c("start", "end"),sep = "-",remove = TRUE)
write.table(DMSO_peaks,"1_diffbind/GSC_CTCF_downSig_peaks.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

stable_peaks = joint_ct[joint_ct$GSC_DPvDMSO_padj>0.05,]
stable_peaks = stable_peaks %>%
  separate(col = region_name,into = c("chr", "coord"),sep = ":",remove = TRUE)
stable_peaks = stable_peaks %>%
  separate(col = coord,into = c("start", "end"),sep = "-",remove = TRUE)
write.table(stable_peaks,"1_diffbind/GSC_CTCF_nonSig_peaks.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



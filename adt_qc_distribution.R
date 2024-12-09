# David Zemmour
# R/4.0.1
# usage: adt_qc_distribution.R [path to FBC raw data with matrix.mtx.gz, genes.tsv.gz, barcodes.tsv.gz] [path to RNA raw data matrix.mtx.gz, genes.tsv.gz, barcodes.tsv.gz]
# usage: adt_qc_distribution.R [FBC/outs/raw_feature_bc_matrix] [GEX/outs/raw_feature_bc_matrix]

args = commandArgs(TRUE)
path_to_hto_adt = args[1]
path_to_rna = args[2]# /Volumes/GoogleDrive/My Drive/ImmgenT/20211126_Pilot_06/ImmgenT_Strategy2/GEX ="/Users/david/googledrive/ImmgenT/20200917_pilot_01/data/Immgen_hash_citeseq_short1/umi_count"

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

message("Load HTO ADT matrix")
htos_adt_reads = Read10X(data.dir = sprintf("%s/read_count", path_to_hto_adt), gene.column = 1, strip.suffix = T)
htos_adt_umi = Read10X(data.dir = sprintf("%s/umi_count", path_to_hto_adt), gene.column = 1, strip.suffix = T)

r = rownames(htos_adt_reads)
r = do.call(rbind, strsplit(x = r, split = "\\-"))[,1]
rownames(htos_adt_reads) = r
rownames(htos_adt_umi) = r
r = sort(r)
htos_adt_reads = htos_adt_reads[r,]
htos_adt_umi = htos_adt_umi[r,]

reads_umapped_percent = sum(htos_adt_reads["unmapped",]) / sum(htos_adt_reads)*100
reads_umapped_sum = sum(htos_adt_reads["unmapped",])

pdf("adt_qc_plot1_barplot_count.pdf", height = 10, width = 40, useDingbats = F)
reads_distrib_percent = rowSums(htos_adt_reads) / sum(htos_adt_reads) *100
reads_distrib = rowSums(htos_adt_reads)
umis_distrib =  rowSums(htos_adt_umi)
umis_distrib_percent =  rowSums(htos_adt_umi)  / sum(htos_adt_umi) *100
par(mar=c(10,3,3,3))
barplot(log10(reads_distrib), las = 2, main = "Number of reads (log10)", cex.names = 0.8)
barplot(log10(umis_distrib), las = 2, main = "Number of umis (log10)", cex.names = 0.8)
barplot(reads_distrib_percent,  las = 2, main = "Percent reads", cex.names = 0.8)
barplot(umis_distrib_percent,  las = 2, main = "Percent umis", cex.names = 0.8)

saturation_mean = sum(htos_adt_umi) / sum(htos_adt_reads) *100
saturation_distrib = rowSums(htos_adt_umi) / rowSums(htos_adt_reads) *100
barplot(saturation_distrib, las= 2, main = "Percent saturation", cex.names = 0.8)

# RNA vs ADT and saturation in cells
rna = Read10X(data.dir = sprintf("%s/raw_feature_bc_matrix", path_to_rna), gene.column = 1, strip.suffix = T)
drops = intersect(colnames(rna_raw), colnames(htos_adt_umi))
cells = Read10X(data.dir = sprintf("%s/filtered_feature_bc_matrix", path_to_rna), gene.column = 1, strip.suffix = T)
cells =  intersect(colnames(cells), colnames(htos_adt_reads))
rna = rna_raw[,drops]
htos_adt  = htos_adt_umi[,drops]

percent_reads_in_cells = sum(htos_adt_reads[,cells]) /  sum(htos_adt_reads) * 100

saturation_mean_cells = sum(htos_adt_umi[,cells]) / sum(htos_adt_reads[,cells]) *100
saturation_distrib_cells = rowSums(htos_adt_umi[,cells]) / rowSums(htos_adt_reads[,cells]) *100
barplot(saturation_distrib_cells, las= 2, main = "Percent saturation in cells", cex.names = 0.8)
dev.off()

pdf("adt_qc_plot1_rna_adt.pdf", height = 8, width = 8, useDingbats = F)
plot(log10(colSums(rna)), log10(colSums(htos_adt)),  pch = ".")
dev.off()

stats = data.frame(Protein.Symbol = names(reads_distrib_percent), reads_distrib_percent = signif(reads_distrib_percent,2), reads_distrib = reads_distrib, umis_distrib = umis_distrib, umis_distrib_percent  = signif(umis_distrib_percent, 2))
write.table(x = stats,file = "adt_qc_stats.txt",quote = F, sep = "\t", row.names = F)

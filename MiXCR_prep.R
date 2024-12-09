library(Seurat)
library(xlsx)

sc <- readRDS("seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_fullTCRinfo_productiveonly.Rds")

meta <- sc@meta.data

hts <- read.csv("feature_refs_2022panel.csv")[1:10,1]
hash_adt_to_remove = read.table("hash_adt_rows_to_remove_from_analysis.txt", header = F)[,1]
hts <- setdiff(hts, hash_adt_to_remove)


for(ht in hts) {
  print(ht)
  dir.create(paste("MiXCR", ht, sep="/"), showWarnings = TRUE, recursive = FALSE)	
  barcodes <- rownames(meta[meta$hash.ID == ht, ])
  df <- data.frame(barcodes)
  rownames(df) <- paste(ht, "cell", rownames(df), sep = "_")
  write.xlsx(x = df, file = "MiXCR/MiXCR_commands.xlsx", sheetName = ht, append = TRUE, col.names = FALSE)
  write.table(df, paste("MiXCR", ht, "barcodes.txt", sep="/"), sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
} 

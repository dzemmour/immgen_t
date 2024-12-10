### Save Doublets - o2 ###


args = commandArgs(TRUE)
path_to_so = args[1]

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))


so <- readRDS(path_to_so)

so$HTO_classification_orig=so$HTO_classification

so@assays$qc_stats_1 <- NULL
for(i in which(so$HTO_classification.simplified=="Doublet")){
  so$HTO_classification[i] <- "Doublet"
  max_hash<-max(so@assays$HTO@counts[,i])
  sum_hash<-sum(so@assays$HTO@counts[,i])
  max2_hash<-max(so@assays$HTO@counts[,i][-which.max(so@assays$HTO@counts[,i])])
  ratio_counts = max_hash/max2_hash
  
  so_max <- so[,so$HTO_classification == names(which.max(so@assays$HTO@counts[,i]))]
  median_max <- median(so_max@assays$HTO@counts[rownames(so_max@assays$HTO@counts) == names(which.max(so@assays$HTO@counts[,i]))])
  norm_max <- max_hash/median_max
  
  so_max2 <-  so[,so$HTO_classification == names(which.max(so@assays$HTO@counts[,i][-which.max(so@assays$HTO@counts[,i])]))]
  median_max2 <- median(so_max2@assays$HTO@counts[rownames(so_max2@assays$HTO@counts) == names(which.max(so@assays$HTO@counts[,i][-which.max(so@assays$HTO@counts[,i])]))])
  norm_max2 <- max2_hash/median_max2
  
  ratio_singlet = median_max/median_max2
  ratio_norm = norm_max/norm_max2
  
  if(max_hash>0.95*sum(so@assays$HTO@counts[,i])){
    so$HTO_classification[i]<-names(which.max(so@assays$HTO@counts[,i]))
    so$HTO_classification.global[i] <- "Rescued_Doublet"
  }
  
  else if(norm_max>15*norm_max2){
    so$HTO_classification[i]<-names(which.max(so@assays$HTO@counts[,i]))
    so$HTO_classification.global[i] <- "Rescued_Doublet_M2"
  }
  
}

so$HTO_classification.simplified = as.character(so$HTO_classification)
so$HTO_classification.simplified[so$HTO_classification.global == "Doublet"] = "Doublet"
Idents(so) <- 'HTO_classification.simplified'

png("make_seurat_rna_adt_hto-1.png",width=1000, height=2500)
RidgePlot(so, assay = "HTO", features = rownames(so[["HTO"]]), ncol = 2)
dev.off()


png("make_seurat_rna_adt_hto-2.png",width=1400, height=1600)

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

pairs(t(so@assays$HTO@data), cex = 4, upper.panel = panel.cor, pch = ".")
dev.off()

Idents(so) <- 'HTO_classification'

png("make_seurat_rna_adt_hto-3.png",width=1400, height=1200)

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}


pairs(t(so@assays$HTO@data), cex = 4, upper.panel = panel.cor, pch = ".")

so = ScaleData(so, assay = "HTO", features = rownames(so@assays$HTO), verbose = FALSE)
so = RunPCA(so, assay ="HTO", features = rownames(so@assays$HTO), reduction.name = "pca_htos", approx = FALSE)
so = RunTSNE(so, assay ="HTO", reduction = "pca_htos", reduction.name = "tsne_htos", dims = 1:ncol(so@reductions$pca_htos@cell.embeddings), perplexity = 100, verbose = T, check_duplicates = FALSE, max_iter = 500)

#so$HTO_classification.simplified = as.character(so$HTO_classification)
#so$HTO_classification.simplified[so$HTO_classification.global == "Doublet"] = "Doublet"
tmp = data.frame(so@meta.data)
tmp$dim1 = so@reductions$tsne_htos@cell.embeddings[,1]
tmp$dim2 = so@reductions$tsne_htos@cell.embeddings[,2]
#tmp$HTO_classification.simplified = as.character(tmp$HTO_classification)
#tmp$HTO_classification.simplified[tmp$HTO_classification.global == "Doublet"] = "Doublet"

#Adding Rescured category in HTO_classification.global
so$HTO_classification.global[so$HTO_classification.global == 'Negative' & so$HTO_classification.simplified != 'Negative'] = "Rescued"

theme_my = theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = HTO_classification.simplified), alpha = I(1), size = 1) + scale_color_manual(values = mypal) + theme_bw() + theme_my + facet_wrap(~HTO_classification.global, nrow = 2)

q = HTOHeatmap(so, assay = "HTO", ncells = 5000)


grid.arrange(
  grobs = list(p,q),
  widths = c(1,1),
  heights = c(1,1,1),
  nrow = 3,
  ncol = 2
  #layout_matrix = rbind(c(1,0),
  #                      c(1,0))
)

dev.off()


message("Saving seuratobject_withHTOADT_alldata.Rda and seuratobject_withHTOADT_singlet.Rda")

tmp = data.frame(so@meta.data)
# tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT))) %>% as.data.frame()
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_alldata_QCstats.txt", quote = F, sep = "\t", row.names = F, col.names = T)
so@assays$qc_stats_1 = tmp

#tmp = so@assays$Experiment_details["Hashtag"]
#tmp = unlist(strsplit(x = tmp, split= ","))
#tmp = unlist(strsplit(x = tmp, split= ";"))
print('hiii')
hts = unique(so$HTO_classification.simplified)
hts = hts[!hts %in% c("Doublet", "Negative")]
so$sample_name = NA
print('bywwww')
# for (i in hts) {
#   so$sample_name[so$HTO_classification.simplified == i] = tmp[grep(pattern = paste(i,"[^0-9]", sep = ""), tmp)]
# }



#so$HTO_classification.global <- factor(so$HTO_classification.global, levels = c('Singlet', 'Doublet','Negative','Rescued'))
Idents(so)<- 'HTO_classification.global' 
table<- t((table(Idents(so))))
#colnames(table)<- c('nSinglets','nDoublets','nNegatives','nRecovered')
write.table(x = table, file = "Demultiplexing_QC_stats.txt", quote = F, sep = "\t", row.names = F, col.names = T)


saveRDS(so, file = "seuratobject_withHTOADT_alldata.Rds")

so = so[,so$HTO_classification.global %in% c("Singlet", "Rescued", "Rescued_Doublet", "Rescued_Doublet_M2")]

Idents(so) <- 'HTO_classification.simplified'


saveRDS(so, file = "seuratobject_withHTOADT_singlet.Rds")
tmp = data.frame(so@meta.data)
# tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT))) %>% as.data.frame()
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_singlet_QCstats.txt", quote = F, sep = "\t", row.names = F, col.names = T)



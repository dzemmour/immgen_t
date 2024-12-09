args = commandArgs(TRUE)
sc_IGT = args[1] 
cutoffs = args[2]
IGT<- args[3]  # Manully set your ncount_ADT cutoff for IGT dataset
spleen<- args[4] #hashtag for spleen control to remove from t cells


libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell","Seurat", "ggplot2", "reshape2", "dplyr", "fitdistrplus", "ggExtra", "ggrastr", "gridExtra", "gatepoints", "rafalib", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))


# Load in seurat_object
sc_IGT = readRDS(sc_IGT)


IGT<- read.table(IGT)
IGT<- as.character(IGT)

spleen<- read.table(spleen)
spleen<- as.character(spleen)


cutoffs<- read.csv('T_cutoffs.csv', header = F)

DotPlotHeatmap = function(data = tmp, title = "ADT_seurat_clusters") {
  tmp = melt(data, variable.name = "SYMBOL")
  tmp$adt_thr = so@assays$ADT@meta.features[tmp$SYMBOL,"cutoff_crl_norm"]
  tmp2 = tmp %>% group_by(cluster, SYMBOL) %>% summarize(mean = mean(value), freq = mean(value > adt_thr)*100, n_cells = n()) %>% as.data.frame()
  ColorRamp = rev(rainbow(10, end = 4/6))
  tmp3 = dcast(data = tmp2, SYMBOL~cluster, value.var = "mean")
  rownames(tmp3) = tmp3$SYMBOL
  tmp3 = tmp3[,-1]
  tmp3 = tmp3[rowSums(tmp3[])>0,]
  hc_col = hclust(dist(1-cor(tmp3, method = "pearson")))
  hc_row = hclust(dist(1-cor(t(tmp3), method = "pearson")))

  tmp2$cluster = factor(tmp2$cluster, levels = hc_col$labels[hc_col$order])
  tmp2$SYMBOL = factor(tmp2$SYMBOL, levels = hc_row$labels[hc_row$order])

  p = ggplot(tmp2) + geom_point(aes(x = SYMBOL, y = cluster, color = mean, size = freq, alpha = freq)) + scale_color_gradient(low = "blue", high = "red") + ggtitle(label = title) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, vjust=0,  hjust=1, size=10), axis.text.y  = element_text(angle=0, size=15), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p)
}




Tcells<- sc_IGT

Seurat_rna_umap_clustering = function(seurat_object = so, tsne.method = c("FIt-SNE", "Rtsne"), print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna") {
  require(Seurat)
  seurat_object = NormalizeData(seurat_object, verbose = T, normalization.method = "LogNormalize", scale.factor = 10000,assay = "RNA")
  seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = T,assay = "RNA")
  DefaultAssay(seurat_object)<- 'RNA'
  genes<- rownames(seurat_object)
  Trav<- genes[grepl('^Trav',genes)]
  Traj<-genes[grepl('^Traj',genes)]
  Trac<-genes[grepl('^Trac',genes)]
  Trbv<-genes[grepl('^Trbv',genes)]
  Trbd<-genes[grepl('^Trbd',genes)]
  Trbj<-genes[grepl('^Trbj',genes)]
  Trd<-genes[grepl('^Trd',genes)]
  Trg<-genes[grepl('^Trg',genes)]
  TCR_genes<- (c(Trav,Traj,Trac,Trbv,Trbd,Trbj,Trd,Trg))
  sex_specific_genes<- c('Ddx3y','Uty','Xist','Eif2s3y','Kdm5d','Tsix')
  remove_var_genes<- c(TCR_genes,sex_specific_genes)
  seurat_object@assays$RNA@var.features<- seurat_object@assays$RNA@var.features[!seurat_object@assays$RNA@var.features %in% remove_var_genes]
  
  seurat_object = ScaleData(seurat_object, assay = "RNA")
  
  seurat_object = RunPCA(seurat_object, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, reduction.name = reduction.name.pca,assay = "RNA", var.features=seurat_object@assays$RNA@var.features)
  
  x = cumsum((seurat_object@reductions[[reduction.name.pca]]@stdev**2 / sum(seurat_object@reductions[[reduction.name.pca]]@stdev**2)))
  ndims = min(which(x >= 0.8))
  ElbowPlot(object = seurat_object, ndims = 100, reduction = reduction.name.pca) + geom_vline(mapping = aes(xintercept = ndims), color = "red", linetype="dashed") + ggtitle(label = sprintf("ndims PCA for further analysis = %s", ndims))
  DimHeatmap(seurat_object, dims = c(1:3, 70:75), cells = 500, balanced = TRUE, reduction = reduction.name.pca)
  
  #message("Running tSNE...")
  #seurat_object = RunTSNE(seurat_object, dims = 1:ndims, tsne.method = tsne.method, nthreads = 8, max_iter = 2000)
  message("Running UMAP")
  seurat_object = RunUMAP(seurat_object, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)
  
  message("RNA clusters...")
  seurat_object = FindNeighbors(seurat_object, reduction = reduction.name.pca, dims = 1:ndims, k.param = 20, verbose = T)
  
  seurat_object = FindClusters(seurat_object, resolution = 0.25, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 0.5, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 1, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 1.5, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 2, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 3, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 4, n.start = 10, algorithm = 1, n.iter = 10)
  #colnames(seurat_object@meta.data) = gsub(x = colnames(seurat_object@meta.data), pattern = "^(.*?)_snn_res.", replacement = "RNACluster_Res")
  
  if (print_clusters == T) {
    message("Printing clusters...")
    tmp = data.frame(seurat_object@meta.data, dim1 = seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[,2] )
    for( i in grep("RNA_snn_res.", colnames(seurat_object@meta.data))) {
      p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = tmp[,i]), size = I(1), alpha = I(1)) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, colnames(seurat_object@meta.data)[i]))
      p2 = LabelClusters(p, id =  colnames(seurat_object@meta.data)[i], size = 5)
      print(p2)
    }
  }
  
  
  return(seurat_object)
}

DefaultAssay(Tcells)<-"RNA"
Final_Tcells= Seurat_rna_umap_clustering(seurat_object = Tcells, tsne.method = "FIt-SNE", print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna")
reduction.name.pca = "pca_adt"
reduction.name.umap = "umap_adt"                                
Seurat_adt_umap_clustering = function(seurat_object = so, VariableFeatures =  rownames(so[["ADT"]])[!grepl(rownames(so[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = "pca_adt", reduction.name.umap = "umap_adt") {
  so = seurat_object
  DefaultAssay(so) = "ADT"
  message("CLR normalization")
  so = NormalizeData(so, assay = "ADT", normalization.method = "CLR") 
  message("VariableFeatures:")
  print(VariableFeatures)
  VariableFeatures(so, assay = "ADT") = VariableFeatures
  so = ScaleData(so, assay = "ADT")
  so = RunPCA(object = so, assay = "ADT", reduction.name = reduction.name.pca)
  x = cumsum((so@reductions[[reduction.name.pca]]@stdev**2 / sum(so@reductions[[reduction.name.pca]]@stdev**2)))
  ndims = min(which(x >= 0.8))
  so = FindNeighbors(so, reduction = "pca_adt", dims = 1:ndims, k.param = 20, verbose = T)
  message("UMAP ADT")
  so = RunUMAP(so, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)
  message("ADT clustering")
  so = FindClusters(so, resolution = 1, n.start = 10, algorithm = 1, n.iter = 10)
  return(so)
}

DefaultAssay(Final_Tcells)<-"ADT"
Final_Tcells= Seurat_adt_umap_clustering(seurat_object = Final_Tcells, VariableFeatures =  rownames(Final_Tcells[["ADT"]])[!grepl(rownames(Final_Tcells[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = reduction.name.pca, reduction.name.umap = reduction.name.umap)


Idents(Final_Tcells)<- 'sample_name'
levels(Final_Tcells)
saveRDS(Final_Tcells,"seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds")


so<- Final_Tcells
cluster <- Final_Tcells$seurat_clusters

tmp = data.frame(t(as.matrix(Final_Tcells@assays$ADT@data)), cluster = Final_Tcells$seurat_clusters)
pdf(file = "adt_qc_postTfilt__Dotplot_ADT_ADTclusters.pdf", 25, 7, useDingbats=FALSE)
DotPlotHeatmap(data = tmp, title = "ADT_seurat_clusters")
dev.off()



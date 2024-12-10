##David Zemmour
##20210304
##Run Seurat analysis (PCA, tSNE, UMAP, Clustering) using Seurat V4
##Usage:
#source("/n/groups/cbdm_lab/dp133/scripts_o2/R_scripts/immgen_t/seurat_rna_umap_clustering.R")
#pdf("PC_heatmap.pdf", 10, 10)
#so = Seurat_rna_umap_clustering(seurat_object = so, tsne.method = "FIt-SNE", print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna")
#dev.off()

Seurat_rna_umap_clustering = function(seurat_object = so, tsne.method = c("FIt-SNE", "Rtsne"), print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna") {
  require(Seurat)
  source("/n/groups/cbdm_lab/dp133/bin/FIt-SNE-master/fast_tsne.R", chdir = T) #path to fast_tsne.R in order to use FIt-SNE
  seurat_object = NormalizeData(seurat_object, verbose = T, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = T)
  seurat_object = ScaleData(seurat_object, assay = "RNA")
  
  seurat_object = RunPCA(seurat_object, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, reduction.name = reduction.name.pca)
  
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
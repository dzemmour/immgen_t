args = commandArgs(TRUE)
path_to_data = args[1]
path_to_data2 = args[2]

# startup #####
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(cowplot)
library(tidyverse)

# load dataset_clean into sc, and dataset into sc.all. The table from sc shows the hashtags you need to include. 

sc<-readRDS(path_to_data)

table(sc@meta.data$HTO_classification.simplified)

sc.all<-readRDS(path_to_data2)

# Tables #####

#outputs cell typeand sample tables that I cut and paste into the PPT alongside the t.filtering and cluster_umap plots
cell_type_table<-table(sc.all@meta.data$cell_type)
sample_table<-table(sc@meta.data$sample_name)

write.csv(cell_type_table, file = "cell_type_table.csv")
write.csv(sample_table, file = "sample_table.csv")

# nonT cell annotation#######
#Shows cell type Ident in both rna and adt reductions. Puts the two side by side and ouputs as a png. The colors should be consistent between data sets: T : grey, B : red, B-T : salmon, NA and n/a : black, ILC : green, MNP : blue.  
DefaultAssay(sc.all)<-"RNA"

Idents(sc.all)<-"cell_type"
t.filtering.rna<-DimPlot(sc.all, order = c("MNP", "ILC","NA", "n/a", "B", "B-T like cells", "T"), cols = c("grey","red", "salmon", "black", "black", "green", "blue"), reduction = "umap_rna", pt.size = 1)
t.filtering.adt<-DimPlot(sc.all, order = c("MNP", "ILC","NA", "n/a", "B", "B-T like cells", "T"), cols = c("grey","red", "salmon", "black", "black", "green", "blue"), reduction = "umap_adt", pt.size = 1)

t.filtering<-ggdraw() +
  draw_plot(t.filtering.rna, x=0, y=0.25, width=0.5, height=0.5) + 
  draw_plot(t.filtering.adt, x=0.5, y=0.25, width=0.5, height=0.5)


# NonT Feature Plots -- expression scales benchmarked to IGT1 Spleen standard#####
#creates a collection of ADT feature plots from filtered cells. Could potentially replace the cell type subsetting with just including any cells in dataset but not dataset_clean.
sc.nonT<-subset(sc.all, idents =c("T"),invert = T )
DefaultAssay(sc.nonT)<-"ADT"

cd19_ADT_plot<-FeaturePlot(sc.nonT, features = "CD19", reduction = "umap_adt", order=TRUE, min.cutoff = 0.5)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 5))

cd138_ADT_plot<-FeaturePlot(sc.nonT, features = "CD138", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 4))

klrbc_ADT_plot<-FeaturePlot(sc.nonT, features = "KLRBC-NK1.1", reduction = "umap_adt", order=TRUE, max.cutoff= 2.5, min.cutoff = 0.5)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 3))

itgam_ADT_plot<-FeaturePlot(sc.nonT, features = "ITAM.CD11B", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 6))

itgax_ADT_plot<-FeaturePlot(sc.nonT, features = "ITAX.CD11C", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 5))

bst2_ADT_plot<-FeaturePlot(sc.nonT, features = "BST2.PDCA1.CD317", reduction = "umap_adt", order=TRUE, min.cutoff = 0.5)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 3))

gr1_ADT_plot<-FeaturePlot(sc.nonT, features = "GR1-LY6G-LY6C1-LY6C2", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 6))

fcera_ADT_plot<-FeaturePlot(sc.nonT, features = "FCERA", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)& scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 4))

ter119_ADT_plot<-FeaturePlot(sc.nonT, features = "TER119", reduction = "umap_adt", order=TRUE, min.cutoff = 0.5) & scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0.25, 3))



# Feature plots #####
#creates a collection RNA and ADT Feature plots. Christophe is pretty particular on the order, scale and arrangement of these plots. max and min cutoffs set off of IGT1 spleen control expression 

DefaultAssay(sc)<-"ADT"

trb_ADT_plot<-FeaturePlot(sc, features = "TCRB", reduction = "umap_adt", order=TRUE, max.cutoff = 2, min.cutoff=0)

trac_RNA_plot<-FeaturePlot(sc, features = "Trac", reduction="umap_rna", order = TRUE, max.cutoff = 3, min.cutoff = -1)

trg_RNA_plot<-FeaturePlot(sc, features = c("Trdc"), order=TRUE, reduction = "umap_rna", max.cutoff = 3, min.cutoff=0)
trg_ADT_plot<-FeaturePlot(sc, features = c("TCRGD"), order=TRUE, reduction = "umap_adt", max.cutoff = 3, min.cutoff=0.75)

cd4_RNA_plot<-FeaturePlot(sc, features = "Cd4", reduction="umap_rna", order = TRUE)
cd4_ADT_plot<-FeaturePlot(sc, features = "CD4", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)

cd8_RNA_plot<-FeaturePlot(sc, features = "Cd8b1", reduction="umap_rna", order = TRUE)
cd8_ADT_plot<-FeaturePlot(sc, features = "CD8B", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)

foxp3_RNA_plot<-FeaturePlot(sc, features = "Foxp3", reduction = "umap_rna", order=TRUE)
il2ra_RNA_plot<-FeaturePlot(sc, features = "Il2ra", reduction = "umap_rna", order=TRUE)
il2ra_ADT_plot<-FeaturePlot(sc, features = "IL2RA.CD25", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)

cd44_RNA_plot <- FeaturePlot(sc, features = "Cd44", reduction = "umap_rna", order=TRUE)
cd44_ADT_plot <- FeaturePlot(sc, features = "CD44", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)

cd62l_RNA_plot <- FeaturePlot(sc, features = "Sell", reduction = "umap_rna", order=TRUE)
cd62l_ADT_plot <- FeaturePlot(sc, features = "CD62L", reduction = "umap_adt", order=TRUE, min.cutoff = 0.25)


# umap by sample ####
#Plots umaps in a 2 x 2 grid for ADT and RNA clusters and reduction. Uses these umaps to set min and max UMAP values for subsequent plots. 

Idents(sc)<-"RNA_clusters"
rna_protein<-DimPlot(sc, reduction = "umap_adt")
rna_rna<-DimPlot(sc, reduction = "umap_rna")

Idents(sc)<-"Protein_clusters"
protein_protein<-DimPlot(sc, reduction = "umap_adt")
protein_rna<-DimPlot(sc, reduction = "umap_rna")

cluster_umaps<- ggdraw() +
  draw_plot(protein_protein, x=0, y=0.5, width = 0.5, height = 0.5)+
  draw_plot(protein_rna, x=0.5, y=0.5, width = 0.5, height = 0.5)+
  draw_plot(rna_protein, x=0, y=0, width = 0.5, height = 0.5)+
  draw_plot(rna_rna, x=0.5, y=0, width = 0.5, height = 0.5)

y.adt<-layer_scales(protein_protein)$y$range$range

x.adt<-layer_scales(protein_protein)$x$range$range

y.rna<-layer_scales(rna_rna)$y$range$range

x.rna<-layer_scales(rna_rna)$x$range$range

# Isotype featureplots ####
# featureplots of isotypes
iso.h.igg<-FeaturePlot(sc, features = "Isotype.Hamster.IgG", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.m.igg1k<-FeaturePlot(sc, features = "Isotype.Mouse.IgG1k", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.m.igg2ak<-FeaturePlot(sc, features = "Isotype.Mouse.IgG2ak", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.m.igg2bk<-FeaturePlot(sc, features = "Isotype.Mouse.IgG2bk", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.r.igg1bk<-FeaturePlot(sc, features = "Isotype.Rat.IgG1bk", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.r.igg1bl<-FeaturePlot(sc, features = "Isotype.Rat.IgG1bl", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.r.igg2ak<-FeaturePlot(sc, features = "Isotype.Rat.IgG2ak", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.r.igg2bk<-FeaturePlot(sc, features = "Isotype.Rat.IgG2bk", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")
iso.r.igg2ck<-FeaturePlot(sc, features = "Isotype.Rat.IgG2ck", order=TRUE, max.cutoff = 2, min.cutoff = 1, keep.scale = "all", reduction = "umap_adt")

# combined plots #####
# arranges the plots and outputs as pngs. Christophe set the layout. 
feature_tcrs<-ggdraw() +
  draw_plot(trb_ADT_plot, x=0, y=0.5, width = 0.5, height=0.464)+
  draw_plot(trac_RNA_plot, x=0, y=0, width = 0.5, height=0.464)+
  draw_plot(trg_ADT_plot, x=0.5, y=0.5, width = 0.5, height=0.464)+
  draw_plot(trg_RNA_plot, x=0.5, y=0, width = 0.5, height=0.464)

feature_cd4.8<-ggdraw() +
  draw_plot(cd4_ADT_plot, x=0, y=0.5, width = 0.5, height=0.464)+
  draw_plot(cd4_RNA_plot, x=0, y=0, width = 0.5, height=0.464)+
  draw_plot(cd8_ADT_plot, x=0.5, y=0.5, width = 0.5, height=0.464)+
  draw_plot(cd8_RNA_plot, x=0.5, y=0, width = 0.5, height=0.464)

feature_treg<-ggdraw() +
  draw_plot(il2ra_RNA_plot, x=0, y=0.5, width = 0.5, height=0.464) +
  draw_plot(il2ra_ADT_plot, x=0.5, y=0.5, width = 0.5, height=0.464)+
  draw_plot(foxp3_RNA_plot, x=0, y=, width = 0.5, height=0.464)

feature_active<-ggdraw() +
  draw_plot(cd44_ADT_plot, x=0, y=0.5, width = 0.5, height=0.464)+
  draw_plot(cd44_RNA_plot, x=0, y=0, width = 0.5, height=0.464)+
  draw_plot(cd62l_ADT_plot, x=0.5, y=0.5, width = 0.5, height=0.464)+
  draw_plot(cd62l_RNA_plot, x=0.5, y=0, width = 0.5, height=0.464)

feature_nonT<-ggarrange(cd19_ADT_plot, cd138_ADT_plot, klrbc_ADT_plot, itgam_ADT_plot, itgax_ADT_plot, bst2_ADT_plot, gr1_ADT_plot, fcera_ADT_plot, ter119_ADT_plot)

iso<-ggarrange(iso.h.igg, iso.m.igg1k, iso.m.igg2ak, iso.m.igg2bk, iso.r.igg1bk, iso.r.igg1bl, iso.r.igg2ak, iso.r.igg2bk, iso.r.igg2ck, ncol = 3, nrow = 3)



# pdf generation #####

pdf(file = "Post_Tcell_analysis.pdf", width = 13.3, height = 7.5)

t.filtering
feature_nonT
cluster_umaps
DimPlot(sc, split.by = "sample_name",group.by = 'RNA_clusters', reduction = "umap_rna", pt.size = 0.1,ncol=4) + NoLegend() + ylim(y.rna) + xlim(x.rna)
DimPlot(sc, split.by = "sample_name",group.by = 'Protein_clusters', reduction = "umap_adt", pt.size = 0.1,ncol=4) + NoLegend() + ylim(y.adt) + xlim(x.adt)
feature_tcrs
feature_cd4.8
feature_treg
feature_active
iso

dev.off()


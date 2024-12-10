# in immgent wrapper, write mkdir -p RNA and mkdir -p ADT 

args = commandArgs(TRUE)
path_to_data = args[1] #path_to_rna = "combined_sample/outs/raw_feature_bc_matrix/"
path_to_data2 = args[2]# path to filtered data = "combined_sample/outs/filtered_feature_bc_matrix"
main_output_dir = args[3] #
IGT = args[4]

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr", "DropletUtils") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))



IGT<- read.table(IGT)
IGT<- as.character(IGT)
sce <- read10xCounts(path_to_data)
barcodes = sce$Barcode
new_barcodes = paste0(IGT,'.',barcodes)
sce$Barcode<-new_barcodes



filtered_data <- read10xCounts(path_to_data2)
barcodes = filtered_data$Barcode
new_barcodes = paste0(IGT,'.',barcodes)
filtered_data$Barcode<-new_barcodes

anti_index<- which(sce@rowRanges@elementMetadata$Type == 'Antibody Capture')    
anti.counts <- (counts(sce[anti_index,]))
rna.counts <- (counts(sce[-anti_index,]))
anti.counts2 <- ((sce[anti_index,]))
rna.counts2 <- ((sce[-anti_index,]))
rna.counts2@rowRanges@elementMetadata
anti.counts2@rowRanges@elementMetadata

# Making adt elbow plot.

pdf(file = paste0(main_output_dir,'/adt_elbow_plot.pdf'), width = 10, height = 10,useDingbats=FALSE)
br.out_anti <- barcodeRanks(anti.counts)
plot(br.out_anti$rank, br.out_anti$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out_anti$rank)
lines(br.out_anti$rank[o], br.out_anti$fitted[o], col="red")
abline(h=metadata(br.out_anti)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out_anti)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
dev.off()

# Making rna elbow plot.

pdf(file = paste0(main_output_dir,'/rna_elbow_plot.pdf'), width = 10, height = 10,useDingbats=FALSE)
br.out_rna <- barcodeRanks(rna.counts)
plot(br.out_rna$rank, br.out_rna$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out_rna$rank)
lines(br.out_rna$rank[o], br.out_rna$fitted[o], col="red")
abline(h=metadata(br.out_rna)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out_rna)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
dev.off()

cells_of_interest_anti<- which(br.out_anti$total > br.out_anti@metadata$inflection)
cells_of_interest_rna<- which(br.out_rna$total > br.out_rna@metadata$inflection)

cells_of_interest_antiM<- anti.counts2[,cells_of_interest_anti]
cells_of_interest_rnaM<- rna.counts2[,cells_of_interest_rna]
cells_of_interest_rnaM <- cells_of_interest_rnaM[,cells_of_interest_rnaM$Barcode %in% cells_of_interest_antiM$Barcode,]
cells_of_interest_antiM <- cells_of_interest_antiM[,cells_of_interest_antiM$Barcode %in% cells_of_interest_rnaM$Barcode,]
cells_of_interest_antiM <-  cells_of_interest_antiM[,cells_of_interest_antiM$Barcode %in% filtered_data$Barcode,]
cells_of_interest_rnaM <-  cells_of_interest_rnaM[,cells_of_interest_rnaM$Barcode %in% filtered_data$Barcode,]


barcodes = cells_of_interest_antiM$Barcode
gene.id = rownames(cells_of_interest_antiM)
gene.symbol = gene.id
write10xCounts(
  "FBC/filtered_feature_bc_matrix",
  counts(cells_of_interest_antiM),
  barcodes = cells_of_interest_antiM$Barcode,
  gene.id = rownames(cells_of_interest_antiM),
  gene.symbol = gene.id,
  overwrite = TRUE
)



barcodes = cells_of_interest_antiM$Barcode
new_barcodes = gsub(paste0(IGT,'.'),'',barcodes)
cells_of_interest_antiM$Barcode<-new_barcodes




barcodes = cells_of_interest_antiM$Barcode
gene.id = rownames(cells_of_interest_antiM)
gene.symbol = gene.id
write10xCounts(
  "FBC/filtered_feature_bc_matrix2",
  counts(cells_of_interest_antiM),
  barcodes = cells_of_interest_antiM$Barcode,
  gene.id = rownames(cells_of_interest_antiM),
  gene.symbol = gene.id,
  overwrite = TRUE
)





barcodes = cells_of_interest_rnaM$Barcode
gene.id = cells_of_interest_rnaM@rowRanges@elementMetadata$Symbol
gene.symbol = gene.id
write10xCounts("GEX/filtered_feature_bc_matrix",
  counts(cells_of_interest_rnaM),
  barcodes = cells_of_interest_rnaM$Barcode,
  gene.id = cells_of_interest_rnaM@rowRanges@elementMetadata$Symbol,
  gene.symbol = gene.id,
  overwrite = TRUE
)


data1 <- Read10X(data.dir = path_to_data)
barcodes = data1$Barcode
data3 <- Read10X(data.dir ="FBC/filtered_feature_bc_matrix2")


seurat_gene<- (data1$`Gene Expression`)
seurat_anti <- (data1$`Antibody Capture`)
seurat_gene2 <- data3

umicount_gene <- (as.matrix(colSums(data1$`Gene Expression`)))
umicount_antibodies <- (as.matrix(colSums(data1$`Antibody Capture`)))
data<- as.data.frame(cbind(umicount_gene,umicount_antibodies))
data[data==0] <- NA
data <-na.omit(data)
data <- log2(data)


add_data <- data[rownames(data) %in% colnames(seurat_gene2),]
print('hello')
head(add_data,10)
subtract_data <- data[!rownames(data) %in% colnames(seurat_gene2),]
head(subtract_data,10)
print('HI')
add_data <- cbind(add_data, 'Cells')
colnames(add_data)[3] <- "Classification"
print('HI2')
subtract_data <- cbind(subtract_data, 'Empty Droplets')
colnames(subtract_data)[3] <- "Classification"
print('HI3')
data_total <- rbind(add_data, subtract_data)
print('HI4')



pdf(file = paste0(main_output_dir,'/ADT_RNA_cutoff_plot.pdf'), width = 20, height = 20,useDingbats=FALSE)
ggplot(data_total,aes(x=V1,y=V2,col=Classification))+geom_point()+
  xlab("Counts (RNA, log2)") + ylab("Counts (antibody, log2)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        legend.text=element_text(size=10), axis.title.x = element_text(size=20) ,
        axis.title.y = element_text(size=20)) +
  scale_color_manual(values=c('Red','Black'))
dev.off()
















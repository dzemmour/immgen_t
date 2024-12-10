# David Zemmour
# R/4.0.1
library(Seurat)
library(Matrix)
library(R.utils)
options(stringsAsFactors = FALSE)

filt.sc.data <- Read10X("combined_sample/outs/filtered_feature_bc_matrix",strip.suffix=T)
filt_gex_mat <- filt.sc.data[[1]]
filt_fbc_mat <- filt.sc.data[[2]]

raw.sc.data <- Read10X("combined_sample/outs/raw_feature_bc_matrix",strip.suffix=T)
raw_gex_mat <- raw.sc.data[[1]]
raw_fbc_mat <- raw.sc.data[[2]]

write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, gene.type, version="2") {
  dir.create(path, showWarnings=FALSE)
  gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)
  
  if (version=="3") {
    gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
    mhandle <- file.path(path, "matrix.mtx")
    bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
    fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
    on.exit({
      close(bhandle)
      close(fhandle)
    })
  } else {
    mhandle <- file.path(path, "matrix.mtx")
    bhandle <- file.path(path, "barcodes.tsv")
    fhandle <- file.path(path, "genes.tsv")
  }
  
  writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  if (version=="3") {
    # Annoyingly, writeMM doesn't take connection objects.
    gzip(mhandle)
  }
  
  return(NULL)
}

type_chooser <- function(path, type) {
  if (type=="auto") {
    type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
  }
  type
}

write10xCounts <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id, gene.type="Gene Expression",
                           overwrite=FALSE, type=c("auto", "sparse", "HDF5"), genome="unknown", version=c("2", "3"),
                           chemistry="Single Cell 3' v3", original.gem.groups=1L, library.ids="custom")
{
  # Doing all the work on a temporary location next to 'path', as we have permissions there.
  # This avoids problems with 'path' already existing.
  temp.path <- tempfile(tmpdir=dirname(path)) 
  on.exit({ 
    if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
  })
  
  # Checking the values.
  if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
    stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
  }
  if (ncol(x)!=length(barcodes)) { 
    stop("'barcodes' must of of the same length as 'ncol(x)'")
  }
  
  # Determining what format to save in.
  version <- match.arg(version)
  type <- type_chooser(path, match.arg(type))
  if (type=="sparse") {
    write_sparse(temp.path, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
  } else {
    .write_hdf5(temp.path, genome, x, barcodes, gene.id, gene.symbol, gene.type, version=version)
  }
  
  # We don't put this at the top as the write functions might fail; 
  # in which case, we would have deleted the existing 'path' for nothing.
  if (overwrite) {
    unlink(path, recursive=TRUE)
  } else if (file.exists(path)) { 
    stop("specified 'path' already exists")
  }
  file.rename(temp.path, path)
  return(invisible(TRUE))
}

write10xCounts(x = filt_gex_mat, path = "GEX/outs/filtered_feature_bc_matrix", version = "3")
write10xCounts(x = filt_fbc_mat, path = "FBC/outs/filtered_feature_bc_matrix", version = "3", gene.type="Antibody Capture")

write10xCounts(x = raw_gex_mat, path = "GEX/outs/raw_feature_bc_matrix", version = "3")
write10xCounts(x = raw_gex_mat, path = "FBC/outs/raw_feature_bc_matrix", version = "3", gene.type="Antibody Capture")

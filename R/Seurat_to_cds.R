#'Seurat to cell data set conversion function.
#'This function takes in a Seurat object and returns a cell data set object.
#'Automatically runs PCA on the new cell data set object.
#'@title seurat_to_cds
#'@description Converts Seurat object to preprocessed cds.
#'@import monocle3
#'@import dplyr
#'@details I've only included this to stop warnings from document()
#'@param seurat The Seurat object to be converted.
#'@param assay Assay within the Seurat object to be imported.
#'@param import_umap Imports UMAP projection into the cell data set object. Default is TRUE.
#'@param preprocess Determines whether the new cds will be preprocessed. Default is TRUE.
#'@param dims Number of principal components returned after preprocessing.
#'@param norm Normalization method used in preprocess_cds. Default is "none".
#'@return A cell data set (cds) object
#'@examples
#' new_cds <- seurat_to_cds(Seurat.object, assay = "SCT")
#'@export
seurat_to_cds <- function(seurat, assay = "SCT", preprocess = TRUE, dims = 50, norm = "none",
                          import_umap = TRUE) {
  #Determine assay to pull data from
  if (assay == "SCT") {
    exp_matrix <- seurat@assays$SCT@data
  } else if (assay == "integrated") {
    exp_matrix <- seurat@assays$integrated@data
  } else {
    exp_matrix <- seurat@assays$RNA@data
  }

  gene_annotation <- as.data.frame(exp_matrix@Dimnames[[1]])
  row.names(gene_annotation) <- exp_matrix@Dimnames[[1]]
  colnames(gene_annotation) <- "gene_short_name"
  metadata <- seurat@meta.data

  #construct the cell_data_set object
  cds <- new_cell_data_set(exp_matrix,
                             cell_metadata = metadata,
                             gene_metadata = gene_annotation)
  if (preprocess) {
    cds <- preprocess_cds(cds, num_dim = dims, norm_method = norm,
                          method = "PCA", scaling = FALSE)
  }
  if (import_umap) {
    cds@int_colData@listData$reducedDims$UMAP <- seurat@reductions$umap@cell.embeddings
  }
  return(cds)
}

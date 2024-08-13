#' Differential Expression Analysis Using Wilcoxon rank-sum Test
#'
#' This function performs differential expression analysis using the Wilcoxon rank-sum test on cell-level data.
#' It is designed to compare two groups within a given cell type.
#'
#' @param i Indicator of cell type.
#' @param object A Seurat object containing the RNA data.
#' @param ident.1 The first group identifier for comparison (default: "disease").
#' @param ident.2 The second group identifier for comparison (default: "control").
#' @param min.cells.feature Minimum number of cells expressing a feature to include (default: 0).
#' @param min.cells.group Minimum number of cells in each group to include a feature (default: 0).
#' @param logfc.threshold The threshold for the log fold change to filter genes (default: 0).
#' @param min.pct Minimum detection percentage required to keep a gene (default: 0).
#' @param test.use Statistical method to use (default: "wilcox").
#' @param verbose Whether to show progress (default: TRUE).
#' @param group_label The label of the group (disease vs. control or treatment vs. control).
#'
#' @return A data frame containing the results of the differential expression analysis.
#'
#' @export
#'

DEs.wilcoxon <- function(i, object, ident.1 = "disease", ident.2 = "control",
                       min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0,
                       min.pct = 0, test.use = "wilcox", verbose = TRUE, group_label){

  object$celltype.status <- paste(Idents(object), object@meta.data[[group_label]], sep = "_")
  object$celltype <- Idents(object)
  Idents(object) <- "celltype.status"

  # Define the group labels for the given cell type
  clst <- i # given a cell type
  clst_status.1 <- paste(clst, ident.1, sep = "_")
  clst_status.2 <- paste(clst, ident.2, sep = "_")

  # Perform the differential expression analysis using FindMarkers
  results <- FindMarkers(object, ident.1 = clst_status.1, ident.2 = clst_status.2,
                         min.cells.feature = min.cells.feature, min.cells.group = min.cells.group,
                         verbose = verbose, logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use)

  results <- results %>%
    mutate(cluster = i, gene = rownames(results),
           FC = 2^(avg_log2FC), FC = ifelse(FC<1, -1/FC, FC)) %>%
    dplyr::select(p_val, FC, pct.1, pct.2, p_val_adj, cluster, gene, avg_log2FC) %>%
    arrange(p_val)

  return(results)
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables(c("avg_log2FC", "FC", "p_val", "pct.1", "pct.2", "p_val_adj", "gene"))


#' Differential Expression Analysis Using Wilcoxon Clustered Test
#'
#' This function performs differential expression analysis using the Wilcoxon clustered test on cell-level data.
#' It is designed to compare two groups within a given cell type, accounting for clustering by subject.
#'
#' @param i Indicator of cell type.
#' @param object A Seurat object containing the RNA data.
#' @param ident.1 The first group identifier for comparison (default: "disease").
#' @param ident.2 The second group identifier for comparison (default: "control").
#' @param subject_ID Subject ID used for identifying subjects.
#' @param logfc.threshold The threshold for the log fold change to filter genes.
#' @param min.pct The minimum detection percentage required to keep a gene.
#' @param min.num_subjects The minimum number of subjects required to keep a gene.
#' @param min.unique_counts The minimum number of unique counts required to keep a gene.
#' @param test.use Statistical test to use (default: "clusWilcox").
#' @param verbose Whether to show progress (default: TRUE).
#' @param parallel Whether to run the function in parallel (default: FALSE).
#' @param group_label Label of the group used for identifying cells.
#'
#' @return A data frame containing the results of the differential expression analysis.
#'
#' @export
#'

DEs.wilcoxon.clustered <- function(i, object, ident.1 = "disease", ident.2 = "control", subject_ID,
                                   logfc.threshold = 0, min.pct = 0, min.num_subjects = 0, min.unique_counts = 0,
                                   test.use = "clusWilcox", verbose = TRUE, parallel = FALSE, group_label){

  # Ensure the Seurat object is valid
  if (!inherits(object, "Seurat")) {
    stop("The input object is not a valid Seurat object.")
  }

  # Filter the data
  gene_filter_list = Gene_filter(object, i, ident.1, ident.2,
                                 logfc.threshold, min.pct, min.num_subjects, min.unique_counts, group_label, subject_ID)
  object_filtered = gene_filter_list$object_filtered
  gene_stats = gene_filter_list$gene_stats

  # Create a combined celltype and group label for each cell
  object$celltype.status <- paste(Idents(object), object@meta.data[[group_label]], sep = "_")
  object_filtered$celltype <- Idents(object_filtered)
  Idents(object_filtered) <- "celltype.status"

  # Define the group labels for the given cell type
  clst <- i # given a cell type
  clst_status.1 <- paste(clst, ident.1, sep = "_")
  clst_status.2 <- paste(clst, ident.2, sep = "_")

  # Subset the data to the relevant groups
  subset_data <- subset(object_filtered, idents = c(clst_status.1, clst_status.2))
  counts_df <- as.data.frame(t(as.matrix(subset_data@assays$RNA@counts)))
  counts_df$cell <- rownames(counts_df)

  # Extract metadata
  metadata <- subset_data@meta.data
  metadata$cell <- rownames(metadata)

  # Combine counts and metadata
  counts_df <- counts_df[order(counts_df$cell), ]
  metadata <- metadata[order(metadata$cell), ]
  data_for_test <- cbind(counts_df, metadata)

  # Initialize results list
  results_list <- list()

  # Loop through each gene
  for (gene in rownames(subset_data)) {
    # Prepare data for clusWilcox.test
    test_data <- data.frame(
      expression = as.numeric(data_for_test[[gene]]),
      group = data_for_test$celltype.status,
      subject = data_for_test[[subject_ID]]  # Assuming 'subject' is the clustering variable in metadata
    )

    # Perform the Wilcoxon clustered test
    test_result <- clusWilcox.test(expression ~ group + clusrank::cluster(subject), test_data, method = "ds",
                                   exact = F)

    # Extract p-value and other statistics
    results_list[[gene]] <- data.frame(
      p_val = test_result$p.value,
      gene = gene,
      test.use = test.use
    )
  }

  # Combine results
  results <- do.call(rbind, results_list)

  # Add fold change and other metrics if needed
  results <- results %>%
    merge(gene_stats, by = "gene") %>%
    mutate(cluster = i,
           p_val_adj = p.adjust(p_val, method = "fdr")) %>%
    arrange(p_val_adj)

  return(results)
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables(c("p_val", "p_val_adj"))

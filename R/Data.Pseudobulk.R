#' Data Preparation for Pseudo-bulk Analysis
#'
#' This function processes a Seurat object to prepare data for differential expression analysis using the Pseudo-bulk method.
#' It filters the data and handles bulk calculation methods (mean or sum).
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
#' @param verbose Whether to show progress (default: TRUE).
#' @param group_label Label of the group used for identifying cells.
#' @param bulk_fun Function to use for bulk calculation ("mean" or "sum").
#' @param fixed_effects Fixed effects (other than "group") in the model, e.g., c("batch", "year").
#'
#' @return A list containing:
#' \item{cellbulk_degs}{A matrix of bulk-calculated differential expression data.}
#' \item{gene_stats}{A data frame containing statistics of the filtered genes.}
#' \item{metadata_subset_distinct}{A data frame with distinct metadata for subjects.}
#'
#' @export


Data.Pseudobulk <- function(i, object, ident.1 = "disease", ident.2 = "control", subject_ID,
                            logfc.threshold = 0, min.pct = 0, min.num_subjects = 0, min.unique_counts = 0,
                            verbose = TRUE, group_label, bulk_fun, fixed_effects){

  # Ensure the Seurat object is valid
  if (!inherits(object, "Seurat")) {
    stop("The input object is not a valid Seurat object.")
  }

  # Filter the data
  gene_filter_list = Gene_filter(object, i, ident.1, ident.2,
                                 logfc.threshold, min.pct, min.num_subjects, min.unique_counts, group_label, subject_ID)
  object_filtered = gene_filter_list$object_filtered
  gene_stats = gene_filter_list$gene_stats

  # Add cell type status to the filtered object
  object_filtered$celltype.status <- paste(Idents(object), object@meta.data[[group_label]], sep = "_")

  # Create a metadata dataframe
  meta <- data.frame(wellKey = colnames(object),
                     DonorID = paste(object@meta.data[[group_label]], object@meta.data[[subject_ID]], sep = "_"),
                     Status = object@meta.data[[group_label]])

  # Calculate average expression for each subject
  average_expression <- AverageExpression(object_filtered, group.by = subject_ID)
  average_expression <- average_expression$RNA

  # Filter metadata for distinct subjects
  subject_ids <- unique(object_filtered[[subject_ID]])[[subject_ID]]
  meta_data <- object_filtered@meta.data[object_filtered@meta.data[[subject_ID]] %in% subject_ids, ]
  metadata_subset <- meta_data[, c(subject_ID, fixed_effects), drop = FALSE]

  # Convert subject_ID to a symbol
  subject_ID_sym <- rlang::sym(subject_ID)

  metadata_subset_distinct <- metadata_subset %>%
    dplyr::group_by(!!subject_ID_sym) %>%
    dplyr::summarize(across(everything(), ~ if(all(is.na(.))) NA else first(na.omit(.)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(all_of(subject_ID), as.character)) %>%
    tidyr::drop_na()

  # Calculate cell bulk differential expression data based on bulk_fun
  if (bulk_fun == "mean"){
    cellbulk_degs <- average_expression
    cellbulk_degs <- t(as.matrix(cellbulk_degs))
  }

  if (bulk_fun == "sum"){
    # Compute the number of cells per subject
    num_cells_per_subject <- table(object_filtered@meta.data[[subject_ID]])
    subject_counts_df <- as.data.frame(num_cells_per_subject)
    subject_counts_numeric <- as.numeric(subject_counts_df$Freq)

    subject_ids <- colnames(average_expression)
    num_cells_per_subject <- num_cells_per_subject[subject_ids]

    # Convert the average expression to a matrix if it's not already
    average_expression_matrix <- as.matrix(average_expression)

    # Multiply the average expression by the number of cells to get the sum
    cellbulk_degs <- average_expression_matrix * subject_counts_numeric
    cellbulk_degs <- t(cellbulk_degs)
  }

  # Filter cellbulk_degs to match metadata_subset_distinct
  cellbulk_degs <- cellbulk_degs[rownames(cellbulk_degs) %in% metadata_subset_distinct[[subject_ID]], ]

  # Extract DonorID and map old names to new names
  DonorID = meta$DonorID

  # Function to extract the numeric part from the new names vector
  extract_suffix <- function(name) {
    sub(".*_", "", name)
  }

  new_suffixes <- sapply(DonorID, extract_suffix)
  names(DonorID) <- new_suffixes
  matched_names <- rownames(cellbulk_degs) %in% names(DonorID)
  rownames(cellbulk_degs)[matched_names] <- DonorID[rownames(cellbulk_degs)[matched_names]]

  return(list(cellbulk_degs = cellbulk_degs,
              gene_stats = gene_stats,
              metadata_subset_distinct = metadata_subset_distinct))
}

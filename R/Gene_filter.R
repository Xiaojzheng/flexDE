#' Gene Filter Function
#'
#' This function filters genes in a Seurat object based on various criteria such as log fold change, detection percentage, minimum number of subjects
#' and unique counts per subject. It returns a filtered Seurat object and a data frame with gene statistics.
#'
#' @param object A Seurat object containing the RNA data to be filtered.
#' @param i The cell type to filter for.
#' @param ident.1 The first group identifier for comparison (disease).
#' @param ident.2 The second group identifier for comparison (control).
#' @param logfc.threshold The threshold for the log fold change to filter genes.
#' @param min.pct The minimum detection percentage required to keep a gene.
#' @param min.num_subjects The minimum number of subjects required to keep a gene.
#' @param min.unique_counts The minimum number of unique counts required to keep a gene.
#' @param group_label The label of the group (disease vs. control or treatment vs. control).
#' @param subject_ID The subject ID used for identifying subjects.
#'
#' @return A list containing:
#' \item{object_filtered}{A Seurat object with the filtered genes.}
#' \item{gene_stats}{A data frame containing statistics of the filtered genes.}
#'
#' @export

Gene_filter <- function(object, i, ident.1, ident.2,
                        logfc.threshold, min.pct, min.num_subjects, min.unique_counts,
                        group_label, subject_ID){

  # Filter the data for cell type 'i'
  cells_to_keep <- Idents(object) == i & object@meta.data[[group_label]] %in% c(ident.1, ident.2)
  sub_grp <- subset(object, cells = Cells(object)[cells_to_keep])

  Idents(sub_grp) <- group_label

  group1_counts <- sub_grp@assays$RNA@counts[, WhichCells(sub_grp, idents = ident.1)]
  group2_counts <- sub_grp@assays$RNA@counts[, WhichCells(sub_grp, idents = ident.2)]

  logFC <- log(x = rowMeans(x = group1_counts) + 0.5, base = 2) - log(x = rowMeans(x = group2_counts) + 0.5, base = 2)
  logFC_sample = logFC; FC_sample = 2^logFC_sample

  # Calculate the detection percentage for each group
  detected_pct_group1 <- rowMeans(group1_counts > 0)
  detected_pct_group2 <- rowMeans(group2_counts > 0)

  # Initialize vectors to store subject counts
  counts = sub_grp@assays$RNA@counts
  num_subjects_total <- numeric(nrow(counts))
  num_subjects_group1 <- numeric(nrow(counts))
  num_subjects_group2 <- numeric(nrow(counts))
  unique_counts_group1 <- numeric(nrow(counts))
  unique_counts_group2 <- numeric(nrow(counts))
  prop_zeros_group1 <- numeric(nrow(counts))
  prop_zeros_group2 <- numeric(nrow(counts))

  # Get subject IDs
  subject_ids <- sub_grp@meta.data[[subject_ID]]

  # Identify cells for each group using logical vectors
  group1_cells <- Idents(sub_grp) == ident.1
  group2_cells <- Idents(sub_grp) == ident.2

  # Ensure subject IDs are correctly aligned
  subject_ids_group1 <- subject_ids[group1_cells]
  subject_ids_group2 <- subject_ids[group2_cells]

  # Calculate mean and standard deviation for each gene by subject and by group
  for (gene_index in 1:nrow(counts)) {
    # Get expression values for the current gene
    gene_expression <- counts[gene_index, ]

    # Calculate mean and standard deviation for group 1
    group1_cells <- which(Idents(sub_grp) == ident.1)
    group1_subjects <- subject_ids[group1_cells]
    group1_expression <- gene_expression[group1_cells]
    unique_counts_group1[gene_index] <- tapply(group1_expression, group1_subjects, function(x) length(unique(x))) %>%
      mean(na.rm = TRUE)
    prop_zeros_group1[gene_index] <- tapply(group1_expression, group1_subjects, function(x) mean(x == 0)) %>%
      mean(na.rm = TRUE)

    # Calculate mean and standard deviation for group 2
    group2_cells <- which(Idents(sub_grp)== ident.2)
    group2_subjects <- subject_ids[group2_cells]
    group2_expression <- gene_expression[group2_cells]
    unique_counts_group2[gene_index] <- tapply(group2_expression, group2_subjects, function(x) length(unique(x))) %>%
      mean(na.rm = TRUE)
    prop_zeros_group2[gene_index] <- tapply(group2_expression, group2_subjects, function(x) mean(x == 0)) %>%
      mean(na.rm = TRUE)
  }

  # Loop through each gene
  for (gene_index in 1:nrow(counts)) {
    # Get expression values for the current gene
    gene_expression <- counts[gene_index, ]

    # Identify non-zero expression cells for the current gene
    non_zero_cells <- which(gene_expression > 0)

    # Get the subjects associated with these cells
    subjects <- subject_ids[non_zero_cells]

    # Count the number of unique subjects
    num_subjects_total[gene_index] <- length(unique(subjects))

    # Count the number of unique control and treatment subjects
    unique_subjects <- unique(subjects)
    num_subjects_group1[gene_index] <- sum(unique_subjects %in% subject_ids_group1)
    num_subjects_group2[gene_index] <- sum(unique_subjects %in% subject_ids_group2)
  }

  # Combine into a data frame
  gene_stats <- data.frame(
    gene = rownames(sub_grp),
    logFC_sample = logFC_sample,
    FC_sample = ifelse(FC_sample<1, -1/FC_sample, FC_sample),
    detected_pct_group1 = detected_pct_group1,
    detected_pct_group2 = detected_pct_group2,
    unique_counts_group1 = unique_counts_group1,
    unique_counts_group2 = unique_counts_group2,
    num_subjects_group1 = num_subjects_group1,
    num_subjects_group2 = num_subjects_group2,
    prop_zeros_group1 = prop_zeros_group1,
    prop_zeros_group2 = prop_zeros_group2
  )

  # Filter genes based on logfc.threshold and min.pct
  filtered_genes <- gene_stats %>%
    dplyr::filter(abs(logFC) >= logfc.threshold) %>%
    dplyr::filter(detected_pct_group1 >= min.pct | detected_pct_group2 >= min.pct) %>%
    dplyr::filter(unique_counts_group1 >= min.unique_counts | unique_counts_group2 >= min.unique_counts) %>%
    dplyr::filter(num_subjects_group1 >= min.num_subjects | num_subjects_group2 >= min.num_subjects)

  # Check if filtered_genes is a data frame and contains the 'gene' column
  if (!is.data.frame(filtered_genes) || !("gene" %in% colnames(filtered_genes))) {
    stop("Filtered genes is not a data frame or does not contain the 'gene' column")
  }

  # Subset the Seurat object to keep only the filtered genes
  object_filtered <- subset(object, features = filtered_genes$gene)

  gene_stats <- gene_stats %>%
    dplyr::filter(gene %in% filtered_genes$gene)

  return(list(object_filtered = object_filtered,
              gene_stats = gene_stats))
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables("gene")


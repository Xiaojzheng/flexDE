#' Differential Expression Analysis Using Nebula
#'
#' This function performs differential expression analysis using the Nebula method to model both within- and between-subject overdispersion.
#' It handles various fixed effects.
#'
#' @param i Indicator of cell type.
#' @param object A Seurat object containing the RNA data.
#' @param ident.1 The first group identifier for comparison (default: "disease").
#' @param ident.2 The second group identifier for comparison (default: "control").
#' @param logfc.threshold The threshold for the log fold change to filter genes.
#' @param min.pct The minimum detection percentage required to keep a gene.
#' @param min.num_subjects The minimum number of subjects required to keep a gene.
#' @param min.unique_counts The minimum number of unique counts required to keep a gene.
#' @param test.use Statistical method to use (default: "NEBULA").
#' @param verbose Whether to show progress (default: TRUE).
#' @param group_label The label of the group (disease vs. control or treatment vs. control).
#' @param subject_ID Subject ID used for identifying subjects.
#' @param fixed_effects Fixed effects (other than "group") in the model, e.g., c("batch", "sex").
#'
#' @return A data frame containing the results of the differential expression analysis.
#'
#' @export
#'

DEs.Nebula <- function(i, object, ident.1 = "disease", ident.2 = "control",
                       logfc.threshold = 0, min.pct = 0, min.num_subjects = 0, min.unique_counts = 0,
                       test.use = "NEBULA", verbose = TRUE, group_label, subject_ID, fixed_effects){

  # Mode function to find the most frequent value
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  # Filter genes using the Gene_filter function
  gene_filter_list = Gene_filter(object, i, ident.1, ident.2,
                                 logfc.threshold, min.pct, min.num_subjects, min.unique_counts, group_label, subject_ID)
  object_filtered = gene_filter_list$object_filtered
  gene_stats = gene_filter_list$gene_stats

  # Extract count matrix and metadata
  counts <- object_filtered[["RNA"]]@counts ### count matrix
  orig.ident <- as.character(object_filtered@meta.data[[subject_ID]]) ### subject ID
  metadata_df <- data.frame(group = object_filtered@meta.data[[group_label]]) ### CKD vs. Control

  # Include fixed effects in metadata
  for (effect in fixed_effects) {
    metadata_df[[effect]] <- object_filtered@meta.data[[effect]]
  }

  # Impute missing values within subjects
  combined_df <- cbind(metadata_df, orig.ident)
  combined_df <- combined_df %>%
    group_by(orig.ident) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), Mode(.[!is.na(.)]), .)))

  # Filter columns with at least two distinct values and no all NA values
  combined_df <- combined_df %>%
    select_if(~ n_distinct(na.omit(.)) > 1 && !all(is.na(.)))

  # Print the number of unique subjects and cells
  message(paste("Number of unique subjects:", n_distinct(combined_df$orig.ident)))
  message(paste("Number of cells:", nrow(combined_df)))

  # Print message about NA values for fixed effects
  for (effect in fixed_effects) {
    num_na <- sum(is.na(combined_df[[effect]]))
    total_obs <- nrow(combined_df)
    message(paste("Fixed effect", effect, "has", num_na, "NA values out of", total_obs, "observations."))
  }

  # Keep only complete cases
  complete_cases <- complete.cases(combined_df)
  combined_df <- combined_df[complete_cases, ]
  sorted_combined_df <- combined_df[order(combined_df$orig.ident), ]

  # Print the number of unique subjects and cells again
  message("After removing the incomplete cases...")
  message(paste("Number of unique subjects:", n_distinct(sorted_combined_df$orig.ident)))
  message(paste("Number of cells:", nrow(sorted_combined_df)))

  # Extract sorted metadata and orig.ident
  metadata_df_sorted <- sorted_combined_df[, -ncol(sorted_combined_df)]
  orig.ident_sorted <- sorted_combined_df$orig.ident

  # Create design matrix
  df <- model.matrix(~ ., data = metadata_df_sorted)

  # Filter and sort the count matrix to match sorted_combined_df
  filtered_counts <- counts[, complete_cases]
  filtered_counts <- filtered_counts[, order(orig.ident[complete_cases])]

  ### Fit the model
  # @model: Default --> 'NBGMM' is for fitting a negative binomial gamma mixed model.
  results <- nebula(count = filtered_counts, id = orig.ident_sorted, pred = df)

  # Extract coefficients and p-values for group
  results_coef <- results$summary
  group_coef <- grep("logFC_group", colnames(results_coef))
  group_coefficients <- results_coef[, group_coef]

  group_pval <- grep("p_group", colnames(results_coef))
  group_pval <- results_coef[, group_pval]

  # Prepare the results data frame
  DE_genes <- results$summary %>%
    mutate(FC = exp(group_coefficients), FC = ifelse(FC<1, -1/FC, FC),
           p_val = group_pval, p_val_adj = p.adjust(p_val, method = "fdr"),
           cluster = i, overdispersion_sub = results$overdispersion[,1],
           overdispersion_cell = results$overdispersion[,2], convergence = results$convergence,
           convergence = ifelse(convergence == 1, "Yes", "No"))

  DE_genes <- DE_genes %>%
    dplyr::select(p_val, FC, overdispersion_sub, overdispersion_cell, convergence, p_val_adj, cluster,
                  gene) %>%
    arrange(p_val)

  DE_genes <- DE_genes %>%
    merge(gene_stats, by = "gene")


  return(DE_genes)
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables(c("FC", "p_val", "convergence", "overdispersion_sub", "overdispersion_cell",
                         "p_val_adj", "gene"))

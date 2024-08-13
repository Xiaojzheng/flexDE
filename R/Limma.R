#' Differential Expression Analysis Using Limma
#'
#' This function performs differential expression analysis using the Limma method on bulk-calculated expression data.
#' It supports normalization and handles both mean and sum aggregation methods.
#'
#' @param cellbulk_degs A matrix of bulk-calculated differential expression data.
#' @param gene_stats_degs A data frame containing statistics of the filtered genes.
#' @param metadata_subset_distinct A data frame with distinct metadata for subjects.
#' @param bulk_fun Function used for bulk calculation ("mean" or "sum").
#' @param subject_ID Subject ID used for identifying subjects.
#' @param normalization Whether to normalize the data (TRUE or FALSE).
#'
#' @return A data frame containing the results of the differential expression analysis.
#'
#' @export
#'

DEs.Limma <- function(cellbulk_degs, gene_stats_degs, metadata_subset_distinct, bulk_fun, subject_ID, normalization){

  # Check if bulk_fun is "mean" and round up values if true
  if (bulk_fun == "mean"){
    cellbulk_degs <- ceiling(cellbulk_degs)
  }

  subject_ID = subject_ID
  gene_stats = gene_stats_degs
  metadata_subset_distinct = metadata_subset_distinct

  ## Input is gene X subject
  # create targets matrix
  x <- t(cellbulk_degs)

  # Create a data frame (targets) with group sample and subject ID
  targets = data.frame(group_sample = colnames(x)) %>%
    mutate(group = gsub("_.+", "", group_sample)) %>%
    mutate(subject_ID_value = sub(".*_(\\d+)$", "\\1", group_sample)) %>%
    rename(!!subject_ID := subject_ID_value)

  # Merge with metadata based on subject_ID
  targets <- targets %>%
    left_join(metadata_subset_distinct, by = subject_ID)

  # Create the design matrix for the model
  design = model.matrix(~ ., data = targets %>% select(-group_sample, -all_of(subject_ID)))

  if (normalization == T){
    y = DGEList(counts = x, group = targets$group) %>%
      calcNormFactors(method = 'TMM')
    x = voom(y, design)
  }

  if (normalization == F){
    y = DGEList(counts = x, group = targets$group)
    x = voom(y, design, normalize.method = "none")
  }

  # Set trend_bool to TRUE to account for the mean-variance trend
  trend_bool <- TRUE

  # Fit the linear model and apply empirical Bayes smoothing
  fit = lmFit(x, design) %>%
    eBayes(trend = trend_bool, robust = trend_bool)

  # Format the results, extracting group coefficients
  res = fit %>%
    topTable(number = Inf, coef = grep("^group", colnames(fit$coefficients), value = TRUE)) %>%
    rownames_to_column('gene') %>%
    # flag metrics in results
    mutate(
      de_family = 'pseudobulk',
      de_method = 'limma',
      de_type = "voom",
      aggregation = bulk_fun)

  DE_genes <- res %>%
    mutate(avg_logFC = logFC, FC = exp(logFC), FC = ifelse(FC<1, -1/FC, FC), p_val = P.Value, p_val_adj = adj.P.Val) %>%
    dplyr::select(gene, p_val, FC, avg_logFC, p_val_adj, de_family, aggregation, de_method, de_type) %>%
    arrange(p_val)

  #  Merge the DE_genes with gene_stats to include additional gene information
  DE_genes <- DE_genes %>%
    merge(gene_stats, by = "gene")

  return(DE_genes)
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables(c("group_sample", "subject_ID_value", "logFC",
                         "FC", "P.Value", "adj.P.Val", "gene", "p_val",
                         "avg_logFC", "p_val_adj", "de_family", "aggregation",
                         "de_method", "de_type"))

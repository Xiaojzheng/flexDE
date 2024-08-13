#' Differential Expression Analysis Using edgeR
#'
#' This function performs differential expression analysis using the edgeR package on bulk-calculated expression data.
#' It supports normalization and handles both LRT and Quasi-Likelihood methods.
#'
#' @param cellbulk_degs A matrix of bulk-calculated differential expression data.
#' @param gene_stats_degs A data frame containing statistics of the filtered genes.
#' @param metadata_subset_distinct A data frame with distinct metadata for subjects.
#' @param subject_ID Subject ID used for identifying subjects.
#' @param bulk_fun Function used for bulk calculation ("mean" or "sum").
#' @param de_type Type of differential expression test to use ("LRT" or "Quasi-Likelihood").
#' @param normalization Whether to normalize the data (TRUE or FALSE).
#'
#' @return A data frame containing the results of the differential expression analysis.
#'
#' @export

DEs.edgeR <- function(cellbulk_degs, gene_stats_degs, metadata_subset_distinct, subject_ID, bulk_fun, de_type = "LRT", normalization){

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
    dplyr::mutate(group = gsub("_.+", "", group_sample)) %>%
    dplyr::mutate(subject_ID_value = sub(".*_(\\d+)$", "\\1", group_sample)) %>%
    dplyr::rename(!!subject_ID := subject_ID_value)

  # Merge with metadata based on subject_ID
  targets <- targets %>%
    dplyr::left_join(metadata_subset_distinct, by = subject_ID)

  # Create the design matrix for the model
  design = model.matrix(~ ., data = targets %>% select(-group_sample, -all_of(subject_ID)))

  if (normalization == T){
    y = DGEList(counts = x, group = targets$group) %>%
      calcNormFactors(method = 'TMM') %>%
      estimateDisp(design)
  }

  if (normalization == F){
    y = DGEList(counts = x, group = targets$group) %>%
      estimateDisp(design)
  }

  ### LRT
  if (de_type == "LRT"){
    fit = glmFit(y, design = design)
    test = glmLRT(fit, coef = grep("^group", colnames(design), value = TRUE))
  }

  ### Quasi-Likelihood
  if (de_type == "Quasi-Likelihood"){
    fit = glmQLFit(y, design = design)
    test = glmQLFTest(fit, coef = grep("^group", colnames(design), value = TRUE))
  }

  ### Results
  res = topTags(test, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    # flag metrics in results
    mutate(de_family = 'pseudobulk',
           de_method = test.use,
           de_type = de_type,
           aggregation = bulk_fun,
           padj = p.adjust(PValue, method = 'BH'))

  DE_genes <- res %>%
    mutate(avg_logFC = logFC, FC = exp(logFC), FC = ifelse(FC<1, -1/FC, FC), p_val = PValue, p_val_adj = padj) %>%
    dplyr::select(gene, p_val, FC, avg_logFC, p_val_adj, de_family, aggregation, de_method, de_type) %>%
    arrange(p_val)

  #  Merge the DE_genes with gene_stats to include additional gene information
  DE_genes <- DE_genes %>%
    merge(gene_stats, by = "gene")

  return(DE_genes)
}

#' @name globalVariables
utils::globalVariables(c("group_sample", "subject_ID_value", "logFC",
                         "FC", "PValue", "padj", "gene", "p_val",
                         "avg_logFC", "p_val_adj", "de_family", "aggregation",
                         "de_method", "test.use"))


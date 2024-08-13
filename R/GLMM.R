#' Differential Expression Analysis using Generalized Linear Mixed Models (GLMM)
#'
#' This function performs differential expression analysis on cell-level data using generalized linear mixed models (GLMM).
#' It supports both negative binomial and tweedie distributions, and offers Wald test or likelihood ratio test (LRT) options
#' for statistical testing. The function can also include zero-inflation models with various configurations.
#'
#' @param object A Seurat object containing the RNA data to be analyzed.
#' @param i Indicator of the cell type to filter for.
#' @param ident.1 The first group identifier for comparison (default is "disease").
#' @param ident.2 The second group identifier for comparison (default is "control").
#' @param logfc.threshold The threshold for the log fold change to filter genes.
#' @param min.pct The minimum detection percentage required to keep a gene.
#' @param min.num_subjects The minimum number of subjects required to keep a gene.
#' @param min.unique_counts The minimum number of unique counts required to keep a gene.
#' @param group_label The label of the group (disease vs. control or treatment vs. control).
#' @param subject_ID The subject ID used for identifying subjects.
#' @param n_threads Number of threads to use for parallel processing (default is the number of cores available).
#' @param de_method The method for differential expression analysis ("negbinom" for negative binomial or "tweedie" for tweedie distribution, default is "negbinom").
#' @param de_type The type of statistical test to use ("Wald" for Wald test or "LRT" for likelihood ratio test, default is "Wald").
#' @param fixed_effects Fixed effects (other than "group") to include in the model, e.g., c("batch", "sex").
#' @param zero_intercept Logical indicating whether to include an intercept in the zero-inflation model (default is FALSE).
#' @param zero_group Logical indicating whether to include "group" as a fixed effect in the zero-inflation model (default is FALSE).
#' @param zero_subject Logical indicating whether to include "subject" as a random effect in the zero-inflation model (default is FALSE).
#' @param zero_fixed Fixed effects (other than "group") to include in the zero-inflation part of the model.
#' @param show_progress Logical indicating whether to show progress of the parallel processing (default is TRUE).
#'
#' @return A data frame containing the differential expression analysis results with columns:
#' \item{gene}{Gene identifier.}
#' \item{p_val}{P-value for the differential expression test.}
#' \item{test_statistic}{Test statistic for the differential expression test.}
#' \item{odds_ratio}{Odds ratio for the differential expression.}
#' \item{zi_odds_ratio}{Zero-inflation model odds ratio (if applicable).}
#' \item{zi_p_val}{P-value for the zero-inflation model (if applicable).}
#' \item{converged}{Logical indicating whether the model converged.}
#' \item{de_method}{Method used for differential expression analysis.}
#' \item{de_type}{Type of statistical test used.}
#' \item{p_val_adj}{Adjusted p-value using the Benjamini-Hochberg method.}
#' \item{logFC_sample}{Log fold change for the sample.}
#' \item{FC_sample}{Fold change for the sample.}
#' \item{detected_pct_group1}{Detection percentage in group 1.}
#' \item{detected_pct_group2}{Detection percentage in group 2.}
#' \item{unique_counts_group1}{Unique counts in group 1.}
#' \item{unique_counts_group2}{Unique counts in group 2.}
#' \item{num_subjects_group1}{Number of subjects in group 1.}
#' \item{num_subjects_group2}{Number of subjects in group 2.}
#' \item{prop_zeros_group1}{Proportion of zeros in group 1.}
#' \item{prop_zeros_group2}{Proportion of zeros in group 2.}
#'
#' @export
#'

DEs.GLMM <- function(object, i, ident.1 = "disease", ident.2 = "nondisease",
                     logfc.threshold, min.pct, min.num_subjects, min.unique_counts, group_label, subject_ID,
                     n_threads = cores, de_method = "negbinom",
                     de_type = 'Wald', fixed_effects, zero_intercept = FALSE,
                     zero_group = FALSE, zero_subject = FALSE, zero_fixed = NULL, show_progress = TRUE){

  # Ensure the Seurat object is valid
  if (!inherits(object, "Seurat")) {
    stop("The input object is not a valid Seurat object.")
  }

  # Filter the data
  gene_filter_list = Gene_filter(object, i, ident.1, ident.2,
                                 logfc.threshold, min.pct, min.num_subjects, min.unique_counts, group_label, subject_ID)
  object_filtered = gene_filter_list$object_filtered
  gene_stats = gene_filter_list$gene_stats

  if (de_method == "negbinom") {
    counts <- as.matrix(object_filtered[["RNA"]]@counts)
  } else if (de_method == "tweedie") {
    counts <- as.matrix(object_filtered[["RNA"]]@data)
  }

  # Prepare the data frame for glmmTMB
  metadata <- object_filtered@meta.data
  metadata$cell <- rownames(metadata) ### Cell name
  metadata$group <- object_filtered@meta.data[[group_label]] ### CKD vs. Control
  metadata$subject <-as.factor(as.character(object_filtered@meta.data[[subject_ID]])) ### subject ID

  # Ensure the fixed effects variables are present in metadata
  required_vars <- unique(c("cell", "group", "subject", fixed_effects, zero_fixed))
  missing_vars <- required_vars[!required_vars %in% colnames(metadata)]
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the metadata:", paste(missing_vars, collapse = ", ")))
  }

  # Subset metadata to include only necessary variables
  metadata_subset <- metadata[, required_vars]

  counts <- t(counts[,rownames(metadata_subset)])
  counts_df <- cbind(metadata_subset, counts)

  #-------------------------------------------------------------------------------------
  # The model formulation consists of two parts: 1. NB counts 2. zero-inflation part
  ## Part 1: NB counts

  # Create the formula string dynamically
  fixed_effects_base <- "group"
  random_effects <- "(1 | subject)"
  all_fixed_effects <- paste(c(fixed_effects_base, fixed_effects), collapse = " + ")
  formula_string <- paste(" ~", all_fixed_effects, "+", random_effects)

  ## Part 2: Zero-inflation model

  # Create the formula string dynamically
  zi_components <- character(0)

  # Add intercept if zero_intercept is TRUE
  if (zero_intercept) {
    zi_components <- c(zi_components, "1")
  }
  # Add 'group' as a fixed effect if zero_group is TRUE
  if (zero_group) {
    zi_components <- c(zi_components, "group")
  }
  # Add other fixed effects
  if (!is.null(zero_fixed)) {
    zi_components <- c(zi_components, zero_fixed)
  }
  # Add 'subject' as a random effect if zero_subject is TRUE
  if (zero_subject) {
    zi_components <- c(zi_components, "(1 | subject)")
  }

  # Combine components into a single formula string
  zi_formula_string <- paste("~", paste(zi_components, collapse = " + "))

  if (!zero_group & !zero_subject & is.null(zero_fixed)){
    zi_formula_string <- "~0"
  }

  # Convert to formula object
  zi_final_formula <- as.formula(zi_formula_string)

  #--------------------------------------------------------------------------------
  # Set the family based on the de_method
  if (de_method == "negbinom") {
    family <- nbinom2
  } else if (de_method == "tweedie") {
    family <- glmmTMB::tweedie(link = "log")
  }

  ### Three model specifications to be considered
  #1. Distribution: a) Negative Binomial b) Tweedie
  #2. Zero_inflation: a) Yes b) No
  #3. Statistical testing: a) Wald test b) Likelihood Ratio test

  # Define the function for parallel processing
  de_analysis <- function(gene) {

    tryCatch({

      final_formula <- as.formula(paste(paste0("`", gene, "`"), formula_string))

      if (de_type == "Wald"){

        # Fit the GLMM

        full_model <- glmmTMB(final_formula, data = counts_df, family = family,
                              ziformula = zi_final_formula, REML = FALSE)

        # Check if the model converged
        converged <- full_model$sdr$pdHess

        # Extract the p-value and test statistic for the Wald test
        tab <- coef(summary(full_model))[[1]]
        coef <- paste0("group", levels(counts_df$group)[2])
        p_val <- tab[coef, "Pr(>|z|)"]
        test_statistic <- tab[coef, "z value"]

      } else if (de_type == 'LRT'){

        # Fit the full model with the treatment group
        full_model <- glmmTMB(final_formula, data = counts_df, family = family,
                              ziformula = zi_final_formula, REML = FALSE)

        # Define the null model without the treatment group

        # Remove the first term "group"
        formula_string_null <- sub("group \\+ ", "", formula_string)
        formula_string_null <- sub("group", "", formula_string_null) # In case "group" is the only term

        zi_formula_string_null <- sub("group \\+ ", "", zi_formula_string)
        zi_formula_string_null <- sub("\\+ group", "", zi_formula_string_null) # In case "group" is the only term

        final_formula_null <-as.formula(paste(paste0("`", gene, "`"), formula_string_null))
        zi_final_formula_null <- as.formula(zi_formula_string_null)

        null_model <- glmmTMB(final_formula_null, data = counts_df, family = family,
                              ziformula = zi_final_formula_null, REML = FALSE)

        # Check if both models converged
        converged <- full_model$sdr$pdHess && null_model$sdr$pdHess

        # Perform the LRT
        lrt <- anova(null_model, full_model)
        p_val <- lrt$`Pr(>Chisq)`[2]
        test_statistic <- lrt$Chisq[2]
      }

      # Calculate log_odds_ratio
      coef_summary = as.data.frame(summary(full_model)$coefficients$cond)
      group_term <- grep("^group", rownames(coef_summary), value = TRUE)
      log_odds_ratio <- coef_summary[group_term, "Estimate"]

      # Calculate the odds ratio (treatment compared to control)
      odds_ratio <- exp(-log_odds_ratio)
      zi_odds_ratio = NA
      zi_p_val = NA

      if(zero_group){
        # Extract coefficients for the zero-inflation part
        zi_coef_summary <- as.data.frame(summary(full_model)$coefficients$zi)
        zi_group_term <- grep("^group", rownames(zi_coef_summary), value = TRUE)
        zi_log_odds_ratio <- zi_coef_summary[zi_group_term, "Estimate"]
        zi_odds_ratio <- exp(-zi_log_odds_ratio)
        zi_p_val <- zi_coef_summary[zi_group_term, "Pr(>|z|)"]
      }

      return(data.frame(
        gene = gene,
        p_val = p_val,
        test_statistic = test_statistic,
        odds_ratio = odds_ratio,
        zi_odds_ratio = zi_odds_ratio,
        zi_p_val = zi_p_val,
        converged = converged,
        de_method = de_method,
        de_type = de_type
      ))
    }, error = function(e) {
      return(data.frame(
        gene = gene,
        p_val = NA,
        test_statistic = NA,
        odds_ratio = NA,
        zi_odds_ratio = NA,
        zi_p_val = NA,
        converged = FALSE,
        de_method = NA,
        de_type = NA
      ))
    })
  }

  # Apply the DE analysis function in parallel
  genes <- unique(colnames(counts))

  if (show_progress) {
    results <- pbmclapply(genes, de_analysis, mc.cores = n_threads)
  } else {
    results <- mclapply(genes, de_analysis, mc.cores = n_threads)
  }

  # Combine the results into a single data frame
  results_df <- do.call(rbind, results)
  results_df$p_val_adj <- p.adjust(results_df$p_val, method = 'BH')
  results_df <- results_df %>%
    merge(gene_stats, by = "gene")

  return(results_df)
}

#' @name globalVariables
#' @keywords internal
utils::globalVariables("cores")

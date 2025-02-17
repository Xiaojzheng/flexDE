% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NEBULA.R
\name{DEs.Nebula}
\alias{DEs.Nebula}
\title{Differential Expression Analysis Using Nebula}
\usage{
DEs.Nebula(
  i,
  object,
  ident.1 = "disease",
  ident.2 = "control",
  logfc.threshold = 0,
  min.pct = 0,
  min.num_subjects = 0,
  min.unique_counts = 0,
  test.use = "NEBULA",
  verbose = TRUE,
  group_label,
  subject_ID,
  fixed_effects
)
}
\arguments{
\item{i}{Indicator of cell type.}

\item{object}{A Seurat object containing the RNA data.}

\item{ident.1}{The first group identifier for comparison (default: "disease").}

\item{ident.2}{The second group identifier for comparison (default: "control").}

\item{logfc.threshold}{The threshold for the log fold change to filter genes.}

\item{min.pct}{The minimum detection percentage required to keep a gene.}

\item{min.num_subjects}{The minimum number of subjects required to keep a gene.}

\item{min.unique_counts}{The minimum number of unique counts required to keep a gene.}

\item{test.use}{Statistical method to use (default: "NEBULA").}

\item{verbose}{Whether to show progress (default: TRUE).}

\item{group_label}{The label of the group (disease vs. control or treatment vs. control).}

\item{subject_ID}{Subject ID used for identifying subjects.}

\item{fixed_effects}{Fixed effects (other than "group") in the model, e.g., c("batch", "sex").}
}
\value{
A data frame containing the results of the differential expression analysis.
}
\description{
This function performs differential expression analysis using the Nebula method to model both within- and between-subject overdispersion.
It handles various fixed effects.
}

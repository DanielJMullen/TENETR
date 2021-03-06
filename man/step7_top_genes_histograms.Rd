% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step7_top_genes_histograms.R
\name{step7_top_genes_histograms}
\alias{step7_top_genes_histograms}
\title{step7_top_genes_histograms}
\usage{
step7_top_genes_histograms(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis
)
}
\arguments{
\item{TENET_directory}{Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_histograms' containing the results of this function.}

\item{hypermeth_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to create histograms of genes/TFs linked to hypermeth probes with G+ links.}

\item{hypermeth_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create histograms of genes/TFs linked to hypermeth probes with G- links.}

\item{hypometh_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to to create histograms of genes/TFs linked to hypometh probes with G+ links.}

\item{hypometh_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create histograms of genes/TFs linked to hypometh probes with G- links.}
}
\value{
Currently returns .pdf files with the histograms showing the number of genes/TFs linked to a given number of enhancer probes of the given analysis type.
}
\description{
This is a step7 function of the TENETR package.
This function takes all genes and only accepted transcription factor
genes from "The Human Transcription Factors" by Lambert et al. 2018.
and generates histograms displaying the number of all genes and TF-only genes
with a links to a given number of enhancer DNA methylation probes for each of
the four hypo or hypermethylated Gplus or Gminus analysis quadrants,
as selected by the user.
}

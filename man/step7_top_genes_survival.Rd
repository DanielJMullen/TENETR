% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step7_top_genes_survival.R
\name{step7_top_genes_survival}
\alias{step7_top_genes_survival}
\title{step7_top_genes_survival}
\usage{
step7_top_genes_survival(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  visualize_survival_plots_genes,
  visualize_survival_plots_probes,
  high_thresh,
  low_thresh,
  core_count
)
}
\arguments{
\item{TENET_directory}{Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_survival' containing the results of this function.}

\item{hypermeth_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to create survival plots for the top genes/TFs by most hypermeth probes with G+ links, and these linked DNA methyation probes if specified.}

\item{hypermeth_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypermeth probes with G- links, and these linked DNA methyation probes if specified.}

\item{hypometh_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G+ links, as well as their linked DNA methyation probes if specified.}

\item{hypometh_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G- links, as well as their linked DNA methyation probes if specified.}

\item{top_gene_number}{Specify a number to generate survival plots for that many of the top genes/TFs, as well as their linked enhancer probes if specified, based on the most linked enhancer probes.}

\item{visualize_survival_plots_genes}{Set to TRUE/FALSE depending on if you want to create .pdfs displaying the survival results for the genes/TFs of interest.}

\item{visualize_survival_plots_probes}{Set to TRUE/FALSE depending on if you want to create .pdfs displaying the survival results for the probes linked to the genes/TFs of interest.}

\item{high_thresh}{Set a number ranging from 0 to 1, as a threshold for proportion of samples above that number to include in the high expression/methylation group, and should be greater than or equal to low_thresh to prevent samples from appearing in both groups.}

\item{low_thresh}{Set a number ranging from 0 to 1, as a threshold for proportion of samples below that number to include in the low expression/methylation group, and should be less than or equal to high_thresh to prevent samples from appearing in both groups.}

\item{core_count}{Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.}
}
\value{
Returns survival information in the form of .tsv files, as well as .pdfs if selected by the user, showing survival information for the expression of the top gene/TFs, as well as the methylation of the enhancer DNA methylation probes linked to them.
}
\description{
This is a step7 function of the TENETR package.
This function takes the top genes/TFs by number of linked probes identified from
the step6_probe_per_gene_tabulation function up to the number as specified by the user
and generates survival plots and information for the expression level of each gene
as well as the DNA methylation of each enhancer probe linked to them,
using percentile cutoffs as specified by the user.
}

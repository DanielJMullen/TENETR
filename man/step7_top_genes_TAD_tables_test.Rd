% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step7_top_genes_TAD_tables_test.R
\name{step7_top_genes_TAD_tables_test}
\alias{step7_top_genes_TAD_tables_test}
\title{step7_top_genes_TAD_tables_test}
\usage{
step7_top_genes_TAD_tables_test(
  TENET_directory,
  TAD_directory,
  DNA_methylation_manifest = "HM450",
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  core_count
)
}
\arguments{
\item{TENET_directory}{Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_TAD_tables' containing the results of this function.}

\item{TAD_directory}{Set a path to a directory which contains tab-delimited bed-like files, see: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 , with information on the TAD compartments of interest. Multiple files can be included in this subdirectory.}

\item{DNA_methylation_manifest}{Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.}

\item{hypermeth_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypermeth probes with G+ links.}

\item{hypermeth_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypermeth probes with G- links.}

\item{hypometh_Gplus_analysis}{Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypometh probes with G+ links.}

\item{hypometh_Gminus_analysis}{Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypometh probes with G- links.}

\item{top_gene_number}{Specify a number of the top genes/TFs based on the most linked enhancer probes to generate TAD tables for the enhancer probes linked to those genes.}

\item{core_count}{Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.}
}
\value{
Exports .tsv files with rows for each of the probes linked to the top genes/TFs as specified by the user, with the probe locations in the hg38 human genome, which of the top genes/TFs they're linked to, and the number of genes as well as their ensembl IDs and names found in the same TAD as each probe for each of the user input TAD files.
}
\description{
This is a step7 function of the TENETR package.
This function takes the top genes/TFs by number of linked probes identified from
the step6_probe_per_gene_tabulation function up to the number as specified by the user
and generates tables of information for each of the enhancer probes linked
to them for each of the four hypo or hypermethylated Gplus or Gminus analysis quadrants,
as selected by the user. These tables note which of the top genes/TFs each probe is linked to,
as well as the the total number of  genes and their names which happen to lie
within the same topologically-associating domain, or TAD, of each probe in each of the
user-supplied TAD files.
}

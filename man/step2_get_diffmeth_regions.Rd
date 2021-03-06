% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step2_get_diffmeth_regions.R
\name{step2_get_diffmeth_regions}
\alias{step2_get_diffmeth_regions}
\title{step2_get_diffmeth_regions}
\usage{
step2_get_diffmeth_regions(
  TENET_directory,
  DNA_methylation_manifest = "HM450",
  TSS_dist = 1500,
  methcutoff,
  hypomethcutoff,
  unmethcutoff,
  hypermethcutoff,
  minExp,
  use_purity_data = FALSE
)
}
\arguments{
\item{TENET_directory}{Set a path to the TENET directory containing the 'step1' subdirectory and results created by the step1_make_external_datasets function as well as a user-supplied .rda file containing methylation and expression data. This function will also create a new 'step2' subdirectory there containing the results of this function.}

\item{DNA_methylation_manifest}{Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.}

\item{TSS_dist}{Set a positive integer to be the buffer in base pairs from GENCODE v22-annotated transcription start sites for DNA methylation probes to not be considered enhancer probes. Defaults to 1500.}

\item{methcutoff}{Set a number from 0 to 1 to be the beta-value cutoff for methylated probes.}

\item{hypomethcutoff}{Set a number from 0 to 1 to be the beta-value cutoff for hypomethylated probes. Should be set lower than the methcutoff.}

\item{unmethcutoff}{Set a number from 0 to 1 to be the beta-value cutoff for unmethylated probes.}

\item{hypermethcutoff}{Set a number from 0 to 1 to be the beta-value cutoff for hypermethylated probes. Should be set higher than the unmethcutoff.}

\item{minExp}{Set a positive integer to be the minimum number of experimental/tumor samples to be considered for the hypo/hypermethylated groups. Should be less than the total number of experimental/tumor groups.}

\item{use_purity_data}{Set TRUE or FALSE to use purity datasets, as .rda files containing DNA methylation values supplied by the user in a 'purity' subdirectory in the TENET_directory, to select for DNA methylation probes not potentially due to differences in cell type purity. Defaults to FALSE.}
}
\value{
Returns two objects, a .rda file with matrices of methylation and expression data for the four types of identified enhancer DNA methylation probes in control/normal and experimental/tumor samples, the clinical data, as well as the user-set parameters for consistency in downstream analyses. Also outputs a .txt file containing metrics on the number of probes found in different categories.
}
\description{
This is the step2 function of the TENETR package.
This function identifies DNA methylation probes that mark putative enhancer
regions. These are probes that lie within at least one enhancer and open chromatin
dataset from the step1_make_external_datasets function and which are located
a user-specified distance from GENCODE v22 transcript transcription start sites.
After identifying enhancer DNA methylation probes, the function classifies
the probes as methylated, unmethylated, hypermethylated, or hypomethylated
based on their differential methylation between the control/normal and
experimental/tumor samples supplied by the user, defined by the cutoff values
also specified by the user.
}

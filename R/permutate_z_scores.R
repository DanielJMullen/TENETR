#' permutate_z_scores
#' #'
#' This is the step4 function of the TENETR package.
#' This function takes the calculated z-scores for the hyper/hypomethylated
#' Gplus or Gminus probe-gene quadrants, and does a pseudo-permutation,
#' calculating individual permutated p-values for the probe-gene z-score
#' based on their relative rank amongst all probes to a given gene.
#'
#'
#' @param TENET_directory Set a path to the directory that contains step3 results from get_analysis_z_score function. This function will also create a new step4 folder there containing the results.
#' @param hypermeth_Gplus_analysis Set TRUE or FALSE if user wants to permutate on hypermeth_Gplus links. Requires hypermeth_analysis from step3 to have been set to TRUE.
#' @param hypermeth_Gminus_analysis Set TRUE or FALSE if user wants to permutate on hypermeth_Gminus links. Requires hypermeth_analysis from step3 to have been set to TRUE.
#' @param hypometh_Gplus_analysis Set TRUE or FALSE if user wants to permutate on hypometh_Gplus links. Requires hypometh_analysis from step3 to have been set to TRUE.
#' @param hypometh_Gminus_analysis Set TRUE or FALSE if user wants to permutate on hypometh_Gminus links. Requires hypometh_analysis from step3 to have been set to TRUE.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns tab-delimited tab-delimited "sig_link_zscores_perm.txt" files for hypo/hyper Gplus/Gminus probe-gene links, similar to step3, but with the permutated p-value for each link.
#' @export

permutate_z_scores <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  core_count
){

  ## If user has not supplied the final '/' in the TENET directory
  ## add it:
  TENET_directory <- ifelse(
    substring(
      TENET_directory,
      nchar(TENET_directory),
      nchar(TENET_directory)
    ) == '/',
    TENET_directory,
    paste(
      TENET_directory,
      '/',
      sep=''
    )
  )

  ## Create a step3 directory to deposit the output paired score files:
  dir.create(
    paste(
      TENET_directory,
      'step4/',
      sep=''
    )
  )

  ## If hypermeth Gplus analysis is selected, run it first:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Check that the significant zscores from step3 exists and load it:
    ## If not, return an error message:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step3/',
          'hypermeth_results/',
          'hyper_Gplus_sig_link_zscores.txt',
          sep=''
        )
      )
    ){

      ## Load the file and change the colnames:
      hypermeth_Gplus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step3/',
          'hypermeth_results/',
          'hyper_Gplus_sig_link_zscores.txt',
          sep=''
        ),
        sep='\t',
        header= TRUE,
        stringsAsFactors = FALSE
      )

      colnames(hypermeth_Gplus_zscores) <- c(
        'geneID',
        'probeID',
        'Zscore'
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hyper_Gplus_sig_link_zscores.txt in step3 of TENET directory was not found. Please check that the file exists and consider rerunning the step3 get_analysis_z_score function.')
    }

    ## Write the function to perform the Z-score permutation:
    getPosEpval_hyperGplus <- function(
      geneID,
      probe,
      Z_score
    ){

      ## Set the relevant folder to find results in:
      step3_folder <- paste(
        TENET_directory,
        'step3/',
        'hypermeth_results',
        sep=''
      )

      ## Load the individual gene file for the given
      ## gene-probe link:
      gene_file <- read.delim(
        list.files(
          step3_folder,
          pattern= geneID,
          full.names = TRUE
        ),
        header= FALSE,
        sep = '\t',
        stringsAsFactors = FALSE
      )

      ## Remove NA and infinite values:
      gene_file <- gene_file[
        !is.infinite(gene_file$V2),
      ]

      gene_file <- gene_file[
        !is.na(gene_file$V2),
      ]

      ## Sort the file by decreasing Z score:
      gene_file <- gene_file[
        order(
          gene_file$V2,
          decreasing = TRUE
        ),
      ]

      ## Add a ranking to each probe link:
      gene_file$rank <- c(
        1:nrow(gene_file)
      )

      ## Calculate the effective permutated p-value
      ## (i.e. which rank among the genes is our given gene):
      perm_p <- (
        gene_file[
          gene_file$V1==probe,
          'rank'
        ]/
        nrow(gene_file)
      )

      ## Return the p-value
      return(perm_p)

      ## Clear the gene_file:
      rm(gene_file)
    }

    ## Execute the function to calculate a permutation p-value
    ## for all significant probe-gene links:
    ## Now let's lapply the function if hypermeth analysis is selected:
    hypermeth_Gplus_zscores$perm_p_value <- parallel::mcmapply(
      getPosEpval_hyperGplus,
      geneID= hypermeth_Gplus_zscores$geneID,
      probe= hypermeth_Gplus_zscores$probeID,
      Z_score= hypermeth_Gplus_zscores$Zscore,
      mc.cores= core_count
    )

    ## If hypermeth_Gplus_zscores are not found, return an error message:
    if(nrow(hypermeth_Gplus_zscores)<1){

      stop("No significant hypermeth Gplus zscores were loaded. Please check the step3 hypermeth to confirm that potentially signficant hyper Gplus zscores were found.")
    }

    ## Save the file:
    write.table(
      hypermeth_Gplus_zscores,
      paste(
        TENET_directory,
        'step4/',
        'hyper_Gplus_sig_link_zscores_perm.txt',
        sep=''
      ),
      sep='\t',
      row.names= FALSE,
      quote= FALSE
    )

    ## Remove the old Zscore file:
    rm(hypermeth_Gplus_zscores)
  }

  ## If hypermeth Gminus analysis is selected, run it:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Check that the significant zscores from step3 exists and load it:
    ## If not, return an error message:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step3/',
          'hypermeth_results/',
          'hyper_Gminus_sig_link_zscores.txt',
          sep=''
        )
      )
    ){

      ## Load the file and change the colnames:
      hypermeth_Gminus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step3/',
          'hypermeth_results/',
          'hyper_Gminus_sig_link_zscores.txt',
          sep=''
        ),
        sep='\t',
        header= TRUE,
        stringsAsFactors = FALSE
      )

      colnames(hypermeth_Gminus_zscores) <- c(
        'geneID',
        'probeID',
        'Zscore'
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hyper_Gminus_sig_link_zscores.txt in step3 of TENET directory was not found. Please check that the file exists and consider rerunning the step3 get_analysis_z_score function.')
    }

    ## Write the function to perform the Z-score permutation:
    getPosEpval_hyperGminus <- function(
      geneID,
      probe,
      Z_score
    ){

      ## Set the relevant folder to find results in:
      step3_folder <- paste(
        TENET_directory,
        'step3/',
        'hypermeth_results',
        sep=''
      )

      ## Load the individual gene file for the given
      ## gene-probe link:
      gene_file <- read.delim(
        list.files(
          step3_folder,
          pattern= geneID,
          full.names = TRUE
        ),
        header= FALSE,
        sep = '\t',
        stringsAsFactors = FALSE
      )

      ## Remove NA and infinite values:
      gene_file <- gene_file[
        !is.infinite(gene_file$V2),
      ]

      gene_file <- gene_file[
        !is.na(gene_file$V2),
      ]

      ## Sort the file by increasing Z score:
      gene_file <- gene_file[
        order(
          gene_file$V2,
          decreasing = FALSE
        ),
      ]

      ## Add a ranking to each probe link:
      gene_file$rank <- c(
        1:nrow(gene_file)
      )

      ## Calculate the effective permutated p-value
      ## (i.e. which rank among the genes is our given gene):
      perm_p <- (
        gene_file[
          gene_file$V1==probe,
          'rank'
        ]/
          nrow(gene_file)
      )

      ## Return the p-value
      return(perm_p)

      ## Clear the gene_file:
      rm(gene_file)
    }

    ## Execute the function to calculate a permutation p-value
    ## for all significant probe-gene links:
    ## Now let's lapply the function if hypermeth analysis is selected:
    hypermeth_Gminus_zscores$perm_p_value <- parallel::mcmapply(
      getPosEpval_hyperGminus,
      geneID= hypermeth_Gminus_zscores$geneID,
      probe= hypermeth_Gminus_zscores$probeID,
      Z_score= hypermeth_Gminus_zscores$Zscore,
      mc.cores= core_count
    )

    ## If hypermeth_Gminus_zscores are not found, return an error message:
    if(nrow(hypermeth_Gminus_zscores)<1){

      stop("No significant hypermeth Gminus zscores were loaded. Please check the step3 hypermeth to confirm that potentially signficant hyper Gminus zscores were found.")
    }

    ## Save the file:
    write.table(
      hypermeth_Gminus_zscores,
      paste(
        TENET_directory,
        'step4/',
        'hyper_Gminus_sig_link_zscores_perm.txt',
        sep=''
      ),
      sep='\t',
      row.names = FALSE,
      quote= FALSE
    )

    ## Remove the old Zscore file:
    rm(hypermeth_Gminus_zscores)
  }

  ## If hypometh Gplus analysis is selected, run it first:
  if(hypometh_Gplus_analysis==TRUE){

    ## Check that the significant zscores from step3 exists and load it:
    ## If not, return an error message:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step3/',
          'hypometh_results/',
          'hypo_Gplus_sig_link_zscores.txt',
          sep=''
        )
      )
    ){

      ## Load the file and change the colnames:
      hypometh_Gplus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step3/',
          'hypometh_results/',
          'hypo_Gplus_sig_link_zscores.txt',
          sep=''
        ),
        sep='\t',
        header= TRUE,
        stringsAsFactors = FALSE
      )

      colnames(hypometh_Gplus_zscores) <- c(
        'geneID',
        'probeID',
        'Zscore'
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hypo_Gplus_sig_link_zscores.txt in step3 of TENET directory was not found. Please check that the file exists and consider rerunning the step3 get_analysis_z_score function.')
    }

    ## Write the function to perform the Z-score permutation:
    getPosEpval_hypoGplus <- function(
      geneID,
      probe,
      Z_score
    ){

      ## Set the relevant folder to find results in:
      step3_folder <- paste(
        TENET_directory,
        'step3/',
        'hypometh_results',
        sep=''
      )

      ## Load the individual gene file for the given
      ## gene-probe link:
      gene_file <- read.delim(
        list.files(
          step3_folder,
          pattern= geneID,
          full.names = TRUE
        ),
        header= FALSE,
        sep = '\t',
        stringsAsFactors = FALSE
      )

      ## Remove NA and infinite values:
      gene_file <- gene_file[
        !is.infinite(gene_file$V2),
      ]

      gene_file <- gene_file[
        !is.na(gene_file$V2),
      ]

      ## Sort the file by decreasing Z score:
      gene_file <- gene_file[
        order(
          gene_file$V2,
          decreasing = FALSE
        ),
      ]

      ## Add a ranking to each probe link:
      gene_file$rank <- c(
        1:nrow(gene_file)
      )

      ## Calculate the effective permutated p-value
      ## (i.e. which rank among the genes is our given gene):
      perm_p <- (
        gene_file[
          gene_file$V1==probe,
          'rank'
        ]/
          nrow(gene_file)
      )

      ## Return the p-value
      return(perm_p)

      ## Clear the gene_file:
      rm(gene_file)
    }

    ## Execute the function to calculate a permutation p-value
    ## for all significant probe-gene links:
    ## Now let's lapply the function if hypometh analysis is selected:
    hypometh_Gplus_zscores$perm_p_value <- parallel::mcmapply(
      getPosEpval_hypoGplus,
      geneID= hypometh_Gplus_zscores$geneID,
      probe= hypometh_Gplus_zscores$probeID,
      Z_score= hypometh_Gplus_zscores$Zscore,
      mc.cores= core_count
    )

    ## If hypometh_Gplus_zscores are not found, return an error message:
    if(nrow(hypometh_Gplus_zscores)<1){

      stop("No significant hypometh Gplus zscores were loaded. Please check the step3 hypometh to confirm that potentially signficant hypo Gplus zscores were found.")
    }

    ## Save the file:
    write.table(
      hypometh_Gplus_zscores,
      paste(
        TENET_directory,
        'step4/',
        'hypo_Gplus_sig_link_zscores_perm.txt',
        sep=''
      ),
      sep='\t',
      row.names= FALSE,
      quote= FALSE
    )

    ## Remove the old Zscore file:
    rm(hypometh_Gplus_zscores)
  }

  ## If hypometh Gminus analysis is selected, run it first:
  if(hypometh_Gminus_analysis==TRUE){

    ## Check that the significant zscores from step3 exists and load it:
    ## If not, return an error message:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step3/',
          'hypometh_results/',
          'hypo_Gminus_sig_link_zscores.txt',
          sep=''
        )
      )
    ){

      ## Load the file and change the colnames:
      hypometh_Gminus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step3/',
          'hypometh_results/',
          'hypo_Gminus_sig_link_zscores.txt',
          sep=''
        ),
        sep='\t',
        header= TRUE,
        stringsAsFactors = FALSE
      )

      colnames(hypometh_Gminus_zscores) <- c(
        'geneID',
        'probeID',
        'Zscore'
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hypo_Gminus_sig_link_zscores.txt in step3 of TENET directory was not found. Please check that the file exists and consider rerunning the step3 get_analysis_z_score function.')
    }

    ## Write the function to perform the Z-score permutation:
    getPosEpval_hypoGminus <- function(
      geneID,
      probe,
      Z_score
    ){

      ## Set the relevant folder to find results in:
      step3_folder <- paste(
        TENET_directory,
        'step3/',
        'hypometh_results',
        sep=''
      )

      ## Load the individual gene file for the given
      ## gene-probe link:
      gene_file <- read.delim(
        list.files(
          step3_folder,
          pattern= geneID,
          full.names = TRUE
        ),
        header= FALSE,
        sep = '\t',
        stringsAsFactors = FALSE
      )

      ## Remove NA and infinite values:
      gene_file <- gene_file[
        !is.infinite(gene_file$V2),
      ]

      gene_file <- gene_file[
        !is.na(gene_file$V2),
      ]

      ## Sort the file by decreasing Z score:
      gene_file <- gene_file[
        order(
          gene_file$V2,
          decreasing = TRUE
        ),
      ]

      ## Add a ranking to each probe link:
      gene_file$rank <- c(
        1:nrow(gene_file)
      )

      ## Calculate the effective permutated p-value
      ## (i.e. which rank among the genes is our given gene):
      perm_p <- (
        gene_file[
          gene_file$V1==probe,
          'rank'
        ]/
          nrow(gene_file)
      )

      ## Return the p-value
      return(perm_p)

      ## Clear the gene_file:
      rm(gene_file)
    }

    ## Execute the function to calculate a permutation p-value
    ## for all significant probe-gene links:
    ## Now let's lapply the function if hypometh analysis is selected:
    hypometh_Gminus_zscores$perm_p_value <- parallel::mcmapply(
      getPosEpval_hypoGminus,
      geneID= hypometh_Gminus_zscores$geneID,
      probe= hypometh_Gminus_zscores$probeID,
      Z_score= hypometh_Gminus_zscores$Zscore,
      mc.cores= core_count
    )

    ## If hypometh_Gminus_zscores are not found, return an error message:
    if(nrow(hypometh_Gminus_zscores)<1){

      stop("No significant hypometh Gminus zscores were loaded. Please check the step3 hypometh to confirm that potentially signficant hypo Gminus zscores were found.")
    }

    ## Save the file:
    write.table(
      hypometh_Gminus_zscores,
      paste(
        TENET_directory,
        'step4/',
        'hypo_Gminus_sig_link_zscores_perm.txt',
        sep=''
      ),
      sep='\t',
      row.names= FALSE,
      quote= FALSE
    )

    ## Remove the old Zscore file:
    rm(hypometh_Gminus_zscores)
  }
}

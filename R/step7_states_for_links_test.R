#' step7_states_for_links_test
#'
#' This is a step7 function of the TENETR package.
#' This function generates .tsv files with information for each of the experimental/tumor,
#' samples and each probe-gene link from each of the four hypo or hypermethylated
#' Gplus or Gminus analysis quadrants, as selected by the user, on if a given sample
#' is said "harbor" each link, depending on if the methylation of the given sample
#' for the probe in the link is above or below the hyper or hypomethcutoff defined
#' from the step2_get_diffmeth_regions function, and the expression of gene in the link
#' is significantly greater than, or less than, the mean expression in the normal/control
#' samples.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_states_for_links' containing the results of this function.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create .tsv files showing which experimental/tumor samples harbor a given hypermeth Gplus probe-gene link.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create .tsv files showing which experimental/tumor samples harbor a given hypermeth Gminus probe-gene link.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create .tsv files showing which experimental/tumor samples harbor a given hypometh Gplus probe-gene link.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create .tsv files showing which experimental/tumor samples harbor a given hypometh Gminus probe-gene link.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Exports .tsv files for each of the specified analysis types with experimental/tumor samples in the columns and each probe-gene link for that analysis type in the rows, with a 1 indicating the sample is positive for that link, and a 0 if that sample isn't.
#' @export

step7_states_for_links_test <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  core_count
){

  ## Check to make sure at least one analysis type has been selected and return
  ## an error message if at least one hasn't been
  if(
    hypermeth_Gplus_analysis==FALSE &
    hypermeth_Gminus_analysis==FALSE &
    hypometh_Gplus_analysis==FALSE &
    hypometh_Gminus_analysis==FALSE
  ){

    stop(
      "All analysis types have been set to false. Set at least one analysis type to TRUE"
    )
  }

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

  ## Create a step7 directory to deposit the output the pdfs and other
  ## downstream analyses if it has not been created yet:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        sep=''
      )
    )

  }

  ## Create a step7 directory to deposit the output the pdfs and other
  ## downstream analyses if it has not been created yet:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7',
        'states_for_links',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        'states_for_links',
        sep=''
      )
    )

  }


  ## List files containing results created by get_diffmeth_regions
  ## in step2:
  step2_files_list <- list.files(
    paste(
      TENET_directory,
      'step2/',
      sep=''
    ),
    full.names = TRUE
  )

  ## Get just the .rda files from the directory:
  ## Seperate the directories into ones with ENH and NDR files:
  step2_rda_files_list <- grep(
    'diff.methylated.datasets.rda',
    step2_files_list,
    value=TRUE
  )

  ## If exactly one .rda file was found, load it. Otherwise return errors
  ## to the user:
  if(length(step2_rda_files_list)==1){

    load(step2_rda_files_list)

  } else if(length(step2_rda_files_list)>1){

    stop(
      "More than one diff.methylated.datasets.rda file was found was found in TENET step2 directory. Please ensure only one file is present in the folder, output from get_diffmeth_regions function."
    )

  } else if(length(step2_rda_files_list)==0){

    stop(
      "No diff.methylated.datasets.rda file was found in TENET step2 directory. Please run get_diffmeth_regions function to create the diff.methylated.datasets.rda file"
    )
  }

  ## Check to ensure that correct objects are found in the loaded .rda file.
  ## If not, return an error message:

  ## Index a vector of objects:
  objects_should_be_present <- c(
    'enhancer_probes',
    'expDataN',
    'expDataT',
    'hypermeth_probes',
    'hypometh_probes',
    'metDataN',
    'metDataT',
    'unmeth_probes',
    'hypermethcutoff',
    'hypomethcutoff',
    'min_experimental_count'
  )

  ## Index an empty vector with the positions of objects that should be present
  ## but are not:
  positions_not_present <- numeric()

  ## Loop through each of the positions of the items that should be present
  ## and return the positions of those that are not:
  for(i in 1:length(objects_should_be_present)){

    if(!exists(objects_should_be_present[i])){

      positions_not_present <- c(
        positions_not_present,
        i
      )
    } else{

      positions_not_present <- positions_not_present
    }
  }

  ## Return an error message noting which objects are missing:
  if(length(positions_not_present)>0){

    stop(
      paste(
        "objects",
        paste(
          c(objects_should_be_present[positions_not_present]),
          collapse=', '
        ),
        "were not imported in the .rda file. Please examine imported data file and/or rerun get_diffmeth_regions function.",
        collapse=' '
      )
    )
  }

  ## Now that the data is loaded, let's reconstitute the datasets of interest:
  'enhmetDataN' <- metDataN[
    enhancer_probes,
  ]
  'enhmetDataT' <- metDataT[
    enhancer_probes,
  ]
  'hypermethDataN' <- metDataN[
    hypermeth_probes,
  ]
  'hypermethDataT' <- metDataT[
    hypermeth_probes,
  ]
  'hypomethDataN' <- metDataN[
    hypometh_probes,
  ]
  'hypomethDataT' <- metDataT[
    hypometh_probes,
  ]
  'unmethDataN' <- metDataN[
    unmeth_probes,
  ]
  'unmethDataT' <- metDataT[
    unmeth_probes,
  ]

  ## Get the dataset of gencode v22 genes:
  gencode_v22_gtf <- TENETR.data::gencode_v22_annotations

  ## Get info for genes
  gencode_v22_genes <- gencode_v22_gtf[
    gencode_v22_gtf$type=='gene',
  ]

  ## Determine whether the start or the end coordinate is the
  ## TSS based on the strand:
  gencode_v22_genes$TSS <- ifelse(
    gencode_v22_genes$strand=='+',
    gencode_v22_genes$start,
    ifelse(
      gencode_v22_genes$strand=='-',
      gencode_v22_genes$end,
      ''
    )
  )

  ## remove the stuff after the period in the ENSG ids:
  gencode_v22_genes$ensembl_ID <- sub(
    '\\..*',
    '',
    gencode_v22_genes$gene_id
  )

  ## Set the rownames to be the ensembl ids without periods:
  rownames(gencode_v22_genes) <- gencode_v22_genes$ensembl_ID

  ## Remove the unneeded gtf file:
  rm(gencode_v22_gtf)

  ## Write the function to calculate "personal links"
  ## This determines if a given tumor sample has a given probe-gene link
  ## Based on if it has bonferroni 1-sided p-value expression greater
  ## Than normal sample mean, and has less than the hypometh cutoff:
  link_evaluator <- function(
    methylation_probe_link,
    gene_link,
    Ggreater_or_Gless,
    hypermeth_or_hypometh
  ){

    ## Get actual gene name depending on input:
    if(substring(gene_link,1,4)=='ENSG' & nchar(gene_link)==15){

      ## Input is in ENSG, get the gene name:
      gene_link_ENSG <- gene_link

      gene_name <- gencode_v22_genes[
        gene_link_ENSG,
        'gene_name'
      ]

      ## Get gene ENSG assuming a name is plugged in:
    } else{

      ## Assume what was given was the gene name, get the ENSG:
      gene_link_name <- gene_link

      gene_link_ENSG <- rownames(
        gencode_v22_genes_df[
          gencode_v22_genes_df$gene_name==gene_link_name,
        ]
      )
    }

    ## extract methylation values from each tumor sample
    tumor_methylation_values_for_probe <- unlist(
      metDataT[
        methylation_probe_link,
      ]
    )

    ## extract gene expression values for each tumor sample
    tumor_expression_values_for_gene <- unlist(
      expDataT[
        gene_link_ENSG,
      ]
    )

    ## Get expression values from the normal samples:
    normal_expression_values_for_gene <- unlist(
      expDataN[
        gene_link_ENSG,
      ]
    )

    ## Use a for loop to determine for which samples the mean in the normal samples
    ## is significantly lower:
    p_values_compared_to_normal <- numeric()

    for(i in tumor_expression_values_for_gene){

      if(Ggreater_or_Gless=='Ggreater'){

        p_values_compared_to_normal <- c(
          p_values_compared_to_normal,
          t.test(
            normal_expression_values_for_gene,
            mu=i,
            alternative="greater"
          )$p.value
        )

      } else if(Ggreater_or_Gless=='Gless'){

        p_values_compared_to_normal <- c(
          p_values_compared_to_normal,
          t.test(
            normal_expression_values_for_gene,
            mu=i,
            alternative="less"
          )$p.value
        )

      }
    }

    ## Bonferroni correct the p-values:
    p_values_compared_to_normal_bonf <- p.adjust(
      p_values_compared_to_normal,
      method= 'bonferroni',
      n=length(p_values_compared_to_normal)
    )

    ## Create a vector that is 1 if the given sample has expression for the gene
    ## significantly different from the normal samples
    significant_from_normal <- ifelse(
      p_values_compared_to_normal_bonf<0.05,
      1,
      0
    )

    ## For each methylation sample, check if it is hyper or hypomethylated for the given probe:
    if(hypermeth_or_hypometh=='hypometh'){

      cutoff_values <- ifelse(
        tumor_methylation_values_for_probe < hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_or_hypometh=='hypermeth'){

      cutoff_values <- ifelse(
        tumor_methylation_values_for_probe > hypermethcutoff,
        1,
        0
      )
    }

    ## Create a final vector with 1 for where both conditions are met, and 0 when
    ## they arent both met:
    final_output <- ifelse(
      significant_from_normal==1,
      ifelse(
        cutoff_values==1,
        1,
        0
      ),
      0
    )

    ## Assign sample names to the output values:
    names(final_output) <- colnames(expDataT)

    ## Return the final output:
    return(final_output)
  }

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus TAD tables:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hyper_Gplus_states_for_links',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hyper_Gplus_states_for_links',
          sep=''
        )
      )
    }

    ## Check that the hyper_Gplus_sig_link_zscores_perm_optimized.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step5/',
          'hyper_Gplus_sig_link_zscores_perm_optimized.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gplus_sig_link_zscores <- read.delim(
        paste(
          TENET_directory,
          'step5/',
          'hyper_Gplus_sig_link_zscores_perm_optimized.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gplus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function.')

    }

    ## Calculate hypermeth Gplus links:
    hyper_Gplus_links_dataset <- parallel::mcmapply(
      link_evaluator,
      methylation_probe_link=hyper_Gplus_sig_link_zscores$probeID,
      gene_link=hyper_Gplus_sig_link_zscores$geneID,
      Ggreater_or_Gless='Ggreater',
      hypermeth_or_hypometh='hypermeth',
      mc.cores = core_count
    )

    ## Change the column names to those of the probe + gene link:
    colnames(hyper_Gplus_links_dataset) <- paste(
      hyper_Gplus_sig_link_zscores$probeID,
      hyper_Gplus_sig_link_zscores$geneID,
      sep='_'
    )

    ## Transpose the dataset to make the tumor samples the columns:
    hyper_Gplus_links_dataset_transposed <- t(hyper_Gplus_links_dataset)

    ## Convert it to a dataframe
    hyper_Gplus_links_dataset_transposed_df <- as.data.frame(
      hyper_Gplus_links_dataset_transposed,
      stringsAsFactors= FALSE
    )

    ## Write the table out to the states subdirectory:
    write.table(
      hyper_Gplus_links_dataset_transposed_df,
      file= paste(
        TENET_directory,
        'step7/',
        'states_for_links/',
        'hyper_Gplus_states_for_links/',
        'hyper_Gplus_links_states_table.tsv',
        sep=''
      ),
      quote= FALSE,
      sep='\t'
    )
  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gminus TAD tables:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hyper_Gminus_states_for_links',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hyper_Gminus_states_for_links',
          sep=''
        )
      )
    }

    ## Check that the hyper_Gminus_sig_link_zscores_perm_optimized.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step5/',
          'hyper_Gminus_sig_link_zscores_perm_optimized.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gminus_sig_link_zscores <- read.delim(
        paste(
          TENET_directory,
          'step5/',
          'hyper_Gminus_sig_link_zscores_perm_optimized.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gminus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function.')

    }

    ## Calculate hypermeth Gminus links:
    hyper_Gminus_links_dataset <- parallel::mcmapply(
      link_evaluator,
      methylation_probe_link=hyper_Gminus_sig_link_zscores$probeID,
      gene_link=hyper_Gminus_sig_link_zscores$geneID,
      Ggreater_or_Gless='Gless',
      hypermeth_or_hypometh='hypermeth',
      mc.cores = core_count
    )

    ## Change the column names to those of the probe + gene link:
    colnames(hyper_Gminus_links_dataset) <- paste(
      hyper_Gminus_sig_link_zscores$probeID,
      hyper_Gminus_sig_link_zscores$geneID,
      sep='_'
    )

    ## Transpose the dataset to make the tumor samples the columns:
    hyper_Gminus_links_dataset_transposed <- t(hyper_Gminus_links_dataset)

    ## Convert it to a dataframe
    hyper_Gminus_links_dataset_transposed_df <- as.data.frame(
      hyper_Gminus_links_dataset_transposed,
      stringsAsFactors= FALSE
    )

    ## Write the table out to the states subdirectory:
    write.table(
      hyper_Gminus_links_dataset_transposed_df,
      file= paste(
        TENET_directory,
        'step7/',
        'states_for_links/',
        'hyper_Gminus_states_for_links/',
        'hyper_Gminus_links_states_table.tsv',
        sep=''
      ),
      quote= FALSE,
      sep='\t'
    )
  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gplus TAD tables:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hypo_Gplus_states_for_links',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hypo_Gplus_states_for_links',
          sep=''
        )
      )
    }

    ## Check that the hypo_Gplus_sig_link_zscores_perm_optimized.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step5/',
          'hypo_Gplus_sig_link_zscores_perm_optimized.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gplus_sig_link_zscores <- read.delim(
        paste(
          TENET_directory,
          'step5/',
          'hypo_Gplus_sig_link_zscores_perm_optimized.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gplus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function.')

    }

    ## Calculate hypometh Gplus links:
    hypo_Gplus_links_dataset <- parallel::mcmapply(
      link_evaluator,
      methylation_probe_link=hypo_Gplus_sig_link_zscores$probeID,
      gene_link=hypo_Gplus_sig_link_zscores$geneID,
      Ggreater_or_Gless='Gless',
      hypermeth_or_hypometh='hypometh',
      mc.cores = core_count
    )

    ## Change the column names to those of the probe + gene link:
    colnames(hypo_Gplus_links_dataset) <- paste(
      hypo_Gplus_sig_link_zscores$probeID,
      hypo_Gplus_sig_link_zscores$geneID,
      sep='_'
    )

    ## Transpose the dataset to make the tumor samples the columns:
    hypo_Gplus_links_dataset_transposed <- t(hypo_Gplus_links_dataset)

    ## Convert it to a dataframe
    hypo_Gplus_links_dataset_transposed_df <- as.data.frame(
      hypo_Gplus_links_dataset_transposed,
      stringsAsFactors= FALSE
    )

    ## Write the table out to the states subdirectory:
    write.table(
      hypo_Gplus_links_dataset_transposed_df,
      file= paste(
        TENET_directory,
        'step7/',
        'states_for_links/',
        'hypo_Gplus_states_for_links/',
        'hypo_Gplus_links_states_table.tsv',
        sep=''
      ),
      quote= FALSE,
      sep='\t'
    )
  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gminus TAD tables:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hypo_Gminus_states_for_links',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'states_for_links/',
          'hypo_Gminus_states_for_links',
          sep=''
        )
      )
    }

    ## Check that the hypo_Gminus_sig_link_zscores_perm_optimized.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step5/',
          'hypo_Gminus_sig_link_zscores_perm_optimized.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gminus_sig_link_zscores <- read.delim(
        paste(
          TENET_directory,
          'step5/',
          'hypo_Gminus_sig_link_zscores_perm_optimized.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gminus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function.')

    }

    ## Calculate hypometh Gminus links:
    hypo_Gminus_links_dataset <- parallel::mcmapply(
      link_evaluator,
      methylation_probe_link=hypo_Gminus_sig_link_zscores$probeID,
      gene_link=hypo_Gminus_sig_link_zscores$geneID,
      Ggreater_or_Gless='Ggreater',
      hypermeth_or_hypometh='hypometh',
      mc.cores = core_count
    )

    ## Change the column names to those of the probe + gene link:
    colnames(hypo_Gminus_links_dataset) <- paste(
      hypo_Gminus_sig_link_zscores$probeID,
      hypo_Gminus_sig_link_zscores$geneID,
      sep='_'
    )

    ## Transpose the dataset to make the tumor samples the columns:
    hypo_Gminus_links_dataset_transposed <- t(hypo_Gminus_links_dataset)

    ## Convert it to a dataframe
    hypo_Gminus_links_dataset_transposed_df <- as.data.frame(
      hypo_Gminus_links_dataset_transposed,
      stringsAsFactors= FALSE
    )

    ## Write the table out to the states subdirectory:
    write.table(
      hypo_Gminus_links_dataset_transposed_df,
      file= paste(
        TENET_directory,
        'step7/',
        'states_for_links/',
        'hypo_Gminus_states_for_links/',
        'hypo_Gminus_links_states_table.tsv',
        sep=''
      ),
      quote= FALSE,
      sep='\t'
    )
  }
}

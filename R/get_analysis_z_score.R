#' get_analysis_z_score
#'
#' This is the step3 function of the TENETR package.
#' This function takes the identified hyper and hypomethylated probes
#' from step2 get_diffmeth_regions function and calculates a zscore comparing
#' the mean expression of each gene in groups that are hyper/hypomethylated for
#' each probe, according to hypermeth and hypometh and those that are not, and permutates this across all
#' hyper/hypomethylated probes
#'
#'
#' @param TENET_directory Set a path to the directory that contains step2 results from get_diffmeth_regions function. This function will also create a new step3 folder there containing the results.
#' @param hypermeth_analysis Set to TRUE/FALSE depending on if you want to calculate z-scores for hypermethylated probes.
#' @param hypometh_analysis Set to TRUE/FALSE depending on if you want to calculate z-scores for hypomethylated probes.
#' @param usecaseonly Set to TRUE/FALSE depending on if you want to include the control/normal samples with the experimental/tumor samples when calculating hyper/hypomethylated groups and z-scores.
#' @param TF_only Set to TRUE/FALSE to determine if you only want to consider genes that are accepted transcription factors in The Human Transcription Factors by Lambert et al (2018) when calculating z-scores.
#' @param significant_p_value Set p-value to identify significant Z-scores for gene expression values selected between hyper/hypomethylated tumor/experimental samples and those that are not.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns tab-delimited "sig_link_zscores.txt" files for hypo/hyper Gplus/Gminus probe-gene links, as well as individual "zscores.txt" files named after each gene in the hypo/hyper analysis with zscores for that gene to all probes from that analysis type.
#' @export

get_analysis_z_score <- function(
  TENET_directory,
  hypermeth_analysis,
  hypometh_analysis,
  usecaseonly,
  TF_only,
  significant_p_value,
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
      'step3/',
      sep=''
    )
  )

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
    'enhmetDataN',
    'enhmetDataT',
    'expDataN',
    'expDataT',
    'hypermethDataN',
    'hypermethDataT',
    'hypomethDataN',
    'hypomethDataT',
    'metDataN',
    'metDataT',
    'unmethDataN',
    'unmethDataT',
    'hypermethcutoff',
    'hypomethcutoff',
    'minExp'
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

  ## Now that the datasets are loaded, let's make tables containing Z scores
  ## for hypo and hypermeth analyses:
  if (usecaseonly==T){

    expData <- expDataT

    if(hypermeth_analysis==TRUE & hypometh_analysis==TRUE){

      metDataHyper <- ifelse(
        hypermethDataT<hypermethcutoff,
        1,
        0
      )

      metDataHypo <- ifelse(
        hypomethDataT<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==TRUE & hypometh_analysis==FALSE){

      metDataHyper <- ifelse(
        hypermethDataT<hypermethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==FALSE & hypometh_analysis==TRUE){

      metDataHypo <- ifelse(
        hypomethDataT<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==FALSE & hypometh_analysis==FALSE){

      ## Return an error because no analyses are selected:
      stop(
        "hypermeth and hypometh analyses are both not set to run. Please set at least one analysis type to TRUE."
      )

    }

  } else if(usecaseonly==F){

    expData <- cbind(
      expDataT,
      expDataN
    )

    if(hypermeth_analysis==TRUE & hypometh_analysis==TRUE){

      metDataHyper <- ifelse(
        cbind(
          hypermethDataT,
          hypermethDataN
        )<hypermethcutoff,
        1,
        0
      )

      metDataHypo <- ifelse(
        cbind(
          hypomethDataT,
          hypomethDataN
        )<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==TRUE & hypometh_analysis==FALSE){

      metDataHyper <- ifelse(
        cbind(
          hypermethDataT,
          hypermethDataN
        )<hypermethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==FALSE & hypometh_analysis==TRUE){

      metDataHypo <- ifelse(
        cbind(
          hypomethDataT,
          hypomethDataN
        )<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_analysis==FALSE & hypometh_analysis==FALSE){

      ## Return an error because no analyses are selected:
      stop(
        "hypermeth and hypometh analyses are both not set to run. Please set at least one analysis type to TRUE."
      )

    }

  }

  ## Now let's import the gtf file again:
  ## Get the dataset of gencode v22 genes:
  gencode_v22_gtf <- TENETR.data::gencode_v22_annotations

  gencode_v22_genes <- gencode_v22_gtf[
    gencode_v22_gtf$type=='gene',
  ]

  rownames(gencode_v22_genes) <- sub(
    '\\..*',
    '',
    gencode_v22_genes$gene_id
  )

  ## Clear the gtf file:
  rm(gencode_v22_gtf)

  ## Now lets sort the expData alpanumerically:
  expData <- expData[
    order(
      rownames(expData)
    ),
  ]

  ## Now let's generate codex files listing the order of hyper/hypometh probes:
  ## As well as files with significant Z-scores by genes:
  if(hypermeth_analysis==TRUE){

    ## Create a hypometh folder:
    dir.create(
      paste(
        TENET_directory,
        'step3/',
        'hypermeth_results/',
        sep=''
      )
    )

    write(
      c('geneID','probeID','Zscore'),
      file= paste(
        TENET_directory,
        'step3/',
        'hypermeth_results/',
        'hyper_Gplus_sig_link_zscores.txt',
        sep=''
      ),
      ncolumns=3,
      append= FALSE,
      sep='\t'
    )

    write(
      c('geneID','probeID','Zscore'),
      file= paste(
        TENET_directory,
        'step3/',
        'hypermeth_results/',
        'hyper_Gminus_sig_link_zscores.txt',
        sep=''
      ),
      ncolumns=3,
      append= FALSE,
      sep='\t'
    )

  }

  if(hypometh_analysis==TRUE){

    ## Create a hypometh folder:
    dir.create(
      paste(
        TENET_directory,
        'step3/',
        'hypometh_results/',
        sep=''
      )
    )

    write(
      c('geneID','probeID','Zscore'),
      file= paste(
        TENET_directory,
        'step3/',
        'hypometh_results/',
        'hypo_Gplus_sig_link_zscores.txt',
        sep=''
      ),
      ncolumns=3,
      append= FALSE,
      sep='\t'
    )

    write(
      c('geneID','probeID','Zscore'),
      file= paste(
        TENET_directory,
        'step3/',
        'hypometh_results/',
        'hypo_Gminus_sig_link_zscores.txt',
        sep=''
      ),
      ncolumns=3,
      append= FALSE,
      sep='\t'
    )

  }

  ## If looking only at TF genes, subset expData to just those genes:
  if(TF_only==TRUE){

    ## Get the names of the human transcription factors:
    Human_TF <- TENETR.data::human_transcription_factors_dataset[
      TENETR.data::human_transcription_factors_dataset$Is.TF.=='Yes',
      "Ensembl.ID"
    ]

    ## Identify those present in the expression data
    Human_TF_present <- Human_TF[
      Human_TF %in% rownames(expData)
    ]

    ## Subset the expression data to just those genes:
    expData <- expData[
      Human_TF_present[order(Human_TF_present)],
    ]
  }

  ## Determine a significant z-score from the input p-value:
  significant_z_score <- qnorm(
    1-significant_p_value
  )

  ## Now lets write a function to calculate z scores for expression of each gene
  ## between groups that are hyper/hypomethylated for each probe, and those that
  ## aren't:
  z_score_calc <- function(
    gene_ENSG,
    hyper_hypo,
    metData_internal,
    minExp_value
  ){

    ## Write Z-score calculation function:
    internal_probe_z_function <- function(
      internal_probe_ID
    ){

      ## Split the methylation dataset based on  whether or not the probe
      ## was hypomethylated in each sample:
      T_P_probe <- metData_internal[
        internal_probe_ID,
      ]

      T_P_probe <- T_P_probe[
        !is.na(T_P_probe)
      ]

      T_P_probe_1 <- T_P_probe[(T_P_probe==1)]
      T_P_probe_0 <- T_P_probe[(T_P_probe==0)]

      ## If probe is viable based on parameters (a double check):
      if(length(T_P_probe_1)>minExp_value & length(T_P_probe_0)>0){

        P_gene_1 <- expData[
          gene_ENSG,
          names(T_P_probe_1)
        ]

        P_gene_0 <- expData[
          gene_ENSG,
          names(T_P_probe_0)
        ]

        z_score_internal <- (
          (mean(P_gene_0, na.rm=T) - mean(P_gene_1, na.rm=T))/
            (sd(P_gene_1, na.rm=T))
        )

        ## Round the Z_score to 6 digits to conserve memory
        z_score_internal <- round(
          z_score_internal,
          6
        )

        return(z_score_internal)

      } else{

        ## For that probe, return an NA value:
        return(NA)
      }
    }

    if(hyper_hypo=='hyper'){

      ## Create dataframe with probe IDs:
      TESTSR=data.frame(
        geneID= rep(
          gene_ENSG,
          nrow(metDataHyper)
        ),
        probeID= rownames(metDataHyper),
        stringsAsFactors = FALSE
      )
      rownames(TESTSR) <- rownames(metDataHyper)

      ## Use the function to calculate z-scores for all probes to the gene:
      TESTSR$Z.real <- sapply(
        TESTSR$probeID,
        internal_probe_z_function
      )

      ## Isolate genes with only significant Zscores:
      TESTSR_sig_neg <- TESTSR[
        (TESTSR$Z.real > significant_z_score),
      ]

      TESTSR_sig_pos <- TESTSR[
        (TESTSR$Z.real < -significant_z_score),
      ]

      ## Remove values with infinite Z scores;
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.infinite(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.infinite(TESTSR_sig_neg$Z.real),
      ]

      ## Remove values with NaN Z scores:
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.nan(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.nan(TESTSR_sig_neg$Z.real),
      ]

      ## Remove values with NA Z scores:
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.na(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.na(TESTSR_sig_neg$Z.real),
      ]

      ## Append the significant Z_scores to the appropriate
      ## sig_link_zscores datasets:
      if(nrow(TESTSR_sig_pos)>0){

        ## Sort the tables by alphanumeric probe:
        TESTSR_sig_pos <- TESTSR_sig_pos[
          order(TESTSR_sig_pos$probeID),
        ]

        ## Add to existing df:
        write.table(
          TESTSR_sig_pos,
          file= paste(
            TENET_directory,
            'step3/',
            'hypermeth_results/',
            'hyper_Gplus_sig_link_zscores.txt',
            sep=''
          ),
          append= TRUE,
          quote= FALSE,
          col.names=FALSE,
          row.names= FALSE,
          sep='\t'
        )

      }

      if(nrow(TESTSR_sig_neg)>0){

        ## Sort the tables by alphanumeric probe:
        TESTSR_sig_neg <- TESTSR_sig_neg[
          order(TESTSR_sig_neg$probeID),
        ]

        ## Add to existing df:
        write.table(
          TESTSR_sig_neg,
          file= paste(
            TENET_directory,
            'step3/',
            'hypermeth_results/',
            'hyper_Gminus_sig_link_zscores.txt',
            sep=''
          ),
          append= TRUE,
          quote= FALSE,
          col.names=FALSE,
          row.names= FALSE,
          sep='\t'
        )

      }

    } else if(hyper_hypo=='hypo'){

      ## Create dataframe with probe IDs:
      TESTSR=data.frame(
        geneID= rep(
          gene_ENSG,
          nrow(metDataHypo)
        ),
        probeID= rownames(metDataHypo),
        stringsAsFactors = FALSE
      )
      rownames(TESTSR) <- rownames(metDataHypo)

      ## Use the function to calculate z-scores for all probes to the gene:
      TESTSR$Z.real <- sapply(
        TESTSR$probeID,
        internal_probe_z_function
      )

      ## Isolate genes with only significant Zscores:
      TESTSR_sig_pos <- TESTSR[
        (TESTSR$Z.real > significant_z_score),
      ]

      TESTSR_sig_neg <- TESTSR[
        (TESTSR$Z.real < -significant_z_score),
      ]

      ## Remove values with infinite Z scores;
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.infinite(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.infinite(TESTSR_sig_neg$Z.real),
      ]

      ## Remove values with NaN Z scores:
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.nan(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.nan(TESTSR_sig_neg$Z.real),
      ]

      ## Remove values with NA Z scores:
      TESTSR_sig_pos <- TESTSR_sig_pos[
        !is.na(TESTSR_sig_pos$Z.real),
      ]

      TESTSR_sig_neg <- TESTSR_sig_neg[
        !is.na(TESTSR_sig_neg$Z.real),
      ]

      ## Append the significant Z_scores to the appropriate
      ## sig_link_zscores datasets:
      if(nrow(TESTSR_sig_pos)>0){

        ## Sort the tables by alphanumeric probe:
        TESTSR_sig_pos <- TESTSR_sig_pos[
          order(TESTSR_sig_pos$probeID),
        ]

        ## Add to existing df:
        write.table(
          TESTSR_sig_pos,
          file= paste(
            TENET_directory,
            'step3/',
            'hypometh_results/',
            'hypo_Gminus_sig_link_zscores.txt',
            sep=''
          ),
          append= TRUE,
          quote= FALSE,
          col.names=FALSE,
          row.names= FALSE,
          sep='\t'
        )

      }

      if(nrow(TESTSR_sig_neg)>0){

        ## Sort the tables by alphanumeric probe:
        TESTSR_sig_neg <- TESTSR_sig_neg[
          order(TESTSR_sig_neg$probeID),
        ]

        ## Add to existing df:
        write.table(
          TESTSR_sig_neg,
          file= paste(
            TENET_directory,
            'step3/',
            'hypometh_results/',
            'hypo_Gplus_sig_link_zscores.txt',
            sep=''
          ),
          append= TRUE,
          quote= FALSE,
          col.names=FALSE,
          row.names= FALSE,
          sep='\t'
        )

      }
    }

    ## Remove the gene name from the TESTSR dataframe:
    TESTSR$geneID <- NULL

    ## Set an analysis folder:
    if(hyper_hypo=='hyper'){

      output_folder <- 'hypermeth_results/'

    } else if(hyper_hypo=='hypo'){

      output_folder <- 'hypometh_results/'

    }

    ## Now write the complete results file with the name of the gene:
    ## and the analysis type:
    write.table(
      TESTSR,
      file= paste(
        TENET_directory,
        'step3/',
        output_folder,
        paste(
          hyper_hypo,
          gene_ENSG,
          'zscores.txt',
          sep='_'
        ),
        sep=''
      ),
      row.names= FALSE,
      col.names= FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Clear memory after each run:
    rm(TESTSR)
    rm(TESTSR_sig_pos)
    rm(TESTSR_sig_neg)
  }

  if(hypermeth_analysis==TRUE){

    ## Now let's lapply the function if hypermeth analysis is selected:
    parallel::mclapply(
      rownames(expData),
      z_score_calc,
      metData_internal= metDataHyper,
      hyper_hypo= 'hyper',
      minExp_value= minExp,
      mc.cores= core_count
    )

  }

  if(hypometh_analysis==TRUE){

    ## Now let's lapply the function if hypometh analysis is selected:
    parallel::mclapply(
      rownames(expData),
      z_score_calc,
      metData_internal= metDataHypo,
      hyper_hypo= 'hypo',
      minExp_value= minExp,
      mc.cores= core_count
    )

  }
}

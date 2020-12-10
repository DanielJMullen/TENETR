#' get_pairs_z_score
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
#' @param hypometh_pairs Set to TRUE/FALSE depending on if you want to calculate z-scores for hypomethylated probes.
#' @param hypermeth_pairs Set to TRUE/FALSE depending on if you want to calculate z-scores for hypermethylated probes.
#' @param usecaseonly Set to TRUE/FALSE depending on if you want to include the control/normal samples with the experimental/tumor samples when calculating hyper/hypomethylated groups and z-scores.
#' @param TF_only Set to TRUE/FALSE to determine if you only want to consider genes that are accepted transcription factors in The Human Transcription Factors by Lambert et al (2018) when calculating z-scores.
#' @param significant_p_value Set p-value to identify significant Z-scores for gene expression values selected between hyper/hypomethylated tumor/experimental samples and those that are not.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns tab-delimited ".txt"zscore_all_genes_rda.txt" for each probe of the selected types analyzed, containing the zscore
#' @export

get_pairs_z_score <- function(
  TENET_directory,
  hypermeth_pairs,
  hypometh_pairs,
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

    if(hypermeth_pairs==TRUE & hypometh_pairs==TRUE){

      metDataHyper <- ifelse(
        hypermethDataT>hypermethcutoff,
        1,
        0
      )

      metDataHypo <- ifelse(
        hypomethDataT<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_pairs==TRUE & hypometh_pairs==FALSE){

      metDataHyper <- ifelse(
        hypermethDataT>hypermethcutoff,
        1,
        0
      )

    } else if(hypermeth_pairs==FALSE & hypometh_pairs==TRUE){

      metDataHypo <- ifelse(
        hypomethDataT<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_pairs==FALSE & hypometh_pairs==FALSE){

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

    if(hypermeth_pairs==TRUE & hypometh_pairs==TRUE){

      metDataHyper <- ifelse(
        cbind(
          hypermethDataT,
          hypermethDataN
        )>hypermethcutoff,
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

    } else if(hypermeth_pairs==TRUE & hypometh_pairs==FALSE){

      metDataHyper <- ifelse(
        cbind(
          hypermethDataT,
          hypermethDataN
        )>hypermethcutoff,
        1,
        0
      )

    } else if(hypermeth_pairs==FALSE & hypometh_pairs==TRUE){

      metDataHypo <- ifelse(
        cbind(
          hypomethDataT,
          hypomethDataN
        )<hypomethcutoff,
        1,
        0
      )

    } else if(hypermeth_pairs==FALSE & hypometh_pairs==FALSE){

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
  if(hypermeth_pairs==TRUE){

    write(
      rownames(metDataHyper),
      file= paste(
        TENET_directory,
        'step3/',
        'hypermeth_probe_names_ordered.txt',
        sep=''
      ),
      ncolumns= 1,
      append= FALSE,
      sep= "\t"
    )

  }

  if(hypometh_pairs==TRUE){

    write(
      rownames(metDataHypo),
      file= paste(
        TENET_directory,
        'step3/',
        'hypometh_probe_names_ordered.txt',
        sep=''
      ),
      ncolumns= 1,
      append= FALSE,
      sep= "\t"
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

  ## Now print out a list of the genes remaining in expData:
  write(
    rownames(expData),
    file= paste(
      TENET_directory,
      'step3/',
      'gene_names_ordered.txt',
      sep=''
    ),
    ncolumns= 1,
    append= FALSE,
    sep= "\t"
  )

  ## Determine a significant z-score from the input p-value:
  significant_z_score <- qnorm(
    1-significant_p_value
  )

  ## Now lets write a function to calculate z scores for expression of each gene
  ## between groups that are hyper/hypomethylated for each probe, and those that
  ## aren't:
  z_score_calc <- function(
    DNA_methylation_probe,
    metData_internal,
    hyper_hypo,
    minExp_value,
    TF_only_internal
  ){

    ## Create dataframe with gene IDs and names:
    TESTSR=data.frame(
      geneID= rownames(expData),
      stringsAsFactors = FALSE
    )
    rownames(TESTSR) <- TESTSR$geneID

    ## Split the methylation dataset based on  whether or not the probe
    ## was hypomethylated in each sample:
    T_P_probe <- metData_internal[
      DNA_methylation_probe,
    ]

    T_P_probe <- T_P_probe[
      !is.na(T_P_probe)
    ]

    T_P_probe_1 <- T_P_probe[(T_P_probe==1)]
    T_P_probe_0 <- T_P_probe[(T_P_probe==0)]

    if(length(T_P_probe_1)>minExp_value & length(T_P_probe_0)>0){

      ## Write an internal function to perform the zscore calculations
      ## on gene expression in the two groups:
      internal_gene_function <- function(
        internal_gene_ID
      ){

        P_gene_1 <- expData[
          internal_gene_ID,
          names(T_P_probe_1)
        ]

        P_gene_0 <- expData[
          internal_gene_ID,
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
      }

      ## Use the function to calculate z-scores for all genes to the probe:
      TESTSR$Z.real <- sapply(
        TESTSR$geneID,
        internal_gene_function
      )

      ## Isolate genes with only significant Zscores:
      TESTSR <- TESTSR[
        (
          (TESTSR$Z.real > significant_z_score) |
          (TESTSR$Z.real < -significant_z_score)
        ),
      ]

      ## Remove values with infinite Z scores;
      TESTSR <- TESTSR[
        !is.infinite(TESTSR$Z.real),
      ]

      ## Remove values with NaN Z scores:
      TESTSR <- TESTSR[
        !is.nan(TESTSR$Z.real),
      ]

      ## Write out the file:
      write.table(
        TESTSR,
        file= paste(
          TENET_directory,
          'step3/',
          DNA_methylation_probe,
          '_zscores.txt',
          sep=''
        ),
        row.names=F,
        col.names=F,
        quote=F,
        sep="\t"
      )

    } else{

      ## Return an error message:
      stop("Improperly labelled hyper/hypomethylated probes detected. Please reexamine your parameters (especially the minExp value) and rerun get_diffmeth_regions")
    }

  }

  if(hypermeth_pairs==TRUE){

    ## Now let's lapply the function if hypermeth analysis is selected:
    parallel::mclapply(
      rownames(metDataHyper),
      z_score_calc,
      metData_internal= metDataHyper,
      hyper_hypo= 'hyper',
      TF_only_internal= TRUE,
      minExp_value= minExp,
      mc.cores= core_count
    )

  }

  if(hypometh_pairs==TRUE){

    ## Now let's lapply the function if hypometh analysis is selected:
    parallel::mclapply(
      rownames(metDataHypo),
      z_score_calc,
      metData_internal= metDataHypo,
      hyper_hypo= 'hypo',
      TF_only_internal= TRUE,
      minExp_value= minExp,
      mc.cores= core_count
    )

  }

  ## Once the files have been generated, clear the workspace
  ## and dump memory:
  rm(
    list=setdiff(
      ls(),
      c(
        "TENET_directory",
        "hypermeth_pairs",
        "hypometh_pairs",
        "significant_z_score"
      )
    )
  )

  gc()

  ## Create a list of the files generated previously
  list_files <- list.files(
    paste(
      TENET_directory,
      'step3/',
      sep=''
    ),
    full.names = TRUE
  )

  ## Find the file containing the list of hypermeth probes
  ## and load it:
  hypermeth_probes_file <- grep(
    'hypermeth_probe_names_ordered',
    list_files,
    value=TRUE
  )

  hypermeth_probes <- read.delim(
    hypermeth_probes_file,
    sep='\t',
    header= FALSE,
    stringsAsFactors = FALSE
  )

  ## Create a listing of all the hypermeth files with
  ## z scores that were created:
  hypermeth_z_scores_files <- paste(
    TENET_directory,
    'step3/',
    hypermeth_probes$V1,
    '_zscores.txt',
    sep=''
  )

  ## Create an empty data frame containing gene IDs
  ## and Z scores (initially)
  hyper_z_scores <- data.frame(
    'geneID'= character(),
    'Zscore'= numeric(),
    'probeID'= character(),
    stringsAsFactors = FALSE
  )

  ## Find the file containing the list of hypometh probes
  ## and load it:
  hypometh_probes_file <- grep(
    'hypometh_probe_names_ordered',
    list_files,
    value=TRUE
  )

  hypometh_probes <- read.delim(
    hypometh_probes_file,
    sep='\t',
    header= FALSE,
    stringsAsFactors = FALSE
  )

  ## Create a listing of all the hypometh files with
  ## z scores that were created:
  hypometh_z_scores_files <- paste(
    TENET_directory,
    'step3/',
    hypometh_probes$V1,
    '_zscores.txt',
    sep=''
  )

  ## Create an empty data frame containing gene IDs
  ## and Z scores (initially)
  hypo_z_scores <- data.frame(
    'geneID'= character(),
    'Zscore'= numeric(),
    'probeID'= character(),
    stringsAsFactors = FALSE
  )

  ## Find the file containing the list of genes of interest
  ## and load it:
  genes_file <- grep(
    'gene_names_ordered',
    list_files,
    value=TRUE
  )

  genes_present <- read.delim(
    genes_file,
    sep='\t',
    header= FALSE,
    stringsAsFactors = FALSE
  )

  ## If hypermeth analysis is selected, load the significant
  ## hypermeth z scores
  if(hypermeth_pairs==TRUE){

    hyper_file_combiner <- function(file_path){

      ## If the file isn't empty, load it:
      if(
        (file.info(file_path)$size)>0
      ){

        ## Load the file:
        file_placeholder <- read.delim(
          file_path,
          header= FALSE,
          sep='\t',
          stringsAsFactors = FALSE
        )

        ## Add a probe ID column:
        file_placeholder$probeID <- rep(
          sub(
            '_zscores.*',
            '',
            basename(
              file_path
            )
          ),
          nrow(
            file_placeholder
          )
        )

        ## Combine the file with the hyper dataset:
        hyper_z_scores <<- rbind(
          hyper_z_scores,
          file_placeholder
        )

        rm(file_placeholder)

      }
    }

    ## Now load the files and data into the previously
    ## created hyper_z_scores dataframe:
    parallel::mclapply(
      hypermeth_z_scores_files,
      hyper_file_combiner,
      mc.cores= core_count
    )

    ## Rename the columns for the time being:
    colnames(hyper_z_scores) <- c(
      'geneID',
      'Zscore',
      'probeID'
    )

    ## Reorganize the columns:
    hyper_z_scores <- hyper_z_scores[
      c(
        'probeID',
        'geneID',
        'Zscore'
      )
    ]

    ## Separate into two dataframes based on Z score:
    hyper.Gplus_z_scores <- hyper_z_scores[
      hyper_z_scores$Zscore<0,
    ]

    hyper.Gminus_z_scores <- hyper_z_scores[
      hyper_z_scores$Zscore>0,
    ]

    ## Save the two dataframes to the output folder:
    ## Write out the file:
    write.table(
      hyper.Gplus_z_scores,
      file= paste(
        TENET_directory,
        'step3/',
        'hyper_Gplus_link_zscores.txt',
        sep=''
      ),
      row.names=F,
      col.names=F,
      quote=F,
      sep="\t"
    )

    write.table(
      hyper.Gminus_z_scores,
      file= paste(
        TENET_directory,
        'step3/',
        'hyper_Gminus_link_zscores.txt',
        sep=''
      ),
      row.names=F,
      col.names=F,
      quote=F,
      sep="\t"
    )

  }

  ## If hypometh analysis is selected, load the significant
  ## hypometh z scores
  if(hypometh_pairs==TRUE){

    hypo_file_combiner <- function(file_path){

      ## If the file isn't empty, load it:
      if(
        (file.info(file_path)$size)>0
      ){

        ## Load the file:
        file_placeholder <- read.delim(
          file_path,
          header= FALSE,
          sep='\t',
          stringsAsFactors = FALSE
        )

        ## Add a probe ID column:
        file_placeholder$probeID <- rep(
          sub(
            '_zscores.*',
            '',
            basename(
              file_path
            )
          ),
          nrow(
            file_placeholder
          )
        )

        ## Combine the file with the hypo dataset:
        hypo_z_scores <<- rbind(
          hypo_z_scores,
          file_placeholder
        )

        rm(file_placeholder)

      }
    }

    ## Now load the files and data into the previously
    ## created hypo_z_scores dataframe:
    parallel::mclapply(
      hypometh_z_scores_files,
      hypo_file_combiner,
      mc.cores= core_count
    )

    ## Rename the columns for the time being:
    colnames(hypo_z_scores) <- c(
      'geneID',
      'Zscore',
      'probeID'
    )

    ## Reorganize the columns:
    hypo_z_scores <- hypo_z_scores[
      c(
        'probeID',
        'geneID',
        'Zscore'
      )
    ]

    ## Separate into two dataframes based on Z score:
    hypo.Gplus_z_scores <- hypo_z_scores[
      hypo_z_scores$Zscore<0,
    ]

    hypo.Gminus_z_scores <- hypo_z_scores[
      hypo_z_scores$Zscore>0,
    ]

    ## Save the two dataframes to the output folder:
    ## Write out the file:
    write.table(
      hypo.Gplus_z_scores,
      file= paste(
        TENET_directory,
        'step3/',
        'hypo_Gplus_link_zscores.txt',
        sep=''
      ),
      row.names=F,
      col.names=F,
      quote=F,
      sep="\t"
    )

    write.table(
      hypo.Gminus_z_scores,
      file= paste(
        TENET_directory,
        'step3/',
        'hypo_Gminus_link_zscores.txt',
        sep=''
      ),
      row.names=F,
      col.names=F,
      quote=F,
      sep="\t"
    )

  }

}

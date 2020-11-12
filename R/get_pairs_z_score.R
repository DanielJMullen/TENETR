#' get_pairs_z_score
#'
#' This is the step3 function of the TENETR package.
#' This function takes the identified hyper and hypomethylated probes
#' from step2 get_diffmeth_regions function and calculates a zscore comparing
#' the mean expression of each gene in groups that are hyper/hypomethylated for
#' each probe and those that are not, and permutates this across all
#' hyper/hypomethylated probes
#'
#'
#' @param TENET_directory Set a path to the directory that contains step2 results from get_diffmeth_regions function. This function will also create a new step3 folder there containing the results.
#' @param hypometh_pairs Set to TRUE/FALSE depending on if you want to calculate z-scores for hypomethylated probes.
#' @param hypermeth_pairs Set to TRUE/FALSE depending on if you want to calculate z-scores for hypermethylated probes.
#' @param usecaseonly Set to TRUE/FALSE depending on if you want to include the control/normal samples with the experimental/tumor samples when calculating hyper/hypomethylated groups and z-scores.
#' @param TF_only Set to TRUE/FALSE to determine if you only want to consider genes that are accepted transcription factors in The Human Transcription Factors by Lambert et al (2018) when calculating z-scores.
#' @param hypomethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypomethylated probes.
#' @param hypermethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypermethylated probes.
#' @param minExp Sets the minimum number of experimental/tumor samples to be considered for the hypo/hypermethylated groups.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns tab-delimited ".txt"zscore_all_genes_rda.txt" for each probe of the selected types analyzed, containing the zscore
#' @export

get_pairs_z_score <- function(
  TENET_directory,
  hypometh_pairs,
  hypermeth_pairs,
  usecaseonly,
  TF_only,
  hypomethcutoff,
  hypermethcutoff,
  minExp,
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
    'unmethDataT'
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

  ## Now lets write a function to calculate z scores for expression of each gene
  ## between groups that are hyper/hypomethylated for each probe, and those that
  ## aren't:
  z_score_calc <- function(
    DNA_methylation_probe,
    expData_internal,
    metData_internal,
    hyper_hypo,
    minExp_value,
    TF_only_internal
  ){

    ## Create dataframe with gene IDs and names:
    TESTSR=data.frame(
      probe= DNA_methylation_probe,
      geneID= rownames(expData_internal),
      stringsAsFactors = FALSE
    )
    rownames(TESTSR) <- TESTSR$geneID

    ## Subset to only TFs if analysis is selected:
    if(TF_only_internal==TRUE){

      Human_TF <- TENETR.data::human_transcription_factors_dataset[
        TENETR.data::human_transcription_factors_dataset$Is.TF.=='Yes',
        "Ensembl.ID"
      ]

      Human_TF_present <- Human_TF[
        Human_TF %in% rownames(TESTSR)
      ]

      TESTSR <- TESTSR[
        Human_TF_present,
      ]
    }

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

        P_gene_1 <- expData_internal[
          internal_gene_ID,
          names(T_P_probe_1)
        ]

        P_gene_0 <- expData_internal[
          internal_gene_ID,
          names(T_P_probe_0)
        ]

        z_score_internal <- (
          (mean(P_gene_0, na.rm=T) - mean(P_gene_1, na.rm=T))/
          (sd(P_gene_1, na.rm=T))
        )

        return(z_score_internal)
      }

      ## Use the function to calculate z-scores for all genes to the probe:
      TESTSR$Z.real <- sapply(
        TESTSR$geneID,
        internal_gene_function
      )

      ## Write out the file:
      write.table(
        TESTSR,
        file= paste(
          paste(
            TENET_directory,
            'step3/',
            hyper_hypo,
            sep=''
          ),
          unique(TESTSR$probe),
          "zscore_all_genes_rda.txt",
          sep ="_"
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

  ## Now let's sapply the function:
  parallel::mclapply(
    rownames(metDataHypo),
    z_score_calc,
    expData_internal= expData,
    metData_internal= metDataHypo,
    hyper_hypo= 'hypo',
    TF_only= TRUE,
    minExp_value= minExp,
    mc.cores= core_count
  )
}

#' step7_linked_probe_motif_searching
#'
#' This is a step7 function of the TENETR package.
#' This function takes a binding motif, supplied by the user, for a
#' TF with linked probes also specified by the user, and for each of those
#' linked probes, it identifies if a motif for that TF is found in
#' a user-specified vicinity of each probe linked to the TF.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_survival' containing the results of this function.
#' @param TF_gene Specify a transcription factor gene by its name or ENSG ID to perform motif searching for.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to do motif searching in the vicintiy of hypermeth probes with G+ links to the TF of interest, if specified.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to do motif searching in the vicintiy of hyometh probes with G- links to the TF of interest, if specified.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to do motif searching in the vicintiy of hypometh probes with G+ links to the TF of interest, if specified.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to do motif searching in the vicintiy of hypermeth probes with G- links to the TF of interest, if specified.
#' @param motif_PCM_PWM Specify a PCM or PWM in the form of a 4xN matrix with the motif for the specified TF_gene. If a PCM or PWM is not available, the user can select one by querying the MotifDb package database. See examples below for more info.
#' @param distance_from_probes Set a positive integer in base pairs to be the distance added to the beginning and from end of the the linked enhancer probes to perform motif searching on. Defaults to 100.
#' @param matchPWM_min_score Set a positive integer, from 1 to 100, to be the min.score passed to the matchPWM argument for motif searching. See ?matchPWM from the Biostrings package for more info. Defaults to 75%.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Returns a .pdf file containing a seqLogo for the specified PCM or PWM, as well as two tables. The first contains information about each motif found in the vicintiy of the linked DNA methylation probes for the specified TF_gene at the specified matchPWM_min_score. The second lists each probe linked to the TF_gene for the analysis types specified, and the number of motif occurrences found in the vicinity of each of those probes.
#' @export
#'
#' @examples
#' # Show available motifs for example TF FOXM1:
#' names(query(MotifDb,"FOXM1")
#' # Once you've selected a PWM for use, you can specify it using the following:
#' motif_PCM_PWM= query(MotifDb, motif_name)[[3]]

step7_linked_probe_motif_searching <- function(
  TENET_directory,
  TF_gene,
  DNA_methylation_manifest="HM450",
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  motif_PCM_PWM,
  distance_from_probes=100,
  matchPWM_min_score="75%",
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

  ## Create a seperate subdirectory to contain the results from the different
  ## motif gene analyses:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7/',
        'probe_motif_search',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        'probe_motif_search',
        sep=''
      )
    )

  }

  ## Now create a seperate subdirectory in the step7/probe_motif_search directory
  ## to contain the motif searching results as well as seqLogos generated for the motif
  ## for the gene of interest:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7/',
        'probe_motif_search/',
        TF_gene,
        '_',
        'linked_probe_motif_search',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        'probe_motif_search/',
        TF_gene,
        '_',
        'linked_probe_motif_search',
        sep=''
      )
    )
  }

  ## Get the dataset of gencode v22 genes:
  gencode_v22_gtf <- TENETR.data::gencode_v22_annotations

  ## Get info for genes
  gencode_v22_genes <- gencode_v22_gtf[
    gencode_v22_gtf$type=='gene',
  ]

  ## Remove the stuff after the period in the ENSG notations:
  gencode_v22_genes$ENSG_ID <- sub(
    '\\..*',
    '',
    gencode_v22_genes$gene_id
  )

  ## Load the hg38 DNA methylation annotations:
  ## Written in by Zexun Wu
  if (DNA_methylation_manifest == "HM450") {

    hg38_manifest_df <- TENETR.data::hm450_hg38_annotations

  } else if (DNA_methylation_manifest == "HM27") {

    hg38_manifest_df <- TENETR.data::hm27_hg38_annotations

  } else if (DNA_methylation_manifest == "EPIC") {

    hg38_manifest_df <- TENETR.data::epic_hg38_annotations

  } else {

    stop("The input for DNA_methylation_manifest is incorrect. Please select one from \"HM450\",\"HM27\",\"EPIC\"!")

  }

  ## Set the rownames of the loaded manifest to be the probeIDs:
  rownames(hg38_manifest_df) <- hg38_manifest_df$probeID

  ## Get the human genome:
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  ## Convert the user input for matchPWM_min_score to a character with a percentage
  ## sign after the number:
  if(
    substring(
      matchPWM_min_score,
      nchar(matchPWM_min_score),
      nchar(matchPWM_min_score)
    )!='%'
  ){

    matchPWM_min_score <- as.character(
      paste(
        matchPWM_min_score,
        '%',
        sep=''
      )
    )
  }

  ## Selecting a PWM to use in the analysis:
  ## First, check to see if the user has happened to supply their own PCM or PWM:
  if(
    'matrix' %in% class(motif_PCM_PWM)
  ){

    ## If they have, check if it's a PWM or PCM. If its the former, it will need
    ## to be converted to a PCM by multiplying the values by 100 and rounding to
    ## convert it to integer form:
    if(
      all(
        motif_PCM_PWM==floor(motif_PCM_PWM)
      )
    ){

      ## The supplied motif_preference is matrix with only integers, so we assume
      ## it is already in the form of a supplied PCM
      motif_PCM_PWM <- motif_PCM_PWM

    } else{

      ## Not all the values are integers, so multiple the values in the matrix by
      ## 100 and round them to convert it to a PCM to be used:

      ## However, we want to preserve the column sums as 100, so we need a special
      ## rounding function:
      round_preserve_sum <- function(
        x,
        digits
      ){

        exp <- 10^digits
        x <- x*exp
        y <- floor(x)
        index <- tail(
          order(x-y),
          round(
            sum(x)
          ) - sum(y)
        )
        y[index] <- y[index]+1

        y/exp

      }

      motif_PCM_PWM <- apply(
        motif_PCM_PWM,
        2,
        round_preserve_sum,
        digits= 2
      ) *100
    }

  } else{

    ## The input PCM or PWM is not in the format of a matrix. Cancel the function and alert the user:
    stop('Input PCM or PWM for the motif_PCM_PWM argument is not given as a matrix. Please specify a 4xN matrix for this argument.')
  }

  ## Now let's create a seqLogo of the motif PCM and save it for the user:

  ## Open a pdf:
  pdf(
    paste(
      TENET_directory,
      'step7/',
      'probe_motif_search/',
      TF_gene,
      '_',
      'linked_probe_motif_search/',
      paste(
        TF_gene,
        'seqLogo.pdf',
        sep='_'
      ),
      sep=''
    ),
    height= 3.5,
    width= 6
  )

  ## Plot the seqLogo:
  seqLogo::seqLogo(
    (motif_PCM_PWM/100)
  )

  ## Close the plot:
  dev.off()

  ## Ensure that the user input is capitalized:
  TF_gene_capitalized <- toupper(
    TF_gene
  )

  ## If the user has input a gene name, get the ENSG ID for it:
  if(
    substring(
      TF_gene_capitalized,
      1,
      4
    )=='ENSG'
  ){

    gene_ENSG <- TF_gene_capitalized

    TF_gene_capitalized <- gencode_v22_genes[
      which(gencode_v22_genes$ENSG_ID%in%TF_gene_capitalized),
      'gene_name'
    ]

  } else{

    gene_ENSG <- gencode_v22_genes[
      which(gencode_v22_genes$gene_name%in%TF_gene_capitalized),
      'ENSG_ID'
    ]
  }

  ## Now let's check the step5 results and find all the probes linked to
  ## the gene in each of the analysis quadrants:

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the TF specified:
    hyper_Gplus_probes_linked_to_TF <- hyper_Gplus_sig_link_zscores[
      hyper_Gplus_sig_link_zscores$geneID==gene_ENSG,
      'probeID'
    ]

    ## Create a data frame with the probes as the first (and only column for now):
    hyper_Gplus_probe_dataset_linked_to_TF <- data.frame(
      'probe_ID'= sort(hyper_Gplus_probes_linked_to_TF),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info:
    hyper_Gplus_probe_dataset_linked_to_TF$seqnames <- hg38_manifest_df[
      hyper_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_chrm'
    ]

    hyper_Gplus_probe_dataset_linked_to_TF$start <- hg38_manifest_df[
      hyper_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_beg'
    ]

    hyper_Gplus_probe_dataset_linked_to_TF$end <- hg38_manifest_df[
      hyper_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_end'
    ]

    ## Set the rownames to be the probe IDs:
    rownames(hyper_Gplus_probe_dataset_linked_to_TF) <- hyper_Gplus_probe_dataset_linked_to_TF$probe_ID

    ## TESTING ONLY, REMOVE LATER:
    # linked_probes <- rownames(hyper_Gplus_probe_dataset_linked_to_TF)[2]

    ## Write a function to find the motif occurences in the specified vicinity
    ## of each linked probe
    motif_finder_in_peaks <- function(
      linked_probes,
      df_or_count_output
    ){

      ## Get the probe name, chromsome, start, and end locations:
      probe <- linked_probes

      chr <- as.character(
        hyper_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          2
        ]
      )

      start <- as.numeric(
        hyper_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          3
        ]
      )

      end  <- as.numeric(
        hyper_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          4
        ]
      )

      ## Subtract and add the buffer specified by the user with the distance_from_probe
      ## argument from the start and end values:
      start_buffer <- start-distance_from_probe+1
      end_buffer <- end+distance_from_probe

      ## Get the DNA string sequence for that segment from the human genome:
      DNA_string <- genome[[chr]][
        start_buffer:end_buffer
      ]

      ## Do the motif search:
      search <- Biostrings::matchPWM(
        motif_PCM_PWM,
        DNA_string,
        matchPWM_min_score
      )

      if(df_or_count_output=='df'){

        ## Return a dataframe with the information about the found motifs:

        ## Get the start and end locations where the motif was found:
        starts <- start(search@ranges)
        ends <- end(search@ranges)

        ## Write an internal function to get the DNA string of the motif matches:
        internal_motif_grabber <- function(
          istart,
          iend
        ){

          return(
            as.character(
              DNA_string[istart:iend]
            )
          )
        }

        ## Return the motif matches:
        motif_matches <- mapply(
          internal_motif_grabber,
          starts,
          ends
        )

        ## Create a dataframe with results:
        return_df <- data.frame(
          'probe'= rep(
            probe,
            length(starts)
          ),
          'motif_chromosome'= rep(
            chr,
            length(starts)
          ),
          'motif_start_position'= (starts+start_buffer-1),
          'motif_end_position'= (ends+start_buffer-1),
          'motif_DNA_sequences'= motif_matches,
          stringsAsFactors = FALSE
        )

        ## Return the data frame:
        return(return_df)

        ## Clear the workspace:
        rm(return_df)

      } else if(df_or_count_output=='count'){

        ## Return the count of motifs found in the vicinity of the probe:
        ## Get the number of matches found
        if(length(as.character(search))==0){

          match_count <- 0

        } else{

          match_count <- length(search)
        }

        ## Return the number of matches
        return(match_count)

        ## Clear the memory:
        rm(match_count)

      }
    }

    ## Get the list of the probe info dfs:
    list_of_probe_info <- parallel::mclapply(
      rownames(hyper_Gplus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='df',
      mc.cores= core_count
    )

    ## Rebind them into a dataframe:
    hyper_Gplus_probe_motifs_df <- do.call(
      rbind,
      list_of_probe_info
    )

    ## Add a column noting the analysis type it came from:
    hyper_Gplus_probe_motifs_df$analysis_type <- rep(
      'hyper_Gplus',
      nrow(hyper_Gplus_probe_motifs_df)
    )

    ## Get a list of the counts of the number of motifs found
    ## per probe:
    count_of_motifs_per_probe_list  <- parallel::mclapply(
      rownames(hyper_Gplus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='count',
      mc.cores= core_count
    )

    ## Convert the list into a vector:
    count_of_motifs_per_probe <- unlist(
      c(
        count_of_motifs_per_probe_list
      )
    )

    ## Now let's convert the list into a data frame with the probes
    ## as the rownames and the analysis type:
    hyper_Gplus_motif_count_per_probe_df <- data.frame(
      'probe_ID'= rownames(
        hyper_Gplus_probe_dataset_linked_to_TF
      ),
      'motif_count'= count_of_motifs_per_probe,
      'analysis_type'= rep(
        'hyper_Gplus',
        length(count_of_motifs_per_probe)
      ),
      stringsAsFactors = FALSE
    )
  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the TF specified:
    hyper_Gminus_probes_linked_to_TF <- hyper_Gminus_sig_link_zscores[
      hyper_Gminus_sig_link_zscores$geneID==gene_ENSG,
      'probeID'
    ]

    ## Create a data frame with the probes as the first (and only column for now):
    hyper_Gminus_probe_dataset_linked_to_TF <- data.frame(
      'probe_ID'= sort(hyper_Gminus_probes_linked_to_TF),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info:
    hyper_Gminus_probe_dataset_linked_to_TF$seqnames <- hg38_manifest_df[
      hyper_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_chrm'
    ]

    hyper_Gminus_probe_dataset_linked_to_TF$start <- hg38_manifest_df[
      hyper_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_beg'
    ]

    hyper_Gminus_probe_dataset_linked_to_TF$end <- hg38_manifest_df[
      hyper_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_end'
    ]

    ## Set the rownames to be the probe IDs:
    rownames(hyper_Gminus_probe_dataset_linked_to_TF) <- hyper_Gminus_probe_dataset_linked_to_TF$probe_ID

    ## TESTING ONLY, REMOVE LATER:
    # linked_probes <- rownames(hyper_Gminus_probe_dataset_linked_to_TF)[2]

    ## Write a function to find the motif occurences in the specified vicinity
    ## of each linked probe
    motif_finder_in_peaks <- function(
      linked_probes,
      df_or_count_output
    ){

      ## Get the probe name, chromsome, start, and end locations:
      probe <- linked_probes

      chr <- as.character(
        hyper_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          2
        ]
      )

      start <- as.numeric(
        hyper_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          3
        ]
      )

      end  <- as.numeric(
        hyper_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          4
        ]
      )

      ## Subtract and add the buffer specified by the user with the distance_from_probe
      ## argument from the start and end values:
      start_buffer <- start-distance_from_probe+1
      end_buffer <- end+distance_from_probe

      ## Get the DNA string sequence for that segment from the human genome:
      DNA_string <- genome[[chr]][
        start_buffer:end_buffer
      ]

      ## Do the motif search:
      search <- Biostrings::matchPWM(
        motif_PCM_PWM,
        DNA_string,
        matchPWM_min_score
      )

      if(df_or_count_output=='df'){

        ## Return a dataframe with the information about the found motifs:

        ## Get the start and end locations where the motif was found:
        starts <- start(search@ranges)
        ends <- end(search@ranges)

        ## Write an internal function to get the DNA string of the motif matches:
        internal_motif_grabber <- function(
          istart,
          iend
        ){

          return(
            as.character(
              DNA_string[istart:iend]
            )
          )
        }

        ## Return the motif matches:
        motif_matches <- mapply(
          internal_motif_grabber,
          starts,
          ends
        )

        ## Create a dataframe with results:
        return_df <- data.frame(
          'probe'= rep(
            probe,
            length(starts)
          ),
          'motif_chromosome'= rep(
            chr,
            length(starts)
          ),
          'motif_start_position'= (starts+start_buffer-1),
          'motif_end_position'= (ends+start_buffer-1),
          'motif_DNA_sequences'= motif_matches,
          stringsAsFactors = FALSE
        )

        ## Return the data frame:
        return(return_df)

        ## Clear the workspace:
        rm(return_df)

      } else if(df_or_count_output=='count'){

        ## Return the count of motifs found in the vicinity of the probe:
        ## Get the number of matches found
        if(length(as.character(search))==0){

          match_count <- 0

        } else{

          match_count <- length(search)
        }

        ## Return the number of matches
        return(match_count)

        ## Clear the memory:
        rm(match_count)

      }
    }

    ## Get the list of the probe info dfs:
    list_of_probe_info <- parallel::mclapply(
      rownames(hyper_Gminus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='df',
      mc.cores= core_count
    )

    ## Rebind them into a dataframe:
    hyper_Gminus_probe_motifs_df <- do.call(
      rbind,
      list_of_probe_info
    )

    ## Add a column noting the analysis type it came from:
    hyper_Gminus_probe_motifs_df$analysis_type <- rep(
      'hyper_Gminus',
      nrow(hyper_Gminus_probe_motifs_df)
    )

    ## Get a list of the counts of the number of motifs found
    ## per probe:
    count_of_motifs_per_probe_list  <- parallel::mclapply(
      rownames(hyper_Gminus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='count',
      mc.cores= core_count
    )

    ## Convert the list into a vector:
    count_of_motifs_per_probe <- unlist(
      c(
        count_of_motifs_per_probe_list
      )
    )

    ## Now let's convert the list into a data frame with the probes
    ## as the rownames and the analysis type:
    hyper_Gminus_motif_count_per_probe_df <- data.frame(
      'probe_ID'= rownames(
        hyper_Gminus_probe_dataset_linked_to_TF
      ),
      'motif_count'= count_of_motifs_per_probe,
      'analysis_type'= rep(
        'hyper_Gminus',
        length(count_of_motifs_per_probe)
      ),
      stringsAsFactors = FALSE
    )
  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the TF specified:
    hypo_Gplus_probes_linked_to_TF <- hypo_Gplus_sig_link_zscores[
      hypo_Gplus_sig_link_zscores$geneID==gene_ENSG,
      'probeID'
    ]

    ## Create a data frame with the probes as the first (and only column for now):
    hypo_Gplus_probe_dataset_linked_to_TF <- data.frame(
      'probe_ID'= sort(hypo_Gplus_probes_linked_to_TF),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info:
    hypo_Gplus_probe_dataset_linked_to_TF$seqnames <- hg38_manifest_df[
      hypo_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_chrm'
    ]

    hypo_Gplus_probe_dataset_linked_to_TF$start <- hg38_manifest_df[
      hypo_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_beg'
    ]

    hypo_Gplus_probe_dataset_linked_to_TF$end <- hg38_manifest_df[
      hypo_Gplus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_end'
    ]

    ## Set the rownames to be the probe IDs:
    rownames(hypo_Gplus_probe_dataset_linked_to_TF) <- hypo_Gplus_probe_dataset_linked_to_TF$probe_ID

    ## TESTING ONLY, REMOVE LATER:
    # linked_probes <- rownames(hypo_Gplus_probe_dataset_linked_to_TF)[2]

    ## Write a function to find the motif occurences in the specified vicinity
    ## of each linked probe
    motif_finder_in_peaks <- function(
      linked_probes,
      df_or_count_output
    ){

      ## Get the probe name, chromsome, start, and end locations:
      probe <- linked_probes

      chr <- as.character(
        hypo_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          2
        ]
      )

      start <- as.numeric(
        hypo_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          3
        ]
      )

      end  <- as.numeric(
        hypo_Gplus_probe_dataset_linked_to_TF[
          linked_probes,
          4
        ]
      )

      ## Subtract and add the buffer specified by the user with the distance_from_probe
      ## argument from the start and end values:
      start_buffer <- start-distance_from_probe+1
      end_buffer <- end+distance_from_probe

      ## Get the DNA string sequence for that segment from the human genome:
      DNA_string <- genome[[chr]][
        start_buffer:end_buffer
      ]

      ## Do the motif search:
      search <- Biostrings::matchPWM(
        motif_PCM_PWM,
        DNA_string,
        matchPWM_min_score
      )

      if(df_or_count_output=='df'){

        ## Return a dataframe with the information about the found motifs:

        ## Get the start and end locations where the motif was found:
        starts <- start(search@ranges)
        ends <- end(search@ranges)

        ## Write an internal function to get the DNA string of the motif matches:
        internal_motif_grabber <- function(
          istart,
          iend
        ){

          return(
            as.character(
              DNA_string[istart:iend]
            )
          )
        }

        ## Return the motif matches:
        motif_matches <- mapply(
          internal_motif_grabber,
          starts,
          ends
        )

        ## Create a dataframe with results:
        return_df <- data.frame(
          'probe'= rep(
            probe,
            length(starts)
          ),
          'motif_chromosome'= rep(
            chr,
            length(starts)
          ),
          'motif_start_position'= (starts+start_buffer-1),
          'motif_end_position'= (ends+start_buffer-1),
          'motif_DNA_sequences'= motif_matches,
          stringsAsFactors = FALSE
        )

        ## Return the data frame:
        return(return_df)

        ## Clear the workspace:
        rm(return_df)

      } else if(df_or_count_output=='count'){

        ## Return the count of motifs found in the vicinity of the probe:
        ## Get the number of matches found
        if(length(as.character(search))==0){

          match_count <- 0

        } else{

          match_count <- length(search)
        }

        ## Return the number of matches
        return(match_count)

        ## Clear the memory:
        rm(match_count)

      }
    }

    ## Get the list of the probe info dfs:
    list_of_probe_info <- parallel::mclapply(
      rownames(hypo_Gplus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='df',
      mc.cores= core_count
    )

    ## Rebind them into a dataframe:
    hypo_Gplus_probe_motifs_df <- do.call(
      rbind,
      list_of_probe_info
    )

    ## Add a column noting the analysis type it came from:
    hypo_Gplus_probe_motifs_df$analysis_type <- rep(
      'hypo_Gplus',
      nrow(hypo_Gplus_probe_motifs_df)
    )

    ## Get a list of the counts of the number of motifs found
    ## per probe:
    count_of_motifs_per_probe_list  <- parallel::mclapply(
      rownames(hypo_Gplus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='count',
      mc.cores= core_count
    )

    ## Convert the list into a vector:
    count_of_motifs_per_probe <- unlist(
      c(
        count_of_motifs_per_probe_list
      )
    )

    ## Now let's convert the list into a data frame with the probes
    ## as the rownames and the analysis type:
    hypo_Gplus_motif_count_per_probe_df <- data.frame(
      'probe_ID'= rownames(
        hypo_Gplus_probe_dataset_linked_to_TF
      ),
      'motif_count'= count_of_motifs_per_probe,
      'analysis_type'= rep(
        'hypo_Gplus',
        length(count_of_motifs_per_probe)
      ),
      stringsAsFactors = FALSE
    )
  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the TF specified:
    hypo_Gminus_probes_linked_to_TF <- hypo_Gminus_sig_link_zscores[
      hypo_Gminus_sig_link_zscores$geneID==gene_ENSG,
      'probeID'
    ]

    ## Create a data frame with the probes as the first (and only column for now):
    hypo_Gminus_probe_dataset_linked_to_TF <- data.frame(
      'probe_ID'= sort(hypo_Gminus_probes_linked_to_TF),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info:
    hypo_Gminus_probe_dataset_linked_to_TF$seqnames <- hg38_manifest_df[
      hypo_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_chrm'
    ]

    hypo_Gminus_probe_dataset_linked_to_TF$start <- hg38_manifest_df[
      hypo_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_beg'
    ]

    hypo_Gminus_probe_dataset_linked_to_TF$end <- hg38_manifest_df[
      hypo_Gminus_probe_dataset_linked_to_TF$probe_ID,
      'CpG_end'
    ]

    ## Set the rownames to be the probe IDs:
    rownames(hypo_Gminus_probe_dataset_linked_to_TF) <- hypo_Gminus_probe_dataset_linked_to_TF$probe_ID

    ## TESTING ONLY, REMOVE LATER:
    # linked_probes <- rownames(hypo_Gminus_probe_dataset_linked_to_TF)[2]

    ## Write a function to find the motif occurences in the specified vicinity
    ## of each linked probe
    motif_finder_in_peaks <- function(
      linked_probes,
      df_or_count_output
    ){

      ## Get the probe name, chromsome, start, and end locations:
      probe <- linked_probes

      chr <- as.character(
        hypo_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          2
        ]
      )

      start <- as.numeric(
        hypo_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          3
        ]
      )

      end  <- as.numeric(
        hypo_Gminus_probe_dataset_linked_to_TF[
          linked_probes,
          4
        ]
      )

      ## Subtract and add the buffer specified by the user with the distance_from_probe
      ## argument from the start and end values:
      start_buffer <- start-distance_from_probe+1
      end_buffer <- end+distance_from_probe

      ## Get the DNA string sequence for that segment from the human genome:
      DNA_string <- genome[[chr]][
        start_buffer:end_buffer
      ]

      ## Do the motif search:
      search <- Biostrings::matchPWM(
        motif_PCM_PWM,
        DNA_string,
        matchPWM_min_score
      )

      if(df_or_count_output=='df'){

        ## Return a dataframe with the information about the found motifs:

        ## Get the start and end locations where the motif was found:
        starts <- start(search@ranges)
        ends <- end(search@ranges)

        ## Write an internal function to get the DNA string of the motif matches:
        internal_motif_grabber <- function(
          istart,
          iend
        ){

          return(
            as.character(
              DNA_string[istart:iend]
            )
          )
        }

        ## Return the motif matches:
        motif_matches <- mapply(
          internal_motif_grabber,
          starts,
          ends
        )

        ## Create a dataframe with results:
        return_df <- data.frame(
          'probe'= rep(
            probe,
            length(starts)
          ),
          'motif_chromosome'= rep(
            chr,
            length(starts)
          ),
          'motif_start_position'= (starts+start_buffer-1),
          'motif_end_position'= (ends+start_buffer-1),
          'motif_DNA_sequences'= motif_matches,
          stringsAsFactors = FALSE
        )

        ## Return the data frame:
        return(return_df)

        ## Clear the workspace:
        rm(return_df)

      } else if(df_or_count_output=='count'){

        ## Return the count of motifs found in the vicinity of the probe:
        ## Get the number of matches found
        if(length(as.character(search))==0){

          match_count <- 0

        } else{

          match_count <- length(search)
        }

        ## Return the number of matches
        return(match_count)

        ## Clear the memory:
        rm(match_count)

      }
    }

    ## Get the list of the probe info dfs:
    list_of_probe_info <- parallel::mclapply(
      rownames(hypo_Gminus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='df',
      mc.cores= core_count
    )

    ## Rebind them into a dataframe:
    hypo_Gminus_probe_motifs_df <- do.call(
      rbind,
      list_of_probe_info
    )

    ## Add a column noting the analysis type it came from:
    hypo_Gminus_probe_motifs_df$analysis_type <- rep(
      'hypo_Gminus',
      nrow(hypo_Gminus_probe_motifs_df)
    )

    ## Get a list of the counts of the number of motifs found
    ## per probe:
    count_of_motifs_per_probe_list  <- parallel::mclapply(
      rownames(hypo_Gminus_probe_dataset_linked_to_TF),
      motif_finder_in_peaks,
      df_or_count_output='count',
      mc.cores= core_count
    )

    ## Convert the list into a vector:
    count_of_motifs_per_probe <- unlist(
      c(
        count_of_motifs_per_probe_list
      )
    )

    ## Now let's convert the list into a data frame:
    hypo_Gminus_motif_count_per_probe_df <- data.frame(
      'probe_ID'= rownames(
        hypo_Gminus_probe_dataset_linked_to_TF
      ),
      'motif_count'= count_of_motifs_per_probe,
      'analysis_type'= rep(
        'hypo_Gminus',
        length(count_of_motifs_per_probe)
      ),
      stringsAsFactors = FALSE
    )
  }

  ## Now for each of the analysis types that is selected, let's assemble the
  ## data frames of the complete motif occurrences and counts:

  ## First let's index empty dataframes with the same columns to rbind to:
  probe_motifs_df <- data.frame(
    'probe'= character(),
    'motif_chromosome'= character(),
    'motif_start_position'= numeric(),
    'motif_end_position'= numeric(),
    'motif_DNA_sequences'= character(),
    'analysis_type'= character(),
    stringsAsFactors = FALSE
  )

  motif_counts_df <- data.frame(
    'probe_name'= character(),
    'motif_count'= numeric(),
    'analysis_type'= character(),
    stringsAsFactors = FALSE
  )

  ## Add the hypermeth Gplus results if specified:
  if(
    hypermeth_Gplus_analysis==TRUE
  ){

    probe_motifs_df <- rbind(
      probe_motifs_df,
      hyper_Gplus_probe_motifs_df
    )

    motif_counts_df <- rbind(
      motif_counts_df,
      hyper_Gplus_motif_count_per_probe_df
    )
  }

  ## Add the hypermeth Gminus results if specified:
  if(
    hypermeth_Gminus_analysis==TRUE
  ){

    probe_motifs_df <- rbind(
      probe_motifs_df,
      hyper_Gminus_probe_motifs_df
    )

    motif_counts_df <- rbind(
      motif_counts_df,
      hyper_Gminus_motif_count_per_probe_df
    )
  }

  ## Add the hypometh Gplus results if specified:
  if(
    hypometh_Gplus_analysis==TRUE
  ){

    probe_motifs_df <- rbind(
      probe_motifs_df,
      hypo_Gplus_probe_motifs_df
    )

    motif_counts_df <- rbind(
      motif_counts_df,
      hypo_Gplus_motif_count_per_probe_df
    )
  }

  ## Add the hypometh Gminus results if specified:
  if(
    hypometh_Gminus_analysis==TRUE
  ){

    probe_motifs_df <- rbind(
      probe_motifs_df,
      hypo_Gminus_probe_motifs_df
    )

    motif_counts_df <- rbind(
      motif_counts_df,
      hypo_Gminus_motif_count_per_probe_df
    )
  }

  ## Now let's save those analyses as well:
  write.table(
    probe_motifs_df,
    file= paste(
      TENET_directory,
      'step7/',
      'probe_motif_search/',
      TF_gene,
      '_',
      'linked_probe_motif_search/',
      paste(
        TF_gene,
        'probe_motif_occurences_table.tsv',
        sep='_'
      ),
      sep=''
    ),
    quote= FALSE,
    sep= '\t',
    row.names= FALSE
  )

  write.table(
    motif_counts_df,
    file= paste(
      TENET_directory,
      'step7/',
      'probe_motif_search/',
      TF_gene,
      '_',
      'linked_probe_motif_search/',
      paste(
        TF_gene,
        'total_motif_occurences_per_probe_count.tsv',
        sep='_'
      ),
      sep=''
    ),
    quote= FALSE,
    sep= '\t',
    row.names= FALSE
  )
}

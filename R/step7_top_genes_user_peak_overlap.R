#' step7_top_genes_user_peak_overlap
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and identifies if the probes linked to those TFs are found in peaks for factors of
#' interest, supplied by the user in the form of .bed, .narrowPeak, .gappedPeak,
#' or .broadPeak files in a user-specified directory.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_user_peak_overlap' containing the results of this function.
#' @param ext_peaks_directory If the user has their own datasets with peaks for factors interest, set this as a path to a directory containing either .bed, .narrowPeak, .broadPeak, or .gappedPeak files with these peaks.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create survival plots for the top genes/TFs by most hypermeth probes with G+ links, and these linked DNA methyation probes if specified.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypermeth probes with G- links, and these linked DNA methyation probes if specified.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G+ links, as well as their linked DNA methyation probes if specified.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G- links, as well as their linked DNA methyation probes if specified.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes to generate files showing overlap with user peak files for the enhancer probes linked to those genes.
#' @param distance_from_probes Set a positive integer in base pairs to be the distance added to the beginning and from end of the the linked enhancer probes to perform peak searching on. Defaults to 100.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Returns a .txt file for each of the top genes/TFs showing which of the user-specified peak files the enhancer DNA methylation probes linked to each of the genes overlap with.
#' @export

step7_top_genes_user_peak_overlap <- function(
  TENET_directory,
  ext_peaks_directory,
  DNA_methylation_manifest="HM450",
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  distance_from_probes=100,
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

  ## Create a subdirectory in the new step7 directory to contain the
  ## peak overlap files:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7/',
        'user_peak_overlap',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        'user_peak_overlap',
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

  ## Create a modified dataframe of the hg38 DNA methylation probe annotations
  ## to later convert to granges:
  hg38_manifest_granges_df <- data.frame(
    'chr'= hg38_manifest_df$CpG_chrm,
    'start'= hg38_manifest_df$CpG_beg,
    'end'= hg38_manifest_df$CpG_end,
    'strand'= rep(
      '*',
      nrow(hg38_manifest_df)
    ),
    'names' = hg38_manifest_df$probeID,
    stringsAsFactors = FALSE
  )

  ## Set the rownames of the hg38_manifest_granges_df
  rownames(hg38_manifest_granges_df) <- rownames(hg38_manifest_df)

  ## Remove the big manifest dataframe:
  rm(hg38_manifest_df)

  ## Remove the probes that have NA values:
  hg38_manifest_no_NA_granges_df <- hg38_manifest_granges_df[
    !is.na(hg38_manifest_granges_df$chr),
  ]

  ## Remove the dataset with NA probes:
  rm(hg38_manifest_granges_df)

  ## Create a granges object from the new
  ## hg38 annotations with no NA df:
  hg38_manifest_annotations_granges <- GenomicRanges::makeGRangesFromDataFrame(
    df= hg38_manifest_no_NA_granges_df,
    keep.extra.columns = FALSE,
    starts.in.df.are.0based = TRUE
  )
  names(hg38_manifest_annotations_granges) <- hg38_manifest_no_NA_granges_df$names

  ## Check to see if the user has supplied ext_peaks_directory files:

  ## First check if the user has specified a external peaks directory:
  if(
    !dir.exists(
      ext_peaks_directory
    )
  ){

    ## No ext_peaks_directory found. Quit the function and note that to the user:
    stop("The specified ext_peaks_directory was not found. Please create such a folder containing bed-like files with .bed, .narrowPeak, .gappedPeak, or .broadPeak extensions containing peaks of interest and provide a path to that directory as the argument")
  }

  ## List all the ENH bed files found:
  peak_bed_files <- list.files(
    path= ext_peaks_directory,
    pattern = ".bed",
    full.names = TRUE
  )

  ## List all the ENH narrowPeak files:
  peak_narrowPeak_files <- list.files(
    path= ext_peaks_directory,
    pattern = ".narrowPeak",
    full.names = TRUE
  )

  ## List all the ENH broadPeak files:
  peak_broadPeak_files <- list.files(
    path= ext_peaks_directory,
    pattern = ".broadPeak",
    full.names = TRUE
  )

  ## List all the ENH gappedPeak files:
  peak_gappedPeak_files <- list.files(
    path= ext_peaks_directory,
    pattern = ".gappedPeak",
    full.names = TRUE
  )

  ## Combine the different file types into
  ## an ENH master list to load:
  peak_file_list <- c(
    peak_bed_files,
    peak_narrowPeak_files,
    peak_broadPeak_files,
    peak_gappedPeak_files
  )

  ## For each of the loaded files, create a granges object for it:
  peak_file_list_peaks_granges <- list()

  for(i in 1:length(peak_file_list)){

    ## Load the peaks in the file:
    peaks <- read.delim(
      peak_file_list[[i]],
      header= FALSE,
      sep='\t',
      stringsAsFactors = FALSE
    )

    ## change the column names for the first three columns in the peak file:
    colnames(peaks)[1:3] <- c(
      'chr',
      'start',
      'end'
    )

    ## Create a granges object:
    ## Asssume starts are 0-based
    peaks_granges <- GenomicRanges::makeGRangesFromDataFrame(
      df= peaks,
      keep.extra.columns = FALSE,
      starts.in.df.are.0based = TRUE
    )

    ## Return that granges to the new list:
    peak_file_list_peaks_granges[[i]] <- peaks_granges

    ## Remove the uneeded datasets:
    rm(peaks)
    rm(peaks_granges)
  }

  ## Add names to the peaks granges list as the names of the files:
  names(peak_file_list_peaks_granges) <- basename(
    unlist(
      c(
        peak_file_list
      )
    )
  )

  ## Make sure that at least one peak file was found:
  if(
    length(peak_file_list)==0
  ){

    ## No .bed, .narrowPeak, .gappedPeak, or .broadPeak files found. Quit the function and note that to the user:
    stop("No .bed, .narrowPeak, .gappedPeak, or .broadPeak files were found in the specified ext_peaks_directory. Please populate the specified ext_peaks_directory with .bed, .narrowPeak, .gappedPeak, or .broadPeak files containing peaks for factors of interest to overlap with linked enhancer DNA methylation probes and rerun the function.")
  }

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus user peak overlap files:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files',
          sep=''
        )
      )
    }

    ## Create further subdirectories to contain the top_gene and top_TF files
    ## separately

    ## Top genes:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    }

    ## Top TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gplus_user_peak_overlap_files/',
          'top_TFs',
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

    ## Check that the hyper_Gplus_links_all_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_gene_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gplus_all_gene_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_gene_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Check that the hyper_Gplus_links_all_TF_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_TF_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gplus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_TF_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Now that the files have been loaded, get the names of the top genes:
    if(
      nrow(hyper_Gplus_all_gene_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hyper_Gplus_all_gene_ENSG <- hyper_Gplus_all_gene_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hyper_Gplus_all_gene_ENSG <- hyper_Gplus_all_gene_freq$geneID

    }

    ## Do the same for the top TF genes only:
    if(
      nrow(hyper_Gplus_all_TF_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hyper_Gplus_all_TF_ENSG <- hyper_Gplus_all_TF_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hyper_Gplus_all_TF_ENSG <- hyper_Gplus_all_TF_freq$geneID

    }

    ## For each of the genes, let's get the probes linked to it, overlap them
    ## with each of the files, then create files noting which datasets the
    ## probes overlapped with:
    internal_dataset_output_function <- function(
      gene_ID,
      gene_or_TF
    ){

      ## Get the probes linked to the gene ID input:
      gene_CpGs_linked <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID==gene_ID,
          'probeID'
        ]
      )

      ## Create dataframes with a row for each of the CpGs in the datasets:
      gene_CpGs_linked_df <- data.frame(
        'chr'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'chr'
        ],
        'start'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'start'
        ]-distance_from_probes,
        'end'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'end'
        ]+distance_from_probes,
        stringsAsFactors = FALSE
      )
      rownames(gene_CpGs_linked_df) <- gene_CpGs_linked


      ## Let's create a granges copy of this dataset:
      gene_CpGs_linked_df_granges <- GenomicRanges::makeGRangesFromDataFrame(
        df= gene_CpGs_linked_df,
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = TRUE
      )

      ## Now for each peak dataset, find the probes that overlap with it:
      for(i in 1:length(peak_file_list_peaks_granges)){

        probe_overlap <- rownames(
          gene_CpGs_linked_df[
            unique(
              S4Vectors::queryHits(
                GenomicRanges::findOverlaps(
                  gene_CpGs_linked_df_granges,
                  peak_file_list_peaks_granges[[i]]
                )
              )
            ),
          ]
        )

        ## Now create a new column in gene_CpGs_linked_df that indicates for the
        ## dataset if each probe overlapped with it or not:
        gene_CpGs_linked_df[[i+3]] <- rownames(
          gene_CpGs_linked_df
        ) %in% probe_overlap

        ## Clear the workspace:
        rm(probe_overlap)

      }

      ## Now let's rename the overlap columns after the dataset that called them:
      colnames(gene_CpGs_linked_df) <- c(
        colnames(gene_CpGs_linked_df)[1:3],
        names(peak_file_list_peaks_granges)
      )

      ## Now let's save the file under the name of the gene:

      ## Get the name of the gene:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ENSG_ID==gene_ID,
        'gene_name'
      ]

      ## Save the data_frame:
      if(gene_or_TF=='gene'){

        ## Save the file as the top genes:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hyper_Gplus_user_peak_overlap_files/',
            'top_genes/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )

      } else if(gene_or_TF=='TF'){

        ## Save the file as the top TFs:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hyper_Gplus_user_peak_overlap_files/',
            'top_TFs/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )
      }
    }

    ## Now sapply the function on the hyper Gplus top genes:
    sapply(
      top_hyper_Gplus_all_gene_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='gene'
    )

    ## Now sapply the function on the hyper Gplus top TFs:
    sapply(
      top_hyper_Gplus_all_TF_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='TF'
    )
  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gminus user peak overlap files:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files',
          sep=''
        )
      )
    }

    ## Create further subdirectories to contain the top_gene and top_TF files
    ## separately

    ## Top genes:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    }

    ## Top TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hyper_Gminus_user_peak_overlap_files/',
          'top_TFs',
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

    ## Check that the hyper_Gminus_links_all_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_gene_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gminus_all_gene_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_gene_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Check that the hyper_Gminus_links_all_TF_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_TF_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gminus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_TF_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Now that the files have been loaded, get the names of the top genes:
    if(
      nrow(hyper_Gminus_all_gene_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hyper_Gminus_all_gene_ENSG <- hyper_Gminus_all_gene_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hyper_Gminus_all_gene_ENSG <- hyper_Gminus_all_gene_freq$geneID

    }

    ## Do the same for the top TF genes only:
    if(
      nrow(hyper_Gminus_all_TF_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hyper_Gminus_all_TF_ENSG <- hyper_Gminus_all_TF_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hyper_Gminus_all_TF_ENSG <- hyper_Gminus_all_TF_freq$geneID

    }

    ## For each of the genes, let's get the probes linked to it, overlap them
    ## with each of the files, then create files noting which datasets the
    ## probes overlapped with:
    internal_dataset_output_function <- function(
      gene_ID,
      gene_or_TF
    ){

      ## Get the probes linked to the gene ID input:
      gene_CpGs_linked <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID==gene_ID,
          'probeID'
        ]
      )

      ## Create dataframes with a row for each of the CpGs in the datasets:
      gene_CpGs_linked_df <- data.frame(
        'chr'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'chr'
        ],
        'start'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'start'
        ]-distance_from_probes,
        'end'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'end'
        ]+distance_from_probes,
        stringsAsFactors = FALSE
      )
      rownames(gene_CpGs_linked_df) <- gene_CpGs_linked


      ## Let's create a granges copy of this dataset:
      gene_CpGs_linked_df_granges <- GenomicRanges::makeGRangesFromDataFrame(
        df= gene_CpGs_linked_df,
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = TRUE
      )

      ## Now for each peak dataset, find the probes that overlap with it:
      for(i in 1:length(peak_file_list_peaks_granges)){

        probe_overlap <- rownames(
          gene_CpGs_linked_df[
            unique(
              S4Vectors::queryHits(
                GenomicRanges::findOverlaps(
                  gene_CpGs_linked_df_granges,
                  peak_file_list_peaks_granges[[i]]
                )
              )
            ),
          ]
        )

        ## Now create a new column in gene_CpGs_linked_df that indicates for the
        ## dataset if each probe overlapped with it or not:
        gene_CpGs_linked_df[[i+3]] <- rownames(
          gene_CpGs_linked_df
        ) %in% probe_overlap

        ## Clear the workspace:
        rm(probe_overlap)

      }

      ## Now let's rename the overlap columns after the dataset that called them:
      colnames(gene_CpGs_linked_df) <- c(
        colnames(gene_CpGs_linked_df)[1:3],
        names(peak_file_list_peaks_granges)
      )

      ## Now let's save the file under the name of the gene:

      ## Get the name of the gene:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ENSG_ID==gene_ID,
        'gene_name'
      ]

      ## Save the data_frame:
      if(gene_or_TF=='gene'){

        ## Save the file as the top genes:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hyper_Gminus_user_peak_overlap_files/',
            'top_genes/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )

      } else if(gene_or_TF=='TF'){

        ## Save the file as the top TFs:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hyper_Gminus_user_peak_overlap_files/',
            'top_TFs/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )
      }
    }

    ## Now sapply the function on the hyper Gminus top genes:
    sapply(
      top_hyper_Gminus_all_gene_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='gene'
    )

    ## Now sapply the function on the hyper Gminus top TFs:
    sapply(
      top_hyper_Gminus_all_TF_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='TF'
    )
  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gplus user peak overlap files:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files',
          sep=''
        )
      )
    }

    ## Create further subdirectories to contain the top_gene and top_TF files
    ## separately

    ## Top genes:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    }

    ## Top TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gplus_user_peak_overlap_files/',
          'top_TFs',
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

    ## Check that the hypo_Gplus_links_all_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_gene_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gplus_all_gene_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_gene_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Check that the hypo_Gplus_links_all_TF_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_TF_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gplus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_TF_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Now that the files have been loaded, get the names of the top genes:
    if(
      nrow(hypo_Gplus_all_gene_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hypo_Gplus_all_gene_ENSG <- hypo_Gplus_all_gene_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hypo_Gplus_all_gene_ENSG <- hypo_Gplus_all_gene_freq$geneID

    }

    ## Do the same for the top TF genes only:
    if(
      nrow(hypo_Gplus_all_TF_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hypo_Gplus_all_TF_ENSG <- hypo_Gplus_all_TF_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hypo_Gplus_all_TF_ENSG <- hypo_Gplus_all_TF_freq$geneID

    }

    ## For each of the genes, let's get the probes linked to it, overlap them
    ## with each of the files, then create files noting which datasets the
    ## probes overlapped with:
    internal_dataset_output_function <- function(
      gene_ID,
      gene_or_TF
    ){

      ## Get the probes linked to the gene ID input:
      gene_CpGs_linked <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID==gene_ID,
          'probeID'
        ]
      )

      ## Create dataframes with a row for each of the CpGs in the datasets:
      gene_CpGs_linked_df <- data.frame(
        'chr'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'chr'
        ],
        'start'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'start'
        ]-distance_from_probes,
        'end'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'end'
        ]+distance_from_probes,
        stringsAsFactors = FALSE
      )
      rownames(gene_CpGs_linked_df) <- gene_CpGs_linked


      ## Let's create a granges copy of this dataset:
      gene_CpGs_linked_df_granges <- GenomicRanges::makeGRangesFromDataFrame(
        df= gene_CpGs_linked_df,
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = TRUE
      )

      ## Now for each peak dataset, find the probes that overlap with it:
      for(i in 1:length(peak_file_list_peaks_granges)){

        probe_overlap <- rownames(
          gene_CpGs_linked_df[
            unique(
              S4Vectors::queryHits(
                GenomicRanges::findOverlaps(
                  gene_CpGs_linked_df_granges,
                  peak_file_list_peaks_granges[[i]]
                )
              )
            ),
          ]
        )

        ## Now create a new column in gene_CpGs_linked_df that indicates for the
        ## dataset if each probe overlapped with it or not:
        gene_CpGs_linked_df[[i+3]] <- rownames(
          gene_CpGs_linked_df
        ) %in% probe_overlap

        ## Clear the workspace:
        rm(probe_overlap)

      }

      ## Now let's rename the overlap columns after the dataset that called them:
      colnames(gene_CpGs_linked_df) <- c(
        colnames(gene_CpGs_linked_df)[1:3],
        names(peak_file_list_peaks_granges)
      )

      ## Now let's save the file under the name of the gene:

      ## Get the name of the gene:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ENSG_ID==gene_ID,
        'gene_name'
      ]

      ## Save the data_frame:
      if(gene_or_TF=='gene'){

        ## Save the file as the top genes:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hypo_Gplus_user_peak_overlap_files/',
            'top_genes/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )

      } else if(gene_or_TF=='TF'){

        ## Save the file as the top TFs:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hypo_Gplus_user_peak_overlap_files/',
            'top_TFs/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )
      }
    }

    ## Now sapply the function on the hypo Gplus top genes:
    sapply(
      top_hypo_Gplus_all_gene_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='gene'
    )

    ## Now sapply the function on the hypo Gplus top TFs:
    sapply(
      top_hypo_Gplus_all_TF_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='TF'
    )
  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gminus user peak overlap files:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files',
          sep=''
        )
      )
    }

    ## Create further subdirectories to contain the top_gene and top_TF files
    ## separately

    ## Top genes:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files/',
          'top_genes',
          sep=''
        )
      )
    }

    ## Top TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'user_peak_overlap/',
          'hypo_Gminus_user_peak_overlap_files/',
          'top_TFs',
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

    ## Check that the hypo_Gminus_links_all_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_gene_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gminus_all_gene_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_gene_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Check that the hypo_Gminus_links_all_TF_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_TF_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gminus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_TF_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation function.')

    }

    ## Now that the files have been loaded, get the names of the top genes:
    if(
      nrow(hypo_Gminus_all_gene_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hypo_Gminus_all_gene_ENSG <- hypo_Gminus_all_gene_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hypo_Gminus_all_gene_ENSG <- hypo_Gminus_all_gene_freq$geneID

    }

    ## Do the same for the top TF genes only:
    if(
      nrow(hypo_Gminus_all_TF_freq)>=top_gene_number
    ){

      ## If there are more genes than the number specified by the user
      ## with top_gene_number, get just those top n genes:
      top_hypo_Gminus_all_TF_ENSG <- hypo_Gminus_all_TF_freq[
        c(1:top_gene_number),
        'geneID'
      ]

    } else{

      ## If there are less genes than the number specified by the user
      ## with top_gene_number, get just all the genes available:
      top_hypo_Gminus_all_TF_ENSG <- hypo_Gminus_all_TF_freq$geneID

    }

    ## For each of the genes, let's get the probes linked to it, overlap them
    ## with each of the files, then create files noting which datasets the
    ## probes overlapped with:
    internal_dataset_output_function <- function(
      gene_ID,
      gene_or_TF
    ){

      ## Get the probes linked to the gene ID input:
      gene_CpGs_linked <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID==gene_ID,
          'probeID'
        ]
      )

      ## Create dataframes with a row for each of the CpGs in the datasets:
      gene_CpGs_linked_df <- data.frame(
        'chr'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'chr'
        ],
        'start'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'start'
        ]-distance_from_probes,
        'end'= hg38_manifest_no_NA_granges_df[
          gene_CpGs_linked,
          'end'
        ]+distance_from_probes,
        stringsAsFactors = FALSE
      )
      rownames(gene_CpGs_linked_df) <- gene_CpGs_linked


      ## Let's create a granges copy of this dataset:
      gene_CpGs_linked_df_granges <- GenomicRanges::makeGRangesFromDataFrame(
        df= gene_CpGs_linked_df,
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = TRUE
      )

      ## Now for each peak dataset, find the probes that overlap with it:
      for(i in 1:length(peak_file_list_peaks_granges)){

        probe_overlap <- rownames(
          gene_CpGs_linked_df[
            unique(
              S4Vectors::queryHits(
                GenomicRanges::findOverlaps(
                  gene_CpGs_linked_df_granges,
                  peak_file_list_peaks_granges[[i]]
                )
              )
            ),
          ]
        )

        ## Now create a new column in gene_CpGs_linked_df that indicates for the
        ## dataset if each probe overlapped with it or not:
        gene_CpGs_linked_df[[i+3]] <- rownames(
          gene_CpGs_linked_df
        ) %in% probe_overlap

        ## Clear the workspace:
        rm(probe_overlap)

      }

      ## Now let's rename the overlap columns after the dataset that called them:
      colnames(gene_CpGs_linked_df) <- c(
        colnames(gene_CpGs_linked_df)[1:3],
        names(peak_file_list_peaks_granges)
      )

      ## Now let's save the file under the name of the gene:

      ## Get the name of the gene:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ENSG_ID==gene_ID,
        'gene_name'
      ]

      ## Save the data_frame:
      if(gene_or_TF=='gene'){

        ## Save the file as the top genes:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hypo_Gminus_user_peak_overlap_files/',
            'top_genes/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )

      } else if(gene_or_TF=='TF'){

        ## Save the file as the top TFs:
        write.table(
          gene_CpGs_linked_df,
          file= paste(
            TENET_directory,
            'step7/',
            'user_peak_overlap/',
            'hypo_Gminus_user_peak_overlap_files/',
            'top_TFs/',
            gene_name,
            '_linked_probes_peak_overlap.tsv',
            sep=''
          ),
          quote= FALSE
        )
      }
    }

    ## Now sapply the function on the hypo Gminus top genes:
    sapply(
      top_hypo_Gminus_all_gene_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='gene'
    )

    ## Now sapply the function on the hypo Gminus top TFs:
    sapply(
      top_hypo_Gminus_all_TF_ENSG,
      internal_dataset_output_function,
      'gene_or_TF'='TF'
    )
  }
}

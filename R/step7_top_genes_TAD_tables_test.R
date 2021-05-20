#' step7_top_genes_TAD_tables_test
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates tables of information for each of the enhancer probes linked
#' to them for each of the four hypo or hypermethylated Gplus or Gminus analysis quadrants,
#' as selected by the user. These tables note which of the top genes/TFs each probe is linked to,
#' as well as the the total number of  genes and their names which happen to lie
#' within the same topologically-associating domain, or TAD, of each probe in each of the
#' user-supplied TAD files.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_TAD_tables' containing the results of this function.
#' @param TAD_directory Set a path to a directory which contains tab-delimited bed-like files, see: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 , with information on the TAD compartments of interest. Multiple files can be included in this subdirectory.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypermeth probes with G+ links.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypermeth probes with G- links.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypometh probes with G+ links.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypometh probes with G- links.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes to generate TAD tables for the enhancer probes linked to those genes.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Exports .tsv files with rows for each of the probes linked to the top genes/TFs as specified by the user, with the probe locations in the hg38 human genome, which of the top genes/TFs they're linked to, and the number of genes as well as their ensembl IDs and names found in the same TAD as each probe for each of the user input TAD files.
#' @export

step7_top_genes_TAD_tables_test <- function(
  TENET_directory,
  TAD_directory,
  DNA_methylation_manifest="HM450",
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
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

  ## Check to see if the user has supplied TAD files:

  ## First check if the user has specified a "TAD" or "tad" subdirectory
  ## in the main level of the TENET folder which will contain the TAD files:
  if(
    !dir.exists(
      TAD_directory
    )
  ){

    ## No TAD subdirectory found. Quit the function and note that to the user:
    stop("The specified TAD_directory was not found. Please create such a folder containing bed-like files with topologically associating domains information and provide a path to that directory as the argument")
  }

  ## Next check to see that there are TAD files included in the folder:
  TAD_file_paths <- list.files(
    path= TAD_directory,
    full.names=TRUE
  )

  ## If no TAD files are found, return an error:
  if(
    length(TAD_file_paths)==0
  ){

    ## No TAD files found. Quit the function and note that to the user:
    stop("No TAD files were found in the TAD subdirectory of the specified TENET directory. Please place bed-delimited files containing TAD information in the 'TAD' subdirectory of the specified TENET directory")
  }

  ## Create an empty list to load each of the files into:
  TAD_file_list <- vector(
    mode = "list",
    length = length(TAD_file_paths)
  )

  ## Prepare and load the TAD files:
  ## This assumes the files are at least bed3 formatted with the first 3 columns:
  ## Being the chromosome, start, and end coordinates:
  ## And are 0-indexed:
  for(i in c(1:length(TAD_file_paths))){

    ## Load the file:
    placeholder_file <- read.delim(
      TAD_file_paths[i],
      header= FALSE,
      stringsAsFactors = FALSE
    )

    ## Remove extra columns and rename the first three columns:
    ## Which should contain chr, start, end locations:
    if(ncol(placeholder_file)>=3){

      placeholder_file <- placeholder_file[1:3]

      colnames(placeholder_file) <- c(
        'seqnames','start','end'
      )

    } else{

      print("File is not at least bed3 file formatted. Please check that the TAD file contains at least 3 columns with chromosome, start, and end coordinates respectively, denoting TAD regions.")
      break()
    }

    ## Add 1 to the starts to make them 1-indexed:
    placeholder_file$start <- placeholder_file$start+1

    ## Add the final processed TAD file to the list:
    TAD_file_list[[i]] <- placeholder_file
  }

  ## Assign names of the files to the TAD file list:
  TAD_file_names <- sub(
    '\\..*',
    '',
    basename(TAD_file_paths)
  )

  names(TAD_file_list) <- TAD_file_names

  ## Write function to identify probes that overlap each TAD in a given file:
  TAD_probe_overlapper <- function(chromosome, start, stop, reference_peaks, buffer){
    paste(
      rownames(
        reference_peaks[
          (
            (reference_peaks[['seqnames']] == chromosome) &
              ((reference_peaks[['start']]-buffer) <= stop) &
              ((reference_peaks[['end']]+buffer) >= start)
          ),
        ]
      ), collapse=","
    )
  }

  ## Write functions to identify genes that overlap each TAD in a given file:
  TAD_gene_overlapper <- function(chromosome, TSS, reference_peaks, buffer){
    paste(
      rownames(
        reference_peaks[
          (
            (reference_peaks[['seqnames']] == chromosome) &
              (((reference_peaks[['start']]-buffer) <= TSS) &
                 ((reference_peaks[['end']]+buffer) >= TSS))
          ),
        ]
      ), collapse=","
    )
  }

  ## Now let's find genes in the same TAD(s):
  TAD_gene_finder <- function(TAD_numbers_comma_delimited, target_list){

    TAD_numbers <- unlist(strsplit(TAD_numbers_comma_delimited,','))

    gene_numbers <- sapply(
      TAD_numbers,
      grep,
      target_list
    )

    unique_gene_numbers <- unique(unlist(gene_numbers))

    unique_gene_numbers_combined <- paste(unique_gene_numbers,collapse=',')

    return(unique_gene_numbers_combined)
  }

  ## Write a function to count the genes in a TAD with a given probe:
  gene_count_in_TAD <- function(gene_numbers){

    gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

    return(length(gene_numbers_separated))
  }

  ## Now let's write a function to convert the gene numbers to the gene name:
  gene_name_from_number_lister <- function(gene_numbers, downregulated_gene_df, return_type){

    gene_numbers_separated <- as.numeric(unlist(strsplit(gene_numbers,',')))

    if(return_type=='gene_name'){

      gene_names <- downregulated_gene_df[
        gene_numbers_separated,
        'gene_name'
      ]

    } else if(return_type=='gene_id'){

      gene_names <- downregulated_gene_df[
        gene_numbers_separated,
        'ensembl_ID'
      ]
    }

    unique_gene_names_combined <- paste(unique(gene_names),collapse=',')

    return(unique_gene_names_combined)
  }

  ## Now let's annotate each gene with the TADs it is found in
  ## for each file, which will be overlapped with the information from the probes
  ## from each analysis later:
  for(TAD_file_index in c(1:length(TAD_file_list))){

    ## Create name for the TAD file overlap columns:
    TAD_file_overlap_column_name <- paste(
      names(TAD_file_list[TAD_file_index]),
      '_overlap_number',
      sep=''
    )

    ## Get overlaps for each of the TADs with the genes of interest:
    gencode_v22_genes[TAD_file_overlap_column_name] <- unname(
      parallel::mcmapply(
        TAD_gene_overlapper,
        chromosome= gencode_v22_genes[,1],
        TSS= gencode_v22_genes[,'TSS'],
        MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0),
        mc.cores = core_count
      )
    )
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
          'hyper_Gplus_tad_tables',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_tad_tables',
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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes and tfs selected:
    hyper_Gplus_probes_linked_to_significant_genes <- unique(
      hyper_Gplus_sig_link_zscores[
        hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_gene_ENSG,
        'probeID'
      ]
    )

    hyper_Gplus_probes_linked_to_significant_TFs <- unique(
      hyper_Gplus_sig_link_zscores[
        hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_TF_ENSG,
        'probeID'
      ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    hyper_Gplus_probe_dataset_linked_to_significant_genes <- data.frame(
      'probe_ID'=sort(hyper_Gplus_probes_linked_to_significant_genes),
      stringsAsFactors = FALSE
    )

    hyper_Gplus_probe_dataset_linked_to_significant_TFs <- data.frame(
      'probe_ID'=sort(hyper_Gplus_probes_linked_to_significant_TFs),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    hyper_Gplus_probe_dataset_linked_to_significant_genes$seqnames <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'chr'
    ]

    hyper_Gplus_probe_dataset_linked_to_significant_TFs$seqnames <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'chr'
    ]

    hyper_Gplus_probe_dataset_linked_to_significant_genes$start <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'start'
    ]

    hyper_Gplus_probe_dataset_linked_to_significant_TFs$start <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'start'
    ]

    hyper_Gplus_probe_dataset_linked_to_significant_genes$end <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'end'
    ]

    hyper_Gplus_probe_dataset_linked_to_significant_TFs$end <-  hg38_manifest_no_NA_granges_df[
      hyper_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'end'
    ]

    ## For each of the genes/TFs of interest, create a new column in the probe dataframes
    ## indicating if each probe is linked to that gene or not:
    for(gene_ENSG in top_hyper_Gplus_all_gene_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_gene <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID == gene_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[column_name] <- hyper_Gplus_probe_dataset_linked_to_significant_genes$probe_ID %in% probes_linked_to_gene

    }

    for(TF_ENSG in top_hyper_Gplus_all_TF_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_TF <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID == TF_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hyper_Gplus_probe_dataset_linked_to_significant_TFs[column_name] <- hyper_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID %in% probes_linked_to_TF

    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

      ## Create name for the TAD file overlap columns:
      TAD_file_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_overlap_number',
        sep=''
      )

      ## Create name for the TAD file gene numbers found in the TAD:
      gene_row_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_numbers',
        sep=''
      )

      ## Create a name for the total number of genes in a TAD of the probe:
      gene_count_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_count_in_TAD',
        sep=''
      )

      ## Create names for the gene ENSGs and gene names found in TADs:
      gene_ENSG_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_ENSGs',
        sep=''
      )

      gene_name_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_names',
        sep=''
      )

      ## Get overlaps for each of the TADs with the probes of interest:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hyper_Gplus_probe_dataset_linked_to_significant_genes[,2],
          start= hyper_Gplus_probe_dataset_linked_to_significant_genes[,3],
          stop= hyper_Gplus_probe_dataset_linked_to_significant_genes[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hyper_Gplus_probe_dataset_linked_to_significant_TFs[,2],
          start= hyper_Gplus_probe_dataset_linked_to_significant_TFs[,3],
          stop= hyper_Gplus_probe_dataset_linked_to_significant_TFs[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      ## Get the numbers of the genes that have overlap TADs with
      ## The probes:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_genes[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_TFs[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      ## Count the number of total genes within a TAD of the probe:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[gene_count_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[gene_count_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      ## Convert the numbers into gene ENSGs:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[gene_ENSG_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[gene_ENSG_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      ## Convert the numbers into gene names:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[gene_name_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[gene_name_column_name] <- unname(
        sapply(
          hyper_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      ## If no TAD is found, make a note of that in the gene names/ENSG columns:
      for(i in c(1:nrow(hyper_Gplus_probe_dataset_linked_to_significant_genes))){

        if(hyper_Gplus_probe_dataset_linked_to_significant_genes[i,TAD_file_overlap_column_name]==''){

          hyper_Gplus_probe_dataset_linked_to_significant_genes[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hyper_Gplus_probe_dataset_linked_to_significant_genes[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      for(i in c(1:nrow(hyper_Gplus_probe_dataset_linked_to_significant_TFs))){

        if(hyper_Gplus_probe_dataset_linked_to_significant_TFs[i,TAD_file_overlap_column_name]==''){

          hyper_Gplus_probe_dataset_linked_to_significant_TFs[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hyper_Gplus_probe_dataset_linked_to_significant_TFs[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      ## Let's remove unneeded columns now:
      hyper_Gplus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- NULL
      hyper_Gplus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- NULL

      hyper_Gplus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- NULL
      hyper_Gplus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- NULL

      ## Write the tables to files:
      write.table(
        hyper_Gplus_probe_dataset_linked_to_significant_genes,
        file= paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_tad_tables/',
          'hyper_Gplus_top_genes_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

      write.table(
        hyper_Gplus_probe_dataset_linked_to_significant_TFs,
        file= paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_tad_tables/',
          'hyper_Gplus_top_TFs_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

    }
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
          'hyper_Gminus_tad_tables',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_tad_tables',
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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes and tfs selected:
    hyper_Gminus_probes_linked_to_significant_genes <- unique(
      hyper_Gminus_sig_link_zscores[
        hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_gene_ENSG,
        'probeID'
      ]
    )

    hyper_Gminus_probes_linked_to_significant_TFs <- unique(
      hyper_Gminus_sig_link_zscores[
        hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_TF_ENSG,
        'probeID'
      ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    hyper_Gminus_probe_dataset_linked_to_significant_genes <- data.frame(
      'probe_ID'=sort(hyper_Gminus_probes_linked_to_significant_genes),
      stringsAsFactors = FALSE
    )

    hyper_Gminus_probe_dataset_linked_to_significant_TFs <- data.frame(
      'probe_ID'=sort(hyper_Gminus_probes_linked_to_significant_TFs),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    hyper_Gminus_probe_dataset_linked_to_significant_genes$seqnames <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'chr'
    ]

    hyper_Gminus_probe_dataset_linked_to_significant_TFs$seqnames <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'chr'
    ]

    hyper_Gminus_probe_dataset_linked_to_significant_genes$start <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'start'
    ]

    hyper_Gminus_probe_dataset_linked_to_significant_TFs$start <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'start'
    ]

    hyper_Gminus_probe_dataset_linked_to_significant_genes$end <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'end'
    ]

    hyper_Gminus_probe_dataset_linked_to_significant_TFs$end <-  hg38_manifest_no_NA_granges_df[
      hyper_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'end'
    ]

    ## For each of the genes/TFs of interest, create a new column in the probe dataframes
    ## indicating if each probe is linked to that gene or not:
    for(gene_ENSG in top_hyper_Gminus_all_gene_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_gene <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID == gene_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[column_name] <- hyper_Gminus_probe_dataset_linked_to_significant_genes$probe_ID %in% probes_linked_to_gene

    }

    for(TF_ENSG in top_hyper_Gminus_all_TF_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_TF <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID == TF_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hyper_Gminus_probe_dataset_linked_to_significant_TFs[column_name] <- hyper_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID %in% probes_linked_to_TF

    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

      ## Create name for the TAD file overlap columns:
      TAD_file_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_overlap_number',
        sep=''
      )

      ## Create name for the TAD file gene numbers found in the TAD:
      gene_row_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_numbers',
        sep=''
      )

      ## Create a name for the total number of genes in a TAD of the probe:
      gene_count_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_count_in_TAD',
        sep=''
      )

      ## Create names for the gene ENSGs and gene names found in TADs:
      gene_ENSG_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_ENSGs',
        sep=''
      )

      gene_name_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_names',
        sep=''
      )

      ## Get overlaps for each of the TADs with the probes of interest:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hyper_Gminus_probe_dataset_linked_to_significant_genes[,2],
          start= hyper_Gminus_probe_dataset_linked_to_significant_genes[,3],
          stop= hyper_Gminus_probe_dataset_linked_to_significant_genes[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hyper_Gminus_probe_dataset_linked_to_significant_TFs[,2],
          start= hyper_Gminus_probe_dataset_linked_to_significant_TFs[,3],
          stop= hyper_Gminus_probe_dataset_linked_to_significant_TFs[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      ## Get the numbers of the genes that have overlap TADs with
      ## The probes:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_genes[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_TFs[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      ## Count the number of total genes within a TAD of the probe:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[gene_count_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[gene_count_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      ## Convert the numbers into gene ENSGs:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[gene_ENSG_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[gene_ENSG_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      ## Convert the numbers into gene names:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[gene_name_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[gene_name_column_name] <- unname(
        sapply(
          hyper_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      ## If no TAD is found, make a note of that in the gene names/ENSG columns:
      for(i in c(1:nrow(hyper_Gminus_probe_dataset_linked_to_significant_genes))){

        if(hyper_Gminus_probe_dataset_linked_to_significant_genes[i,TAD_file_overlap_column_name]==''){

          hyper_Gminus_probe_dataset_linked_to_significant_genes[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hyper_Gminus_probe_dataset_linked_to_significant_genes[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      for(i in c(1:nrow(hyper_Gminus_probe_dataset_linked_to_significant_TFs))){

        if(hyper_Gminus_probe_dataset_linked_to_significant_TFs[i,TAD_file_overlap_column_name]==''){

          hyper_Gminus_probe_dataset_linked_to_significant_TFs[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hyper_Gminus_probe_dataset_linked_to_significant_TFs[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      ## Let's remove unneeded columns now:
      hyper_Gminus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- NULL
      hyper_Gminus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- NULL

      hyper_Gminus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- NULL
      hyper_Gminus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- NULL

      ## Write the tables to files:
      write.table(
        hyper_Gminus_probe_dataset_linked_to_significant_genes,
        file= paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_tad_tables/',
          'hyper_Gminus_top_genes_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

      write.table(
        hyper_Gminus_probe_dataset_linked_to_significant_TFs,
        file= paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_tad_tables/',
          'hyper_Gminus_top_TFs_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

    }
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
          'hypo_Gplus_tad_tables',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_tad_tables',
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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes and tfs selected:
    hypo_Gplus_probes_linked_to_significant_genes <- unique(
      hypo_Gplus_sig_link_zscores[
        hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_gene_ENSG,
        'probeID'
      ]
    )

    hypo_Gplus_probes_linked_to_significant_TFs <- unique(
      hypo_Gplus_sig_link_zscores[
        hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_TF_ENSG,
        'probeID'
      ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    hypo_Gplus_probe_dataset_linked_to_significant_genes <- data.frame(
      'probe_ID'=sort(hypo_Gplus_probes_linked_to_significant_genes),
      stringsAsFactors = FALSE
    )

    hypo_Gplus_probe_dataset_linked_to_significant_TFs <- data.frame(
      'probe_ID'=sort(hypo_Gplus_probes_linked_to_significant_TFs),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    hypo_Gplus_probe_dataset_linked_to_significant_genes$seqnames <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'chr'
    ]

    hypo_Gplus_probe_dataset_linked_to_significant_TFs$seqnames <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'chr'
    ]

    hypo_Gplus_probe_dataset_linked_to_significant_genes$start <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'start'
    ]

    hypo_Gplus_probe_dataset_linked_to_significant_TFs$start <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'start'
    ]

    hypo_Gplus_probe_dataset_linked_to_significant_genes$end <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_genes$probe_ID,
      'end'
    ]

    hypo_Gplus_probe_dataset_linked_to_significant_TFs$end <-  hg38_manifest_no_NA_granges_df[
      hypo_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'end'
    ]

    ## For each of the genes/TFs of interest, create a new column in the probe dataframes
    ## indicating if each probe is linked to that gene or not:
    for(gene_ENSG in top_hypo_Gplus_all_gene_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_gene <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID == gene_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[column_name] <- hypo_Gplus_probe_dataset_linked_to_significant_genes$probe_ID %in% probes_linked_to_gene

    }

    for(TF_ENSG in top_hypo_Gplus_all_TF_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_TF <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID == TF_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hypo_Gplus_probe_dataset_linked_to_significant_TFs[column_name] <- hypo_Gplus_probe_dataset_linked_to_significant_TFs$probe_ID %in% probes_linked_to_TF

    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

      ## Create name for the TAD file overlap columns:
      TAD_file_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_overlap_number',
        sep=''
      )

      ## Create name for the TAD file gene numbers found in the TAD:
      gene_row_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_numbers',
        sep=''
      )

      ## Create a name for the total number of genes in a TAD of the probe:
      gene_count_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_count_in_TAD',
        sep=''
      )

      ## Create names for the gene ENSGs and gene names found in TADs:
      gene_ENSG_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_ENSGs',
        sep=''
      )

      gene_name_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_names',
        sep=''
      )

      ## Get overlaps for each of the TADs with the probes of interest:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hypo_Gplus_probe_dataset_linked_to_significant_genes[,2],
          start= hypo_Gplus_probe_dataset_linked_to_significant_genes[,3],
          stop= hypo_Gplus_probe_dataset_linked_to_significant_genes[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hypo_Gplus_probe_dataset_linked_to_significant_TFs[,2],
          start= hypo_Gplus_probe_dataset_linked_to_significant_TFs[,3],
          stop= hypo_Gplus_probe_dataset_linked_to_significant_TFs[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      ## Get the numbers of the genes that have overlap TADs with
      ## The probes:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_genes[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_TFs[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      ## Count the number of total genes within a TAD of the probe:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[gene_count_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[gene_count_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      ## Convert the numbers into gene ENSGs:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[gene_ENSG_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[gene_ENSG_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      ## Convert the numbers into gene names:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[gene_name_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[gene_name_column_name] <- unname(
        sapply(
          hypo_Gplus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      ## If no TAD is found, make a note of that in the gene names/ENSG columns:
      for(i in c(1:nrow(hypo_Gplus_probe_dataset_linked_to_significant_genes))){

        if(hypo_Gplus_probe_dataset_linked_to_significant_genes[i,TAD_file_overlap_column_name]==''){

          hypo_Gplus_probe_dataset_linked_to_significant_genes[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hypo_Gplus_probe_dataset_linked_to_significant_genes[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      for(i in c(1:nrow(hypo_Gplus_probe_dataset_linked_to_significant_TFs))){

        if(hypo_Gplus_probe_dataset_linked_to_significant_TFs[i,TAD_file_overlap_column_name]==''){

          hypo_Gplus_probe_dataset_linked_to_significant_TFs[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hypo_Gplus_probe_dataset_linked_to_significant_TFs[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      ## Let's remove unneeded columns now:
      hypo_Gplus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- NULL
      hypo_Gplus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- NULL

      hypo_Gplus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- NULL
      hypo_Gplus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- NULL

      ## Write the tables to files:
      write.table(
        hypo_Gplus_probe_dataset_linked_to_significant_genes,
        file= paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_tad_tables/',
          'hypo_Gplus_top_genes_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

      write.table(
        hypo_Gplus_probe_dataset_linked_to_significant_TFs,
        file= paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_tad_tables/',
          'hypo_Gplus_top_TFs_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

    }
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
          'hypo_Gminus_tad_tables',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_tad_tables',
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

    ## Let's index a data frame with some basic info about the probes of interest
    ## linked to the top genes and tfs selected:
    hypo_Gminus_probes_linked_to_significant_genes <- unique(
      hypo_Gminus_sig_link_zscores[
        hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_gene_ENSG,
        'probeID'
      ]
    )

    hypo_Gminus_probes_linked_to_significant_TFs <- unique(
      hypo_Gminus_sig_link_zscores[
        hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_TF_ENSG,
        'probeID'
      ]
    )

    ## Create a data frame with the probes as the first (and only column for now):
    hypo_Gminus_probe_dataset_linked_to_significant_genes <- data.frame(
      'probe_ID'=sort(hypo_Gminus_probes_linked_to_significant_genes),
      stringsAsFactors = FALSE
    )

    hypo_Gminus_probe_dataset_linked_to_significant_TFs <- data.frame(
      'probe_ID'=sort(hypo_Gminus_probes_linked_to_significant_TFs),
      stringsAsFactors = FALSE
    )

    ## Now let's add the probe location info for hg38:
    hypo_Gminus_probe_dataset_linked_to_significant_genes$seqnames <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'chr'
    ]

    hypo_Gminus_probe_dataset_linked_to_significant_TFs$seqnames <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'chr'
    ]

    hypo_Gminus_probe_dataset_linked_to_significant_genes$start <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'start'
    ]

    hypo_Gminus_probe_dataset_linked_to_significant_TFs$start <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'start'
    ]

    hypo_Gminus_probe_dataset_linked_to_significant_genes$end <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_genes$probe_ID,
      'end'
    ]

    hypo_Gminus_probe_dataset_linked_to_significant_TFs$end <-  hg38_manifest_no_NA_granges_df[
      hypo_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID,
      'end'
    ]

    ## For each of the genes/TFs of interest, create a new column in the probe dataframes
    ## indicating if each probe is linked to that gene or not:
    for(gene_ENSG in top_hypo_Gminus_all_gene_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_gene <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID == gene_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[column_name] <- hypo_Gminus_probe_dataset_linked_to_significant_genes$probe_ID %in% probes_linked_to_gene

    }

    for(TF_ENSG in top_hypo_Gminus_all_TF_ENSG){

      ## Get the gene name corresponding to the ENSG:
      corresponding_gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      column_name <- paste(
        corresponding_gene_name,
        '_linked',
        sep=''
      )

      ## Get the probes linked to that gene specifically:
      probes_linked_to_TF <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID == TF_ENSG,
          'probeID'
        ]
      )

      ## Add a TRUE/FALSE vector to tell if each of the unique probes linked to
      ## at least one of the top n TFs is linked to this specific one:
      hypo_Gminus_probe_dataset_linked_to_significant_TFs[column_name] <- hypo_Gminus_probe_dataset_linked_to_significant_TFs$probe_ID %in% probes_linked_to_TF

    }

    ## Now for each of the probes, let's see what genes are in the same TAD as it:
    ## To do this, we will need to identify the TAD number(s) each probe is in and each gene is in:
    for(TAD_file_index in c(1:length(TAD_file_list))){

      ## Create name for the TAD file overlap columns:
      TAD_file_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_overlap_number',
        sep=''
      )

      ## Create name for the TAD file gene numbers found in the TAD:
      gene_row_overlap_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_numbers',
        sep=''
      )

      ## Create a name for the total number of genes in a TAD of the probe:
      gene_count_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_gene_count_in_TAD',
        sep=''
      )

      ## Create names for the gene ENSGs and gene names found in TADs:
      gene_ENSG_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_ENSGs',
        sep=''
      )

      gene_name_column_name <- paste(
        names(TAD_file_list[TAD_file_index]),
        '_TAD_gene_names',
        sep=''
      )

      ## Get overlaps for each of the TADs with the probes of interest:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hypo_Gminus_probe_dataset_linked_to_significant_genes[,2],
          start= hypo_Gminus_probe_dataset_linked_to_significant_genes[,3],
          stop= hypo_Gminus_probe_dataset_linked_to_significant_genes[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- unname(
        mapply(
          TAD_probe_overlapper,
          chromosome= hypo_Gminus_probe_dataset_linked_to_significant_TFs[,2],
          start= hypo_Gminus_probe_dataset_linked_to_significant_TFs[,3],
          stop= hypo_Gminus_probe_dataset_linked_to_significant_TFs[,4],
          MoreArgs = list(reference_peaks=TAD_file_list[[TAD_file_index]], buffer=0)
        )
      )

      ## Get the numbers of the genes that have overlap TADs with
      ## The probes:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_genes[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_TFs[,TAD_file_overlap_column_name],
          TAD_gene_finder,
          target_list=gencode_v22_genes[,TAD_file_overlap_column_name]
        )
      )

      ## Count the number of total genes within a TAD of the probe:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[gene_count_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[gene_count_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_count_in_TAD
        )
      )

      ## Convert the numbers into gene ENSGs:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[gene_ENSG_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[gene_ENSG_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_id'
        )
      )

      ## Convert the numbers into gene names:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[gene_name_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_genes[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[gene_name_column_name] <- unname(
        sapply(
          hypo_Gminus_probe_dataset_linked_to_significant_TFs[,gene_row_overlap_column_name],
          gene_name_from_number_lister,
          downregulated_gene_df=gencode_v22_genes,
          return_type='gene_name'
        )
      )

      ## If no TAD is found, make a note of that in the gene names/ENSG columns:
      for(i in c(1:nrow(hypo_Gminus_probe_dataset_linked_to_significant_genes))){

        if(hypo_Gminus_probe_dataset_linked_to_significant_genes[i,TAD_file_overlap_column_name]==''){

          hypo_Gminus_probe_dataset_linked_to_significant_genes[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hypo_Gminus_probe_dataset_linked_to_significant_genes[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      for(i in c(1:nrow(hypo_Gminus_probe_dataset_linked_to_significant_TFs))){

        if(hypo_Gminus_probe_dataset_linked_to_significant_TFs[i,TAD_file_overlap_column_name]==''){

          hypo_Gminus_probe_dataset_linked_to_significant_TFs[i,gene_ENSG_column_name] <- "No_TAD_indentified"
          hypo_Gminus_probe_dataset_linked_to_significant_TFs[i,gene_name_column_name] <- "No_TAD_indentified"
        }

      }

      ## Let's remove unneeded columns now:
      hypo_Gminus_probe_dataset_linked_to_significant_genes[TAD_file_overlap_column_name] <- NULL
      hypo_Gminus_probe_dataset_linked_to_significant_genes[gene_row_overlap_column_name] <- NULL

      hypo_Gminus_probe_dataset_linked_to_significant_TFs[TAD_file_overlap_column_name] <- NULL
      hypo_Gminus_probe_dataset_linked_to_significant_TFs[gene_row_overlap_column_name] <- NULL

      ## Write the tables to files:
      write.table(
        hypo_Gminus_probe_dataset_linked_to_significant_genes,
        file= paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_tad_tables/',
          'hypo_Gminus_top_genes_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

      write.table(
        hypo_Gminus_probe_dataset_linked_to_significant_TFs,
        file= paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_tad_tables/',
          'hypo_Gminus_top_TFs_TAD_analysis.tsv',
          sep=''
        ),
        quote= FALSE,
        sep='\t',
        row.names= FALSE
      )

    }
  }

}

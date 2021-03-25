#' step7_top_genes_UCSC_bed_files
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates interact files see: https://genome.ucsc.edu/goldenPath/help/interact.html
#' that can be uploaded to the UCSC genome browser and visualize the
#' links between each of the top specified genes/TFs and the enhancer probes linked
#' to them for each of the four hypo or hypermethylated Gplus or Gminus analysis quadrants,
#' as selected by the user.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_UCSC_bed_files' containing the results of this function.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create TAD tables for the enhancer probes linked to the top genes/TFs by most hypermeth probes with G+ links.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create interact .bed files that can be uploaded to the UCSC genome browser showing links between the top genes/TFs by most hypermeth probes with G+ links, and said linked probes.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create interact .bed files that can be uploaded to the UCSC genome browser showing links between the top genes/TFs by most hypermeth probes with G- links, and said linked probes.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create interact .bed files that can be uploaded to the UCSC genome browser showing links between the top genes/TFs by most hypometh probes with G+ links, and said linked probes.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create interact .bed files that can be uploaded to the UCSC genome browser showing links between the top genes/TFs by most hypometh probes with G- links, and said linked probes.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes to generate interact bed files for showing the links between the genes and each of their linked enhancer probes.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Exports .bed formatted interact files to upload to the UCSC genome browser. These files display the interactions between the top genes/TFs and their linked enhancer probes for the given analysis types.
#' @export

step7_top_genes_UCSC_bed_files <- function(
  TENET_directory,
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

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus TAD tables:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_ucsc_bed_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_ucsc_bed_files',
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

    ## Now create dataframes listing the top genes/TFs in one column (repeated)
    ## And all the probes linked to them in a second column:

    ## Index empty dataframes:
    top_hyper_Gplus_all_genes_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    top_hyper_Gplus_all_TFs_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    ## Add the probes linked to each of the top genes to the dataframe:
    for(gene_ENSG in top_hyper_Gplus_all_gene_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          gene_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hyper_Gplus_all_genes_intersect_df <- rbind(
        top_hyper_Gplus_all_genes_intersect_df,
        temp_df
      )
    }

    ## Add the probes linked to each of the top TFs to the dataframe:
    for(TF_ENSG in top_hyper_Gplus_all_TF_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID==TF_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          TF_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hyper_Gplus_all_TFs_intersect_df <- rbind(
        top_hyper_Gplus_all_TFs_intersect_df,
        temp_df
      )
    }

    ## For each of the genes in the intersect dfs
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    top_hyper_Gplus_all_genes_intersect_df$gene_chr <- gencode_v22_genes[
      top_hyper_Gplus_all_genes_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hyper_Gplus_all_genes_intersect_df$gene_start <- gencode_v22_genes[
      top_hyper_Gplus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gplus_all_genes_intersect_df$gene_end <- gencode_v22_genes[
      top_hyper_Gplus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$gene_chr <- gencode_v22_genes[
      top_hyper_Gplus_all_TFs_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$gene_start <- gencode_v22_genes[
      top_hyper_Gplus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$gene_end <- gencode_v22_genes[
      top_hyper_Gplus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    ## For each of the probes in the intersect dfs
    ## get the chromosome, "start", and "end" of the probe:
    top_hyper_Gplus_all_genes_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_genes_intersect_df$CpG_ID,
      'chr'
    ]

    top_hyper_Gplus_all_genes_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_genes_intersect_df$CpG_ID,
      'start'
    ]

    top_hyper_Gplus_all_genes_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_genes_intersect_df$CpG_ID,
      'end'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_TFs_intersect_df$CpG_ID,
      'chr'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_TFs_intersect_df$CpG_ID,
      'start'
    ]

    top_hyper_Gplus_all_TFs_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gplus_all_TFs_intersect_df$CpG_ID,
      'end'
    ]

    ## Do some numeric conversions:
    top_hyper_Gplus_all_genes_intersect_df$gene_start <- as.integer(top_hyper_Gplus_all_genes_intersect_df$gene_start)
    top_hyper_Gplus_all_genes_intersect_df$gene_end <- as.integer(top_hyper_Gplus_all_genes_intersect_df$gene_end)

    top_hyper_Gplus_all_TFs_intersect_df$gene_start <- as.integer(top_hyper_Gplus_all_TFs_intersect_df$gene_start)
    top_hyper_Gplus_all_TFs_intersect_df$gene_end <- as.integer(top_hyper_Gplus_all_TFs_intersect_df$gene_end)

    ## Now let's use the rainbow color function to set up a gradient of colors equal to
    ## the number of genes/TFs analyzed and assign each gene/TF a color:
    jet_color_numeric_color_grad_genes <- rainbow(top_gene_number)
    jet_color_numeric_color_grad_TFs <- rainbow(top_gene_number)

    names(jet_color_numeric_color_grad_genes) <- top_hyper_Gplus_all_gene_ENSG
    names(jet_color_numeric_color_grad_TFs) <- top_hyper_Gplus_all_TF_ENSG

    ## Now that we have the information, let's assemble the output file
    ## for the top genes:
    top_hyper_Gplus_all_genes_output_df <- data.frame(
      'chrom'= top_hyper_Gplus_all_genes_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hyper_Gplus_all_genes_intersect_df$gene_chr==top_hyper_Gplus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gplus_all_genes_intersect_df$CpG_start < top_hyper_Gplus_all_genes_intersect_df$gene_start),
          top_hyper_Gplus_all_genes_intersect_df$CpG_start,
          top_hyper_Gplus_all_genes_intersect_df$gene_start
        ),
        top_hyper_Gplus_all_genes_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hyper_Gplus_all_genes_intersect_df$gene_chr==top_hyper_Gplus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gplus_all_genes_intersect_df$CpG_start < top_hyper_Gplus_all_genes_intersect_df$gene_start),
          top_hyper_Gplus_all_genes_intersect_df$gene_start-1,
          top_hyper_Gplus_all_genes_intersect_df$CpG_start
        ),
        top_hyper_Gplus_all_genes_intersect_df$gene_end
      ),
      'name'= paste(
        top_hyper_Gplus_all_genes_intersect_df$gene_ID,
        top_hyper_Gplus_all_genes_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hyper_Gplus_all_genes_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hyper_Gplus_all_genes_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hyper_Gplus_all_genes_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hyper_Gplus_all_genes_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hyper_Gplus_all_genes_intersect_df$gene_chr,
      'sourceStart'= top_hyper_Gplus_all_genes_intersect_df$gene_start-1,
      'sourceEnd'= top_hyper_Gplus_all_genes_intersect_df$gene_end,
      'sourceName'= top_hyper_Gplus_all_genes_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hyper_Gplus_all_genes_intersect_df)
      ),
      'targetChrom'= top_hyper_Gplus_all_genes_intersect_df$CpG_chr,
      'targetStart'= top_hyper_Gplus_all_genes_intersect_df$CpG_start,
      'targetEnd'= top_hyper_Gplus_all_genes_intersect_df$CpG_end,
      'targetName'= top_hyper_Gplus_all_genes_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hyper_Gplus_all_genes_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Now that we have the information, let's assemble the output file
    ## for the top TFs:
    top_hyper_Gplus_all_TFs_output_df <- data.frame(
      'chrom'= top_hyper_Gplus_all_TFs_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hyper_Gplus_all_TFs_intersect_df$gene_chr==top_hyper_Gplus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gplus_all_TFs_intersect_df$CpG_start < top_hyper_Gplus_all_TFs_intersect_df$gene_start),
          top_hyper_Gplus_all_TFs_intersect_df$CpG_start,
          top_hyper_Gplus_all_TFs_intersect_df$gene_start
        ),
        top_hyper_Gplus_all_TFs_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hyper_Gplus_all_TFs_intersect_df$gene_chr==top_hyper_Gplus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gplus_all_TFs_intersect_df$CpG_start < top_hyper_Gplus_all_TFs_intersect_df$gene_start),
          top_hyper_Gplus_all_TFs_intersect_df$gene_start-1,
          top_hyper_Gplus_all_TFs_intersect_df$CpG_start
        ),
        top_hyper_Gplus_all_TFs_intersect_df$gene_end
      ),
      'name'= paste(
        top_hyper_Gplus_all_TFs_intersect_df$gene_ID,
        top_hyper_Gplus_all_TFs_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hyper_Gplus_all_TFs_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hyper_Gplus_all_TFs_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hyper_Gplus_all_TFs_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hyper_Gplus_all_TFs_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hyper_Gplus_all_TFs_intersect_df$gene_chr,
      'sourceStart'= top_hyper_Gplus_all_TFs_intersect_df$gene_start-1,
      'sourceEnd'= top_hyper_Gplus_all_TFs_intersect_df$gene_end,
      'sourceName'= top_hyper_Gplus_all_TFs_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hyper_Gplus_all_TFs_intersect_df)
      ),
      'targetChrom'= top_hyper_Gplus_all_TFs_intersect_df$CpG_chr,
      'targetStart'= top_hyper_Gplus_all_TFs_intersect_df$CpG_start,
      'targetEnd'= top_hyper_Gplus_all_TFs_intersect_df$CpG_end,
      'targetName'= top_hyper_Gplus_all_TFs_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hyper_Gplus_all_TFs_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Create text for the header lines:
    top_hyper_Gplus_all_genes_header_text <- "track type=interact name=\"TENETR_hyper.G+_interactions\" description=\"TENETR top genes to enhancer DNA methylation probe links\""

    top_hyper_Gplus_all_TFs_header_text <- "track type=interact name=\"TENETR_hyper.G+_interactions\" description=\"TENETR top TFs to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    top_hyper_Gplus_all_genes_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gplus_ucsc_bed_files/',
      'top_hyper_Gplus_genes_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    top_hyper_Gplus_all_TFs_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gplus_ucsc_bed_files/',
      'top_hyper_Gplus_TFs_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    ## Add the header line to the new bed files:
    cat(
      top_hyper_Gplus_all_genes_header_text,
      "\n",
      file = top_hyper_Gplus_all_genes_bed_file_name
    )

    cat(
      top_hyper_Gplus_all_TFs_header_text,
      "\n",
      file = top_hyper_Gplus_all_TFs_bed_file_name
    )

    ## Write the info to the files:
    write.table(
      top_hyper_Gplus_all_genes_output_df,
      file = top_hyper_Gplus_all_genes_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
    )

    write.table(
      top_hyper_Gplus_all_TFs_output_df,
      file = top_hyper_Gplus_all_TFs_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
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
          'hyper_Gminus_ucsc_bed_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_ucsc_bed_files',
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

    ## Now create dataframes listing the top genes/TFs in one column (repeated)
    ## And all the probes linked to them in a second column:

    ## Index empty dataframes:
    top_hyper_Gminus_all_genes_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    top_hyper_Gminus_all_TFs_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    ## Add the probes linked to each of the top genes to the dataframe:
    for(gene_ENSG in top_hyper_Gminus_all_gene_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          gene_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hyper_Gminus_all_genes_intersect_df <- rbind(
        top_hyper_Gminus_all_genes_intersect_df,
        temp_df
      )
    }

    ## Add the probes linked to each of the top TFs to the dataframe:
    for(TF_ENSG in top_hyper_Gminus_all_TF_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID==TF_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          TF_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hyper_Gminus_all_TFs_intersect_df <- rbind(
        top_hyper_Gminus_all_TFs_intersect_df,
        temp_df
      )
    }

    ## For each of the genes in the intersect dfs
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    top_hyper_Gminus_all_genes_intersect_df$gene_chr <- gencode_v22_genes[
      top_hyper_Gminus_all_genes_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hyper_Gminus_all_genes_intersect_df$gene_start <- gencode_v22_genes[
      top_hyper_Gminus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gminus_all_genes_intersect_df$gene_end <- gencode_v22_genes[
      top_hyper_Gminus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$gene_chr <- gencode_v22_genes[
      top_hyper_Gminus_all_TFs_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$gene_start <- gencode_v22_genes[
      top_hyper_Gminus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$gene_end <- gencode_v22_genes[
      top_hyper_Gminus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    ## For each of the probes in the intersect dfs
    ## get the chromosome, "start", and "end" of the probe:
    top_hyper_Gminus_all_genes_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_genes_intersect_df$CpG_ID,
      'chr'
    ]

    top_hyper_Gminus_all_genes_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_genes_intersect_df$CpG_ID,
      'start'
    ]

    top_hyper_Gminus_all_genes_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_genes_intersect_df$CpG_ID,
      'end'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_TFs_intersect_df$CpG_ID,
      'chr'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_TFs_intersect_df$CpG_ID,
      'start'
    ]

    top_hyper_Gminus_all_TFs_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hyper_Gminus_all_TFs_intersect_df$CpG_ID,
      'end'
    ]

    ## Do some numeric conversions:
    top_hyper_Gminus_all_genes_intersect_df$gene_start <- as.integer(top_hyper_Gminus_all_genes_intersect_df$gene_start)
    top_hyper_Gminus_all_genes_intersect_df$gene_end <- as.integer(top_hyper_Gminus_all_genes_intersect_df$gene_end)

    top_hyper_Gminus_all_TFs_intersect_df$gene_start <- as.integer(top_hyper_Gminus_all_TFs_intersect_df$gene_start)
    top_hyper_Gminus_all_TFs_intersect_df$gene_end <- as.integer(top_hyper_Gminus_all_TFs_intersect_df$gene_end)

    ## Now let's use the rainbow color function to set up a gradient of colors equal to
    ## the number of genes/TFs analyzed and assign each gene/TF a color:
    jet_color_numeric_color_grad_genes <- rainbow(top_gene_number)
    jet_color_numeric_color_grad_TFs <- rainbow(top_gene_number)

    names(jet_color_numeric_color_grad_genes) <- top_hyper_Gminus_all_gene_ENSG
    names(jet_color_numeric_color_grad_TFs) <- top_hyper_Gminus_all_TF_ENSG

    ## Now that we have the information, let's assemble the output file
    ## for the top genes:
    top_hyper_Gminus_all_genes_output_df <- data.frame(
      'chrom'= top_hyper_Gminus_all_genes_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hyper_Gminus_all_genes_intersect_df$gene_chr==top_hyper_Gminus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gminus_all_genes_intersect_df$CpG_start < top_hyper_Gminus_all_genes_intersect_df$gene_start),
          top_hyper_Gminus_all_genes_intersect_df$CpG_start,
          top_hyper_Gminus_all_genes_intersect_df$gene_start
        ),
        top_hyper_Gminus_all_genes_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hyper_Gminus_all_genes_intersect_df$gene_chr==top_hyper_Gminus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gminus_all_genes_intersect_df$CpG_start < top_hyper_Gminus_all_genes_intersect_df$gene_start),
          top_hyper_Gminus_all_genes_intersect_df$gene_start-1,
          top_hyper_Gminus_all_genes_intersect_df$CpG_start
        ),
        top_hyper_Gminus_all_genes_intersect_df$gene_end
      ),
      'name'= paste(
        top_hyper_Gminus_all_genes_intersect_df$gene_ID,
        top_hyper_Gminus_all_genes_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hyper_Gminus_all_genes_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hyper_Gminus_all_genes_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hyper_Gminus_all_genes_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hyper_Gminus_all_genes_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hyper_Gminus_all_genes_intersect_df$gene_chr,
      'sourceStart'= top_hyper_Gminus_all_genes_intersect_df$gene_start-1,
      'sourceEnd'= top_hyper_Gminus_all_genes_intersect_df$gene_end,
      'sourceName'= top_hyper_Gminus_all_genes_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hyper_Gminus_all_genes_intersect_df)
      ),
      'targetChrom'= top_hyper_Gminus_all_genes_intersect_df$CpG_chr,
      'targetStart'= top_hyper_Gminus_all_genes_intersect_df$CpG_start,
      'targetEnd'= top_hyper_Gminus_all_genes_intersect_df$CpG_end,
      'targetName'= top_hyper_Gminus_all_genes_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hyper_Gminus_all_genes_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Now that we have the information, let's assemble the output file
    ## for the top TFs:
    top_hyper_Gminus_all_TFs_output_df <- data.frame(
      'chrom'= top_hyper_Gminus_all_TFs_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hyper_Gminus_all_TFs_intersect_df$gene_chr==top_hyper_Gminus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gminus_all_TFs_intersect_df$CpG_start < top_hyper_Gminus_all_TFs_intersect_df$gene_start),
          top_hyper_Gminus_all_TFs_intersect_df$CpG_start,
          top_hyper_Gminus_all_TFs_intersect_df$gene_start
        ),
        top_hyper_Gminus_all_TFs_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hyper_Gminus_all_TFs_intersect_df$gene_chr==top_hyper_Gminus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hyper_Gminus_all_TFs_intersect_df$CpG_start < top_hyper_Gminus_all_TFs_intersect_df$gene_start),
          top_hyper_Gminus_all_TFs_intersect_df$gene_start-1,
          top_hyper_Gminus_all_TFs_intersect_df$CpG_start
        ),
        top_hyper_Gminus_all_TFs_intersect_df$gene_end
      ),
      'name'= paste(
        top_hyper_Gminus_all_TFs_intersect_df$gene_ID,
        top_hyper_Gminus_all_TFs_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hyper_Gminus_all_TFs_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hyper_Gminus_all_TFs_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hyper_Gminus_all_TFs_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hyper_Gminus_all_TFs_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hyper_Gminus_all_TFs_intersect_df$gene_chr,
      'sourceStart'= top_hyper_Gminus_all_TFs_intersect_df$gene_start-1,
      'sourceEnd'= top_hyper_Gminus_all_TFs_intersect_df$gene_end,
      'sourceName'= top_hyper_Gminus_all_TFs_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hyper_Gminus_all_TFs_intersect_df)
      ),
      'targetChrom'= top_hyper_Gminus_all_TFs_intersect_df$CpG_chr,
      'targetStart'= top_hyper_Gminus_all_TFs_intersect_df$CpG_start,
      'targetEnd'= top_hyper_Gminus_all_TFs_intersect_df$CpG_end,
      'targetName'= top_hyper_Gminus_all_TFs_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hyper_Gminus_all_TFs_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Create text for the header lines:
    top_hyper_Gminus_all_genes_header_text <- "track type=interact name=\"TENETR_hyper.G+_interactions\" description=\"TENETR top genes to enhancer DNA methylation probe links\""

    top_hyper_Gminus_all_TFs_header_text <- "track type=interact name=\"TENETR_hyper.G+_interactions\" description=\"TENETR top TFs to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    top_hyper_Gminus_all_genes_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gminus_ucsc_bed_files/',
      'top_hyper_Gminus_genes_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    top_hyper_Gminus_all_TFs_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gminus_ucsc_bed_files/',
      'top_hyper_Gminus_TFs_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    ## Add the header line to the new bed files:
    cat(
      top_hyper_Gminus_all_genes_header_text,
      "\n",
      file = top_hyper_Gminus_all_genes_bed_file_name
    )

    cat(
      top_hyper_Gminus_all_TFs_header_text,
      "\n",
      file = top_hyper_Gminus_all_TFs_bed_file_name
    )

    ## Write the info to the files:
    write.table(
      top_hyper_Gminus_all_genes_output_df,
      file = top_hyper_Gminus_all_genes_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
    )

    write.table(
      top_hyper_Gminus_all_TFs_output_df,
      file = top_hyper_Gminus_all_TFs_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
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
          'hypo_Gplus_ucsc_bed_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_ucsc_bed_files',
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

    ## Now create dataframes listing the top genes/TFs in one column (repeated)
    ## And all the probes linked to them in a second column:

    ## Index empty dataframes:
    top_hypo_Gplus_all_genes_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    top_hypo_Gplus_all_TFs_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    ## Add the probes linked to each of the top genes to the dataframe:
    for(gene_ENSG in top_hypo_Gplus_all_gene_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          gene_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hypo_Gplus_all_genes_intersect_df <- rbind(
        top_hypo_Gplus_all_genes_intersect_df,
        temp_df
      )
    }

    ## Add the probes linked to each of the top TFs to the dataframe:
    for(TF_ENSG in top_hypo_Gplus_all_TF_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID==TF_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          TF_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hypo_Gplus_all_TFs_intersect_df <- rbind(
        top_hypo_Gplus_all_TFs_intersect_df,
        temp_df
      )
    }

    ## For each of the genes in the intersect dfs
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    top_hypo_Gplus_all_genes_intersect_df$gene_chr <- gencode_v22_genes[
      top_hypo_Gplus_all_genes_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hypo_Gplus_all_genes_intersect_df$gene_start <- gencode_v22_genes[
      top_hypo_Gplus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gplus_all_genes_intersect_df$gene_end <- gencode_v22_genes[
      top_hypo_Gplus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$gene_chr <- gencode_v22_genes[
      top_hypo_Gplus_all_TFs_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$gene_start <- gencode_v22_genes[
      top_hypo_Gplus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$gene_end <- gencode_v22_genes[
      top_hypo_Gplus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    ## For each of the probes in the intersect dfs
    ## get the chromosome, "start", and "end" of the probe:
    top_hypo_Gplus_all_genes_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_genes_intersect_df$CpG_ID,
      'chr'
    ]

    top_hypo_Gplus_all_genes_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_genes_intersect_df$CpG_ID,
      'start'
    ]

    top_hypo_Gplus_all_genes_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_genes_intersect_df$CpG_ID,
      'end'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_TFs_intersect_df$CpG_ID,
      'chr'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_TFs_intersect_df$CpG_ID,
      'start'
    ]

    top_hypo_Gplus_all_TFs_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gplus_all_TFs_intersect_df$CpG_ID,
      'end'
    ]

    ## Do some numeric conversions:
    top_hypo_Gplus_all_genes_intersect_df$gene_start <- as.integer(top_hypo_Gplus_all_genes_intersect_df$gene_start)
    top_hypo_Gplus_all_genes_intersect_df$gene_end <- as.integer(top_hypo_Gplus_all_genes_intersect_df$gene_end)

    top_hypo_Gplus_all_TFs_intersect_df$gene_start <- as.integer(top_hypo_Gplus_all_TFs_intersect_df$gene_start)
    top_hypo_Gplus_all_TFs_intersect_df$gene_end <- as.integer(top_hypo_Gplus_all_TFs_intersect_df$gene_end)

    ## Now let's use the rainbow color function to set up a gradient of colors equal to
    ## the number of genes/TFs analyzed and assign each gene/TF a color:
    jet_color_numeric_color_grad_genes <- rainbow(top_gene_number)
    jet_color_numeric_color_grad_TFs <- rainbow(top_gene_number)

    names(jet_color_numeric_color_grad_genes) <- top_hypo_Gplus_all_gene_ENSG
    names(jet_color_numeric_color_grad_TFs) <- top_hypo_Gplus_all_TF_ENSG

    ## Now that we have the information, let's assemble the output file
    ## for the top genes:
    top_hypo_Gplus_all_genes_output_df <- data.frame(
      'chrom'= top_hypo_Gplus_all_genes_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hypo_Gplus_all_genes_intersect_df$gene_chr==top_hypo_Gplus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gplus_all_genes_intersect_df$CpG_start < top_hypo_Gplus_all_genes_intersect_df$gene_start),
          top_hypo_Gplus_all_genes_intersect_df$CpG_start,
          top_hypo_Gplus_all_genes_intersect_df$gene_start
        ),
        top_hypo_Gplus_all_genes_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hypo_Gplus_all_genes_intersect_df$gene_chr==top_hypo_Gplus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gplus_all_genes_intersect_df$CpG_start < top_hypo_Gplus_all_genes_intersect_df$gene_start),
          top_hypo_Gplus_all_genes_intersect_df$gene_start-1,
          top_hypo_Gplus_all_genes_intersect_df$CpG_start
        ),
        top_hypo_Gplus_all_genes_intersect_df$gene_end
      ),
      'name'= paste(
        top_hypo_Gplus_all_genes_intersect_df$gene_ID,
        top_hypo_Gplus_all_genes_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hypo_Gplus_all_genes_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hypo_Gplus_all_genes_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hypo_Gplus_all_genes_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hypo_Gplus_all_genes_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hypo_Gplus_all_genes_intersect_df$gene_chr,
      'sourceStart'= top_hypo_Gplus_all_genes_intersect_df$gene_start-1,
      'sourceEnd'= top_hypo_Gplus_all_genes_intersect_df$gene_end,
      'sourceName'= top_hypo_Gplus_all_genes_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hypo_Gplus_all_genes_intersect_df)
      ),
      'targetChrom'= top_hypo_Gplus_all_genes_intersect_df$CpG_chr,
      'targetStart'= top_hypo_Gplus_all_genes_intersect_df$CpG_start,
      'targetEnd'= top_hypo_Gplus_all_genes_intersect_df$CpG_end,
      'targetName'= top_hypo_Gplus_all_genes_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hypo_Gplus_all_genes_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Now that we have the information, let's assemble the output file
    ## for the top TFs:
    top_hypo_Gplus_all_TFs_output_df <- data.frame(
      'chrom'= top_hypo_Gplus_all_TFs_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hypo_Gplus_all_TFs_intersect_df$gene_chr==top_hypo_Gplus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gplus_all_TFs_intersect_df$CpG_start < top_hypo_Gplus_all_TFs_intersect_df$gene_start),
          top_hypo_Gplus_all_TFs_intersect_df$CpG_start,
          top_hypo_Gplus_all_TFs_intersect_df$gene_start
        ),
        top_hypo_Gplus_all_TFs_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hypo_Gplus_all_TFs_intersect_df$gene_chr==top_hypo_Gplus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gplus_all_TFs_intersect_df$CpG_start < top_hypo_Gplus_all_TFs_intersect_df$gene_start),
          top_hypo_Gplus_all_TFs_intersect_df$gene_start-1,
          top_hypo_Gplus_all_TFs_intersect_df$CpG_start
        ),
        top_hypo_Gplus_all_TFs_intersect_df$gene_end
      ),
      'name'= paste(
        top_hypo_Gplus_all_TFs_intersect_df$gene_ID,
        top_hypo_Gplus_all_TFs_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hypo_Gplus_all_TFs_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hypo_Gplus_all_TFs_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hypo_Gplus_all_TFs_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hypo_Gplus_all_TFs_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hypo_Gplus_all_TFs_intersect_df$gene_chr,
      'sourceStart'= top_hypo_Gplus_all_TFs_intersect_df$gene_start-1,
      'sourceEnd'= top_hypo_Gplus_all_TFs_intersect_df$gene_end,
      'sourceName'= top_hypo_Gplus_all_TFs_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hypo_Gplus_all_TFs_intersect_df)
      ),
      'targetChrom'= top_hypo_Gplus_all_TFs_intersect_df$CpG_chr,
      'targetStart'= top_hypo_Gplus_all_TFs_intersect_df$CpG_start,
      'targetEnd'= top_hypo_Gplus_all_TFs_intersect_df$CpG_end,
      'targetName'= top_hypo_Gplus_all_TFs_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hypo_Gplus_all_TFs_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Create text for the header lines:
    top_hypo_Gplus_all_genes_header_text <- "track type=interact name=\"TENETR_hypo.G+_interactions\" description=\"TENETR top genes to enhancer DNA methylation probe links\""

    top_hypo_Gplus_all_TFs_header_text <- "track type=interact name=\"TENETR_hypo.G+_interactions\" description=\"TENETR top TFs to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    top_hypo_Gplus_all_genes_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gplus_ucsc_bed_files/',
      'top_hypo_Gplus_genes_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    top_hypo_Gplus_all_TFs_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gplus_ucsc_bed_files/',
      'top_hypo_Gplus_TFs_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    ## Add the header line to the new bed files:
    cat(
      top_hypo_Gplus_all_genes_header_text,
      "\n",
      file = top_hypo_Gplus_all_genes_bed_file_name
    )

    cat(
      top_hypo_Gplus_all_TFs_header_text,
      "\n",
      file = top_hypo_Gplus_all_TFs_bed_file_name
    )

    ## Write the info to the files:
    write.table(
      top_hypo_Gplus_all_genes_output_df,
      file = top_hypo_Gplus_all_genes_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
    )

    write.table(
      top_hypo_Gplus_all_TFs_output_df,
      file = top_hypo_Gplus_all_TFs_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
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
          'hypo_Gminus_ucsc_bed_files',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_ucsc_bed_files',
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

    ## Now create dataframes listing the top genes/TFs in one column (repeated)
    ## And all the probes linked to them in a second column:

    ## Index empty dataframes:
    top_hypo_Gminus_all_genes_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    top_hypo_Gminus_all_TFs_intersect_df <- data.frame(
      'gene_ID'=character(),
      'gene_name'=character(),
      'CpG_ID'=character(),
      stringsAsFactors = FALSE
    )

    ## Add the probes linked to each of the top genes to the dataframe:
    for(gene_ENSG in top_hypo_Gminus_all_gene_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          gene_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hypo_Gminus_all_genes_intersect_df <- rbind(
        top_hypo_Gminus_all_genes_intersect_df,
        temp_df
      )
    }

    ## Add the probes linked to each of the top TFs to the dataframe:
    for(TF_ENSG in top_hypo_Gminus_all_TF_ENSG){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        TF_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given TF:
      linked_CpGs <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID==TF_ENSG,
          'probeID'
        ]
      )

      ## Create temporary dataframe with the 3 vectors of info:
      temp_df <- data.frame(
        'gene_ID'= rep(
          TF_ENSG,
          length(linked_CpGs)
        ),
        'gene_name'= rep(
          gene_name,
          length(linked_CpGs)
        ),
        'CpG_ID'= linked_CpGs,
        stringsAsFactors = FALSE
      )

      ## Bind the temp_df to the intersect_df:
      top_hypo_Gminus_all_TFs_intersect_df <- rbind(
        top_hypo_Gminus_all_TFs_intersect_df,
        temp_df
      )
    }

    ## For each of the genes in the intersect dfs
    ## get the chromosome, "start", and "end" of the gene:
    ## Start and end will both be the TSS
    top_hypo_Gminus_all_genes_intersect_df$gene_chr <- gencode_v22_genes[
      top_hypo_Gminus_all_genes_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hypo_Gminus_all_genes_intersect_df$gene_start <- gencode_v22_genes[
      top_hypo_Gminus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gminus_all_genes_intersect_df$gene_end <- gencode_v22_genes[
      top_hypo_Gminus_all_genes_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$gene_chr <- gencode_v22_genes[
      top_hypo_Gminus_all_TFs_intersect_df$gene_ID,
      'seqnames'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$gene_start <- gencode_v22_genes[
      top_hypo_Gminus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$gene_end <- gencode_v22_genes[
      top_hypo_Gminus_all_TFs_intersect_df$gene_ID,
      'TSS'
    ]

    ## For each of the probes in the intersect dfs
    ## get the chromosome, "start", and "end" of the probe:
    top_hypo_Gminus_all_genes_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_genes_intersect_df$CpG_ID,
      'chr'
    ]

    top_hypo_Gminus_all_genes_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_genes_intersect_df$CpG_ID,
      'start'
    ]

    top_hypo_Gminus_all_genes_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_genes_intersect_df$CpG_ID,
      'end'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$CpG_chr <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_TFs_intersect_df$CpG_ID,
      'chr'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$CpG_start <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_TFs_intersect_df$CpG_ID,
      'start'
    ]

    top_hypo_Gminus_all_TFs_intersect_df$CpG_end <- hg38_manifest_no_NA_granges_df[
      top_hypo_Gminus_all_TFs_intersect_df$CpG_ID,
      'end'
    ]

    ## Do some numeric conversions:
    top_hypo_Gminus_all_genes_intersect_df$gene_start <- as.integer(top_hypo_Gminus_all_genes_intersect_df$gene_start)
    top_hypo_Gminus_all_genes_intersect_df$gene_end <- as.integer(top_hypo_Gminus_all_genes_intersect_df$gene_end)

    top_hypo_Gminus_all_TFs_intersect_df$gene_start <- as.integer(top_hypo_Gminus_all_TFs_intersect_df$gene_start)
    top_hypo_Gminus_all_TFs_intersect_df$gene_end <- as.integer(top_hypo_Gminus_all_TFs_intersect_df$gene_end)

    ## Now let's use the rainbow color function to set up a gradient of colors equal to
    ## the number of genes/TFs analyzed and assign each gene/TF a color:
    jet_color_numeric_color_grad_genes <- rainbow(top_gene_number)
    jet_color_numeric_color_grad_TFs <- rainbow(top_gene_number)

    names(jet_color_numeric_color_grad_genes) <- top_hypo_Gminus_all_gene_ENSG
    names(jet_color_numeric_color_grad_TFs) <- top_hypo_Gminus_all_TF_ENSG

    ## Now that we have the information, let's assemble the output file
    ## for the top genes:
    top_hypo_Gminus_all_genes_output_df <- data.frame(
      'chrom'= top_hypo_Gminus_all_genes_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hypo_Gminus_all_genes_intersect_df$gene_chr==top_hypo_Gminus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gminus_all_genes_intersect_df$CpG_start < top_hypo_Gminus_all_genes_intersect_df$gene_start),
          top_hypo_Gminus_all_genes_intersect_df$CpG_start,
          top_hypo_Gminus_all_genes_intersect_df$gene_start
        ),
        top_hypo_Gminus_all_genes_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hypo_Gminus_all_genes_intersect_df$gene_chr==top_hypo_Gminus_all_genes_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gminus_all_genes_intersect_df$CpG_start < top_hypo_Gminus_all_genes_intersect_df$gene_start),
          top_hypo_Gminus_all_genes_intersect_df$gene_start-1,
          top_hypo_Gminus_all_genes_intersect_df$CpG_start
        ),
        top_hypo_Gminus_all_genes_intersect_df$gene_end
      ),
      'name'= paste(
        top_hypo_Gminus_all_genes_intersect_df$gene_ID,
        top_hypo_Gminus_all_genes_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hypo_Gminus_all_genes_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hypo_Gminus_all_genes_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hypo_Gminus_all_genes_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hypo_Gminus_all_genes_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hypo_Gminus_all_genes_intersect_df$gene_chr,
      'sourceStart'= top_hypo_Gminus_all_genes_intersect_df$gene_start-1,
      'sourceEnd'= top_hypo_Gminus_all_genes_intersect_df$gene_end,
      'sourceName'= top_hypo_Gminus_all_genes_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hypo_Gminus_all_genes_intersect_df)
      ),
      'targetChrom'= top_hypo_Gminus_all_genes_intersect_df$CpG_chr,
      'targetStart'= top_hypo_Gminus_all_genes_intersect_df$CpG_start,
      'targetEnd'= top_hypo_Gminus_all_genes_intersect_df$CpG_end,
      'targetName'= top_hypo_Gminus_all_genes_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hypo_Gminus_all_genes_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Now that we have the information, let's assemble the output file
    ## for the top TFs:
    top_hypo_Gminus_all_TFs_output_df <- data.frame(
      'chrom'= top_hypo_Gminus_all_TFs_intersect_df$gene_chr,
      'chromStart'= ifelse(
        top_hypo_Gminus_all_TFs_intersect_df$gene_chr==top_hypo_Gminus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gminus_all_TFs_intersect_df$CpG_start < top_hypo_Gminus_all_TFs_intersect_df$gene_start),
          top_hypo_Gminus_all_TFs_intersect_df$CpG_start,
          top_hypo_Gminus_all_TFs_intersect_df$gene_start
        ),
        top_hypo_Gminus_all_TFs_intersect_df$gene_start-1
      ),
      'chromEnd'= ifelse(
        top_hypo_Gminus_all_TFs_intersect_df$gene_chr==top_hypo_Gminus_all_TFs_intersect_df$CpG_chr,
        ifelse(
          (top_hypo_Gminus_all_TFs_intersect_df$CpG_start < top_hypo_Gminus_all_TFs_intersect_df$gene_start),
          top_hypo_Gminus_all_TFs_intersect_df$gene_start-1,
          top_hypo_Gminus_all_TFs_intersect_df$CpG_start
        ),
        top_hypo_Gminus_all_TFs_intersect_df$gene_end
      ),
      'name'= paste(
        top_hypo_Gminus_all_TFs_intersect_df$gene_ID,
        top_hypo_Gminus_all_TFs_intersect_df$CpG_ID,
        'link',
        sep='_'
      ),
      'score'= rep(
        0,
        nrow(top_hypo_Gminus_all_TFs_intersect_df)
      ),
      'value'= rep(
        0,
        nrow(top_hypo_Gminus_all_TFs_intersect_df)
      ),
      'exp'= rep(
        '.',
        nrow(top_hypo_Gminus_all_TFs_intersect_df)
      ),
      'color'= substring(
        jet_color_numeric_color_grad_genes[
          top_hypo_Gminus_all_TFs_intersect_df$gene_ID
        ],
        1,
        7
      ),
      'sourceChrom'= top_hypo_Gminus_all_TFs_intersect_df$gene_chr,
      'sourceStart'= top_hypo_Gminus_all_TFs_intersect_df$gene_start-1,
      'sourceEnd'= top_hypo_Gminus_all_TFs_intersect_df$gene_end,
      'sourceName'= top_hypo_Gminus_all_TFs_intersect_df$gene_name,
      'sourceStrand'= rep(
        '.',
        nrow(top_hypo_Gminus_all_TFs_intersect_df)
      ),
      'targetChrom'= top_hypo_Gminus_all_TFs_intersect_df$CpG_chr,
      'targetStart'= top_hypo_Gminus_all_TFs_intersect_df$CpG_start,
      'targetEnd'= top_hypo_Gminus_all_TFs_intersect_df$CpG_end,
      'targetName'= top_hypo_Gminus_all_TFs_intersect_df$CpG_ID,
      'targetStrand'= rep(
        '.',
        nrow(top_hypo_Gminus_all_TFs_intersect_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Create text for the header lines:
    top_hypo_Gminus_all_genes_header_text <- "track type=interact name=\"TENETR_hypo.G+_interactions\" description=\"TENETR top genes to enhancer DNA methylation probe links\""

    top_hypo_Gminus_all_TFs_header_text <- "track type=interact name=\"TENETR_hypo.G+_interactions\" description=\"TENETR top TFs to enhancer DNA methylation probe links\""

    ## Create a file name for the output bed file:
    top_hypo_Gminus_all_genes_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gminus_ucsc_bed_files/',
      'top_hypo_Gminus_genes_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    top_hypo_Gminus_all_TFs_bed_file_name <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gminus_ucsc_bed_files/',
      'top_hypo_Gminus_TFs_to_enhancer_probe_links_hg38.bed',
      sep=''
    )

    ## Add the header line to the new bed files:
    cat(
      top_hypo_Gminus_all_genes_header_text,
      "\n",
      file = top_hypo_Gminus_all_genes_bed_file_name
    )

    cat(
      top_hypo_Gminus_all_TFs_header_text,
      "\n",
      file = top_hypo_Gminus_all_TFs_bed_file_name
    )

    ## Write the info to the files:
    write.table(
      top_hypo_Gminus_all_genes_output_df,
      file = top_hypo_Gminus_all_genes_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
    )

    write.table(
      top_hypo_Gminus_all_TFs_output_df,
      file = top_hypo_Gminus_all_TFs_bed_file_name,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE,
      quote= FALSE
    )

  }

}

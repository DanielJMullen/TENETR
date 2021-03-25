#' step7_top_genes_circos
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates circos plots for each gene showing the genomic links between
#' each gene and each enhancer probe linked to them for the analysis types specified.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_circos' containing the results of this function.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create circos plots displaying links between the top genes/TFs by most hypermeth probes with G+ links and their linked enhancer probes of that type.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create circos plots displaying links between the top genes/TFs by most hypermeth probes with G- links and their linked enhancer probes of that type.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create circos plots displaying links between the top genes/TFs by most hypometh probes with G+ links and their linked enhancer probes of that type.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create circos plots displaying links between the top genes/TFs by most hypometh probes with G- links and their linked enhancer probes of that type.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes of a given analysis type to generate circos plots for.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns .html files with the circos plots visualizing the links between the top genes/TFs and their linked enhancer probes for the selected analysis types.
#' @export

step7_top_genes_circos <- function(
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

  ## Load the chromosome sizes file:
  chrom_sizes <- TENETR.data::hg38_chrom_sizes

  ## Ensure chrom sizes are ordered from largest to smallest:
  chrom_sizes <- chrom_sizes[
    order(chrom_sizes$size, decreasing = TRUE),
  ]

  ## Get only the first 24 entries (extras are weird constructs):
  chrom_sizes <- chrom_sizes[1:24,]

  ## Get just the chromosome names:
  chrom_names <- chrom_sizes$chromosome

  ## Remove the 'chr' from the start of the chrom_names
  ## to match the annotation used by BioCircos
  chrom_names <- sub('.*\\chr', '', chrom_names)

  ## Get the chromosome sizes only as a vector now:
  chrom_sizes <- chrom_sizes$size

  ## Set chrom_names as the names for the chrom_sizes:
  names(chrom_sizes) <- chrom_names

  ## Remove X and Y chromosomes :
  chrom_sizes <- chrom_sizes[as.character(c(1:22))]

  ## Create a genome build for later use with BioCircos:
  circos_genome <- setNames(
    as.list(chrom_sizes),
    names(chrom_sizes)
  )

  ## remove uneeded datasets:
  rm(gencode_v22_gtf)
  rm(chrom_names)

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus circos plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_circos',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_circos',
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
      stop('hyper_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hyper_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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

    ## Write a function to generate and save a circos plot
    ## for a given Tgene to each of its linked enhancer probes:
    hyper_Gplus_circos_function <- function(
      gene_ENSG
    ){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ensembl_ID==gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given gene:
      linked_CpGs <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## List the promoter chromosome and start for gene
      ## a number of times equal to linked CpGs:
      circos_start_chromosome <- rep(
        as.character(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'seqnames'
          ]
        ),
        length(linked_CpGs)
      )

      circos_start_position <- rep(
        as.numeric(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'TSS'
          ]
        ),
        length(linked_CpGs)
      )

      ## List the chromosome and start for probes
      ## linked to given gene as end points on circos:
      circos_end_chromosome <- hg38_manifest_no_NA_granges_df[
        linked_CpGs,
        'chr'
      ]

      circos_end_position <- (
        hg38_manifest_no_NA_granges_df[
          linked_CpGs,
          'start'
        ]+1
      )

      ## Convert chromosomes to individual chromsome names
      ## without 'chr'
      circos_start_chromosome <- substring(
        circos_start_chromosome,
        4,
        nchar(circos_start_chromosome)
      )

      circos_end_chromosome <- substring(
        circos_end_chromosome,
        4,
        nchar(circos_end_chromosome)
      )

      ## Get the path to the step 7 folder to save
      ## the pdf in:
      path_to_folder <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_circos',
        sep=''
      )

      ## Create a title for the circos plot html:
      circos_plot_html_title <- paste(
        path_to_folder,
        '/',
        gene_name,
        '_circos_plot.html',
        sep=''
      )

      ## Create a BioCircos tracklist:
      tracklist_gene <- BioCircos::BioCircosBackgroundTrack(
        "myBackgroundTrack",
        minRadius = 0,
        maxRadius = 0.95,
        borderSize = 0,
        fillColors = "#FFFFFF"
      )

      tracklist_gene <- tracklist_gene +
        BioCircos::BioCircosLinkTrack(
          'myLinkTrack',
          circos_start_chromosome,
          circos_start_position,
          circos_start_position + 10000,
          circos_end_chromosome,
          circos_end_position,
          circos_end_position + 10000,
          maxRadius = 0.95,
          color='red'
        )

      ## Plot the BioCircos track:
      x <- BioCircos::BioCircos(
        tracklist_gene,
        genomeFillColor = "PuOr",
        chrPad = 0.02,
        displayGenomeBorder = FALSE,
        yChr = TRUE,
        genomeTicksDisplay = FALSE,
        genomeLabelTextSize = "8pt",
        genomeLabelDy = 0,
        genome = circos_genome
      )

      ## Save the circos plot html file:
      htmlwidgets::saveWidget(
        x,
        circos_plot_html_title
      )

      ## Clear files created by saveWidget:
      unlink(
        paste(
          path_to_folder,
          '/',
          gene_name,
          '_circos_plot_files',
          sep=''
        ),
        recursive=TRUE
      )

      ## clear plots:
      rm(x)
    }

    ## lapply the function on the genes of interest:
    parallel::mclapply(
      top_hyper_Gplus_all_gene_ENSG,
      hyper_Gplus_circos_function,
      mc.cores= core_count
    )

    parallel::mclapply(
      top_hyper_Gplus_all_TF_ENSG,
      hyper_Gplus_circos_function,
      mc.cores= core_count
    )

  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gminus circos plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_circos',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_circos',
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
      stop('hyper_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

    }

    ## Check that the hyper_Gminus_links_TF_gene_freq.txt file exists:
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
      stop('hyper_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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

    ## Write a function to generate and save a circos plot
    ## for a given Tgene to each of its linked enhancer probes:
    hyper_Gminus_circos_function <- function(
      gene_ENSG
    ){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ensembl_ID==gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given gene:
      linked_CpGs <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## List the promoter chromosome and start for gene
      ## a number of times equal to linked CpGs:
      circos_start_chromosome <- rep(
        as.character(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'seqnames'
          ]
        ),
        length(linked_CpGs)
      )

      circos_start_position <- rep(
        as.numeric(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'TSS'
          ]
        ),
        length(linked_CpGs)
      )

      ## List the chromosome and start for probes
      ## linked to given gene as end points on circos:
      circos_end_chromosome <- hg38_manifest_no_NA_granges_df[
        linked_CpGs,
        'chr'
      ]

      circos_end_position <- (
        hg38_manifest_no_NA_granges_df[
          linked_CpGs,
          'start'
        ]+1
      )

      ## Convert chromosomes to individual chromsome names
      ## without 'chr'
      circos_start_chromosome <- substring(
        circos_start_chromosome,
        4,
        nchar(circos_start_chromosome)
      )

      circos_end_chromosome <- substring(
        circos_end_chromosome,
        4,
        nchar(circos_end_chromosome)
      )

      ## Get the path to the step 7 folder to save
      ## the pdf in:
      path_to_folder <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_circos',
        sep=''
      )

      ## Create a title for the circos plot html:
      circos_plot_html_title <- paste(
        path_to_folder,
        '/',
        gene_name,
        '_circos_plot.html',
        sep=''
      )

      ## Create a BioCircos tracklist:
      tracklist_gene <- BioCircos::BioCircosBackgroundTrack(
        "myBackgroundTrack",
        minRadius = 0,
        maxRadius = 0.95,
        borderSize = 0,
        fillColors = "#FFFFFF"
      )

      tracklist_gene <- tracklist_gene +
        BioCircos::BioCircosLinkTrack(
          'myLinkTrack',
          circos_start_chromosome,
          circos_start_position,
          circos_start_position + 10000,
          circos_end_chromosome,
          circos_end_position,
          circos_end_position + 10000,
          maxRadius = 0.95,
          color='red'
        )

      ## Plot the BioCircos track:
      x <- BioCircos::BioCircos(
        tracklist_gene,
        genomeFillColor = "PuOr",
        chrPad = 0.02,
        displayGenomeBorder = FALSE,
        yChr = TRUE,
        genomeTicksDisplay = FALSE,
        genomeLabelTextSize = "8pt",
        genomeLabelDy = 0,
        genome = circos_genome
      )

      ## Save the circos plot html file:
      htmlwidgets::saveWidget(
        x,
        circos_plot_html_title
      )

      ## Clear files created by saveWidget:
      unlink(
        paste(
          path_to_folder,
          '/',
          gene_name,
          '_circos_plot_files',
          sep=''
        ),
        recursive=TRUE
      )

      ## clear plots:
      rm(x)
    }

    ## lapply the function on the genes of interest:
    parallel::mclapply(
      top_hyper_Gminus_all_gene_ENSG,
      hyper_Gminus_circos_function,
      mc.cores= core_count
    )

    parallel::mclapply(
      top_hyper_Gminus_all_TF_ENSG,
      hyper_Gminus_circos_function,
      mc.cores= core_count
    )

  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gplus circos plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_circos',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_circos',
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
      stop('hypo_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

    }

    ## Check that the hypo_Gminus_links_TF_gene_freq.txt file exists:
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
      stop('hypo_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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

    ## Write a function to generate and save a circos plot
    ## for a given Tgene to each of its linked enhancer probes:
    hypo_Gplus_circos_function <- function(
      gene_ENSG
    ){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ensembl_ID==gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given gene:
      linked_CpGs <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## List the promoter chromosome and start for gene
      ## a number of times equal to linked CpGs:
      circos_start_chromosome <- rep(
        as.character(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'seqnames'
          ]
        ),
        length(linked_CpGs)
      )

      circos_start_position <- rep(
        as.numeric(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'TSS'
          ]
        ),
        length(linked_CpGs)
      )

      ## List the chromosome and start for probes
      ## linked to given gene as end points on circos:
      circos_end_chromosome <- hg38_manifest_no_NA_granges_df[
        linked_CpGs,
        'chr'
      ]

      circos_end_position <- (
        hg38_manifest_no_NA_granges_df[
          linked_CpGs,
          'start'
        ]+1
      )

      ## Convert chromosomes to individual chromsome names
      ## without 'chr'
      circos_start_chromosome <- substring(
        circos_start_chromosome,
        4,
        nchar(circos_start_chromosome)
      )

      circos_end_chromosome <- substring(
        circos_end_chromosome,
        4,
        nchar(circos_end_chromosome)
      )

      ## Get the path to the step 7 folder to save
      ## the pdf in:
      path_to_folder <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_circos',
        sep=''
      )

      ## Create a title for the circos plot html:
      circos_plot_html_title <- paste(
        path_to_folder,
        '/',
        gene_name,
        '_circos_plot.html',
        sep=''
      )

      ## Create a BioCircos tracklist:
      tracklist_gene <- BioCircos::BioCircosBackgroundTrack(
        "myBackgroundTrack",
        minRadius = 0,
        maxRadius = 0.95,
        borderSize = 0,
        fillColors = "#FFFFFF"
      )

      tracklist_gene <- tracklist_gene +
        BioCircos::BioCircosLinkTrack(
          'myLinkTrack',
          circos_start_chromosome,
          circos_start_position,
          circos_start_position + 10000,
          circos_end_chromosome,
          circos_end_position,
          circos_end_position + 10000,
          maxRadius = 0.95,
          color='red'
        )

      ## Plot the BioCircos track:
      x <- BioCircos::BioCircos(
        tracklist_gene,
        genomeFillColor = "PuOr",
        chrPad = 0.02,
        displayGenomeBorder = FALSE,
        yChr = TRUE,
        genomeTicksDisplay = FALSE,
        genomeLabelTextSize = "8pt",
        genomeLabelDy = 0,
        genome = circos_genome
      )

      ## Save the circos plot html file:
      htmlwidgets::saveWidget(
        x,
        circos_plot_html_title
      )

      ## Clear files created by saveWidget:
      unlink(
        paste(
          path_to_folder,
          '/',
          gene_name,
          '_circos_plot_files',
          sep=''
        ),
        recursive=TRUE
      )

      ## clear plots:
      rm(x)
    }

    ## lapply the function on the genes of interest:
    parallel::mclapply(
      top_hypo_Gplus_all_gene_ENSG,
      hypo_Gplus_circos_function,
      mc.cores= core_count
    )

    parallel::mclapply(
      top_hypo_Gplus_all_TF_ENSG,
      hypo_Gplus_circos_function,
      mc.cores= core_count
    )

  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gminus circos plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_circos',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_circos',
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
      stop('hypo_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

    }

    ## Check that the hypo_Gminus_links_TF_gene_freq.txt file exists:
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
      stop('hypo_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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

    ## Write a function to generate and save a circos plot
    ## for a given Tgene to each of its linked enhancer probes:
    hypo_Gminus_circos_function <- function(
      gene_ENSG
    ){

      ## Get the name of the gene of interest from
      ## the ENSG:
      gene_name <- gencode_v22_genes[
        gencode_v22_genes$ensembl_ID==gene_ENSG,
        'gene_name'
      ]

      ## Get a list of all the CpGs linked to the given gene:
      linked_CpGs <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID==gene_ENSG,
          'probeID'
        ]
      )

      ## List the promoter chromosome and start for gene
      ## a number of times equal to linked CpGs:
      circos_start_chromosome <- rep(
        as.character(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'seqnames'
          ]
        ),
        length(linked_CpGs)
      )

      circos_start_position <- rep(
        as.numeric(
          gencode_v22_genes[
            gencode_v22_genes$ensembl_ID==gene_ENSG,
            'TSS'
          ]
        ),
        length(linked_CpGs)
      )

      ## List the chromosome and start for probes
      ## linked to given gene as end points on circos:
      circos_end_chromosome <- hg38_manifest_no_NA_granges_df[
        linked_CpGs,
        'chr'
      ]

      circos_end_position <- (
        hg38_manifest_no_NA_granges_df[
          linked_CpGs,
          'start'
        ]+1
      )

      ## Convert chromosomes to individual chromsome names
      ## without 'chr'
      circos_start_chromosome <- substring(
        circos_start_chromosome,
        4,
        nchar(circos_start_chromosome)
      )

      circos_end_chromosome <- substring(
        circos_end_chromosome,
        4,
        nchar(circos_end_chromosome)
      )

      ## Get the path to the step 7 folder to save
      ## the pdf in:
      path_to_folder <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_circos',
        sep=''
      )

      ## Create a title for the circos plot html:
      circos_plot_html_title <- paste(
        path_to_folder,
        '/',
        gene_name,
        '_circos_plot.html',
        sep=''
      )

      ## Create a BioCircos tracklist:
      tracklist_gene <- BioCircos::BioCircosBackgroundTrack(
        "myBackgroundTrack",
        minRadius = 0,
        maxRadius = 0.95,
        borderSize = 0,
        fillColors = "#FFFFFF"
      )

      tracklist_gene <- tracklist_gene +
        BioCircos::BioCircosLinkTrack(
          'myLinkTrack',
          circos_start_chromosome,
          circos_start_position,
          circos_start_position + 10000,
          circos_end_chromosome,
          circos_end_position,
          circos_end_position + 10000,
          maxRadius = 0.95,
          color='red'
        )

      ## Plot the BioCircos track:
      x <- BioCircos::BioCircos(
        tracklist_gene,
        genomeFillColor = "PuOr",
        chrPad = 0.02,
        displayGenomeBorder = FALSE,
        yChr = TRUE,
        genomeTicksDisplay = FALSE,
        genomeLabelTextSize = "8pt",
        genomeLabelDy = 0,
        genome = circos_genome
      )

      ## Save the circos plot html file:
      htmlwidgets::saveWidget(
        x,
        circos_plot_html_title
      )

      ## Clear files created by saveWidget:
      unlink(
        paste(
          path_to_folder,
          '/',
          gene_name,
          '_circos_plot_files',
          sep=''
        ),
        recursive=TRUE
      )

      ## clear plots:
      rm(x)
    }

    ## lapply the function on the genes of interest:
    parallel::mclapply(
      top_hypo_Gminus_all_gene_ENSG,
      hypo_Gminus_circos_function,
      mc.cores= core_count
    )

    parallel::mclapply(
      top_hypo_Gminus_all_TF_ENSG,
      hypo_Gminus_circos_function,
      mc.cores= core_count
    )

  }
}

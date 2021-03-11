#' top_tr_met_heatmaps
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TRs by number of linked probes identified from
#' the step6 top_tr_tabulation function up to the number as specified by the user
#' and generates heatmaps showing the methylation level of the enhancer
#' DNA methylation probes linked to those genes.
#'
#'
#' @param TENET_directory Set a path to the directory that contains step6 results from the top_tr_tabulation function. This function will also create a new step7 folder there if it has not been created, with a subdirectory with 'met_heatmaps' containing the results.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create heatmaps showing DNA methylation levels of probes linked to the top TRs with the most hypermeth probes with G+ links.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of probes linked the top TRs with the most hypermeth probes with G- links.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of probes linked the top TRs with the most hypometh probes with G+ links.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of probes linked the top TRs with the most hypometh probes with G- links.
#' @param top_gene_number Specify a number to generate heatmaps for the probes linked to that number of the top genes/TFs based on the most linked enhancer probes.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns .pdf files with the heatmaps showing the DNA methylation levels of all probes of a specific analysis type linked to those genes, as well as the expression of said genes across the samples in the column labels
#' @export

top_tr_met_heatmaps <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  core_count
){

  ## Load latest version of heatmap.3 function
  devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

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

  ## Load the methylation and expression file from step2:
  if(
    file.exists(
      paste(
        TENET_directory,
        'step2/',
        'diff.methylated.datasets.rda',
        sep=''
      )
    )
  ){

    ## Load the file if it exists:
    load(
      paste(
        TENET_directory,
        'step2/',
        'diff.methylated.datasets.rda',
        sep=''
      )
    )

  } else{

    # Return an error message that the file wasn't found:
    stop('diff.methylated.datasets.rda in step2 of TENET directory was not found. Please check that the file exists and consider rerunning the step2 get_diffmeth_regions function.')
  }

  ## Create function to convert numbers to jet color vectors:
  jet_color_function <- function(numeric_values){

    jet_color_numeric_color_grad <- matlab::jet.colors(200)

    color_values <- jet_color_numeric_color_grad[numeric_values]

    return(color_values)
  }

  ## Write a function to scale expression values across samples
  ## on a proportional scale from 0 to 200, ignoring 0 values:
  rescale_func_zero_ignored <- function(x){

    ## get minimum non-zero, non-NA value:
    non_zero_vec <- x[x!=0]

    non_zero_vec <- non_zero_vec[is.na(non_zero_vec)!=TRUE]

    ## get minimum of non-zero values:
    used_min <- min(non_zero_vec)
    used_max <- max(non_zero_vec)

    place <- ((x-used_min)/(used_max-used_min))*200

    return_value <- ifelse(
      place<=0,
      0.001,
      place
    )

    return(
      ceiling(return_value)
    )
  }

  ## Write function for later to get expression values for
  ## the genes of interest:
  tumor_expression_grabber <- function(gene_id){

    unlist(
      c(
        expDataT[
          c(gene_id),
        ]
      )
    )

  }

  ## Define clustering functions:
  distf=function(d){
    dist(
      d,
      method="euclidean"
    )
  }

  clustf=function(e){
    hclust(
      e,
      method="ward.D2"
    )
  }

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus survival plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_met_heatmaps',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_met_heatmaps',
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
      stop('hyper_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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
      stop('hyper_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

    }

    ## Check that the hyper_Gplus_links_all_TR_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_TR_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gplus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gplus_links_all_TR_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gplus_links_all_TR_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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

    ## Get the names for each of the top genes/TFs:
    top_hyper_Gplus_all_gene_name <- gencode_v22_genes[
      top_hyper_Gplus_all_gene_ENSG,
      'gene_name'
    ]

    top_hyper_Gplus_all_TF_name <- gencode_v22_genes[
      top_hyper_Gplus_all_TF_ENSG,
      'gene_name'
    ]

    ## Get the list of probes linked to the significant genes/TFs:
    top_gene_CpGs_linked <- hyper_Gplus_sig_link_zscores[
      hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_gene_ENSG,
      'probeID'
    ]

    top_TF_CpGs_linked <- hyper_Gplus_sig_link_zscores[
      hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_TF_ENSG,
      'probeID'
    ]

    ## Create methylation datasets from the tumor samples of interest:
    top_gene_CpGs_methylation <- as.matrix(
      metDataT[
        top_gene_CpGs_linked,
      ]
    )

    top_TF_CpGs_methylation <- as.matrix(
      metDataT[
        top_TF_CpGs_linked,
      ]
    )

    ## Get expression from tumor samples for all the genes of interest:
    top_hyper_Gplus_all_gene_expression <- sapply(
      top_hyper_Gplus_all_gene_ENSG,
      tumor_expression_grabber
    )

    top_hyper_Gplus_all_TF_expression <- sapply(
      top_hyper_Gplus_all_TF_ENSG,
      tumor_expression_grabber
    )

    ## Rescale the expression values of the genes of interest:
    top_hyper_Gplus_all_gene_expression_rescaled <- apply(
      top_hyper_Gplus_all_gene_expression,
      2,
      rescale_func_zero_ignored
    )

    top_hyper_Gplus_all_TF_expression_rescaled <- apply(
      top_hyper_Gplus_all_TF_expression,
      2,
      rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    top_hyper_Gplus_all_gene_expression_rescaled_jet_color <- apply(
      top_hyper_Gplus_all_gene_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hyper_Gplus_all_gene_expression_rescaled_jet_color) <- rownames(top_hyper_Gplus_all_gene_expression_rescaled)

    top_hyper_Gplus_all_TF_expression_rescaled_jet_color <- apply(
      top_hyper_Gplus_all_TF_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hyper_Gplus_all_TF_expression_rescaled_jet_color) <- rownames(top_hyper_Gplus_all_TF_expression_rescaled)

    ## Create row labels for the histograms
    ## (blank for now):
    top_hyper_Gplus_all_gene_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_gene_CpGs_methylation)
      )
    )
    rownames(top_hyper_Gplus_all_gene_row_color_labels) <- c(
      ""
    )

    top_hyper_Gplus_all_TF_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_TF_CpGs_methylation)
      )
    )
    rownames(top_hyper_Gplus_all_TF_row_color_labels) <- c(
      ""
    )

    ## Create column labels for the histograms:
    top_hyper_Gplus_all_gene_col_color_labels <- top_hyper_Gplus_all_gene_expression_rescaled_jet_color[
      colnames(top_gene_CpGs_methylation),
      rev(
        colnames(top_hyper_Gplus_all_gene_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hyper_Gplus_all_gene_col_color_labels) <- rev(top_hyper_Gplus_all_gene_name)

    top_hyper_Gplus_all_TF_col_color_labels <- top_hyper_Gplus_all_TF_expression_rescaled_jet_color[
      colnames(top_TF_CpGs_methylation),
      rev(
        colnames(top_hyper_Gplus_all_TF_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hyper_Gplus_all_TF_col_color_labels) <- rev(top_hyper_Gplus_all_TF_name)

    ## Create a natural clustering of the heatmaps:
    top_hyper_Gplus_all_gene_col_dist <- distf(
      t(top_gene_CpGs_methylation)
    )
    top_hyper_Gplus_all_gene_col_clust <- clustf(top_hyper_Gplus_all_gene_col_dist)
    top_hyper_Gplus_all_gene_col_dend <- as.dendrogram(top_hyper_Gplus_all_gene_col_clust)

    top_hyper_Gplus_all_TF_col_dist <- distf(
      t(top_TF_CpGs_methylation)
    )
    top_hyper_Gplus_all_TF_col_clust <- clustf(top_hyper_Gplus_all_TF_col_dist)
    top_hyper_Gplus_all_TF_col_dend <- as.dendrogram(top_hyper_Gplus_all_TF_col_clust)

    ## Create a row clustering for the heatmaps:
    top_hyper_Gplus_all_gene_row_dist <- distf(
      top_gene_CpGs_methylation
    )
    top_hyper_Gplus_all_gene_row_clust <- clustf(top_hyper_Gplus_all_gene_row_dist)
    top_hyper_Gplus_all_gene_row_dend <- as.dendrogram(top_hyper_Gplus_all_gene_row_clust)

    top_hyper_Gplus_all_TF_row_dist <- distf(
      top_TF_CpGs_methylation
    )
    top_hyper_Gplus_all_TF_row_clust <- clustf(top_hyper_Gplus_all_TF_row_dist)
    top_hyper_Gplus_all_TF_row_dend <- as.dendrogram(top_hyper_Gplus_all_TF_row_clust)

    ## Create titles for the heatmap plots:
    top_hyper_Gplus_all_gene_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gplus_met_heatmaps/',
      'hyper_Gplus_top_genes_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    top_hyper_Gplus_all_TF_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gplus_met_heatmaps/',
      'hyper_Gplus_top_TFs_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    ## Open a pdf for saving the all genes plot:
    pdf(
      top_hyper_Gplus_all_gene_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_gene_CpGs_methylation,
      Rowv= top_hyper_Gplus_all_gene_row_dend,
      Colv= top_hyper_Gplus_all_gene_col_dend,
      RowSideColors= top_hyper_Gplus_all_gene_row_color_labels,
      ColSideColors= top_hyper_Gplus_all_gene_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top genes",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()

    ## Open a pdf for saving the all TFs plot:
    pdf(
      top_hyper_Gplus_all_TF_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_TF_CpGs_methylation,
      Rowv= top_hyper_Gplus_all_TF_row_dend,
      Colv= top_hyper_Gplus_all_TF_col_dend,
      RowSideColors= top_hyper_Gplus_all_TF_row_color_labels,
      ColSideColors= top_hyper_Gplus_all_TF_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top transcriptional regulators",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gminus survival plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_met_heatmaps',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_met_heatmaps',
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
      stop('hyper_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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
      stop('hyper_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

    }

    ## Check that the hyper_Gminus_links_TF_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_TR_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hyper_Gminus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hyper_Gminus_links_all_TR_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hyper_Gminus_links_all_TR_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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

    ## Get the names for each of the top genes/TFs:
    top_hyper_Gminus_all_gene_name <- gencode_v22_genes[
      top_hyper_Gminus_all_gene_ENSG,
      'gene_name'
    ]

    top_hyper_Gminus_all_TF_name <- gencode_v22_genes[
      top_hyper_Gminus_all_TF_ENSG,
      'gene_name'
    ]

    ## Get the list of probes linked to the significant genes/TFs:
    top_gene_CpGs_linked <- hyper_Gminus_sig_link_zscores[
      hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_gene_ENSG,
      'probeID'
    ]

    top_TF_CpGs_linked <- hyper_Gminus_sig_link_zscores[
      hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_TF_ENSG,
      'probeID'
    ]

    ## Create methylation datasets from the tumor samples of interest:
    top_gene_CpGs_methylation <- as.matrix(
      metDataT[
        top_gene_CpGs_linked,
      ]
    )

    top_TF_CpGs_methylation <- as.matrix(
      metDataT[
        top_TF_CpGs_linked,
      ]
    )

    ## Get expression from tumor samples for all the genes of interest:
    top_hyper_Gminus_all_gene_expression <- sapply(
      top_hyper_Gminus_all_gene_ENSG,
      tumor_expression_grabber
    )

    top_hyper_Gminus_all_TF_expression <- sapply(
      top_hyper_Gminus_all_TF_ENSG,
      tumor_expression_grabber
    )

    ## Rescale the expression values of the genes of interest:
    top_hyper_Gminus_all_gene_expression_rescaled <- apply(
      top_hyper_Gminus_all_gene_expression,
      2,
      rescale_func_zero_ignored
    )

    top_hyper_Gminus_all_TF_expression_rescaled <- apply(
      top_hyper_Gminus_all_TF_expression,
      2,
      rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    top_hyper_Gminus_all_gene_expression_rescaled_jet_color <- apply(
      top_hyper_Gminus_all_gene_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hyper_Gminus_all_gene_expression_rescaled_jet_color) <- rownames(top_hyper_Gminus_all_gene_expression_rescaled)

    top_hyper_Gminus_all_TF_expression_rescaled_jet_color <- apply(
      top_hyper_Gminus_all_TF_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hyper_Gminus_all_TF_expression_rescaled_jet_color) <- rownames(top_hyper_Gminus_all_TF_expression_rescaled)

    ## Create row labels for the histograms
    ## (blank for now):
    top_hyper_Gminus_all_gene_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_gene_CpGs_methylation)
      )
    )
    rownames(top_hyper_Gminus_all_gene_row_color_labels) <- c(
      ""
    )

    top_hyper_Gminus_all_TF_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_TF_CpGs_methylation)
      )
    )
    rownames(top_hyper_Gminus_all_TF_row_color_labels) <- c(
      ""
    )

    ## Create column labels for the histograms:
    top_hyper_Gminus_all_gene_col_color_labels <- top_hyper_Gminus_all_gene_expression_rescaled_jet_color[
      colnames(top_gene_CpGs_methylation),
      rev(
        colnames(top_hyper_Gminus_all_gene_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hyper_Gminus_all_gene_col_color_labels) <- rev(top_hyper_Gminus_all_gene_name)

    top_hyper_Gminus_all_TF_col_color_labels <- top_hyper_Gminus_all_TF_expression_rescaled_jet_color[
      colnames(top_TF_CpGs_methylation),
      rev(
        colnames(top_hyper_Gminus_all_TF_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hyper_Gminus_all_TF_col_color_labels) <- rev(top_hyper_Gminus_all_TF_name)

    ## Create a natural clustering of the heatmaps:
    top_hyper_Gminus_all_gene_col_dist <- distf(
      t(top_gene_CpGs_methylation)
    )
    top_hyper_Gminus_all_gene_col_clust <- clustf(top_hyper_Gminus_all_gene_col_dist)
    top_hyper_Gminus_all_gene_col_dend <- as.dendrogram(top_hyper_Gminus_all_gene_col_clust)

    top_hyper_Gminus_all_TF_col_dist <- distf(
      t(top_TF_CpGs_methylation)
    )
    top_hyper_Gminus_all_TF_col_clust <- clustf(top_hyper_Gminus_all_TF_col_dist)
    top_hyper_Gminus_all_TF_col_dend <- as.dendrogram(top_hyper_Gminus_all_TF_col_clust)

    ## Create a row clustering for the heatmaps:
    top_hyper_Gminus_all_gene_row_dist <- distf(
      top_gene_CpGs_methylation
    )
    top_hyper_Gminus_all_gene_row_clust <- clustf(top_hyper_Gminus_all_gene_row_dist)
    top_hyper_Gminus_all_gene_row_dend <- as.dendrogram(top_hyper_Gminus_all_gene_row_clust)

    top_hyper_Gminus_all_TF_row_dist <- distf(
      top_TF_CpGs_methylation
    )
    top_hyper_Gminus_all_TF_row_clust <- clustf(top_hyper_Gminus_all_TF_row_dist)
    top_hyper_Gminus_all_TF_row_dend <- as.dendrogram(top_hyper_Gminus_all_TF_row_clust)

    ## Create titles for the heatmap plots:
    top_hyper_Gminus_all_gene_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gminus_met_heatmaps/',
      'hyper_Gminus_top_genes_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    top_hyper_Gminus_all_TF_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hyper_Gminus_met_heatmaps/',
      'hyper_Gminus_top_TFs_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    ## Open a pdf for saving the all genes plot:
    pdf(
      top_hyper_Gminus_all_gene_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_gene_CpGs_methylation,
      Rowv= top_hyper_Gminus_all_gene_row_dend,
      Colv= top_hyper_Gminus_all_gene_col_dend,
      RowSideColors= top_hyper_Gminus_all_gene_row_color_labels,
      ColSideColors= top_hyper_Gminus_all_gene_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top genes",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()

    ## Open a pdf for saving the all TFs plot:
    pdf(
      top_hyper_Gminus_all_TF_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_TF_CpGs_methylation,
      Rowv= top_hyper_Gminus_all_TF_row_dend,
      Colv= top_hyper_Gminus_all_TF_col_dend,
      RowSideColors= top_hyper_Gminus_all_TF_row_color_labels,
      ColSideColors= top_hyper_Gminus_all_TF_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top transcriptional regulators",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gplus survival plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_met_heatmaps',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_met_heatmaps',
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
      stop('hypo_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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
      stop('hypo_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

    }

    ## Check that the hypo_Gminus_links_TF_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_TR_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gplus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gplus_links_all_TR_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gplus_links_all_TR_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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

    ## Get the names for each of the top genes/TFs:
    top_hypo_Gplus_all_gene_name <- gencode_v22_genes[
      top_hypo_Gplus_all_gene_ENSG,
      'gene_name'
    ]

    top_hypo_Gplus_all_TF_name <- gencode_v22_genes[
      top_hypo_Gplus_all_TF_ENSG,
      'gene_name'
    ]

    ## Get the list of probes linked to the significant genes/TFs:
    top_gene_CpGs_linked <- hypo_Gplus_sig_link_zscores[
      hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_gene_ENSG,
      'probeID'
    ]

    top_TF_CpGs_linked <- hypo_Gplus_sig_link_zscores[
      hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_TF_ENSG,
      'probeID'
    ]

    ## Create methylation datasets from the tumor samples of interest:
    top_gene_CpGs_methylation <- as.matrix(
      metDataT[
        top_gene_CpGs_linked,
      ]
    )

    top_TF_CpGs_methylation <- as.matrix(
      metDataT[
        top_TF_CpGs_linked,
      ]
    )

    ## Get expression from tumor samples for all the genes of interest:
    top_hypo_Gplus_all_gene_expression <- sapply(
      top_hypo_Gplus_all_gene_ENSG,
      tumor_expression_grabber
    )

    top_hypo_Gplus_all_TF_expression <- sapply(
      top_hypo_Gplus_all_TF_ENSG,
      tumor_expression_grabber
    )

    ## Rescale the expression values of the genes of interest:
    top_hypo_Gplus_all_gene_expression_rescaled <- apply(
      top_hypo_Gplus_all_gene_expression,
      2,
      rescale_func_zero_ignored
    )

    top_hypo_Gplus_all_TF_expression_rescaled <- apply(
      top_hypo_Gplus_all_TF_expression,
      2,
      rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    top_hypo_Gplus_all_gene_expression_rescaled_jet_color <- apply(
      top_hypo_Gplus_all_gene_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hypo_Gplus_all_gene_expression_rescaled_jet_color) <- rownames(top_hypo_Gplus_all_gene_expression_rescaled)

    top_hypo_Gplus_all_TF_expression_rescaled_jet_color <- apply(
      top_hypo_Gplus_all_TF_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hypo_Gplus_all_TF_expression_rescaled_jet_color) <- rownames(top_hypo_Gplus_all_TF_expression_rescaled)

    ## Create row labels for the histograms
    ## (blank for now):
    top_hypo_Gplus_all_gene_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_gene_CpGs_methylation)
      )
    )
    rownames(top_hypo_Gplus_all_gene_row_color_labels) <- c(
      ""
    )

    top_hypo_Gplus_all_TF_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_TF_CpGs_methylation)
      )
    )
    rownames(top_hypo_Gplus_all_TF_row_color_labels) <- c(
      ""
    )

    ## Create column labels for the histograms:
    top_hypo_Gplus_all_gene_col_color_labels <- top_hypo_Gplus_all_gene_expression_rescaled_jet_color[
      colnames(top_gene_CpGs_methylation),
      rev(
        colnames(top_hypo_Gplus_all_gene_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hypo_Gplus_all_gene_col_color_labels) <- rev(top_hypo_Gplus_all_gene_name)

    top_hypo_Gplus_all_TF_col_color_labels <- top_hypo_Gplus_all_TF_expression_rescaled_jet_color[
      colnames(top_TF_CpGs_methylation),
      rev(
        colnames(top_hypo_Gplus_all_TF_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hypo_Gplus_all_TF_col_color_labels) <- rev(top_hypo_Gplus_all_TF_name)

    ## Create a natural clustering of the heatmaps:
    top_hypo_Gplus_all_gene_col_dist <- distf(
      t(top_gene_CpGs_methylation)
    )
    top_hypo_Gplus_all_gene_col_clust <- clustf(top_hypo_Gplus_all_gene_col_dist)
    top_hypo_Gplus_all_gene_col_dend <- as.dendrogram(top_hypo_Gplus_all_gene_col_clust)

    top_hypo_Gplus_all_TF_col_dist <- distf(
      t(top_TF_CpGs_methylation)
    )
    top_hypo_Gplus_all_TF_col_clust <- clustf(top_hypo_Gplus_all_TF_col_dist)
    top_hypo_Gplus_all_TF_col_dend <- as.dendrogram(top_hypo_Gplus_all_TF_col_clust)

    ## Create a row clustering for the heatmaps:
    top_hypo_Gplus_all_gene_row_dist <- distf(
      top_gene_CpGs_methylation
    )
    top_hypo_Gplus_all_gene_row_clust <- clustf(top_hypo_Gplus_all_gene_row_dist)
    top_hypo_Gplus_all_gene_row_dend <- as.dendrogram(top_hypo_Gplus_all_gene_row_clust)

    top_hypo_Gplus_all_TF_row_dist <- distf(
      top_TF_CpGs_methylation
    )
    top_hypo_Gplus_all_TF_row_clust <- clustf(top_hypo_Gplus_all_TF_row_dist)
    top_hypo_Gplus_all_TF_row_dend <- as.dendrogram(top_hypo_Gplus_all_TF_row_clust)

    ## Create titles for the heatmap plots:
    top_hypo_Gplus_all_gene_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gplus_met_heatmaps/',
      'hypo_Gplus_top_genes_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    top_hypo_Gplus_all_TF_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gplus_met_heatmaps/',
      'hypo_Gplus_top_TFs_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    ## Open a pdf for saving the all genes plot:
    pdf(
      top_hypo_Gplus_all_gene_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_gene_CpGs_methylation,
      Rowv= top_hypo_Gplus_all_gene_row_dend,
      Colv= top_hypo_Gplus_all_gene_col_dend,
      RowSideColors= top_hypo_Gplus_all_gene_row_color_labels,
      ColSideColors= top_hypo_Gplus_all_gene_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top genes",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()

    ## Open a pdf for saving the all TFs plot:
    pdf(
      top_hypo_Gplus_all_TF_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_TF_CpGs_methylation,
      Rowv= top_hypo_Gplus_all_TF_row_dend,
      Colv= top_hypo_Gplus_all_TF_col_dend,
      RowSideColors= top_hypo_Gplus_all_TF_row_color_labels,
      ColSideColors= top_hypo_Gplus_all_TF_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top transcriptional regulators",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gminus survival plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_met_heatmaps',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_met_heatmaps',
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
      stop('hypo_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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
      stop('hypo_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

    }

    ## Check that the hypo_Gminus_links_TF_gene_freq.txt file exists:
    if(
      file.exists(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_TR_freq.txt',
          sep=''
        )
      )
    ){

      ## Load the file if it does:
      hypo_Gminus_all_TF_freq <- read.delim(
        paste(
          TENET_directory,
          'step6/',
          'hypo_Gminus_links_all_TR_freq.txt',
          sep=''
        ),
        header= TRUE,
        stringsAsFactors = FALSE
      )

    } else{

      ## Return an error message that the file wasn't found:
      stop('hypo_Gminus_links_all_TR_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6 top_tr_tabulation function.')

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

    ## Get the names for each of the top genes/TFs:
    top_hypo_Gminus_all_gene_name <- gencode_v22_genes[
      top_hypo_Gminus_all_gene_ENSG,
      'gene_name'
    ]

    top_hypo_Gminus_all_TF_name <- gencode_v22_genes[
      top_hypo_Gminus_all_TF_ENSG,
      'gene_name'
    ]

    ## Get the list of probes linked to the significant genes/TFs:
    top_gene_CpGs_linked <- hypo_Gminus_sig_link_zscores[
      hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_gene_ENSG,
      'probeID'
    ]

    top_TF_CpGs_linked <- hypo_Gminus_sig_link_zscores[
      hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_TF_ENSG,
      'probeID'
    ]

    ## Create methylation datasets from the tumor samples of interest:
    top_gene_CpGs_methylation <- as.matrix(
      metDataT[
        top_gene_CpGs_linked,
      ]
    )

    top_TF_CpGs_methylation <- as.matrix(
      metDataT[
        top_TF_CpGs_linked,
      ]
    )

    ## Get expression from tumor samples for all the genes of interest:
    top_hypo_Gminus_all_gene_expression <- sapply(
      top_hypo_Gminus_all_gene_ENSG,
      tumor_expression_grabber
    )

    top_hypo_Gminus_all_TF_expression <- sapply(
      top_hypo_Gminus_all_TF_ENSG,
      tumor_expression_grabber
    )

    ## Rescale the expression values of the genes of interest:
    top_hypo_Gminus_all_gene_expression_rescaled <- apply(
      top_hypo_Gminus_all_gene_expression,
      2,
      rescale_func_zero_ignored
    )

    top_hypo_Gminus_all_TF_expression_rescaled <- apply(
      top_hypo_Gminus_all_TF_expression,
      2,
      rescale_func_zero_ignored
    )

    ## Convert the rescaled values to jetcolor:
    top_hypo_Gminus_all_gene_expression_rescaled_jet_color <- apply(
      top_hypo_Gminus_all_gene_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hypo_Gminus_all_gene_expression_rescaled_jet_color) <- rownames(top_hypo_Gminus_all_gene_expression_rescaled)

    top_hypo_Gminus_all_TF_expression_rescaled_jet_color <- apply(
      top_hypo_Gminus_all_TF_expression_rescaled,
      2,
      jet_color_function
    )
    rownames(top_hypo_Gminus_all_TF_expression_rescaled_jet_color) <- rownames(top_hypo_Gminus_all_TF_expression_rescaled)

    ## Create row labels for the histograms
    ## (blank for now):
    top_hypo_Gminus_all_gene_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_gene_CpGs_methylation)
      )
    )
    rownames(top_hypo_Gminus_all_gene_row_color_labels) <- c(
      ""
    )

    top_hypo_Gminus_all_TF_row_color_labels <- rbind(
      rep(
        'white',
        nrow(top_TF_CpGs_methylation)
      )
    )
    rownames(top_hypo_Gminus_all_TF_row_color_labels) <- c(
      ""
    )

    ## Create column labels for the histograms:
    top_hypo_Gminus_all_gene_col_color_labels <- top_hypo_Gminus_all_gene_expression_rescaled_jet_color[
      colnames(top_gene_CpGs_methylation),
      rev(
        colnames(top_hypo_Gminus_all_gene_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hypo_Gminus_all_gene_col_color_labels) <- rev(top_hypo_Gminus_all_gene_name)

    top_hypo_Gminus_all_TF_col_color_labels <- top_hypo_Gminus_all_TF_expression_rescaled_jet_color[
      colnames(top_TF_CpGs_methylation),
      rev(
        colnames(top_hypo_Gminus_all_TF_expression_rescaled_jet_color)
      )
    ]
    colnames(top_hypo_Gminus_all_TF_col_color_labels) <- rev(top_hypo_Gminus_all_TF_name)

    ## Create a natural clustering of the heatmaps:
    top_hypo_Gminus_all_gene_col_dist <- distf(
      t(top_gene_CpGs_methylation)
    )
    top_hypo_Gminus_all_gene_col_clust <- clustf(top_hypo_Gminus_all_gene_col_dist)
    top_hypo_Gminus_all_gene_col_dend <- as.dendrogram(top_hypo_Gminus_all_gene_col_clust)

    top_hypo_Gminus_all_TF_col_dist <- distf(
      t(top_TF_CpGs_methylation)
    )
    top_hypo_Gminus_all_TF_col_clust <- clustf(top_hypo_Gminus_all_TF_col_dist)
    top_hypo_Gminus_all_TF_col_dend <- as.dendrogram(top_hypo_Gminus_all_TF_col_clust)

    ## Create a row clustering for the heatmaps:
    top_hypo_Gminus_all_gene_row_dist <- distf(
      top_gene_CpGs_methylation
    )
    top_hypo_Gminus_all_gene_row_clust <- clustf(top_hypo_Gminus_all_gene_row_dist)
    top_hypo_Gminus_all_gene_row_dend <- as.dendrogram(top_hypo_Gminus_all_gene_row_clust)

    top_hypo_Gminus_all_TF_row_dist <- distf(
      top_TF_CpGs_methylation
    )
    top_hypo_Gminus_all_TF_row_clust <- clustf(top_hypo_Gminus_all_TF_row_dist)
    top_hypo_Gminus_all_TF_row_dend <- as.dendrogram(top_hypo_Gminus_all_TF_row_clust)

    ## Create titles for the heatmap plots:
    top_hypo_Gminus_all_gene_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gminus_met_heatmaps/',
      'hypo_Gminus_top_genes_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    top_hypo_Gminus_all_TF_heatmap_title <- paste(
      TENET_directory,
      'step7/',
      'hypo_Gminus_met_heatmaps/',
      'hypo_Gminus_top_TFs_linked_probe_methylation_heatmap.pdf',
      sep=''
    )

    ## Open a pdf for saving the all genes plot:
    pdf(
      top_hypo_Gminus_all_gene_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_gene_CpGs_methylation,
      Rowv= top_hypo_Gminus_all_gene_row_dend,
      Colv= top_hypo_Gminus_all_gene_col_dend,
      RowSideColors= top_hypo_Gminus_all_gene_row_color_labels,
      ColSideColors= top_hypo_Gminus_all_gene_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top genes",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()

    ## Open a pdf for saving the all TFs plot:
    pdf(
      top_hypo_Gminus_all_TF_heatmap_title,
      height= 7,
      width= 10
    )

    ## Create the plot for the top genes of interest:
    heatmap.3(
      x= top_TF_CpGs_methylation,
      Rowv= top_hypo_Gminus_all_TF_row_dend,
      Colv= top_hypo_Gminus_all_TF_col_dend,
      RowSideColors= top_hypo_Gminus_all_TF_row_color_labels,
      ColSideColors= top_hypo_Gminus_all_TF_col_color_labels,
      dendrogram= "col",
      labCol= NA,
      labRow= NA,
      lmat= rbind( c(0,0,5), c(0,0,2), c(4,1,3), c(0,0,6)),
      lwid= c(0.25,0.02,2),
      lhei= c(0.4,0.6,2,0.001),
      margins= c(2,2),
      col= matlab::jet.colors(200),
      trace= "none",
      key= FALSE,
      main= NULL,
      ylab= "Enhancer probes linked to top transcriptional regulators",
      xlab= "Samples"
    )

    ## Close the plot:
    dev.off()
  }
}

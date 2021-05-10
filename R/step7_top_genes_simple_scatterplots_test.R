#' step7_top_genes_simple_scatterplots_test
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates scatterplots displaying the expression level of each of these genes
#' and the methylation level of each enhancer probe linked to them for each of the
#' four hypo or hypermethylated Gplus or Gminus analysis quadrants, as selected by the user.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_simple_scatterplots' containing the results of this function.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create scatterplots for the top genes/TRs by most hypermeth probes with G+ links and each of their linked probes.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create scatterplots for the top genes/TRs by most hypermeth probes with G- links and each of their linked probes.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create scatterplots for the top genes/TRs by most hypometh probes with G+ links and each of their linked probes.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create scatterplots for the top genes/TRs by most hypometh probes with G- links and each of their linked probes.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes of a given analysis type to generate scatterplots for showing expression of the genes and methylation of each of their linked enhancer probes.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns .pdf files with scatterplots showing the expression of the genes of interest on the x-axis and the methylation of the linked probes on the y-axis.
#' @export

step7_top_genes_simple_scatterplots_test <- function(
  TENET_directory,
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

  ## Combine the methylation and expression data from both tumor and normal samples:
  metDataF_subC <- cbind(metDataT, metDataN)
  expDataF_subC <- cbind(expDataT, expDataN)

  ## Create dataframe with color info for the tumor and normal samples
  ## blue for normal, red for tumor data points:
  DichF <- data.frame(
    group=c(
      colnames(metDataT), colnames(metDataN)
    ),
    cluster=c(
      rep(
        "Tumor",
        dim(metDataT)[2]
      ),
      rep(
        "Normal",
        dim(metDataN)[2]
      )
    ),
    stringsAsFactors = FALSE
  )

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus circos plots:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hyper_Gplus_simple_scatterplots
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    }

    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_simple_scatterplots/',
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

    ## Index an empty list for all the top genes
    ## and the TF genes as well:
    top_hyper_Gplus_all_genes_linked_cpgs_list <- list()

    top_hyper_Gplus_all_TF_linked_cpgs_list <- list()

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hyper_Gplus_all_gene_ENSG)
    )){

      ## Get the genes' ENSG:
      gene_ENSG_placeholder <- top_hyper_Gplus_all_gene_ENSG[i]

      ## Get the probes linked to each gene
      probes_linked_to_significant_gene <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_gene_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hyper_Gplus_all_genes_linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hyper_Gplus_all_TF_ENSG)
    )){

      ## Get the TFs ENSG:
      TF_ENSG_placeholder <- top_hyper_Gplus_all_TF_ENSG[i]

      ## Get the probes linked to each TF
      probes_linked_to_significant_TF <- unique(
        hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_TF_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hyper_Gplus_all_TF_linked_cpgs_list[[i]] <- probes_linked_to_significant_TF
    }

    ## Add the names of the genes to the lists:
    names(top_hyper_Gplus_all_genes_linked_cpgs_list) <- top_hyper_Gplus_all_gene_ENSG

    names(top_hyper_Gplus_all_TF_linked_cpgs_list) <- top_hyper_Gplus_all_TF_ENSG

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    hyper_Gplus_scatterplot_function <- function(
      CpGs_list,
      gene_ENSG,
      gene_or_TF
    ){

      ## Unlist the CpGs linked to each probe:
      unlisted_CpGs <- c(
        unlist(
          CpGs_list
        )
      )

      ## Convert the gene ENSG into the gene name:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Get the expression of the gene:
      TF_expression <- c(
        unlist(
          expDataF_subC[
            gene_ENSG,
            DichF$group
          ]
        )
      )

      ## Now write an internal function that will get each linked CpGs methylation,
      ## and use the methylation and gene expression to create a ggplot2 scatterplot
      ## and save it:
      internal_scatterplot_function <- function(
        CpGs_linked_to_TF,
        internal_gene_or_TF
      ){

        ## Save the CpG name:
        CpG_name_placeholder <- CpGs_linked_to_TF

        ## Get DNA methylation values:
        unlisted_CpG_methylation <- c(
          unlist(
            metDataF_subC[
              CpGs_linked_to_TF,
              DichF$group
            ]
          )
        )

        ## Manually coloring samples:
        t_v_n_group_colors <- c(
          'Normal'='dodgerblue3',
          'Tumor'='red3'
        )

        ## Creating scatter with ggplot2:
        scatter_plot <- ggplot2::qplot(
          x=TF_expression,
          y=unlisted_CpG_methylation,
          geom=c("point"),
          colour=DichF$cluster
        )

        ## Create the plot:
        scatter_plot_updated <- scatter_plot +
          ggplot2::ggtitle(
            paste(
              gene_name,
              ' gene expression vs.\n',
              CpG_name_placeholder,
              ' DNA methylation',
              sep=''
            )
          ) +
          ggplot2::ylab(
            paste(
              CpG_name_placeholder,
              'DNA methylation',
              sep=' '
            )
          ) +
          ggplot2::xlab(
            paste(
              gene_name,
              'gene expression',
              sep=' '
            )
          ) +
          ggplot2::theme_bw() +
          ggplot2::scale_color_manual(
            values=t_v_n_group_colors,
            name="Sample type"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust=0.5, size=18),
            legend.title = ggplot2::element_text(hjust=0.5, size=12),
            legend.text = ggplot2::element_text(size=10, colour='black'),
            panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
            axis.title.x = ggplot2::element_text(size=16),
            axis.title.y = ggplot2::element_text(size=16),
            axis.text.x = ggplot2::element_text(size=14, colour='black'),
            axis.text.y = ggplot2::element_text(size=14, colour='black'),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
          )

        ## Determine where to save the results:
        if(internal_gene_or_TF=='gene'){

          subdirectory <- 'top_genes/'

        } else if(internal_gene_or_TF=='TF'){

          subdirectory <- 'top_TFs/'
        }

        ## Create a title for the scatter plot pdf:
        ## This is a comination of the probe name with the linked gene:
        scatterplot_pdf_title <- paste(
          paste(
            TENET_directory,
            'step7/',
            'hyper_Gplus_simple_scatterplots/',
            subdirectory,
            sep=''
          ),
          gene_name,
          '_',
          CpG_name_placeholder,
          '_scatterplot.pdf',
          sep=''
        )

        ## Open a pdf for saving the plot:
        pdf(
          scatterplot_pdf_title,
          height= 7,
          width= 10
        )

        plot(scatter_plot_updated)

        ## Close the plot:
        dev.off()
      }

      lapply(
        X= unlisted_CpGs,
        FUN= internal_scatterplot_function,
        internal_gene_or_TF= gene_or_TF
      )
    }

    ## Generate the plots for all TFs of interest:
    suppressWarnings(
      parallel::mcmapply(
        FUN= hyper_Gplus_scatterplot_function,
        CpGs_list= top_hyper_Gplus_all_genes_linked_cpgs_list,
        gene_ENSG= names(top_hyper_Gplus_all_genes_linked_cpgs_list),
        gene_or_TF= 'gene',
        mc.cores= core_count
      )
    )

    suppressWarnings(
      parallel::mcmapply(
        FUN= hyper_Gplus_scatterplot_function,
        CpGs_list= top_hyper_Gplus_all_TF_linked_cpgs_list,
        gene_ENSG= names(top_hyper_Gplus_all_TF_linked_cpgs_list),
        gene_or_TF= 'TF',
        mc.cores= core_count
      )
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
          'hyper_Gminus_simple_scatterplots',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_simple_scatterplots',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hyper_Gminus_simple_scatterplots
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    }

    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_simple_scatterplots/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_simple_scatterplots/',
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
      stop('hyper_Gminus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function')

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

    ## Index an empty list for all the top genes
    ## and the TF genes as well:
    top_hyper_Gminus_all_genes_linked_cpgs_list <- list()

    top_hyper_Gminus_all_TF_linked_cpgs_list <- list()

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hyper_Gminus_all_gene_ENSG)
    )){

      ## Get the genes' ENSG:
      gene_ENSG_placeholder <- top_hyper_Gminus_all_gene_ENSG[i]

      ## Get the probes linked to each gene
      probes_linked_to_significant_gene <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_gene_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hyper_Gminus_all_genes_linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hyper_Gminus_all_TF_ENSG)
    )){

      ## Get the TFs ENSG:
      TF_ENSG_placeholder <- top_hyper_Gminus_all_TF_ENSG[i]

      ## Get the probes linked to each TF
      probes_linked_to_significant_TF <- unique(
        hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_TF_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hyper_Gminus_all_TF_linked_cpgs_list[[i]] <- probes_linked_to_significant_TF
    }

    ## Add the names of the genes to the lists:
    names(top_hyper_Gminus_all_genes_linked_cpgs_list) <- top_hyper_Gminus_all_gene_ENSG

    names(top_hyper_Gminus_all_TF_linked_cpgs_list) <- top_hyper_Gminus_all_TF_ENSG

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    hyper_Gminus_scatterplot_function <- function(
      CpGs_list,
      gene_ENSG,
      gene_or_TF
    ){

      ## Unlist the CpGs linked to each probe:
      unlisted_CpGs <- c(
        unlist(
          CpGs_list
        )
      )

      ## Convert the gene ENSG into the gene name:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Add the expression of the gene of interest to DichF:
      TF_expression <- c(
        unlist(
          expDataF_subC[
            gene_ENSG,
            DichF$group
          ]
        )
      )

      ## Now write an internal function that will get each linked CpGs methylation,
      ## and use the methylation and gene expression to create a ggplot2 scatterplot
      ## and save it:
      internal_scatterplot_function <- function(
        CpGs_linked_to_TF,
        internal_gene_or_TF
      ){

        ## Save the CpG name:
        CpG_name_placeholder <- CpGs_linked_to_TF

        ## Get DNA methylation values:
        unlisted_CpG_methylation <- c(
          unlist(
            metDataF_subC[
              CpGs_linked_to_TF,
              DichF$group
            ]
          )
        )

        ## Manually coloring samples:
        t_v_n_group_colors <- c(
          'Normal'='dodgerblue3',
          'Tumor'='red3'
        )

        ## Creating scatter with ggplot2:
        scatter_plot <- ggplot2::qplot(
          x=TF_expression,
          y=unlisted_CpG_methylation,
          geom=c("point"),
          colour=DichF$cluster
        )

        ## Create the plot:
        scatter_plot_updated <- scatter_plot +
          ggplot2::ggtitle(
            paste(
              gene_name,
              ' gene expression vs.\n',
              CpG_name_placeholder,
              ' DNA methylation',
              sep=''
            )
          ) +
          ggplot2::ylab(
            paste(
              CpG_name_placeholder,
              'DNA methylation',
              sep=' '
            )
          ) +
          ggplot2::xlab(
            paste(
              gene_name,
              'gene expression',
              sep=' '
            )
          ) +
          ggplot2::theme_bw() +
          ggplot2::scale_color_manual(
            values=t_v_n_group_colors,
            name="Sample type"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust=0.5, size=18),
            legend.title = ggplot2::element_text(hjust=0.5, size=12),
            legend.text = ggplot2::element_text(size=10, colour='black'),
            panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
            axis.title.x = ggplot2::element_text(size=16),
            axis.title.y = ggplot2::element_text(size=16),
            axis.text.x = ggplot2::element_text(size=14, colour='black'),
            axis.text.y = ggplot2::element_text(size=14, colour='black'),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
          )

        ## Determine where to save the results:
        if(internal_gene_or_TF=='gene'){

          subdirectory <- 'top_genes/'

        } else if(internal_gene_or_TF=='TF'){

          subdirectory <- 'top_TFs/'
        }

        ## Create a title for the scatter plot pdf:
        ## This is a comination of the probe name with the linked gene:
        scatterplot_pdf_title <- paste(
          paste(
            TENET_directory,
            'step7/',
            'hyper_Gminus_simple_scatterplots/',
            subdirectory,
            sep=''
          ),
          gene_name,
          '_',
          CpG_name_placeholder,
          '_scatterplot.pdf',
          sep=''
        )

        ## Open a pdf for saving the plot:
        pdf(
          scatterplot_pdf_title,
          height= 7,
          width= 10
        )

        plot(scatter_plot_updated)

        ## Close the plot:
        dev.off()
      }

      lapply(
        X= unlisted_CpGs,
        FUN= internal_scatterplot_function,
        internal_gene_or_TF= gene_or_TF
      )
    }

    ## Generate the plots for all genes/TFs of interest:
    suppressWarnings(
      parallel::mcmapply(
        FUN= hyper_Gminus_scatterplot_function,
        CpGs_list= top_hyper_Gminus_all_genes_linked_cpgs_list,
        gene_ENSG= names(top_hyper_Gminus_all_genes_linked_cpgs_list),
        gene_or_TF= 'gene',
        mc.cores= core_count
      )
    )

    suppressWarnings(
      parallel::mcmapply(
        FUN= hyper_Gminus_scatterplot_function,
        CpGs_list= top_hyper_Gminus_all_TF_linked_cpgs_list,
        gene_ENSG= names(top_hyper_Gminus_all_TF_linked_cpgs_list),
        gene_or_TF= 'TF',
        mc.cores= core_count
      )
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
          'hypo_Gplus_simple_scatterplots',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_simple_scatterplots',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hypo_Gplus_simple_scatterplots
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    }

    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_simple_scatterplots/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_simple_scatterplots/',
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
      stop('hypo_Gplus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function')

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

    ## Index an empty list for all the top genes
    ## and the TF genes as well:
    top_hypo_Gplus_all_genes_linked_cpgs_list <- list()

    top_hypo_Gplus_all_TF_linked_cpgs_list <- list()

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hypo_Gplus_all_gene_ENSG)
    )){

      ## Get the genes' ENSG:
      gene_ENSG_placeholder <- top_hypo_Gplus_all_gene_ENSG[i]

      ## Get the probes linked to each gene
      probes_linked_to_significant_gene <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_gene_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hypo_Gplus_all_genes_linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hypo_Gplus_all_TF_ENSG)
    )){

      ## Get the TFs ENSG:
      TF_ENSG_placeholder <- top_hypo_Gplus_all_TF_ENSG[i]

      ## Get the probes linked to each TF
      probes_linked_to_significant_TF <- unique(
        hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_TF_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hypo_Gplus_all_TF_linked_cpgs_list[[i]] <- probes_linked_to_significant_TF
    }

    ## Add the names of the genes to the lists:
    names(top_hypo_Gplus_all_genes_linked_cpgs_list) <- top_hypo_Gplus_all_gene_ENSG

    names(top_hypo_Gplus_all_TF_linked_cpgs_list) <- top_hypo_Gplus_all_TF_ENSG

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    hypo_Gplus_scatterplot_function <- function(
      CpGs_list,
      gene_ENSG,
      gene_or_TF
    ){

      ## Unlist the CpGs linked to each probe:
      unlisted_CpGs <- c(
        unlist(
          CpGs_list
        )
      )

      ## Convert the gene ENSG into the gene name:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Add the expression of the gene of interest to DichF:
      TF_expression <- c(
        unlist(
          expDataF_subC[
            gene_ENSG,
            DichF$group
          ]
        )
      )

      ## Now write an internal function that will get each linked CpGs methylation,
      ## and use the methylation and gene expression to create a ggplot2 scatterplot
      ## and save it:
      internal_scatterplot_function <- function(
        CpGs_linked_to_TF,
        internal_gene_or_TF
      ){

        ## Save the CpG name:
        CpG_name_placeholder <- CpGs_linked_to_TF

        ## Get DNA methylation values:
        unlisted_CpG_methylation <- c(
          unlist(
            metDataF_subC[
              CpGs_linked_to_TF,
              DichF$group
            ]
          )
        )

        ## Manually coloring samples:
        t_v_n_group_colors <- c(
          'Normal'='dodgerblue3',
          'Tumor'='red3'
        )

        ## Creating scatter with ggplot2:
        scatter_plot <- ggplot2::qplot(
          x=TF_expression,
          y=unlisted_CpG_methylation,
          geom=c("point"),
          colour=DichF$cluster
        )

        ## Create the plot:
        scatter_plot_updated <- scatter_plot +
          ggplot2::ggtitle(
            paste(
              gene_name,
              ' gene expression vs.\n',
              CpG_name_placeholder,
              ' DNA methylation',
              sep=''
            )
          ) +
          ggplot2::ylab(
            paste(
              CpG_name_placeholder,
              'DNA methylation',
              sep=' '
            )
          ) +
          ggplot2::xlab(
            paste(
              gene_name,
              'gene expression',
              sep=' '
            )
          ) +
          ggplot2::theme_bw() +
          ggplot2::scale_color_manual(
            values=t_v_n_group_colors,
            name="Sample type"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust=0.5, size=18),
            legend.title = ggplot2::element_text(hjust=0.5, size=12),
            legend.text = ggplot2::element_text(size=10, colour='black'),
            panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
            axis.title.x = ggplot2::element_text(size=16),
            axis.title.y = ggplot2::element_text(size=16),
            axis.text.x = ggplot2::element_text(size=14, colour='black'),
            axis.text.y = ggplot2::element_text(size=14, colour='black'),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
          )

        ## Determine where to save the results:
        if(internal_gene_or_TF=='gene'){

          subdirectory <- 'top_genes/'

        } else if(internal_gene_or_TF=='TF'){

          subdirectory <- 'top_TFs/'
        }

        ## Create a title for the scatter plot pdf:
        ## This is a comination of the probe name with the linked gene:
        scatterplot_pdf_title <- paste(
          paste(
            TENET_directory,
            'step7/',
            'hypo_Gplus_simple_scatterplots/',
            subdirectory,
            sep=''
          ),
          gene_name,
          '_',
          CpG_name_placeholder,
          '_scatterplot.pdf',
          sep=''
        )

        ## Open a pdf for saving the plot:
        pdf(
          scatterplot_pdf_title,
          height= 7,
          width= 10
        )

        plot(scatter_plot_updated)

        ## Close the plot:
        dev.off()
      }

      lapply(
        X= unlisted_CpGs,
        FUN= internal_scatterplot_function,
        internal_gene_or_TF= gene_or_TF
      )
    }

    ## Generate the plots for all genes/TFs of interest:
    suppressWarnings(
      parallel::mcmapply(
        FUN= hypo_Gplus_scatterplot_function,
        CpGs_list= top_hypo_Gplus_all_genes_linked_cpgs_list,
        gene_ENSG= names(top_hypo_Gplus_all_genes_linked_cpgs_list),
        gene_or_TF= 'gene',
        mc.cores= core_count
      )
    )

    suppressWarnings(
      parallel::mcmapply(
        FUN= hypo_Gplus_scatterplot_function,
        CpGs_list= top_hypo_Gplus_all_TF_linked_cpgs_list,
        gene_ENSG= names(top_hypo_Gplus_all_TF_linked_cpgs_list),
        gene_or_TF= 'TF',
        mc.cores= core_count
      )
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
          'hypo_Gminus_simple_scatterplots',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_simple_scatterplots',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hypo_Gminus_simple_scatterplots
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_simple_scatterplots/',
          'top_genes',
          sep=''
        )
      )
    }

    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_simple_scatterplots/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_simple_scatterplots/',
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
      stop('hypo_Gminus_sig_link_zscores_perm_optimized.txt in step5 of TENET directory was not found. Please check that the file exists and consider rerunning the step5_optimize_links function')

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

    ## Index an empty list for all the top genes
    ## and the TF genes as well:
    top_hypo_Gminus_all_genes_linked_cpgs_list <- list()

    top_hypo_Gminus_all_TF_linked_cpgs_list <- list()

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hypo_Gminus_all_gene_ENSG)
    )){

      ## Get the genes' ENSG:
      gene_ENSG_placeholder <- top_hypo_Gminus_all_gene_ENSG[i]

      ## Get the probes linked to each gene:
      probes_linked_to_significant_gene <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_gene_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hypo_Gminus_all_genes_linked_cpgs_list[[i]] <- probes_linked_to_significant_gene
    }

    ## For each gene and TF of interest, get the list of CpGs associated with it:
    for(i in c(
      1:length(top_hypo_Gminus_all_TF_ENSG)
    )){

      ## Get the TFs ENSG:
      TF_ENSG_placeholder <- top_hypo_Gminus_all_TF_ENSG[i]

      ## Get the probes linked to each TF
      probes_linked_to_significant_TF <- unique(
        hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_TF_ENSG[i],
          'probeID'
        ]
      )

      ## Add the probes to the list:
      top_hypo_Gminus_all_TF_linked_cpgs_list[[i]] <- probes_linked_to_significant_TF
    }

    ## Add the names of the genes to the lists:
    names(top_hypo_Gminus_all_genes_linked_cpgs_list) <- top_hypo_Gminus_all_gene_ENSG

    names(top_hypo_Gminus_all_TF_linked_cpgs_list) <- top_hypo_Gminus_all_TF_ENSG

    ## Write a function that, when given a list of CpGs linked to a TF,
    ## will plot scatterplots showing the methylation of the probe on the X
    ## and expression of that TF on the Y across all the tumor and normal smaples:
    hypo_Gminus_scatterplot_function <- function(
      CpGs_list,
      gene_ENSG,
      gene_or_TF
    ){

      ## Unlist the CpGs linked to each probe:
      unlisted_CpGs <- c(
        unlist(
          CpGs_list
        )
      )

      ## Convert the gene ENSG into the gene name:
      gene_name <- gencode_v22_genes[
        gene_ENSG,
        'gene_name'
      ]

      ## Add the expression of the gene of interest to DichF:
      TF_expression <- c(
        unlist(
          expDataF_subC[
            gene_ENSG,
            DichF$group
          ]
        )
      )

      ## Now write an internal function that will get each linked CpGs methylation,
      ## and use the methylation and gene expression to create a ggplot2 scatterplot
      ## and save it:
      internal_scatterplot_function <- function(
        CpGs_linked_to_TF,
        internal_gene_or_TF
      ){

        ## Save the CpG name:
        CpG_name_placeholder <- CpGs_linked_to_TF

        ## Get DNA methylation values:
        unlisted_CpG_methylation <- c(
          unlist(
            metDataF_subC[
              CpGs_linked_to_TF,
              DichF$group
            ]
          )
        )

        ## Manually coloring samples:
        t_v_n_group_colors <- c(
          'Normal'='dodgerblue3',
          'Tumor'='red3'
        )

        ## Creating scatter with ggplot2:
        scatter_plot <- ggplot2::qplot(
          x=TF_expression,
          y=unlisted_CpG_methylation,
          geom=c("point"),
          colour=DichF$cluster
        )

        ## Create the plot:
        scatter_plot_updated <- scatter_plot +
          ggplot2::ggtitle(
            paste(
              gene_name,
              ' gene expression vs.\n',
              CpG_name_placeholder,
              ' DNA methylation',
              sep=''
            )
          ) +
          ggplot2::ylab(
            paste(
              CpG_name_placeholder,
              'DNA methylation',
              sep=' '
            )
          ) +
          ggplot2::xlab(
            paste(
              gene_name,
              'gene expression',
              sep=' '
            )
          ) +
          ggplot2::theme_bw() +
          ggplot2::scale_color_manual(
            values=t_v_n_group_colors,
            name="Sample type"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust=0.5, size=18),
            legend.title = ggplot2::element_text(hjust=0.5, size=12),
            legend.text = ggplot2::element_text(size=10, colour='black'),
            panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
            axis.title.x = ggplot2::element_text(size=16),
            axis.title.y = ggplot2::element_text(size=16),
            axis.text.x = ggplot2::element_text(size=14, colour='black'),
            axis.text.y = ggplot2::element_text(size=14, colour='black'),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
          )

        ## Determine where to save the results:
        if(internal_gene_or_TF=='gene'){

          subdirectory <- 'top_genes/'

        } else if(internal_gene_or_TF=='TF'){

          subdirectory <- 'top_TFs/'
        }

        ## Create a title for the scatter plot pdf:
        ## This is a comination of the probe name with the linked gene:
        scatterplot_pdf_title <- paste(
          paste(
            TENET_directory,
            'step7/',
            'hypo_Gminus_simple_scatterplots/',
            subdirectory,
            sep=''
          ),
          gene_name,
          '_',
          CpG_name_placeholder,
          '_scatterplot.pdf',
          sep=''
        )

        ## Open a pdf for saving the plot:
        pdf(
          scatterplot_pdf_title,
          height= 7,
          width= 10
        )

        plot(scatter_plot_updated)

        ## Close the plot:
        dev.off()
      }

      lapply(
        X= unlisted_CpGs,
        FUN= internal_scatterplot_function,
        internal_gene_or_TF= gene_or_TF
      )
    }

    ## Generate the plots for all genes/TFs of interest:
    suppressWarnings(
      parallel::mcmapply(
        FUN= hypo_Gminus_scatterplot_function,
        CpGs_list= top_hypo_Gminus_all_genes_linked_cpgs_list,
        gene_ENSG= names(top_hypo_Gminus_all_genes_linked_cpgs_list),
        gene_or_TF= 'gene',
        mc.cores= core_count
      )
    )

    suppressWarnings(
      parallel::mcmapply(
        FUN= hypo_Gminus_scatterplot_function,
        CpGs_list= top_hypo_Gminus_all_TF_linked_cpgs_list,
        gene_ENSG= names(top_hypo_Gminus_all_TF_linked_cpgs_list),
        gene_or_TF= 'TF',
        mc.cores= core_count
      )
    )
  }
}

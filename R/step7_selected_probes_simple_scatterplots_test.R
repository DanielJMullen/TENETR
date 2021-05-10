#' step7_selected_probes_simple_scatterplots
#'
#' This is a step7 function of the TENETR package.
#' This function takes a custom list of DNA methylation probes specified by the user
#' and generates scatterplots displaying the expression level of the genes
#' linked to these probes for each of the four hypo or hypermethylated
#' Gplus or Gminus analysis quadrants, as selected by the user.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with an additional a subdirectory 'selected_probes_simple_scatterplots' containing the results of this function.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create scatterplots for the user-supplied probes with all genes the probes may have hypermeth Gplus links to.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to create scatterplots for the user-supplied probes with all genes the probes may have hypermeth Gminus links to.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to create scatterplots for the user-supplied probes with all genes the probes may have hypometh Gplus links to.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to create scatterplots for the user-supplied probes with all genes the probes may have hypometh Gminus links to.
#' @param probe_list Supply a vector of DNA methylation probe names for which scatterplots of those probes with expression of any linked genes of the specified analysis types above will be generated.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns .pdf files with scatterplots showing the expression of any genes linked to the specified probes on the x-axis and the methylation level of the probes on the y-axis.
#' @export

step7_selected_probes_simple_scatterplots <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  probe_list,
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
  ## selected probes simple scatterplots:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step7/',
        'selected_probes_simple_scatterplots',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step7/',
        'selected_probes_simple_scatterplots',
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
      colnames(metDataT),
      colnames(metDataN)
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

  ## For each of the analysis types selected, check to see if the
  ## sig_link_zscores_perm_optimized.txt from step 5 exists:

  ## Check results for hypermeth Gplus analysis:
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
  }

  ## Check results for hypermeth Gminus analysis:
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
  }

  ## Check results for hypometh Gplus analysis:
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
  }

  ## Check results for hypometh Gminus analysis:
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
  }

  ## For each probe, write a function to create a subdirectory for it,
  ## identify which probes it might be linked to, the analysis type of those links,
  ## and generate a scatterplot for each link (and note if no links are found):
  probe_scatterplot_function <- function(
    probe_of_interest
  ){

    ## If a subdirectory hasn't been created for the probe, create it:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'selected_probes_simple_scatterplots/',
          probe_of_interest,
          '_scatterplots',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'selected_probes_simple_scatterplots/',
          probe_of_interest,
          '_scatterplots',
          sep=''
        )
      )
    }

    ## Get the methylation of the given probe:
    ## if the probe is in the analysis.
    ## If it isn't return an error message:
    if(
      probe_of_interest %in% rownames(metDataF_subC)
    ){

      probe_methylation <- c(
        unlist(
          metDataF_subC[
            probe_of_interest,
            DichF$group
          ]
        )
      )

      ## Index empty vectors for the genes the probe is linked to
      ## as well as the analysis type the each gene came from:
      linked_genes <- c()
      analysis_type <- c()

      ## For each analysis type, identify genes linked to the probe:
      ## and indicate the analysis type it came from:
      if(
        hypermeth_Gplus_analysis==TRUE
      ){

        hyper_Gplus_genes <- hyper_Gplus_sig_link_zscores[
          hyper_Gplus_sig_link_zscores$probeID==probe_of_interest,
          'geneID'
        ]

        linked_genes <- c(
          linked_genes,
          hyper_Gplus_genes
        )

        analysis_type <- c(
          analysis_type,
          rep(
            'hyper_Gplus',
            times= length(hyper_Gplus_genes)
          )
        )
      }

      if(
        hypermeth_Gminus_analysis==TRUE
      ){

        hyper_Gminus_genes <- hyper_Gminus_sig_link_zscores[
          hyper_Gminus_sig_link_zscores$probeID==probe_of_interest,
          'geneID'
        ]

        linked_genes <- c(
          linked_genes,
          hyper_Gminus_genes
        )

        analysis_type <- c(
          analysis_type,
          rep(
            'hyper_Gminus',
            times= length(hyper_Gminus_genes)
          )
        )
      }

      if(
        hypometh_Gplus_analysis==TRUE
      ){

        hypo_Gplus_genes <- hypo_Gplus_sig_link_zscores[
          hypo_Gplus_sig_link_zscores$probeID==probe_of_interest,
          'geneID'
        ]

        linked_genes <- c(
          linked_genes,
          hypo_Gplus_genes
        )

        analysis_type <- c(
          analysis_type,
          rep(
            'hypo_Gplus',
            times= length(hypo_Gplus_genes)
          )
        )
      }

      if(
        hypometh_Gminus_analysis==TRUE
      ){

        hypo_Gminus_genes <- hypo_Gminus_sig_link_zscores[
          hypo_Gminus_sig_link_zscores$probeID==probe_of_interest,
          'geneID'
        ]

        linked_genes <- c(
          linked_genes,
          hypo_Gminus_genes
        )

        analysis_type <- c(
          analysis_type,
          rep(
            'hypo_Gminus',
            times= length(hypo_Gminus_genes)
          )
        )
      }

      ## Now let's write an internal function to generate scatterplots
      ## for each of the probe-gene links:
      internal_scatterplot_function <- function(
        gene_ENSG_linked_to_CpG,
        analysis_type
      ){

        ## Get expression of the gene of interest
        ## in the order of DichF:
        ## Get the expression of the gene:
        gene_expression <- c(
          unlist(
            expDataF_subC[
              gene_ENSG_linked_to_CpG,
              DichF$group
            ]
          )
        )

        ## Get the gene's name:
        gene_name <- gencode_v22_genes[
          gene_ENSG_linked_to_CpG,
          'gene_name'
        ]

        ## Manually coloring samples:
        t_v_n_group_colors <- c(
          'Normal'='dodgerblue3',
          'Tumor'='red3'
        )

        ## Creating scatter with ggplot2:
        scatter_plot <- ggplot2::qplot(
          x=gene_expression,
          y=probe_methylation,
          geom=c("point"),
          colour=DichF$cluster
        )

        ## Create the plot:
        scatter_plot_updated <- scatter_plot +
          ggplot2::ggtitle(
            paste(
              gene_name,
              ' gene expression vs.\n',
              probe_of_interest,
              ' DNA methylation',
              sep=''
            )
          ) +
          ggplot2::ylab(
            paste(
              probe_of_interest,
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

        ## Create a title for the scatter plot pdf:
        ## This is a comination of the probe name with the linked gene:
        scatterplot_pdf_title <- paste(
          paste(
            TENET_directory,
            'step7/',
            'selected_probes_simple_scatterplots/',
            probe_of_interest,
            '_scatterplots/',
            sep=''
          ),
          analysis_type,
          '_',
          gene_name,
          '_',
          probe_of_interest,
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

      ## If genes were found to be linked to the probe, create the scatterplots
      ## otherwise, generate a file noting no probes were found:
      if(
        length(linked_genes)>0
      ){

        mapply(
          FUN= internal_scatterplot_function,
          gene_ENSG_linked_to_CpG= linked_genes,
          analysis_type= analysis_type
        )

      } else{

        fileConn <- file(
          paste(
            TENET_directory,
            'step7/',
            'selected_probes_simple_scatterplots/',
            probe_of_interest,
            '_scatterplots/',
            'no_results_found.txt',
            sep=''
          )
        )

        writeLines(
          c(
            "No genes were found to be linked to the probe in the analysis types specified",
            "Check to see if additional the analysis types should be specified and set to TRUE"
          ),
          fileConn
        )

        close(fileConn)

      }

    } else{

      fileConn <- file(
        paste(
          TENET_directory,
          'step7/',
          'selected_probes_simple_scatterplots/',
          probe_of_interest,
          '_scatterplots/',
          'probe_not_found.txt',
          sep=''
        )
      )

      writeLines(
        c(
          "The probe of interest was not found in the supplied methylation/expression dataset.",
          "Please check the probe name or the dataset if probe should be present."
        ),
        fileConn
      )

      close(fileConn)

    }
  }

  ## Apply the function to each probe:
  parallel::mclapply(
    X= probe_list,
    FUN= probe_scatterplot_function,
    mc.cores= core_count
  )

}

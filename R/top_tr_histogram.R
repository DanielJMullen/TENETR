#' top_tr_histogram
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top TRs by number of linked probes identified from
#' the step6 top_tr_tabulation function and generates heatmaps showing the
#' number of all genes and TF-only genes with a given number of enhancer
#' DNA methylation probes of a given type linked to them.
#'
#'
#' @param TENET_directory Set a path to the directory that contains step6 results from the top_tr_tabulation function. This function will also create a new step7 folder there if it has not been created, with a subdirectory called histogram containing the results.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create histograms for the top TRs by most hypermeth probes with G+ links.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create histograms for the top TRs by most hypermeth probes with G- links.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create histograms for the top TRs by most hypometh probes with G+ links.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create histograms for the top TRs by most hypometh probes with G- links.
#' @return Currently returns .pdf files with the histograms showing the number of TRs linked to a given number of probes of the given analysis type.
#' @export

top_tr_histogram <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis
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

  ## Generate results for hypermeth Gplus probes:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gplus histograms:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_histograms',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_histograms',
          sep=''
        )
      )

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

    ## Check that thehyper_Gplus_links_all_TR_freq.txt file exists:
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

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of genes with a given number of
    ## hyper.G+ links linked to them:
    hyper_Gplus_all_genes_histogram <- ggplot2::ggplot(
      hyper_Gplus_all_gene_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hyper_Gplus_all_gene_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of Hyper.G+ linked enhancer probes per gene") +
      ggplot2::xlab("Number of linked enhancer probes per gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_histograms/',
        'hyper_Gplus_links_all_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hyper_Gplus_all_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of TF-only genes with a given number of
    ## hyper.G+ links linked to them:
    hyper_Gplus_TF_genes_histogram <- ggplot2::ggplot(
      hyper_Gplus_all_TF_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hyper_Gplus_all_TF_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of Hyper.G+ linked enhancer probes per TF gene") +
      ggplot2::xlab("Number of linked enhancer probes per TF gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_histograms/',
        'hyper_Gplus_links_TF_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hyper_Gplus_TF_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )
  }

  ## Generate results for hypermeth Gminus probes:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypermeth Gminus histograms:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_histograms',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_histograms',
          sep=''
        )
      )

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

    ## Check that the hyper_Gminus_links_all_gene_freq.txt file exists:
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

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of genes with a given number of
    ## hyper.G- links linked to them:
    hyper_Gminus_all_genes_histogram <- ggplot2::ggplot(
      hyper_Gminus_all_gene_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hyper_Gminus_all_gene_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of Hyper.G- linked enhancer probes per gene") +
      ggplot2::xlab("Number of linked enhancer probes per gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_histograms/',
        'hyper_Gminus_links_all_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hyper_Gminus_all_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of TF-only genes with a given number of
    ## hyper.G- links linked to them:
    hyper_Gminus_TF_genes_histogram <- ggplot2::ggplot(
      hyper_Gminus_all_TF_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hyper_Gminus_all_TF_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of Hyper.G- linked enhancer probes per TF gene") +
      ggplot2::xlab("Number of linked enhancer probes per TF gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_histograms/',
        'hyper_Gminus_links_TF_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hyper_Gminus_TF_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )
  }

  ## Generate results for hypometh Gplus probes:
  if(hypometh_Gplus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gplus histograms:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_histograms',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_histograms',
          sep=''
        )
      )

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

    ## Check that the hypo_Gminus_links_all_gene_freq.txt file exists:
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

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of genes with a given number of
    ## hypo.G+ links linked to them:
    hypo_Gplus_all_genes_histogram <- ggplot2::ggplot(
      hypo_Gplus_all_gene_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hypo_Gplus_all_gene_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of hypo.G+ linked enhancer probes per gene") +
      ggplot2::xlab("Number of linked enhancer probes per gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_histograms/',
        'hypo_Gplus_links_all_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hypo_Gplus_all_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of TF-only genes with a given number of
    ## hypo.G+ links linked to them:
    hypo_Gplus_TF_genes_histogram <- ggplot2::ggplot(
      hypo_Gplus_all_TF_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hypo_Gplus_all_TF_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of hypo.G+ linked enhancer probes per TF gene") +
      ggplot2::xlab("Number of linked enhancer probes per TF gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_histograms/',
        'hypo_Gplus_links_TF_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hypo_Gplus_TF_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )
  }

  ## Generate results for hypometh Gminus probes:
  if(hypometh_Gminus_analysis==TRUE){

    ## Create a subdirectory in the new step7 directory to contain the
    ## hypometh Gminus histograms:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_histograms',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_histograms',
          sep=''
        )
      )

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

    ## Check that the hypo_Gminus_links_all_gene_freq.txt file exists:
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

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of genes with a given number of
    ## hypo.G- links linked to them:
    hypo_Gminus_all_genes_histogram <- ggplot2::ggplot(
      hypo_Gminus_all_gene_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hypo_Gminus_all_gene_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of hypo.G- linked enhancer probes per gene") +
      ggplot2::xlab("Number of linked enhancer probes per gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_histograms/',
        'hypo_Gminus_links_all_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hypo_Gminus_all_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )

    ## Now that the files have been loaded, create a histogram
    ## that shows the number of TF-only genes with a given number of
    ## hypo.G- links linked to them:
    hypo_Gminus_TF_genes_histogram <- ggplot2::ggplot(
      hypo_Gminus_all_TF_freq,
      ggplot2::aes(x=Freq)
    ) +
      ggplot2::geom_histogram(
        color='black',
        fill='darkgrey',
        binwidth=ceiling(
          max(hypo_Gminus_all_TF_freq$Freq)/200
        )
      ) +
      ggplot2::ggtitle("Histogram of hypo.G- linked enhancer probes per TF gene") +
      ggplot2::xlab("Number of linked enhancer probes per TF gene") +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust=0.5, size=20),
        panel.border = ggplot2::element_rect(colour = 'black', fill=NA, size=1),
        axis.title.x = ggplot2::element_text(size=20),
        axis.title.y = ggplot2::element_text(size=20),
        axis.text.x = ggplot2::element_text(size=18, colour = 'black'),
        axis.text.y = ggplot2::element_text(size=16, colour = 'black'),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    ## Save the plot:
    ggplot2::ggsave(
      filename= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_histograms/',
        'hypo_Gminus_links_TF_gene_freq_histogram.pdf',
        sep=''
      ),
      plot= hypo_Gminus_TF_genes_histogram,
      width = 8,
      height = 8,
      units = c("in")
    )
  }
}

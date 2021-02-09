#' optimize_links
#'
#' This is the step5 function of the TENETR package.
#' This function takes the calculated permutated p-values for the hyper/hypomethylated
#' Gplus or Gminus probe-gene quadrants, and selects optimized links based on relative
#' expression of the gene in hyper/hypomethylated experimental/tumor samples and the control/normal samples
#' based on their means and in a wilcoxon test, checks that the hyper/hypometh samples for that given
#' gene-probe link also show appropriately higher/lower expression of the linked gene in a number of
#' experimental/tumor samples greater than the minExp number specified in step2 get_diffmeth_regions function
#' and have maximum/minimum methylation above/below the hyper_stringency and hypo_stringency cutoffs.
#'
#'
#' @param TENET_directory Set a path to the directory that contains step4 results from permutate_z_scores function. This function will also create a new step5 folder there containing the results.
#' @param hypermeth_Gplus_analysis Set TRUE or FALSE if user wants to optimize hypermeth_Gplus links. Requires hypermeth_analysis from step4 to have been set to TRUE.
#' @param hypermeth_Gminus_analysis Set TRUE or FALSE if user wants to optimize hypermeth_Gminus links. Requires hypermeth_analysis from step4 to have been set to TRUE.
#' @param hypometh_Gplus_analysis Set TRUE or FALSE if user wants to optimize hypometh_Gplus links. Requires hypometh_analysis from step4 to have been set to TRUE.
#' @param hypometh_Gminus_analysis Set TRUE or FALSE if user wants to optimize hypometh_Gminus links. Requires hypometh_analysis from step4 to have been set to TRUE.
#' @param adj_pval_cutoff Set p-value for BH-corrected Wilcoxon p-values in comparison of gene expression values between hyper/hypomethylated tumor/experimental samples and the control/normal samples..
#' @param hyper_stringency Set a number from 0 to 1 to be the beta-value cutoff to optimize for hypermeth links with maximum methylation values above the cutoff
#' @param hypo_stringency Set a number from 0 to 1 to be the beta-value cutoff to optimize for hypometh links with minimum methylation values below the cutoff
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns tab-delimited "sig_link_zscores_perm_optimized.txt" files for hypo/hyper Gplus/Gminus probe-gene links, similar to step4, but only with the optimized links.
#' @export

optimize_links <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  adj_pval_cutoff,
  hyper_stringency,
  hypo_stringency,
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

  ## Create a step5 directory to deposit the output paired score files:
  dir.create(
    paste(
      TENET_directory,
      'step5/',
      sep=''
    )
  )

  ## Check that the diff methylated dataset from step2 exists and load it:
  ## If not, return an error message:
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

    load(
      paste(
        TENET_directory,
        'step2/',
        'diff.methylated.datasets.rda',
        sep=''
      )
    )

  } else{

    ## Return the error if the file wasn't found:
    stop('diff.methylated.datasets.rda in step2 of TENET directory was not found. Please check that the file exists and consider rerunning the step2 get_diffmeth_regions function.')
  }

  ## Create expression datasets noting for each gene, which samples
  ## have 0 expression for it as 0s, and those with expression as 1s:
  expDataT0=ifelse(
    expDataT==0,
    0,
    1
  )
  expDataN0=ifelse(
    expDataN==0,
    0,
    1
  )

  ## Write functions to get the mean expression of genes
  ## in the normal and tumor samples:
  getmeanexpN <- function(geneID){
    mean(
      expDataN[
        as.character(geneID),
      ],
      na.rm=T
    )
  }

  getmeanexpT <- function(geneID){
    mean(
      expDataT[
        as.character(geneID),
      ],
      na.rm=T
    )
  }

  ## If hypermeth_Gplus_analysis is selected, do the optimization:
  if(hypermeth_Gplus_analysis==TRUE){

    ## Check that the hypo_Gminus_sig_link_zscores_perm.txt file exists;
    if(
      file.exists(
        paste(
          TENET_directory,
          'step4/',
          'hyper_Gplus_sig_link_zscores_perm.txt',
          sep=''
        )
      )
    ){

      hyper_Gplus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step4/',
          'hyper_Gplus_sig_link_zscores_perm.txt',
          sep=''
        ),
        header= TRUE,
        sep='\t',
        stringsAsFactors = FALSE
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hyper_Gplus_sig_link_zscores_perm.txt in step4 of TENET directory was not found. Please check that the file exists and consider rerunning the step4 permutate_z_scores function.')
    }

    ## Index empty columns with NA values for now:
    hyper_Gplus_zscores$mean.expN=paste(NA)
    hyper_Gplus_zscores$mean.expT=paste(NA)
    hyper_Gplus_zscores$wilcox.expNcTc=paste(NA)
    hyper_Gplus_zscores$hypermeth.tumor.length=paste(NA)
    hyper_Gplus_zscores$meanhypermethT=paste(NA)
    hyper_Gplus_zscores$hypermeth.lower.tumor.length=paste(NA)
    hyper_Gplus_zscores$max.metTc=paste(NA)
    hyper_Gplus_zscores$mean.expN.high.expT=paste(NA)

    ## For each probe, note which of the samples are above the hypermeth cutoff
    ## with 1s:
    metDataTcat <- ifelse(
      hypermethDataT>hypermethcutoff,
      1,
      0
    )

    ## Now write a function to calculate a Wilcoxon p-value for expression
    ## of a given gene for a probe between the normal and tumor samples:
    getWilcoxexpNcTc <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the normal samples:
      expN <- expDataN[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed normal dataset:
      expN1 <- expDataN0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expN1_names <- names(
        expN1[
          expN1==1
        ]
      )

      ## Get expression values for only the expressed normal samples:
      expNc <- expN[
        expN1_names
      ]

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## If there are samples expressing it in both datasets, return the
      ## Wilcoxon p-value comparing their differential expression.
      ## Otherwise return NA:
      if(length(expNc)>0 & length(expTc)>0){

        wilcoxPval <- suppressWarnings(
          wilcox.test(
            expNc,
            expTc,
          )$p.value
        )

      } else{

        wilcoxPval <- NA

      }

      ## Return the final wilcoxPval:
      return(wilcoxPval)
    }

    ## Now write a function to calculate the number
    ## of hypermethylated samples for that probe:
    gethypermethTlength <- function(probe){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Return the length of these samples:
      return(
        length(
          metT1_names
        )
      )
    }

    ## Get the mean expression of the gene of interest in the
    ## hypermethylated tumor samples:
    getmeanhypermethT <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean expression in these samples:
      return(
        mean(
          expTc,
          na.rm=T
        )
      )
    }

    ## Write a function calculating the number
    ## of hypermeth samples less than the tumor expression mean:
    gethypermethlowerTlength <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean of these samples:
      expT_mean <- mean(
        expT,
        na.rm=T
      )

      ## Calculate the number of hypermeth samples less the mean:
      less_than_mean_samples <- sum(
        expTc < expT_mean,
        na.rm=T
      )

      ## Return that value:
      return(less_than_mean_samples)
    }

    getmaxmetTc <- function(
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the maximum methylation value for these hypermeth probes:
      max_val <-  max(
        hypermethDataT[
          as.character(probe),
          match(
            metT1_names,
            colnames(hypermethDataT)
          )
        ],
        na.rm=TRUE
      )

      return(max_val)
    }

    ## Execute these functions to calculate key
    hyper_Gplus_zscores$mean.expN <- parallel::mcmapply(
      getmeanexpN,
      hyper_Gplus_zscores$geneID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$mean.expT <- parallel::mcmapply(
      getmeanexpT,
      hyper_Gplus_zscores$geneID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$wilcox.expNcTc <- parallel::mcmapply(
      getWilcoxexpNcTc,
      geneID= hyper_Gplus_zscores$geneID,
      probe= hyper_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$hypermeth.tumor.length <- parallel::mcmapply(
      gethypermethTlength,
      hyper_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$meanhypermethT <- parallel::mcmapply(
      getmeanhypermethT,
      geneID= hyper_Gplus_zscores$geneID,
      probe= hyper_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$hypermeth.lower.tumor.length <- parallel::mcmapply(
      gethypermethlowerTlength,
      geneID= hyper_Gplus_zscores$geneID,
      probe= hyper_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$max.metTc <- parallel::mcmapply(
      getmaxmetTc,
      hyper_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gplus_zscores$mean.expN.high.expT <- ifelse(
      as.numeric(
        hyper_Gplus_zscores$mean.expN
      ) >
      as.numeric(
        hyper_Gplus_zscores$mean.expT
      ),
      1,
      0
    )

    ## Now perform BH multiple testing correction on the wilcox p-values:
    hyper_Gplus_zscores$wilcox.expNcTc.adj.pval <- p.adjust(
      as.numeric(
        hyper_Gplus_zscores$wilcox.expNcTc
      ),
      "BH"
    )

    ## Now select the final links to carry over into post-hoc analysis:
    hyper_Gplus_links <- hyper_Gplus_zscores[
      which(
        as.numeric(hyper_Gplus_zscores$mean.expN)!=0 &
        as.numeric(hyper_Gplus_zscores$mean.expT)!=0 &
        as.numeric(hyper_Gplus_zscores$mean.expN.high.expT)==1
      ),
    ]

    hyper_Gplus_links <- hyper_Gplus_links[
      which(
        hyper_Gplus_links$wilcox.expNcTc.adj.pval<adj_pval_cutoff &
        as.numeric(
          hyper_Gplus_links$hypermeth.tumor.length
        )>minExp &
        as.numeric(
          hyper_Gplus_links$hypermeth.lower.tumor.length
        )>minExp &
        as.numeric(
          hyper_Gplus_links$max.metTc
        ) > hyper_stringency
      ),
    ]

    ## Write file to step5:
    write.table(
      hyper_Gplus_links,
      paste(
        TENET_directory,
        'step5/',
        'hyper_Gplus_sig_link_zscores_perm_optimized.txt',
        sep=''
      ),
      quote= FALSE,
      row.names = FALSE,
      sep='\t'
    )
  }

  ## If hypermeth_Gminus_analysis is selected, do the optimization:
  if(hypermeth_Gminus_analysis==TRUE){

    ## Check that the hypo_Gminus_sig_link_zscores_perm.txt
    if(
      file.exists(
        paste(
          TENET_directory,
          'step4/',
          'hyper_Gminus_sig_link_zscores_perm.txt',
          sep=''
        )
      )
    ){

      hyper_Gminus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step4/',
          'hyper_Gminus_sig_link_zscores_perm.txt',
          sep=''
        ),
        header= TRUE,
        sep='\t',
        stringsAsFactors = FALSE
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hyper_Gminus_sig_link_zscores_perm.txt in step4 of TENET directory was not found. Please check that the file exists and consider rerunning the step4 permutate_z_scores function.')
    }

    ## Index empty columns with NA values for now:
    hyper_Gminus_zscores$mean.expN=paste(NA)
    hyper_Gminus_zscores$mean.expT=paste(NA)
    hyper_Gminus_zscores$wilcox.expNcTc=paste(NA)
    hyper_Gminus_zscores$hypermeth.tumor.length=paste(NA)
    hyper_Gminus_zscores$meanhypomethT=paste(NA)
    hyper_Gminus_zscores$hypermeth.higher.tumor.length=paste(NA)
    hyper_Gminus_zscores$max.metTc=paste(NA)
    hyper_Gminus_zscores$mean.expN.low.expT=paste(NA)

    ## For each probe, note which of the samples are above the hypermeth cutoff
    ## with 1s:
    metDataTcat <- ifelse(
      hypermethDataT>hypermethcutoff,
      1,
      0
    )

    ## Now write a function to calculate a Wilcoxon p-value for expression
    ## of a given gene for a probe between the normal and tumor samples:
    getWilcoxexpNcTc <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the normal samples:
      expN <- expDataN[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed normal dataset:
      expN1 <- expDataN0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expN1_names <- names(
        expN1[
          expN1==1
        ]
      )

      ## Get expression values for only the expressed normal samples:
      expNc <- expN[
        expN1_names
      ]

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## If there are samples expressing it in both datasets, return the
      ## Wilcoxon p-value comparing their differential expression.
      ## Otherwise return NA:
      if(length(expNc)>0 & length(expTc)>0){

        wilcoxPval <- suppressWarnings(
          wilcox.test(
            expNc,
            expTc,
          )$p.value
        )

      } else{

        wilcoxPval <- NA

      }

      ## Return the final wilcoxPval:
      return(wilcoxPval)
    }

    ## Now write a function to calculate the number
    ## of hypermethylated samples for that probe:
    gethypermethTlength <- function(probe){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Return the length of these samples:
      return(
        length(
          metT1_names
        )
      )
    }

    ## Get the mean expression of the gene of interest in the
    ## hypermethylated tumor samples:
    getmeanhypermethT <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean expression in these samples:
      return(
        mean(
          expTc,
          na.rm=T
        )
      )
    }

    ## Write a function calculating the number
    ## of hypermeth samples less than the tumor expression mean:
    gethypermethhigherTlength <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypermethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean of these samples:
      expT_mean <- mean(
        expT,
        na.rm=T
      )

      ## Calculate the number of hypermeth samples less the mean:
      less_than_mean_samples <- sum(
        expTc > expT_mean,
        na.rm=T
      )

      ## Return that value:
      return(less_than_mean_samples)
    }

    getmaxmetTc <- function(
      probe
    ){

      ## Get the samples above/not above the hypermeth cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypermeth cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the maximum methylation value for these hypermeth probes:
      max_val <-  max(
        hypermethDataT[
          as.character(probe),
          match(
            metT1_names,
            colnames(hypermethDataT)
          )
        ],
        na.rm=TRUE
      )

      return(max_val)
    }

    ## Execute these functions to calculate key
    hyper_Gminus_zscores$mean.expN <- parallel::mcmapply(
      getmeanexpN,
      hyper_Gminus_zscores$geneID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$mean.expT <- parallel::mcmapply(
      getmeanexpT,
      hyper_Gminus_zscores$geneID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$wilcox.expNcTc <- parallel::mcmapply(
      getWilcoxexpNcTc,
      geneID= hyper_Gminus_zscores$geneID,
      probe= hyper_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$hypermeth.tumor.length <- parallel::mcmapply(
      gethypermethTlength,
      hyper_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$meanhypermethT <- parallel::mcmapply(
      getmeanhypermethT,
      geneID= hyper_Gminus_zscores$geneID,
      probe= hyper_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$hypermeth.higher.tumor.length <- parallel::mcmapply(
      gethypermethhigherTlength,
      geneID= hyper_Gminus_zscores$geneID,
      probe= hyper_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$max.metTc <- parallel::mcmapply(
      getmaxmetTc,
      hyper_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hyper_Gminus_zscores$mean.expN.low.expT <- ifelse(
      as.numeric(
        hyper_Gminus_zscores$mean.expN
      ) <
      as.numeric(
        hyper_Gminus_zscores$mean.expT
      ),
      1,
      0
    )

    ## Now perform BH multiple testing correction on the wilcox p-values:
    hyper_Gminus_zscores$wilcox.expNcTc.adj.pval <- p.adjust(
      as.numeric(
        hyper_Gminus_zscores$wilcox.expNcTc
      ),
      "BH"
    )

    ## Now select the final links to carry over into post-hoc analysis:
    hyper_Gminus_links <- hyper_Gminus_zscores[
      which(
        as.numeric(hyper_Gminus_zscores$mean.expN)!=0 &
          as.numeric(hyper_Gminus_zscores$mean.expT)!=0 &
          as.numeric(hyper_Gminus_zscores$mean.expN.low.expT)==1
      ),
    ]

    hyper_Gminus_links <- hyper_Gminus_links[
      which(
        hyper_Gminus_links$wilcox.expNcTc.adj.pval<adj_pval_cutoff &
          as.numeric(
            hyper_Gminus_links$hypermeth.tumor.length
          )>minExp &
          as.numeric(
            hyper_Gminus_links$hypermeth.higher.tumor.length
          )>minExp &
          as.numeric(
            hyper_Gminus_links$max.metTc
          ) > hyper_stringency
      ),
    ]

    ## Write file to step5:
    write.table(
      hyper_Gminus_links,
      paste(
        TENET_directory,
        'step5/',
        'hyper_Gminus_sig_link_zscores_perm_optimized.txt',
        sep=''
      ),
      quote= FALSE,
      row.names = FALSE,
      sep='\t'
    )
  }

  ## If hypometh_Gplus_analysis is selected, do the optimization:
  if(hypometh_Gplus_analysis==TRUE){

    ## Check that the hypo_Gminus_sig_link_zscores_perm.txt
    if(
      file.exists(
        paste(
          TENET_directory,
          'step4/',
          'hypo_Gplus_sig_link_zscores_perm.txt',
          sep=''
        )
      )
    ){

      hypo_Gplus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step4/',
          'hypo_Gplus_sig_link_zscores_perm.txt',
          sep=''
        ),
        header= TRUE,
        sep='\t',
        stringsAsFactors = FALSE
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hypo_Gplus_sig_link_zscores_perm.txt in step4 of TENET directory was not found. Please check that the file exists and consider rerunning the step4 permutate_z_scores function.')
    }

    ## Index empty columns with NA values for now:
    hypo_Gplus_zscores$mean.expN=paste(NA)
    hypo_Gplus_zscores$mean.expT=paste(NA)
    hypo_Gplus_zscores$wilcox.expNcTc=paste(NA)
    hypo_Gplus_zscores$hypometh.tumor.length=paste(NA)
    hypo_Gplus_zscores$meanhypomethT=paste(NA)
    hypo_Gplus_zscores$hypometh.higher.tumor.length=paste(NA)
    hypo_Gplus_zscores$min.metTc=paste(NA)
    hypo_Gplus_zscores$mean.expN.low.expT=paste(NA)

    ## For each probe, note which of the samples are above the hypometh cutoff
    ## with 1s:
    metDataTcat <- ifelse(
      hypomethDataT<hypomethcutoff,
      1,
      0
    )

    ## Now write a function to calculate a Wilcoxon p-value for expression
    ## of a given gene for a probe between the normal and tumor samples:
    getWilcoxexpNcTc <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the normal samples:
      expN <- expDataN[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed normal dataset:
      expN1 <- expDataN0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expN1_names <- names(
        expN1[
          expN1==1
        ]
      )

      ## Get expression values for only the expressed normal samples:
      expNc <- expN[
        expN1_names
      ]

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## If there are samples expressing it in both datasets, return the
      ## Wilcoxon p-value comparing their differential expression.
      ## Otherwise return NA:
      if(length(expNc)>0 & length(expTc)>0){

        wilcoxPval <- suppressWarnings(
          wilcox.test(
            expNc,
            expTc,
          )$p.value
        )

      } else{

        wilcoxPval <- NA

      }

      ## Return the final wilcoxPval:
      return(wilcoxPval)
    }

    ## Now write a function to calculate the number
    ## of hypomethylated samples for that probe:
    gethypomethTlength <- function(probe){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Return the length of these samples:
      return(
        length(
          metT1_names
        )
      )
    }

    ## Get the mean expression of the gene of interest in the
    ## hypomethylated tumor samples:
    getmeanhypomethT <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean expression in these samples:
      return(
        mean(
          expTc,
          na.rm=T
        )
      )
    }

    ## Write a function calculating the number
    ## of hypometh samples less than the tumor expression mean:
    gethypomethhigherTlength <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean of these samples:
      expT_mean <- mean(
        expT,
        na.rm=T
      )

      ## Calculate the number of hypometh samples less the mean:
      less_than_mean_samples <- sum(
        expTc > expT_mean,
        na.rm=T
      )

      ## Return that value:
      return(less_than_mean_samples)
    }

    getminmetTc <- function(
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the maximum methylation value for these hypometh probes:
      max_val <-  min(
        hypomethDataT[
          as.character(probe),
          match(
            metT1_names,
            colnames(hypomethDataT)
          )
        ],
        na.rm=TRUE
      )

      return(max_val)
    }

    ## Execute these functions to calculate key
    hypo_Gplus_zscores$mean.expN <- parallel::mcmapply(
      getmeanexpN,
      hypo_Gplus_zscores$geneID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$mean.expT <- parallel::mcmapply(
      getmeanexpT,
      hypo_Gplus_zscores$geneID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$wilcox.expNcTc <- parallel::mcmapply(
      getWilcoxexpNcTc,
      geneID= hypo_Gplus_zscores$geneID,
      probe= hypo_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$hypometh.tumor.length <- parallel::mcmapply(
      gethypomethTlength,
      hypo_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$meanhypomethT <- parallel::mcmapply(
      getmeanhypomethT,
      geneID= hypo_Gplus_zscores$geneID,
      probe= hypo_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$hypometh.higher.tumor.length <- parallel::mcmapply(
      gethypomethhigherTlength,
      geneID= hypo_Gplus_zscores$geneID,
      probe= hypo_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$min.metTc <- parallel::mcmapply(
      getminmetTc,
      hypo_Gplus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gplus_zscores$mean.expN.low.expT <- ifelse(
      as.numeric(
        hypo_Gplus_zscores$mean.expN
      ) <
        as.numeric(
          hypo_Gplus_zscores$mean.expT
        ),
      1,
      0
    )

    ## Now perform BH multiple testing correction on the wilcox p-values:
    hypo_Gplus_zscores$wilcox.expNcTc.adj.pval <- p.adjust(
      as.numeric(
        hypo_Gplus_zscores$wilcox.expNcTc
      ),
      "BH"
    )

    ## Now select the final links to carry over into post-hoc analysis:
    hypo_Gplus_links <- hypo_Gplus_zscores[
      which(
        as.numeric(hypo_Gplus_zscores$mean.expN)!=0 &
          as.numeric(hypo_Gplus_zscores$mean.expT)!=0 &
          as.numeric(hypo_Gplus_zscores$mean.expN.low.expT)==1
      ),
    ]

    hypo_Gplus_links <- hypo_Gplus_links[
      which(
        hypo_Gplus_links$wilcox.expNcTc.adj.pval<adj_pval_cutoff &
          as.numeric(
            hypo_Gplus_links$hypometh.tumor.length
          )>minExp &
          as.numeric(
            hypo_Gplus_links$hypometh.higher.tumor.length
          )>minExp &
          as.numeric(
            hypo_Gplus_links$min.metTc
          ) < hypo_stringency
      ),
    ]

    ## Write file to step5:
    write.table(
      hypo_Gplus_links,
      paste(
        TENET_directory,
        'step5/',
        'hypo_Gplus_sig_link_zscores_perm_optimized.txt',
        sep=''
      ),
      quote= FALSE,
      row.names = FALSE,
      sep='\t'
    )
  }

  ## If hypometh_Gminus_analysis is selected, do the optimization:
  if(hypometh_Gminus_analysis==TRUE){

    ## Check that the hypo_Gminus_sig_link_zscores_perm.txt
    if(
      file.exists(
        paste(
          TENET_directory,
          'step4/',
          'hypo_Gminus_sig_link_zscores_perm.txt',
          sep=''
        )
      )
    ){

      hypo_Gminus_zscores <- read.delim(
        paste(
          TENET_directory,
          'step4/',
          'hypo_Gminus_sig_link_zscores_perm.txt',
          sep=''
        ),
        header= TRUE,
        sep='\t',
        stringsAsFactors = FALSE
      )

    } else{

      ## Return the error if the file wasn't found:
      stop('hypo_Gminus_sig_link_zscores_perm.txt in step4 of TENET directory was not found. Please check that the file exists and consider rerunning the step4 permutate_z_scores function.')
    }

    ## Index empty columns with NA values for now:
    hypo_Gminus_zscores$mean.expN=paste(NA)
    hypo_Gminus_zscores$mean.expT=paste(NA)
    hypo_Gminus_zscores$wilcox.expNcTc=paste(NA)
    hypo_Gminus_zscores$hypometh.tumor.length=paste(NA)
    hypo_Gminus_zscores$meanhypomethT=paste(NA)
    hypo_Gminus_zscores$hypometh.lower.tumor.length=paste(NA)
    hypo_Gminus_zscores$min.metTc=paste(NA)
    hypo_Gminus_zscores$mean.expN.high.expT=paste(NA)

    ## For each probe, note which of the samples are above the hypometh cutoff
    ## with 1s:
    metDataTcat <- ifelse(
      hypomethDataT<hypomethcutoff,
      1,
      0
    )

    ## Now write a function to calculate a Wilcoxon p-value for expression
    ## of a given gene for a probe between the normal and tumor samples:
    getWilcoxexpNcTc <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the normal samples:
      expN <- expDataN[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed normal dataset:
      expN1 <- expDataN0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expN1_names <- names(
        expN1[
          expN1==1
        ]
      )

      ## Get expression values for only the expressed normal samples:
      expNc <- expN[
        expN1_names
      ]

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## If there are samples expressing it in both datasets, return the
      ## Wilcoxon p-value comparing their differential expression.
      ## Otherwise return NA:
      if(length(expNc)>0 & length(expTc)>0){

        wilcoxPval <- suppressWarnings(
          wilcox.test(
            expNc,
            expTc,
          )$p.value
        )

      } else{

        wilcoxPval <- NA

      }

      ## Return the final wilcoxPval:
      return(wilcoxPval)
    }

    ## Now write a function to calculate the number
    ## of hypomethylated samples for that probe:
    gethypomethTlength <- function(probe){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Return the length of these samples:
      return(
        length(
          metT1_names
        )
      )
    }

    ## Get the mean expression of the gene of interest in the
    ## hypomethylated tumor samples:
    getmeanhypomethT <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean expression in these samples:
      return(
        mean(
          expTc,
          na.rm=T
        )
      )
    }

    ## Write a function calculating the number
    ## of hypometh samples less than the tumor expression mean:
    gethypomethlowerTlength <- function(
      geneID,
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the expression for the gene in the tumor samples:
      expT <- expDataT[
        as.character(geneID),
      ]

      ## Get the data for the gene from the expressed/no expressed tumor dataset:
      expT1 <- expDataT0[
        as.character(geneID),
      ]

      ## Get names the samples that have expression for the gene:
      expT1_names <- names(
        expT1[
          expT1==1
        ]
      )

      ## Now get samples that have expression of the gene
      ## AND that are hypomethylated:
      expT1_names_matched <- intersect(
        expT1_names,
        metT1_names
      )

      ## Get expression values for only the expressed tumor samples:
      expTc <- expT[
        expT1_names_matched
      ]

      ## Calculate the mean of these samples:
      expT_mean <- mean(
        expT,
        na.rm=T
      )

      ## Calculate the number of hypometh samples less the mean:
      less_than_mean_samples <- sum(
        expTc < expT_mean,
        na.rm=T
      )

      ## Return that value:
      return(less_than_mean_samples)
    }

    getminmetTc <- function(
      probe
    ){

      ## Get the samples above/not above the hypometh cutoff for the probe:
      metT1 <- metDataTcat[
        as.character(probe),
      ]

      ## Get the sample names that are above the hypometh cutoff for
      ## the probe:
      metT1_names <- names(
        metT1[
          metT1==1
        ]
      )

      ## Get the maximum methylation value for these hypometh probes:
      max_val <-  min(
        hypomethDataT[
          as.character(probe),
          match(
            metT1_names,
            colnames(hypomethDataT)
          )
        ],
        na.rm=TRUE
      )

      return(max_val)
    }

    ## Execute these functions to calculate key
    hypo_Gminus_zscores$mean.expN <- parallel::mcmapply(
      getmeanexpN,
      hypo_Gminus_zscores$geneID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$mean.expT <- parallel::mcmapply(
      getmeanexpT,
      hypo_Gminus_zscores$geneID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$wilcox.expNcTc <- parallel::mcmapply(
      getWilcoxexpNcTc,
      geneID= hypo_Gminus_zscores$geneID,
      probe= hypo_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$hypometh.tumor.length <- parallel::mcmapply(
      gethypomethTlength,
      hypo_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$meanhypomethT <- parallel::mcmapply(
      getmeanhypomethT,
      geneID= hypo_Gminus_zscores$geneID,
      probe= hypo_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$hypometh.lower.tumor.length <- parallel::mcmapply(
      gethypomethlowerTlength,
      geneID= hypo_Gminus_zscores$geneID,
      probe= hypo_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$min.metTc <- parallel::mcmapply(
      getminmetTc,
      hypo_Gminus_zscores$probeID,
      mc.cores= core_count
    )

    hypo_Gminus_zscores$mean.expN.high.expT <- ifelse(
      as.numeric(
        hypo_Gminus_zscores$mean.expN
      ) >
        as.numeric(
          hypo_Gminus_zscores$mean.expT
        ),
      1,
      0
    )

    ## Now perform BH multiple testing correction on the wilcox p-values:
    hypo_Gminus_zscores$wilcox.expNcTc.adj.pval <- p.adjust(
      as.numeric(
        hypo_Gminus_zscores$wilcox.expNcTc
      ),
      "BH"
    )

    ## Now select the final links to carry over into post-hoc analysis:
    hypo_Gminus_links <- hypo_Gminus_zscores[
      which(
        as.numeric(hypo_Gminus_zscores$mean.expN)!=0 &
          as.numeric(hypo_Gminus_zscores$mean.expT)!=0 &
          as.numeric(hypo_Gminus_zscores$mean.expN.high.expT)==1
      ),
    ]

    hypo_Gminus_links <- hypo_Gminus_links[
      which(
        hypo_Gminus_links$wilcox.expNcTc.adj.pval<adj_pval_cutoff &
          as.numeric(
            hypo_Gminus_links$hypometh.tumor.length
          )>minExp &
          as.numeric(
            hypo_Gminus_links$hypometh.lower.tumor.length
          )>minExp &
          as.numeric(
            hypo_Gminus_links$min.metTc
          ) < hypo_stringency
      ),
    ]

    ## Write file to step5:
    write.table(
      hypo_Gminus_links,
      paste(
        TENET_directory,
        'step5/',
        'hypo_Gminus_sig_link_zscores_perm_optimized.txt',
        sep=''
      ),
      quote= FALSE,
      row.names = FALSE,
      sep='\t'
    )
  }
}

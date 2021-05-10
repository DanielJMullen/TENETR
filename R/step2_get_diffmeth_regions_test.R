#' step2_get_diffmeth_regions
#'
#' This is the step2 function of the TENETR package.
#' This function identifies DNA methylation probes that mark putative regulatory
#' elements, including enhancer and promoter regions. These are probes that lie within
#' regions with specific histone modifications and open chromatin regions,
#' from the step1_make_external_datasets function, and which are located
#' a user-specified distance relative to GENCODE v22 transcript transcription start sites.
#' After identifying DNA methylation probes representing the the specified regulatory elements,
#' the function classifies the probes as methylated, unmethylated, hypermethylated,
#' or hypomethylated based on their differential methylation between the control/normal
#' and experimental/tumor samples supplied by the user, defined by the cutoff values
#' also specified by the user.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step1' subdirectory and results created by the step1_make_external_datasets function. This function will also create a new 'step2' subdirectory there containing the results of this function.
#' @param methylation_expression_dataset Set a path to the an .rda file containing methylation and expression data, possibly created using the TCGA_downloader function. This dataset should contain 'expDataN' and 'expDataT' expression datasets, as well as 'metDataN' and 'metDataT' methylation datasets, and optionally a 'clinical' dataset with sample clinical information.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param assess_promoter Set TRUE or FALSE depending on if you want to identify DNA methylation probes that mark promoter regions or distal enhancer regions, respectively. Defaults to FALSE.
#' @param TSS_dist Set a positive integer to be the buffer in base pairs from GENCODE v22-annotated transcription start sites for DNA methylation probes to be considered promoter probes. Probes outside of the TSS_dist from annotated transcription start sites will be considered enhancer probes. Defaults to 1500.
#' @param use_ext_HM Set TRUE or FALSE depending on if you want to use DNA methylation probes in user-supplied histone modification datasets from the step1_make_external_datasets function. These probes will be combined with probes found in the consensus enhancer regions if 'use_consensus_ENH' is also set to TRUE. These will be overlapped with probes from nucleosome depleted regions if either 'use_ext_NDR' or 'use_consensus_NDR' are set to TRUE. Setting to TRUE requires that the user have run the step1_make_external_datasets function on user-supplied histone modification datasets previously. Defaults to FALSE.
#' @param use_ext_NDR Set TRUE or FALSE depending on if you want to use DNA methylation probes in user-supplied nucleosome depeleted regions datasets from the step1_make_external_datasets function. These probes will be combined with probes found in the consensus nucleosome depeleted regions if 'use_consensus_NDR' is also set to TRUE. These will be overlapped with probes from regions with histone modifications or enhancer regions if either 'use_ext_HM' or 'use_consensus_ENH' are set to TRUE. Setting to TRUE requires that the user have run the step1_make_external_datasets function on user-supplied nucleosome depleted regions datasets previously. Defaults to FALSE.
#' @param use_consensus_ENH Set TRUE or FALSE depending on if you want to use DNA methylation probes in consensus enhancer regions from the TENETR.data package and step1_make_external_datasets. These probes will be combined with probes found in user-supplied histone modification regions if 'use_ext_HM' is also set to TRUE. These will be overlapped with probes in nucleosome depleted regions if either 'use_ext_NDR' or 'use_consensus_NDR' are set to TRUE. Setting to TRUE requires that the user have run step1_make_external_datasets on consensus enhancer datasets from the TENETR.data package previously. Defaults to TRUE.
#' @param use_consensus_NDR Set TRUE or FALSE depending on if you want to use DNA methylation probes in consensus nucleosome depeleted regions from the TENETR.data package and step1_make_external_datasets. These probes will be combined with probes found in user-supplied nucleosome depeleted regions if 'use_ext_NDR' is also set to TRUE. These will be overlapped with probes from regions with histone modifications or enhancer regions if either 'use_ext_HM' or 'use_consensus_ENH' are set to TRUE. Setting to TRUE requires that the user have run step1_make_external_datasets on consensus nucleosome depleted regions datasets from the TENETR.data package previously. Defaults to TRUE.
#' @param methcutoff Set a number from 0 to 1 to be the beta-value cutoff for methylated probes.
#' @param hypomethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypomethylated probes. Should be set lower than the methcutoff.
#' @param unmethcutoff Set a number from 0 to 1 to be the beta-value cutoff for unmethylated probes.
#' @param hypermethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypermethylated probes. Should be set higher than the unmethcutoff.
#' @param min_experimental_count Set a positive integer to be the minimum number of experimental/tumor samples to be considered for the hypo/hypermethylated groups. Should be less than the total number of experimental/tumor groups.
#' @param purity_directory If this user has their own datasets containing data for possible cell types that may contaminate the data, as .rda files containing a data matrix with DNA methylation values, set this as the path to the directory containing the purity .rda files. Otherwise set to FALSE. Defaults to FALSE.
#' @return Returns two objects, a .rda file with matrices of methylation and expression data for control/normal and experimental/tumor sample, lists of the four types of identified DNA methylation probes marking the regulatory elements of interest, the clinical data if provided, and the user-set parameters for consistency in downstream analyses. Also outputs a .txt file containing metrics on the number of probes found in the different analysis categories.
#' @export

step2_get_diffmeth_regions <- function(
  TENET_directory,
  methylation_expression_dataset,
  DNA_methylation_manifest="HM450",
  assess_promoter=FALSE,
  TSS_dist=1500,
  use_ext_HM=FALSE,
  use_ext_NDR=FALSE,
  use_consensus_ENH=TRUE,
  use_consensus_NDR=TRUE,
  methcutoff,
  hypomethcutoff,
  unmethcutoff,
  hypermethcutoff,
  min_experimental_count,
  purity_directory= FALSE
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

  ## Create a step2 directory to deposit overlapped probe files if it doesn't
  ## already exist:
  if(
    !dir.exists(
      paste(
        TENET_directory,
        'step2/',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        TENET_directory,
        'step2/',
        sep=''
      )
    )
  }

  ## Now let's load the combined methylation/expression supplied by the user:

  ## First check that the file exists, and if it doesn't return an error message:
  ## Otherwise, load it:
  if(!file.exists(methylation_expression_dataset)){

    stop('The .rda dataset specified by methylation_expression_dataset was not found. Please check the path that was provided to the dataset and rerun.')

  } else{

    load(methylation_expression_dataset)
  }


  ## Next check to make sure the methylation and expression datasets
  ## are present:
  if(!exists('expDataN')){

    stop(
      'Control (Normal) condition expression data not loaded. \nPlease load Control (Normal) expression data in the .rda file as a matrix named expDataN.'
    )

  } else if(!exists('expDataT')){

    stop(
      'Experimental (Tumor) condition expression data not loaded. \nPlease load experimental (tumor) expression data in the .rda file as a matrix named expDataT.'
    )

  } else if(!exists('metDataN')){

    stop(
      'Control (Normal) condition methylation data not loaded. \nPlease load Control (Normal) methylation data in the .rda file as a matrix named metDataN.'
    )

  } else if(!exists('metDataT')){

    stop(
      'Experimental (Tumor) condition methylation data not loaded. \nPlease load experimental (tumor) methylation data in the .rda file as a matrix named metDataT.'
    )

  }

  ## Load the hg38 DNA methylation annotations:
  ## Written in by Zexun Wu
  if(DNA_methylation_manifest == "HM450"){

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

  ## Get the dataset of gencode v22 genes:
  gencode_v22_gtf <- TENETR.data::gencode_v22_annotations

  ## Get info for both genes and transcripts
  ## This will cause some near duplicate ranges but that is ok
  ## Since intersecting will be the same no matter how many identical
  ## ranges there are:
  gencode_v22_transcript <- gencode_v22_gtf[
    gencode_v22_gtf$type=='transcript',
  ]

  gencode_v22_genes <- gencode_v22_gtf[
    gencode_v22_gtf$type=='gene',
  ]

  ## Create the dataset of both genes and transcripts:
  gencode_v22_TSS <- rbind(gencode_v22_genes, gencode_v22_transcript)

  ## Remove uneeded datasets:
  rm(gencode_v22_transcript)
  rm(gencode_v22_genes)

  ## Determine whether the start or the end coordinate is the
  ## TSS based on the strand:
  gencode_v22_TSS$TSS <- ifelse(
    gencode_v22_TSS$strand=='+',
    gencode_v22_TSS$start,
    ifelse(
      gencode_v22_TSS$strand=='-',
      gencode_v22_TSS$end,
      ''
    )
  )

  ## Create a dataset with just the TSS and the spacing representing
  ## the TSS_dist value input by user:
  gencode_v22_TSS_df <- data.frame(
    'chr'= gencode_v22_TSS$seqnames,
    'start'= as.numeric(gencode_v22_TSS$TSS)-TSS_dist,
    'end'= as.numeric(gencode_v22_TSS$TSS)+TSS_dist,
    'strand'= rep('*',nrow(gencode_v22_TSS)),
    'name'= paste(
      sub(
        '\\..*',
        '',
        gencode_v22_TSS$gene_id
      ),
      gencode_v22_TSS$type,
      sep='_'
    ),
    stringsAsFactors = FALSE
  )
  rownames(gencode_v22_TSS_df) <- make.unique(
    sub(
      '\\..*',
      '',
      gencode_v22_TSS$gene_id
    )
  )

  ## remove uneeded datasets:
  rm(gencode_v22_TSS)

  ## Create a granges object from the TSS dataframe:
  gencode_v22_TSS_granges <- GenomicRanges::makeGRangesFromDataFrame(
    df= gencode_v22_TSS_df,
    keep.extra.columns = FALSE,
    starts.in.df.are.0based = FALSE
  )

  ## Find probes not overlap with TSS windows

  ## If user is looking for promoter probes, instead set probes_of_interest to
  ## isntead be the probes within the TSS buffer specified by the user.
  if(assess_promoter==FALSE){

    probes_of_interest <- rownames(
      hg38_manifest_no_NA_granges_df[
        setdiff(
          1:nrow(hg38_manifest_no_NA_granges_df),
          unique(
            S4Vectors::queryHits(
              GenomicRanges::findOverlaps(
                hg38_manifest_annotations_granges,
                gencode_v22_TSS_granges
              )
            )
          )
        ),
      ]
    )

  } else if(assess_promoter==TRUE){

    probes_of_interest <- rownames(
      hg38_manifest_no_NA_granges_df[
        unique(
          S4Vectors::queryHits(
            GenomicRanges::findOverlaps(
              hg38_manifest_annotations_granges,
              gencode_v22_TSS_granges
            )
          )
        ),
      ]
    )

  } else{

    ## If assess promoter is not set to TRUE or FALSE, return an error:
    stop('assess_promoter argument is not set to either TRUE or FALSE. Please set it to either TRUE or FALSE depending on if you want to assess DNA methylation probes in either promoter regions or distal enhance regions, respectively.')
  }

  ## Clear the workspace:
  rm(gencode_v22_TSS_df, gencode_v22_TSS_granges, hg38_manifest_annotations_granges)

  ## Now that probes marking regulatory regions of interest have been identified,
  ## let's import the probes in regions of histone modifications, NDR,
  ## or consensus enhancers from step1:

  ## check to see if the user has specified any histone modification, NDR, or
  ## enhancer data to be used:
  if(
    use_ext_HM==TRUE | use_ext_NDR==TRUE | use_consensus_ENH==TRUE | use_consensus_NDR==TRUE
  ){

    ## Check to see that the step1 directory at least exists
    if(
      dir.exists(
        paste(
          TENET_directory,
          "step1/",
          sep=''
        )
      )
    ){

      ## If the directory exists, depending on the analysis types selected, insure
      ## that their specific subdirectory exists in the step1 directory, and if it
      ## does, load the probelist files contained within:

      ## ext_HM:
      if(use_ext_HM==TRUE){

        ## Check to see that the ext_HM directory exists:
        if(
          dir.exists(
            paste(
              TENET_directory,
              "step1/",
              'ext_HM',
              sep=''
            )
          )
        ){

          ## If there directory exists, load the probelist files contained within:
          ## list all files in the directory:
          ext_HM_filelist_placeholder <- list.files(
            paste(
              TENET_directory,
              "step1/",
              'ext_HM',
              sep=''
            ),
            full.names = TRUE
          )

          ## Make sure only probelist.txt files are grabbed:
          ext_HM_filelist <- grep(
            'probelist.txt',
            ext_HM_filelist_placeholder,
            value=TRUE
          )

          ## Check to see that any files were loaded.
          ## If no files were loaded, issue an error
          if(
            length(ext_HM_filelist)==0
          ){

            stop("No ext_HM .probelist files were found for user supplied HM regions in the step1/consensus_NDR subdirectory of the specified 'TENET_directory'. Please insure the step1_make_external_datasets function was run with ext_HM set to TRUE, and that the .probelist files for user supplied HM regions files are still contained within.")
          }

        } else{

          ## Directory doesn't exist so return an error:
          stop('The ext_HM subdirectory with .probelist files does not exist in the step1 subdirectory of the TENET_directory. Please check that the step1_make_external_datasets function was run on the user-supplied histone modification files and rerun if necessary.')
        }
      }

      ## ext_NDR
      if(use_ext_NDR==TRUE){

        ## Check to see that the ext_NDR directory exists:
        if(
          dir.exists(
            paste(
              TENET_directory,
              "step1/",
              'ext_NDR',
              sep=''
            )
          )
        ){

          ## If there directory exists, load the probelist files contained within:
          ## list all files in the directory:
          ext_NDR_filelist_placeholder <- list.files(
            paste(
              TENET_directory,
              "step1/",
              'ext_NDR',
              sep=''
            ),
            full.names = TRUE
          )

          ## Make sure only probelist.txt files are grabbed:
          ext_NDR_filelist <- grep(
            'probelist.txt',
            ext_NDR_filelist_placeholder,
            value=TRUE
          )

          ## Check to see that any files were loaded.
          ## If no files were loaded, issue an error
          if(
            length(ext_NDR_filelist)==0
          ){

            stop("No ext_NDR .probelist files were found for user supplied NDR regions in the step1/consensus_NDR subdirectory of the specified 'TENET_directory'. Please insure the step1_make_external_datasets function was run with ext_NDR set to TRUE, and that the .probelist files for user supplied NDR regions files are still contained within.")
          }

        } else{

          ## Directory doesn't exist so return an error:
          stop('The ext_NDR subdirectory with .probelist files does not exist in the step1 subdirectory of the TENET_directory. Please check that the step1_make_external_datasets function was run on the user-supplied nucleosome depleted regions files and rerun if necessary.')
        }
      }

      ## consensus_ENH:
      if(use_consensus_ENH==TRUE){

        ## Check to see that the ext_HM directory exists:
        if(
          dir.exists(
            paste(
              TENET_directory,
              "step1/",
              'consensus_ENH',
              sep=''
            )
          )
        ){

          ## If there directory exists, load the probelist files contained within:
          ## list all files in the directory:
          consensus_ENH_filelist_placeholder <- list.files(
            paste(
              TENET_directory,
              "step1/",
              'consensus_ENH',
              sep=''
            ),
            full.names = TRUE
          )

          ## Make sure only probelist.txt files are grabbed:
          consensus_ENH_filelist <- grep(
            'probelist.txt',
            consensus_ENH_filelist_placeholder,
            value=TRUE
          )

          ## Check to see that any files were loaded.
          ## If no files were loaded, issue an error
          if(
            length(consensus_ENH_filelist)==0
          ){

            stop("No consensus_ENH.probelist file was found in the step1/consensus_ENH subdirectory of the specified 'TENET_directory'. Please insure the step1_make_external_datasets function was run with consensus_ENH set to TRUE, and that the .probelist file is still contained within.")
          }

        } else{

          ## Directory doesn't exist so return an error:
          stop('The consensus_ENH subdirectory with .probelist files does not exist in the step1 subdirectory of the TENET_directory. Please check that the step1_make_external_datasets function was run on the consensus enhancer dataset and rerun if necessary.')
        }
      }

      ## consensus_NDR:
      if(use_consensus_NDR==TRUE){

        ## Check to see that the ext_HM directory exists:
        if(
          dir.exists(
            paste(
              TENET_directory,
              "step1/",
              'consensus_NDR',
              sep=''
            )
          )
        ){

          ## If there directory exists, load the probelist files contained within:
          ## list all files in the directory:
          consensus_NDR_filelist_placeholder <- list.files(
            paste(
              TENET_directory,
              "step1/",
              'consensus_NDR',
              sep=''
            ),
            full.names = TRUE
          )

          ## Make sure only probelist.txt files are grabbed:
          consensus_NDR_filelist <- grep(
            'probelist.txt',
            consensus_NDR_filelist_placeholder,
            value=TRUE
          )

          ## Check to see that any files were loaded.
          ## If no files were loaded, issue an error
          if(
            length(consensus_NDR_filelist)==0
          ){

            stop("No consensus_NDR.probelist file was found in the step1/consensus_NDR subdirectory of the specified 'TENET_directory'. Please insure the step1_make_external_datasets function was run with consensus_NDR set to TRUE, and that the .probelist file is still contained within.")
          }

        } else{

          ## Directory doesn't exist so return an error:
          stop('The consensus_NDR subdirectory with .probelist files does not exist in the step1 subdirectory of the TENET_directory. Please check that the step1_make_external_datasets function was run on the consensus nucleosome depleted regions dataset and rerun if necessary.')
        }
      }

      ## Combine files into a HM and NDR listing:

      ## Combine HM and ENH files:
      HM_probelist_files <- c()

      if(use_ext_HM==TRUE){

        HM_probelist_files <- c(
          HM_probelist_files,
          ext_HM_filelist
        )
      }

      if(use_consensus_ENH==TRUE){

        HM_probelist_files <- c(
          HM_probelist_files,
          consensus_ENH_filelist
        )
      }

      ## Combine NDR files:
      NDR_probelist_files <- c()

      if(use_ext_NDR==TRUE){

        NDR_probelist_files <- c(
          NDR_probelist_files,
          ext_NDR_filelist
        )
      }

      if(use_consensus_NDR==TRUE){

        NDR_probelist_files <- c(
          NDR_probelist_files,
          consensus_NDR_filelist
        )
      }

      ## Now for each of the probelist files assembled,
      ## get the CpGs from it and add them to master list
      ## if probelist files were found

      ## Get the HM probes listed:
      if(length(HM_probelist_files)>0){

        HM_raw_probes <- character()

        for(i in HM_probelist_files){

          ## Read in the first file:
          probelist_HM_placeholder <- read.delim(
            i,
            header= FALSE,
            stringsAsFactors = FALSE
          )

          ## Add the probes to the probelist:
          HM_raw_probes <- c(
            HM_raw_probes,
            probelist_HM_placeholder$V1
          )

        }
      }

      ## Get the NDR probes listed:
      if(length(NDR_probelist_files)>0){

        NDR_raw_probes <- character()

        for(i in NDR_probelist_files){

          ## Read in the first file:
          probelist_NDR_placeholder <- read.delim(
            i,
            header= FALSE,
            stringsAsFactors = FALSE
          )

          ## Add the probes to the probelist:
          NDR_raw_probes <- c(
            NDR_raw_probes,
            probelist_NDR_placeholder$V1
          )

        }
      }

      ## Now get the unique probes from each dataset:

      ## Get the unique HM probes:
      if(length(HM_raw_probes)>0){

        HM_unique_probes <- unique(HM_raw_probes)

      }

      ## Get the unique NDR probes:
      if(length(NDR_raw_probes)>0){

        NDR_unique_probes <- unique(NDR_raw_probes)

      }

      ## If both sets of probes are available, get the ones that are present
      ## in both. Otherwise set the probes present to only be the ones present
      ## in the sets available:
      if(
        (length(HM_raw_probes)>0) & (length(NDR_raw_probes)>0)
      ){

        HM_NDR_probes <- intersect(
          HM_unique_probes,
          NDR_unique_probes
        )

        probes_of_interest <- intersect(
          HM_NDR_probes,
          probes_of_interest
        )

      } else if(
        length(HM_raw_probes)>0
      ){

        ## If only HM probes are selected, set the probes of interest
        ## to be the ones that overlap with all the probes:
        probes_of_interest <- intersect(
          HM_unique_probes,
          probes_of_interest
        )

      } else if(
        length(NDR_raw_probes)>0
      ){

        ## If only NDR probes are selected, set the probes of interest
        ## to be the ones that overlap with all the probes:
        probes_of_interest <- intersect(
          NDR_unique_probes,
          probes_of_interest
        )
      }

    } else{

      ## One of the HM, NDR, or enhancer datasets is selected, but no step1 folder
      ## is found, so return an error:
      stop(
        "The step1 folder containing results from the step1_make_external_datasets function was not found in the specified 'TENET_directory'. Please ensure the step1_make_external_datasets function was run on the specified hisone modification, nucleosome depleted regions, or enhancer datasets and the step1 subdirectory containing further subdirectories with .probelist files is present in the 'TENET_directory'."
      )
    }

  } else{

    ## No HM, enhancer, or NDR datasets were selected for inclusion
    ## so use all probes available for analysis:
    probes_of_interest <- probes_of_interest
  }

  ## Next match the genes in the expression datasets with those
  ## from the GENCODE v22 annotation:
  matched_genes_control <- intersect(
    sub(
      '\\..*',
      '',
      gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
        'gene_id'
      ]
    ),
    rownames(expDataN)
  )

  matched_genes_experimental <- intersect(
    sub(
      '\\..*',
      '',
      gencode_v22_gtf[
        gencode_v22_gtf$type=='gene',
        'gene_id'
      ]
    ),
    rownames(expDataT)
  )

  matched_genes <- intersect(
    matched_genes_control,
    matched_genes_experimental
  )

  ## Remove uneeded datasets:
  rm(gencode_v22_gtf)

  ## Next match the genes in the methylation datasets with those
  ## from the hm450 array annotation
  matched_probes_control <- intersect(
    hg38_manifest_no_NA_granges_df$names,
    rownames(metDataN)
  )

  matched_probes_experimental <- intersect(
    hg38_manifest_no_NA_granges_df$names,
    rownames(metDataT)
  )

  matched_probes <- intersect(
    matched_probes_control,
    matched_probes_experimental
  )

  ## Now subset each of the datasets to only the
  ## matched genes/probes:
  expDataN <- expDataN[
    matched_genes,
  ]

  expDataT <- expDataT[
    matched_genes,
  ]

  metDataN <- metDataN[
    matched_probes,
  ]

  metDataT <- metDataT[
    matched_probes,
  ]

  ## Now make sure the datasets still have values
  ## and return appropriate error message if they do not
  ## (i.e. data was not found in GENCODE v22 or HM450 array):
  if(length(matched_genes_control)==0 & length(matched_genes_experimental)==0){

    stop(
      "No GENCODE v22-annotated genes were found in the expression datasets. /nPlease reload the expression datasets with GENCODE v22-annotated Ensembl gene IDs in the rownames."
    )

  } else if(length(matched_genes_control)==0){

    stop(
      "No GENCODE v22-annotated genes were found in the control (normal) expression dataset. /nPlease reload the expression dataset with GENCODE v22-annotated Ensembl gene IDs in the rownames."
    )

  } else if(length(matched_genes_experimental)==0){

    stop(
      "No GENCODE v22-annotated genes were found in the experimental (tumor) expression dataset. /nPlease reload the expression dataset with GENCODE v22-annotated Ensembl gene IDs in the rownames."
    )

  } else if(length(matched_genes)==0){

    stop(
      "No matched GENCODE v22-annotated genes were found between the two expression datasets. /nPlease reload the expression datasets with GENCODE v22-annotated Ensembl gene IDs in the rownames."
    )

  }

  if(length(matched_probes_control)==0 & length(matched_probes_experimental)==0){

    stop(
      "No DNA methylation probes of the specified array were found in the methylation datasets. /nPlease reload the methylation datasets with HM450 probe IDs in the rownames."
    )

  } else if(length(matched_probes_control)==0){

    stop(
      "No DNA methylation probes of the specified array were found in the control (normal) methylation dataset. /nPlease reload the methylation dataset with HM450 probe IDs in the rownames."
    )

  } else if(length(matched_probes_experimental)==0){

    stop(
      "No DNA methylation probes of the specified array were found in the experimental (tumor) methylation dataset. /nPlease reload the methylation dataset with HM450 probe IDs in the rownames."
    )

  } else if(length(matched_probes)==0){

    stop(
      "No matched DNA methylation probes of the specified array were found between the two methylation datasets. /nPlease reload the methylation datasets with HM450 probe IDs in the rownames."
    )

  }

  ## First get quadrant probes if purity data isn't specified:
  if(purity_directory==FALSE){

    ## Let's identify which of the probes of interest were found in the
    ## resulting methylation datasets:
    probes_of_interest_matched <- intersect(
      probes_of_interest,
      matched_probes
    )

    ## For each of these probes, calculate the mean methylation
    ## in tumor and normal samples:
    probelist_df <- data.frame(
      probe= probes_of_interest_matched,
      stringsAsFactors = FALSE
    )

    probelist_df$mean_normal <- apply(
      metDataN[
        probes_of_interest_matched,
      ],
      1,
      mean,
      na.rm=T
    )

    probelist_df$mean_tumor <- apply(
      metDataT[
        probes_of_interest_matched,
      ],
      1,
      mean,
      na.rm=T
    )

    ## Remove the probes with NaN values:
    probelist_df_present <- probelist_df[
      complete.cases(probelist_df),
    ]

    ## Create methylation datasets with only the non-NaN enhancer probes:
    metDataN_present <- metDataN[
      probelist_df_present$probe,
    ]

    metDataT_present <- metDataT[
      probelist_df_present$probe,
    ]

    ## Now let's get the methylated quadrant probes:

    ## First get a copy of the experimental (tumor) data for probes
    ## with average methylation above the methcutoff
    ## in normal and tumor samples:
    temp_t <- metDataT_present[
      match(
        probelist_df_present[
          which(
            probelist_df_present$mean_normal>methcutoff & probelist_df_present$mean_tumor>methcutoff
          ),
          'probe'
        ],
        rownames(metDataT_present)
      ),
    ]

    ## Then for each experimental (tumor) sample,
    ## identify if it individually is below the methcutoff
    temp_t_cat <- ifelse(
      temp_t<methcutoff,
      1,
      0
    )

    ## For each of these probes, then calculate the number of individual
    ## experimental (tumor) samples below the methcutoff
    length_temp_t_cat <- apply(
      temp_t_cat,
      1,
      sum,
      na.rm=T
    )

    ## Create a new matrix with the experimental (tumor) methylation data
    ## for probes that have a number of individual samples below the methcutoff
    ## less than the min_experimental_count value.
    ## These represent the methylated probes:
    temp_t1 <- temp_t[
      which(
        length_temp_t_cat<min_experimental_count
      ),
    ]

    ## Now lets create dataframes with this methylated data:
    methDataT <- temp_t1

    methDataN <- metDataN[
      match(
        rownames(methDataT),
        rownames(metDataN)
      ),
    ]

    ## Now let's get the hypomethylated quadrant probes:

    ## First get a copy of the experimental (tumor) data for probes
    ## with average methylation above the methcutoff
    ## in just the normal samples:
    temp_t <- metDataT_present[
      match(
        probelist_df_present[
          which(
            probelist_df_present$mean_normal>methcutoff
          ),
          'probe'
        ],
        rownames(metDataT_present)
      ),
    ]

    ## Then for each experimental (tumor) sample,
    ## identify if it individually is below the hypomethcutoff
    temp_t_cat <- ifelse(
      temp_t<hypomethcutoff,
      1,
      0
    )

    ## For each of these probes, then calculate the number of individual
    ## experimental (tumor) samples below the hypomethcutoff
    length_temp_t_cat <- apply(
      temp_t_cat,
      1,
      sum,
      na.rm=T
    )

    ## Create a new matrix with the experimental (tumor) methylation data
    ## for probes that have a number of individual samples below the hypomethcutoff
    ## greater than the min_experimental_count value.
    ## These represent the hypomethylated probes:
    temp_t1 <- temp_t[
      which(
        length_temp_t_cat>min_experimental_count
      ),
    ]

    ## Now lets create dataframes with this hypomethylated data:
    hypomethDataT <- temp_t1

    hypomethDataN <- metDataN[
      match(
        rownames(hypomethDataT),
        rownames(metDataN)
      ),
    ]

    ## Now let's get the unmethylated quadrant probes:

    ## First get a copy of the experimental (tumor) data for probes
    ## with average methylation blow the unmethcutoff
    ## in normal and tumor samples:
    temp_t <- metDataT_present[
      match(
        probelist_df_present[
          which(
            probelist_df_present$mean_normal<unmethcutoff & probelist_df_present$mean_tumor<unmethcutoff
          ),
          'probe'
        ],
        rownames(metDataT_present)
      ),
    ]

    ## Then for each experimental (tumor) sample,
    ## identify if it individually is above the unmethcutoff
    temp_t_cat <- ifelse(
      temp_t>unmethcutoff,
      1,
      0
    )

    ## For each of these probes, then calculate the number of individual
    ## experimental (tumor) samples above the unmethcutoff
    length_temp_t_cat <- apply(
      temp_t_cat,
      1,
      sum,
      na.rm=T
    )

    ## Create a new matrix with the experimental (tumor) methylation data
    ## for probes that have a number of individual samples below the hypomethcutoff
    ## greater than the min_experimental_count value.
    ## These represent the unmethylated probes:
    temp_t1 <- temp_t[
      which(
        length_temp_t_cat<min_experimental_count
      ),
    ]

    ## Now lets create dataframes with this hypomethylated data:
    unmethDataT <- temp_t1

    unmethDataN <- metDataN[
      match(
        rownames(unmethDataT),
        rownames(metDataN)
      ),
    ]

    ## Now let's get the hypermethylated quadrant probes:

    ## First get a copy of the experimental (tumor) data for probes
    ## with average methylation blow the unmethcutoff
    ## in normal samples:
    temp_t <- metDataT_present[
      match(
        probelist_df_present[
          which(
            probelist_df_present$mean_normal<unmethcutoff
          ),
          'probe'
        ],
        rownames(metDataT_present)
      ),
    ]

    ## Then for each experimental (tumor) sample,
    ## identify if it individually is above the hypermethcutoff
    temp_t_cat <- ifelse(
      temp_t>hypermethcutoff,
      1,
      0
    )

    ## For each of these probes, then calculate the number of individual
    ## experimental (tumor) samples above the hypermethcutoff
    length_temp_t_cat <- apply(
      temp_t_cat,
      1,
      sum,
      na.rm=T
    )

    ## Create a new matrix with the experimental (tumor) methylation data
    ## for probes that have a number of individual samples above the hypermethcutoff
    ## greater than the min_experimental_count value.
    ## These represent the hypermethylated probes:
    temp_t1 <- temp_t[
      which(
        length_temp_t_cat>min_experimental_count
      ),
    ]

    ## Now lets create dataframes with this hypomethylated data:
    hypermethDataT <- temp_t1

    hypermethDataN <- metDataN[
      match(
        rownames(hypermethDataT),
        rownames(metDataN)
      ),
    ]

  } else if(purity_directory!=TRUE){

    ## Check to make sure the user has supplied a 'purity' subdirectory:
    ## If subdirectory isn't found return an error.
    if(!dir.exists(purity_directory)){

      stop(
        'The purity_directory was not found, although purity information was specified. Please put .rda files containing methylation data from other cell types for purity analyses in the specified purity_directory'
      )

    }

    ## Let's list all the files in the purity subdirectory:
    purity_file_list <- list.files(
      purity_directory,
      pattern=".rda",
      full.names = TRUE
    )

    purity_file_list_trunc <- list.files(
      purity_directory,
      pattern=".rda",
      full.names = FALSE
    )

    ## Make sure .rda files have been loaded:
    if(length(purity_file_list)==0){

      stop(
        'No purity .rda files were found in the purity subdirectory, although purity information was specified. Please put .rda files containing methylation data from other cell types for purity analyses in a separate subdirectory in the TENET_directory called "purity"'
      )

    }

    ## For each .rda file, load it and check the contents:

    ## First initialize an empty vector to contain name of objects in each
    ## .rda file:
    purity_rda_content_vector <- c()

    ## For each of the .rdas, load them, and check what was loaded and add it
    ## to the purity_rda_content_vector variable:
    for(purity_rda_file in purity_file_list){

      ## Create a snapshot of the environment to compare as .rdas are loaded:
      environment_snapshot <- ls()

      ## Load the .rda file:
      load(purity_rda_file)

      ## Check the new environment:
      new_environment_snapshot <- ls()

      ## Identify the objects loaded in the .rda file:
      purity_rda_objects_loaded <- setdiff(
        new_environment_snapshot,
        environment_snapshot
      )

      ## Remove the environment_snapshot variable from the objects loaded:
      purity_rda_objects_loaded <- purity_rda_objects_loaded[
        !purity_rda_objects_loaded=='environment_snapshot'
      ]

      ## Add the loaded objects to the purity_rda_content_vector:
      purity_rda_content_vector <- c(
        purity_rda_content_vector,
        purity_rda_objects_loaded
      )

      ## Remove the purity_rda_objects_loaded object:
      rm(purity_rda_objects_loaded)
    }

    ## Make sure that the objects loaded in the rdas is equal to
    ## the number of rdas present (i.e. 1 data object per rda):
    if(length(purity_rda_content_vector)!=length(purity_file_list)){

      stop(
        'The number of data objects loaded from the .rda files placed in the purity subdirectory are not equal to the number of .rda files discovered. Please ensure each of the purity .rda files included have only a single data object in each.'
      )
    }

    ## Initialize empty vectors to contain probe information from the
    ## different analysis quadrants as well as all the probes that overlapped
    ## each purity dataset:
    purity_dataset_probe_overlap_list <- list()
    purity_methylated_probe_vector <- list()
    purity_hypomethylated_probe_vector <- list()
    purity_unmethylated_probe_vector <- list()
    purity_hypermethylated_probe_vector <- list()

    ## Now for each of the loaded purity datasets, let's get a vector of the
    ## probes in them that match with the HM450 probe annotation:
    for(i in 1:length(purity_rda_content_vector)){

      ## Next match the genes in the expression datasets with those
      ## from the GENCODE v22 annotation:
      matched_probes_purity_dataset <- intersect(
        hg38_manifest_no_NA_granges_df$names,
        rownames(get(purity_rda_content_vector[i]))
      )

      ## If no probes are found, return an error message
      if(length(matched_probes_purity_dataset)==0){

        stop(
          paste(
            'No DNA methylation probes of the specified array were found in the',
            purity_file_list_trunc[i],
            'purity dataset. /nPlease reload the purity dataset with HM450 probe IDs in the rownames.',
            sep=' '
          )
        )
      }

      ## Let's identify which of the probes of interest were found in the
      ## resulting methylation datasets:
      probes_of_interest_matched <- intersect(
        probes_of_interest,
        matched_probes
      )

      ## Now check to make sure there are overlapping probes between the
      ## uploaded control/normal and experimental/tumor methylation datasets
      ## and the purity dataset:
      matched_probes_three_datasets <- intersect(
        matched_probes_purity_dataset,
        probes_of_interest_matched
      )

      ## If no probes were found overlapping the purity datasets
      ## and the control/normal and experimental/tumor methylation datasets
      ## return an error.
      if(length(matched_probes_three_datasets)==0){

        stop(
          paste(
            'No DNA methylation probes of the specified array were found overlapping between the',
            purity_file_list_trunc[i],
            'purity dataset, and the uploaded control (normal) experimental (tumor) methylation datasets. /nPlease check these datasets, which should include HM450 probe IDs in the rownames.',
            sep=' '
          )
        )

      }

      ## For each of these probes overlapped between the three datasets,
      ## calculate the mean methylation in tumor and normal samples:
      probelist_df <- data.frame(
        probe= matched_probes_three_datasets,
        stringsAsFactors = FALSE
      )

      probelist_df$mean_normal <- apply(
        metDataN[
          matched_probes_three_datasets,
        ],
        1,
        mean,
        na.rm=T
      )

      probelist_df$mean_tumor <- apply(
        metDataT[
          matched_probes_three_datasets,
        ],
        1,
        mean,
        na.rm=T
      )

      probelist_df$mean_purity <- apply(
        get(purity_rda_content_vector[i])[
          matched_probes_three_datasets,
        ],
        1,
        mean,
        na.rm=T
      )

      ## Remove the probes with NaN values:
      probelist_df_present <- probelist_df[
        complete.cases(probelist_df),
      ]

      ## Check to make sure there are still probes present and return an
      ## error if there aren't. If there are some, add them to the
      ## purity_dataset_probe_overlap_list list.
      if(length(rownames(probelist_df_present))==0){

        stop(
          paste(
            'No non-NA DNA methylation probes of the specified array were found overlapping between the',
            purity_file_list_trunc[i],
            'purity dataset, and the uploaded control (normal) experimental (tumor) methylation datasets. /nPlease check these datasets, making sure the probes have non-NA values and HM450 probe IDs in the rownames.',
            sep=' '
          )
        )

      } else{

        ## Add them to the list if they have overlap:
        purity_dataset_probe_overlap_list[[i]] <- probelist_df_present$probe

      }

      ## Create methylation datasets with only the non-NaN enhancer probes:
      metDataN_present <- metDataN[
        probelist_df_present$probe,
      ]

      metDataT_present <- metDataT[
        probelist_df_present$probe,
      ]

      ## First let's get the methylated quadrant probes:

      ## First get a copy of the experimental (tumor) data for probes
      ## with average methylation above the methcutoff
      ## in normal and tumor samples:
      temp_t <- metDataT_present[
        match(
          probelist_df_present[
            which(
              probelist_df_present$mean_normal>methcutoff & probelist_df_present$mean_tumor>methcutoff
            ),
            'probe'
          ],
          rownames(metDataT_present)
        ),
      ]

      ## Then for each experimental (tumor) sample,
      ## identify if it individually is below the methcutoff
      temp_t_cat <- ifelse(
        temp_t<methcutoff,
        1,
        0
      )

      ## For each of these probes, then calculate the number of individual
      ## experimental (tumor) samples below the methcutoff
      length_temp_t_cat <- apply(
        temp_t_cat,
        1,
        sum,
        na.rm=T
      )

      ## Get the probes which have fewer than min_experimental_count samples and save them
      ## as the methylated probes:
      purity_methylated_probe_vector[[i]] <- rownames(
        temp_t[
          which(
            length_temp_t_cat<min_experimental_count
          ),
        ]
      )

      ## Now let's get the hypomethylated quadrant probes:

      ## First get a copy of the experimental (tumor) data for probes
      ## with average methylation above the methcutoff
      ## in just the normal samples and purity dataset:
      temp_t <- metDataT_present[
        match(
          probelist_df_present[
            which(
              probelist_df_present$mean_normal>methcutoff & probelist_df_present$mean_purity>methcutoff
            ),
            'probe'
          ],
          rownames(metDataT_present)
        ),
      ]

      ## Then for each experimental (tumor) sample,
      ## identify if it individually is below the hypomethcutoff
      temp_t_cat <- ifelse(
        temp_t<hypomethcutoff,
        1,
        0
      )

      ## For each of these probes, then calculate the number of individual
      ## experimental (tumor) samples below the hypomethcutoff
      length_temp_t_cat <- apply(
        temp_t_cat,
        1,
        sum,
        na.rm=T
      )

      ## Get the probes which have more than min_experimental_count samples and save them
      ## as the hypomethylated probes:
      purity_hypomethylated_probe_vector[[i]] <- rownames(
        temp_t[
          which(
            length_temp_t_cat>min_experimental_count
          ),
        ]
      )

      ## Now let's get the unmethylated quadrant probes:

      ## First get a copy of the experimental (tumor) data for probes
      ## with average methylation blow the unmethcutoff
      ## in normal and tumor samples:
      temp_t <- metDataT_present[
        match(
          probelist_df_present[
            which(
              probelist_df_present$mean_normal<unmethcutoff & probelist_df_present$mean_tumor<unmethcutoff
            ),
            'probe'
          ],
          rownames(metDataT_present)
        ),
      ]

      ## Then for each experimental (tumor) sample,
      ## identify if it individually is above the unmethcutoff
      temp_t_cat <- ifelse(
        temp_t>unmethcutoff,
        1,
        0
      )

      ## For each of these probes, then calculate the number of individual
      ## experimental (tumor) samples above the unmethcutoff
      length_temp_t_cat <- apply(
        temp_t_cat,
        1,
        sum,
        na.rm=T
      )

      ## Get the probes which have less than min_experimental_count samples and save them
      ## as the unmethylated probes:
      purity_unmethylated_probe_vector[[i]] <- rownames(
        temp_t[
          which(
            length_temp_t_cat<min_experimental_count
          ),
        ]
      )

      ## Now let's get the hypermethylated quadrant probes:

      ## First get a copy of the experimental (tumor) data for probes
      ## with average methylation blow the unmethcutoff
      ## in normal samples and purity samples:
      temp_t <- metDataT_present[
        match(
          probelist_df_present[
            which(
              probelist_df_present$mean_normal<unmethcutoff & probelist_df_present$mean_purity<unmethcutoff
            ),
            'probe'
          ],
          rownames(metDataT_present)
        ),
      ]

      ## Then for each experimental (tumor) sample,
      ## identify if it individually is above the hypermethcutoff
      temp_t_cat <- ifelse(
        temp_t>hypermethcutoff,
        1,
        0
      )

      ## For each of these probes, then calculate the number of individual
      ## experimental (tumor) samples above the hypermethcutoff
      length_temp_t_cat <- apply(
        temp_t_cat,
        1,
        sum,
        na.rm=T
      )

      ## Get the probes which have more than min_experimental_count samples and save them
      ## as the hypermethylated probes:
      purity_hypermethylated_probe_vector[[i]] <- rownames(
        temp_t[
          which(
            length_temp_t_cat>min_experimental_count
          ),
        ]
      )
    }

    ## Now that we have assembled lists of the quadrant probes,
    ## as well as all the probes analyzed between each purity file
    ## identify those in common between each as the final datasets:
    purity_consensus_probes <- Reduce(
      intersect,
      purity_dataset_probe_overlap_list
    )

    purity_consensus_methylated_probes <- Reduce(
      intersect,
      purity_methylated_probe_vector
    )

    purity_consensus_hypomethylated_probes <- Reduce(
      intersect,
      purity_hypomethylated_probe_vector
    )

    purity_consensus_unmethylated_probes <- Reduce(
      intersect,
      purity_unmethylated_probe_vector
    )

    purity_consensus_hypermethylated_probes <- Reduce(
      intersect,
      purity_hypermethylated_probe_vector
    )

    ## Now using these consensus probes, create methylation dataframes
    ## for control (normal) and experimental (tumor) samples
    ## for just the probes in each quadrant as well as all
    ## enhancer probes that were analyzed:

    ## Start with enhancer probes common in all datasets:
    metDataT_present <- metDataT[
      purity_consensus_probes,
    ]

    metDataN_present <- metDataN[
      match(
        rownames(metDataT_present),
        rownames(metDataN)
      ),
    ]

    ## methylated:
    methDataT <- metDataT[
      purity_consensus_methylated_probes,
    ]

    methDataN <- metDataN[
      match(
        rownames(methDataT),
        rownames(metDataN)
      ),
    ]

    ## hypomethylated:
    hypomethDataT <- metDataT[
      purity_consensus_hypomethylated_probes,
    ]

    hypomethDataN <- metDataN[
      match(
        rownames(hypomethDataT),
        rownames(metDataN)
      ),
    ]

    ## unmethylated:
    unmethDataT <- metDataT[
      purity_consensus_unmethylated_probes,
    ]

    unmethDataN <- metDataN[
      match(
        rownames(unmethDataT),
        rownames(metDataN)
      ),
    ]

    ## hypermethylated:
    hypermethDataT <- metDataT[
      purity_consensus_hypermethylated_probes,
    ]

    hypermethDataN <- metDataN[
      match(
        rownames(hypermethDataT),
        rownames(metDataN)
      ),
    ]
  }

  ## Now let's create the final dataset as well as step2 metadata:

  ## First let's give new names to the final enhancer datasets:
  enhmetDataT <- metDataT_present
  enhmetDataN <- metDataN_present

  ## temporarily rebind expression and methylation data:
  temp_met_enh <- cbind(enhmetDataT, enhmetDataN)
  temp_met <- cbind(metDataT, metDataN)
  temp_exp <- cbind(expDataT, expDataN)

  ## Ensure that rows with all NA values have been removed from
  ## the datasets:
  temp_met_enh_NA_Index <- apply(temp_met_enh, 1, function(x) all(is.na(x)))
  temp_met_NA_Index <- apply(temp_met, 1, function(x) all(is.na(x)))
  temp_exp_NA_Index <- apply(temp_exp, 1, function(x) all(is.na(x)))

  temp_met_enh <- temp_met_enh[ !temp_met_enh_NA_Index, ]
  temp_met <- temp_met[ !temp_met_NA_Index, ]
  temp_exp <- temp_exp[ !temp_exp_NA_Index, ]

  expDataN <- expDataN[rownames(temp_exp),]
  expDataT <- expDataT[rownames(temp_exp),]

  metDataN <- metDataN[rownames(temp_met),]
  metDataT <- metDataT[rownames(temp_met),]

  enhmetDataN <- enhmetDataN[rownames(temp_met_enh),]
  enhmetDataT <- enhmetDataT[rownames(temp_met_enh),]

  ## Remove non-cg probes:
  enhmetDataN <- enhmetDataN[which(substring(rownames(enhmetDataN),1,2)=='cg'),]
  enhmetDataT <- enhmetDataT[which(substring(rownames(enhmetDataT),1,2)=='cg'),]

  hypermethDataN <- hypermethDataN[which(substring(rownames(hypermethDataN),1,2)=='cg'),]
  hypermethDataT <- hypermethDataT[which(substring(rownames(hypermethDataT),1,2)=='cg'),]

  methDataN <- methDataN[which(substring(rownames(methDataN),1,2)=='cg'),]
  methDataT <- methDataT[which(substring(rownames(methDataT),1,2)=='cg'),]

  hypomethDataN <- hypomethDataN[which(substring(rownames(hypomethDataN),1,2)=='cg'),]
  hypomethDataT <- hypomethDataT[which(substring(rownames(hypomethDataT),1,2)=='cg'),]

  unmethDataN <- unmethDataN[which(substring(rownames(unmethDataN),1,2)=='cg'),]
  unmethDataT <- unmethDataT[which(substring(rownames(unmethDataT),1,2)=='cg'),]

  ## Let's get the names of all the non-cg probes and save them as variables:
  unmeth_probes <- rownames(unmethDataN)
  meth_probes <- rownames(methDataN)
  hypometh_probes <- rownames(hypomethDataN)
  hypermeth_probes <- rownames(hypermethDataN)

  ## Note: these probes will represent the promoter probes if that analysis is selected:
  enhancer_probes <- rownames(enhmetDataN)


  ## Create a dataset with the number of different types of probes found
  ## as well as the total number of genes and probes analyzed
  ## and create a metadata table for it:
  names <- c(
    'unmeth_probes',
    'hypermeth_probes',
    'meth_probes',
    'hypometh_probes',
    'enhmet_probes',
    'total_genes',
    'total_probes'
  )

  values <- c(
    nrow(unmethDataN),
    nrow(hypermethDataN),
    nrow(methDataN),
    nrow(hypomethDataN),
    nrow(enhmetDataN),
    nrow(expDataN),
    nrow(metDataN)
  )

  overall_data <- data.frame(
    'data'=names,
    'values'=values,
    stringsAsFactors = FALSE
  )

  ## Now write the metadata to the step2 folder:
  write.table(
    x= overall_data,
    file= paste(
      TENET_directory,
      'step2/',
      "TENET_step2_overall_metadata.txt",
      sep=''
    ),
    quote= FALSE,
    sep= '\t',
    row.names= FALSE,
    col.names = FALSE
  )

  if(exists('clinical')){

    ## Now let's save the final dataset:
    save(
      metDataT, metDataN, expDataT, expDataN, enhancer_probes, unmeth_probes, meth_probes, hypometh_probes, hypermeth_probes, hypermethcutoff, hypomethcutoff, min_experimental_count, clinical,
      file= paste(
        TENET_directory,
        'step2/',
        "diff.methylated.datasets.rda",
        sep=''
      )
    )

  } else{

    ## Now let's save the final dataset without the clinical data:
    save(
      metDataT, metDataN, expDataT, expDataN, enhancer_probes, unmeth_probes, meth_probes, hypometh_probes, hypermeth_probes, hypermethcutoff, hypomethcutoff, min_experimental_count,
      file= paste(
        TENET_directory,
        'step2/',
        "diff.methylated.datasets.rda",
        sep=''
      )
    )
  }
}

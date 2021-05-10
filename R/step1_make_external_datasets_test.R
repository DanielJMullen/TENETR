#' step1_make_external_datasets
#'
#' This is the step1 function of the TENETR package.
#' This function allows users to use bed-like files see: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#' which contain either histone modification or open chromatin data contained within directories
#' specified by the user, or to use pre developed consensus enhancer and open chromatin
#' datasets from many cell/tissue types contained in the TENETR.data package,
#' in order to identify the hg38-annotated DNA methylation probes of the
#' user-specified array that fall in the peaks/regions specified by each file,
#' then saves .probelist.txt files listing the probes to a directory specified
#' by the user, which will also act as the 'TENET_directory' for future functions
#' in the TENETR package.
#'
#' @param output_directory Set a path to a directory, which will be created by this function, if it doesn't already exist, where you want the .probelist.txt files created by this function to be exported to. This directory will be the 'TENET_directory' utilized in later functions and a 'step1' subdirectory will be created in that directory, with further subdirectories containing the output .probelist.txt files made by this function.
#' @param DNA_methylation_manifest Set to 'HM27', 'HM450', or 'EPIC' depending on the DNA methylation array of interest for the user's data. hg38 array annotations come from https://zwdzwd.github.io/InfiniumAnnotation. Defaults to 'HM450'.
#' @param ext_HM If the user has their own histone modification datasets, set this as a path to a directory containing either .bed, .narrowPeak, .broadPeak, or .gappedPeak files with the user's own histone modification datasets. Otherwise, set to FALSE. Defaults to FALSE.
#' @param ext_NDR If the user has their own open chromatin, or NDR, datasets, set this as a path to a directory containing either .bed, .narrowPeak, .broadPeak, or .gappedPeak files with the user's own open chromatin datasets. Otherwise, set to FALSE. Defaults to FALSE.
#' @param consensus_ENH Set TRUE or FALSE depending on if user would like to use the consensus enhancer histone modification data from TENETR.data for analysis. Defaults to TRUE.
#' @param consensus_NDR Set TRUE or FALSE depending on if user would like to use the consensus open chromatin data from TENETR.data for analysis. Defaults to TRUE.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Returns .probelist.txt fils with the same name as the input files listing the hg38-annotated DNA methylation probes that fell within the regions specified by each file. These files will be used in later TENETR analyses.
#' @export

step1_make_external_datasets <- function(
  output_directory,
  DNA_methylation_manifest='HM450',
  ext_HM=FALSE,
  ext_NDR=FALSE,
  consensus_ENH=TRUE,
  consensus_NDR=TRUE,
  core_count
){

  ## If user has not supplied the final '/' in the output directory
  ## add it:
  output_directory <- ifelse(
    substring(
      output_directory,
      nchar(output_directory),
      nchar(output_directory)
    ) == '/',
    output_directory,
    paste(
      output_directory,
      '/',
      sep=''
    )
  )

  ## Create the output_directory if it does not already exist;
  if(
    !dir.exists(output_directory)
  ){

    dir.create(output_directory)
  }

  ## Create a step1 directory to deposit overlapped probe files:
  if(
    !dir.exists(
      paste(
        output_directory,
        'step1',
        sep=''
      )
    )
  ){

    dir.create(
      paste(
        output_directory,
        'step1/',
        sep=''
      )
    )
  }

  ## Load the hg38 DNA methylation annotations:
  ## Written in by Zexun Wu
  if (DNA_methylation_manifest == "HM450" | missing(DNA_methylation_manifest)) {

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

  ## Write a function to determine hg38 HM450 probes that
  ## overlap with bed-like file and write them out
  ## as a .text file:
  user_probe_overlap_writing_function <- function(
    bed_like_file_path,
    output_directory_for_file
  ){

    ## Get the bed file name without the path :
    bedlike_file_basename <- sub('.*/', '', bed_like_file_path)

    ## Load the bed file as a dataframe:
    bedlike_file_df <- read.delim(
      bed_like_file_path,
      header= FALSE,
      stringsAsFactors = FALSE
    )

    ## Create a modified dataframe to later convert to granges:
    bedlike_file_granges_df <- data.frame(
      'chr'= bedlike_file_df$V1,
      'start'= bedlike_file_df$V2,
      'end'= bedlike_file_df$V3,
      'strand'= rep(
        '*',
        nrow(bedlike_file_df)
      ),
      stringsAsFactors = FALSE
    )

    ## Create a granges object from the new bed file df:
    bedlike_file_granges <- GenomicRanges::makeGRangesFromDataFrame(
      df= bedlike_file_granges_df ,
      keep.extra.columns = FALSE,
      starts.in.df.are.0based = TRUE
    )

    ## Name the peaks after the chromosome, start, and end positions
    names(bedlike_file_granges) <- paste(
      bedlike_file_granges_df$chr,
      ':',
      bedlike_file_granges_df$start,
      '-',
      bedlike_file_granges_df$end,
      sep=''
    )

    ## Get the names of the CpGs that overlapped with the bed file:
    CpG_overlap_with_bedlike <- hg38_manifest_no_NA_granges_df[
      unique(
        S4Vectors::subjectHits(
          suppressWarnings(
            GenomicRanges::findOverlaps(
              bedlike_file_granges,
              hg38_manifest_annotations_granges
            )
          )
        )
      ),
      'names'
      ]

    ## Sort the CpGs and write to file:
    CpG_overlap_with_bedlike_sorted <- sort(CpG_overlap_with_bedlike)

    ## Write the unique sorted CpGs out to the file containing the overlap data:
    write(
      CpG_overlap_with_bedlike_sorted,
      file= paste(
        output_directory_for_file,
        bedlike_file_basename,
        '.probelist.txt',
        sep=''
      ),
      ncolumns = 1
    )
  }

  ## If external ENH files are used run this section:
  if(!ext_HM==FALSE){

    ## Create a directory in the output folder to deposit files if it does
    ## not already exist:
    if(
      !dir.exists(
        paste(
          output_directory,
          'step1/',
          'ext_HM',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          output_directory,
          'step1/',
          'ext_HM',
          sep=''
        )
      )
    }

    ## Check that the supplied directory exists.
    ## If it does, load any .bed, .narrowPeak, or .broadPeak files inside
    ## otherwise, return an error:
    if(
      dir.exists(
        ext_HM
      )
    ){

      ## List all the ENH bed files found:
      HM_bed_files <- list.files(
        path= ext_HM,
        pattern = ".bed",
        full.names = TRUE
      )

      ## List all the ENH narrowPeak files:
      HM_narrowPeak_files <- list.files(
        path= ext_HM,
        pattern = ".narrowPeak",
        full.names = TRUE
      )

      ## List all the ENH broadPeak files:
      HM_broadPeak_files <- list.files(
        path= ext_HM,
        pattern = ".broadPeak",
        full.names = TRUE
      )

      ## List all the ENH gappedPeak files:
      HM_gappedPeak_files <- list.files(
        path= ext_HM,
        pattern = ".gappedPeak",
        full.names = TRUE
      )

      ## Combine the different file types into
      ## an ENH master list to load:
      HM_file_list <- c(
        HM_bed_files,
        HM_narrowPeak_files,
        HM_broadPeak_files,
        HM_gappedPeak_files
      )

      ## Sapply the function on the HM files input by the user:
      invisible(
        parallel::mclapply(
          HM_file_list,
          FUN= user_probe_overlap_writing_function,
          mc.cores= core_count,
          output_directory_for_file= paste(
            output_directory,
            'step1/',
            'ext_HM/',
            sep=''
          )
        )
      )

    } else{

      ## Supplied directory not found, return an error:
      stop('The directory path given for ext_HM was not found. Please check the argument and supply a path to a directory with .bed, .narrowPeak, or .gappedPeak files with regions of histone modification information. Otherwise, set this as FALSE')
    }
  }

  ## If external NDR files are used run this section:
  if(ext_NDR!=FALSE){

    ## Create a directory in the output folder to deposit files if it does
    ## not already exist:
    if(
      !dir.exists(
        paste(
          output_directory,
          'step1/',
          'ext_NDR',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          output_directory,
          'step1/',
          'ext_NDR',
          sep=''
        )
      )
    }

    ## Check that the supplied directory exists.
    ## If it does, load any .bed, .narrowPeak, or .broadPeak files inside
    ## otherwise, return an error:
    if(
      dir.exists(
        ext_NDR
      )
    ){

      ## List all the NDR bed files found:
      NDR_bed_files <- list.files(
        path= ext_NDR,
        pattern = ".bed",
        full.names = TRUE
      )

      ## List all the NDR narrowPeak files:
      NDR_narrowPeak_files <- list.files(
        path= ext_NDR,
        pattern = ".narrowPeak",
        full.names = TRUE
      )

      ## List all the NDR broadPeak files:
      NDR_broadPeak_files <- list.files(
        path= ext_NDR,
        pattern = ".broadPeak",
        full.names = TRUE
      )

      ## List all the NDR gappedPeak files:
      NDR_gappedPeak_files <- list.files(
        path= ext_NDR,
        pattern = ".gappedPeak",
        full.names = TRUE
      )

      ## Combine the different file types into
      ## an NDR master list to load:
      NDR_file_list <- c(
        NDR_bed_files,
        NDR_narrowPeak_files,
        NDR_broadPeak_files,
        NDR_gappedPeak_files
      )

      ## Sapply the function on the NDR files input by the user:
      invisible(
        parallel::mclapply(
          NDR_file_list,
          FUN= user_probe_overlap_writing_function,
          mc.cores= core_count,
          output_directory_for_file= paste(
            output_directory,
            'step1/',
            'ext_NDR/',
            sep=''
          )
        )
      )

    } else{

      ## Supplied directory not found, return an error:
      stop('The directory path given for ext_NDR was not found. Please check the argument and supply a path to a directory with .bed, .narrowPeak, or .gappedPeak files with NDR/open chromatin region information. Otherwise, set this as FALSE')
    }
  }

  ## If consensus ENH files are used run this section:
  if(consensus_ENH==TRUE){

    ## Create a directory in the output folder to deposit files if it does
    ## not already exist:
    if(
      !dir.exists(
        paste(
          output_directory,
          'step1/',
          'consensus_ENH',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          output_directory,
          'step1/',
          'consensus_ENH',
          sep=''
        )
      )
    }

    ## Create a granges object from the consensus enhancer regions:
    bedlike_file_granges <- GenomicRanges::makeGRangesFromDataFrame(
      df= TENETR.data::consensus_enhancer_regions,
      keep.extra.columns = FALSE,
      starts.in.df.are.0based = FALSE
    )

    ## Name the peaks after the chromosome, start, and end positions
    names(bedlike_file_granges) <- paste(
      TENETR.data::consensus_enhancer_regions$chr,
      ':',
      TENETR.data::consensus_enhancer_regions$start,
      '-',
      TENETR.data::consensus_enhancer_regions$end,
      sep=''
    )

    ## Get the names of the CpGs that overlapped with the bed file:
    CpG_overlap_with_bedlike <- hg38_manifest_no_NA_granges_df[
      unique(
        S4Vectors::subjectHits(
          suppressWarnings(
            GenomicRanges::findOverlaps(
              bedlike_file_granges,
              hg38_manifest_annotations_granges
            )
          )
        )
      ),
      'names'
      ]

    ## Sort the CpGs and write to file:
    CpG_overlap_with_bedlike_sorted <- sort(CpG_overlap_with_bedlike)

    ## Write the unique sorted CpGs out to the file containing the overlap data:
    if(
      !file.exists(
        paste(
          output_directory,
          'consensus_ENH.probelist.txt',
          sep=''
        )
      )
    ){

      write(
        CpG_overlap_with_bedlike_sorted,
        file= paste(
          output_directory,
          'step1/',
          'consensus_ENH/',
          'consensus_ENH.probelist.txt',
          sep=''
        ),
        ncolumns = 1
      )

    } else{

      write(
        CpG_overlap_with_bedlike_sorted,
        file= paste(
          output_directory,
          'step1/',
          'consensus_ENH/',
          'consensus_ENH.TENETR.data.probelist.txt',
          sep=''
        ),
        ncolumns = 1
      )

    }
  }

  ## If consensus NDR (open chromatin) files are used run this section:
  if(consensus_NDR==TRUE){

    ## Create a directory in the output folder to deposit files if it does
    ## not already exist:
    if(
      !dir.exists(
        paste(
          output_directory,
          'step1/',
          'consensus_NDR',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          output_directory,
          'step1/',
          'consensus_NDR',
          sep=''
        )
      )
    }

    ## Create a granges object from the consensus enhancer regions:
    bedlike_file_granges <- GenomicRanges::makeGRangesFromDataFrame(
      df= TENETR.data::consensus_open_chromatin_regions,
      keep.extra.columns = FALSE,
      starts.in.df.are.0based = FALSE
    )

    ## Name the peaks after the chromosome, start, and end positions
    names(bedlike_file_granges) <- paste(
      TENETR.data::consensus_open_chromatin_regions$chr,
      ':',
      TENETR.data::consensus_open_chromatin_regions$start,
      '-',
      TENETR.data::consensus_open_chromatin_regions$end,
      sep=''
    )

    ## Get the names of the CpGs that overlapped with the bed file:
    CpG_overlap_with_bedlike <- hg38_manifest_no_NA_granges_df[
      unique(
        S4Vectors::subjectHits(
          suppressWarnings(
            GenomicRanges::findOverlaps(
              bedlike_file_granges,
              hg38_manifest_annotations_granges
            )
          )
        )
      ),
      'names'
    ]

    ## Sort the CpGs and write to file:
    CpG_overlap_with_bedlike_sorted <- sort(CpG_overlap_with_bedlike)

    ## Write the unique sorted CpGs out to the file containing the overlap data:
    if(
      !file.exists(
        paste(
          output_directory,
          'consensus_NDR.probelist.txt',
          sep=''
        )
      )
    ){

      write(
        CpG_overlap_with_bedlike_sorted,
        file= paste(
          output_directory,
          'step1/',
          'consensus_NDR/',
          'consensus_NDR.probelist.txt',
          sep=''
        ),
        ncolumns = 1
      )

    } else{

      write(
        CpG_overlap_with_bedlike_sorted,
        file= paste(
          output_directory,
          'step1/',
          'consensus_NDR/',
          'consensus_NDR.TENETR.data.probelist.txt',
          sep=''
        ),
        ncolumns = 1
      )
    }
  }
}

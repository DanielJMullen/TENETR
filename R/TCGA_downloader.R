#' TCGA_downloader
#'
#' This is a function to download and compile TCGA gene expression and DNA
#' methylation datasets, as well as clinical data primarily for the use with
#' the TENETR package. This function takes advantage of the TCGAbiolinks
#' downloading functions in a more simplified manner and identifies samples with
#' matching gene expression and DNA methylation data and can also remove
#' duplicate tumor samples taken from the same patient donor.
#'
#'
#' @param TCGA_directory Set a path to the directory where TCGAbiolinks can download data to. Note that this dataset can be very sizeable.
#' @param TCGA_study_abbreviation Input a four letter code for a TCGA dataset to download data for. See: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations for more info and complete list of options.
#' @param RNA_seq_workflow Select the type of RNA-seq data to download. For TENET purposes, choose either "HTSeq - FPKM" or "HTSeq - FPKM-UQ". "HTSeq - Counts" can also be used but is not recommended for TENET analyses. See: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/ for more info on these.
#' @param RNA_seq_log2_normalization Set TRUE or FALSE to do log2 normalization of RNA_seq expression values.
#' @param matching_exp_met_samples Set TRUE or FALSE to isolate only samples with matching expression and methylation data in both these datasets. Set to TRUE by default.
#' @param remove_dup_tumor Set TRUE or FALSE to remove duplicate tumor samples taking from the same subject, leaving only one per subject in alphanumeric order. Set to TRUE by default.
#' @param TENET_directory Set a path to the directory where you want the downloaded datasets to be exported to. This could be a directory where you want to eventually run TENETR analyses.
#' @return Currently returns a .rda file in the specified TENET_directory containing 5 elements: 'clinical' containing clinical data, 'expDataN' and 'expDataT' containing gene expression data for adjacent normal and tumor samples respectively, and 'metDataN' and 'metDataT' containing methylation data for normal and tumor samples respectively.
#' @export

## Write TCGA_downloader function:
TCGA_downloader <- function(
  TCGA_directory,
  TCGA_study_abbreviation,
  RNA_seq_workflow,
  RNA_seq_log2_normalization,
  matching_exp_met_samples=TRUE,
  remove_dup_tumor=TRUE,
  TENET_directory
){

  ## If user has not supplied the final '/' in the TCGA and TENET directories
  ## add it:
  TCGA_directory <- ifelse(
    substring(
      TCGA_directory,
      nchar(TCGA_directory),
      nchar(TCGA_directory)
    ) == '/',
    TCGA_directory,
    paste(
      TCGA_directory,
      '/',
      sep=''
    )
  )

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

  ## Change to the directory where the user wants files downloaded:
  setwd(TCGA_directory)

  ## Convert the supplied TCGA study abbreviation to an all caps
  ## and all lowercase version:
  TCGA_study_abbreviation_all_upper <- toupper(
    TCGA_study_abbreviation
  )

  TCGA_study_abbreviation_all_lower <- tolower(
    TCGA_study_abbreviation
  )

  ## Create the TCGA abbreviation used for downloads:
  TCGA_study_abbreviation_download <- paste(
    'TCGA-',
    TCGA_study_abbreviation_all_upper,
    sep=''
  )

  ## Set up the expression  query:
  expression_query <- TCGAbiolinks::GDCquery(
    project= TCGA_study_abbreviation_download,
    data.category= "Transcriptome Profiling",
    data.type= "Gene Expression Quantification",
    experimental.strategy= 'RNA-Seq',
    workflow.type= RNA_seq_workflow,
    legacy = FALSE
  )

  ## Download the expression data
  TCGAbiolinks::GDCdownload(expression_query)

  ## Put the expression data in a dataframe:
  expression_data <- TCGAbiolinks::GDCprepare(
    query=expression_query,
    summarizedExperiment= FALSE
  )
  expression_data <- as.data.frame(expression_data)

  ## Remove the stuff after the periods from the gene annotation:
  expression_data$X1 <- sub(
    '\\..*',
    '',
    expression_data$X1
  )

  ## Set the rownames of expression data to be the gene annotations:
  rownames(expression_data) <- expression_data$X1

  ## Delete that column afterwards:
  expression_data$X1 <- NULL

  ## Perform log2 normalization if selected:
  if(RNA_seq_log2_normalization==TRUE){

    expression_data <- log2(expression_data+1)

  }

  ## Download methylation data:
  methylation_query <- TCGAbiolinks::GDCquery(
    project= TCGA_study_abbreviation_download,
    data.category= "DNA Methylation",
    platform= "Illumina Human Methylation 450",
    legacy= FALSE
  )

  ## Download the methylation data
  TCGAbiolinks::GDCdownload(methylation_query)

  ## Put the methylation data in a dataframe:
  ## TCGAbiolinks GDCprepare function doesn't do this well
  ## So I have coded my own:

  ## First navigate to the directory that contains just the DNA methylation
  ## data for the analysis of interest:
  setwd(
    paste(
      '.',
      'GDCdata',
      TCGA_study_abbreviation_download,
      'harmonized',
      'DNA_Methylation',
      sep='/'
    )
  )

  ## Then, list out all the downloaded methylation files
  ## from the TCGA downloaded using the query:
  methylation_files_list_full_path <- list.files(
    pattern='HumanMethylation450',
    full.names = TRUE,
    recursive= TRUE
  )

  ## Write a for loop to load and bind the files to a dataset:
  for(i in 1:length(methylation_files_list_full_path)){

    ## Load the file and set the rownames to be the probeIDs
    file_placeholder <- read.delim(
      methylation_files_list_full_path[i],
      header= TRUE,
      sep='\t',
      stringsAsFactors = FALSE
    )

    ## Get the base name of the file
    file_name <- basename(methylation_files_list_full_path[i])

    ## Get part of file name starting with 'TCGA.'
    file_name_TCGA <- substr(
      file_name,
      regexpr(
        'TCGA',
        file_name
      ),
      nchar(file_name)
    )

    ## Get part of file name before the period and remaining info
    ## This should return the full sample name:
    file_name_sample <- sub(
      '\\..*',
      '',
      file_name_TCGA
    )

    if(!exists('methylation_data')){

      ## If a methylation file to compile info hasn't been created, create
      ## it now and add the methylation data from the first file to it.
      methylation_data <- data.frame(
        'placeholder'=file_placeholder$Beta_value
      )

      ## Then rename the first column to the sample name and the rownames
      ## to be the probe names:
      colnames(methylation_data) <- file_name_sample
      rownames(methylation_data) <- file_placeholder$Composite.Element.REF

    } else{

      methylation_data[[i]] <- file_placeholder$Beta_value

      ## Change the column names to reflect the sample name of the newest file:
      colnames(methylation_data) <- c(
        colnames(methylation_data)[
          1:(
            length(
              colnames(methylation_data)
            )-1
          )
        ],
        file_name_sample
      )
    }

    rm(file_placeholder)
    rm(file_name)
    rm(file_name_TCGA)
    rm(file_name_sample)
  }

  ## Set up the clinical query:
  clinical_query <- TCGAbiolinks::GDCquery(
    project= TCGA_study_abbreviation_download,
    data.category= "Clinical",
    file.type= "xml")

  ## Download the clinical data
  TCGAbiolinks::GDCdownload(
    clinical_query
  )

  ## Put the clinical data in a dataframe:
  ## Put a supress warning on this due to warning
  ## that pops up from TCGAbiolinks/dplyr usage:
  clinical_data <- suppressWarnings(
    TCGAbiolinks::GDCprepare_clinic(
      clinical_query,
      clinical.info= "patient"
    )
  )

  ## Get the unique samples from the data:
  clinical <- unique(clinical_data)

  ## Regardless of matching, cut down names of expression
  ## and methylation data to the first 19 characters:
  colnames(expression_data) <- substring(
    colnames(expression_data),
    1,
    19
  )

  colnames(methylation_data) <- substring(
    colnames(methylation_data),
    1,
    19
  )

  ## If matching_exp_met_samples is TRUE
  ## Find the samples with just matched tumor and expression data:
  if(matching_exp_met_samples==TRUE){

    ## Find the sample names that are present in both datasets:
    matched_exp_met_names <- intersect(
      colnames(expression_data),
      colnames(methylation_data)
    )

    ## Limit the expression and methylation data to just those samples:
    expression_data <- expression_data[matched_exp_met_names]

    methylation_data <- methylation_data[matched_exp_met_names]

  }

  ## Now let's isolate the tumor and normal samples:
  ## Tumor samples have 01-09 in their sample ID (characters 14 and 15)
  ## while normal samples have 10 -19.
  ## See: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
  expDataT <- expression_data[
    as.numeric(
      substring(
        colnames(expression_data),
        14,
        15
      )
    ) < 10
  ]

  expDataN <- expression_data[
    as.numeric(
      substring(
        colnames(expression_data),
        14,
        15
      )
    ) >= 10
  ]

  metDataT <- methylation_data[
    as.numeric(
      substring(
        colnames(methylation_data),
        14,
        15
      )
    ) < 10
  ]

  metDataN <- methylation_data[
    as.numeric(
      substring(
        colnames(methylation_data),
        14,
        15
      )
    ) >= 10
  ]

  ## Let's sort the samples of each dataset alphanumerically:
  expDataN <- expDataN[
    sort(
      colnames(expDataN),
      decreasing = FALSE
    )
  ]

  expDataT <- expDataT[
    sort(
      colnames(expDataT),
      decreasing = FALSE
    )
  ]

  metDataN <- metDataN[
    sort(
      colnames(metDataN),
      decreasing = FALSE
    )
  ]

  metDataT <- metDataT[
    sort(
      colnames(metDataT),
      decreasing = FALSE
    )
  ]

  ## Now if remove_dup_tumor is set to TRUE
  ## Remove the tumor samples that are duplicates from the same patient
  ## (leaving one)
  if(remove_dup_tumor==TRUE){

    ## Get the substring of the tumor sample names equal to the
    ## part of the barcode through the participant ID
    ## See: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    expDataT_colnames_trunc <- substring(
      colnames(expDataT),
      1,
      12
    )

    names(expDataT_colnames_trunc) <- colnames(expDataT)

    metDataT_colnames_trunc <- substring(
      colnames(metDataT),
      1,
      12
    )

    names(metDataT_colnames_trunc) <- colnames(metDataT)

    ## Now remove the extra duplicate samples, leaving the first of each
    expDataT_colnames_trunc_dup_removed <- expDataT_colnames_trunc[
      !duplicated(expDataT_colnames_trunc)
    ]

    metDataT_colnames_trunc_dup_removed <- metDataT_colnames_trunc[
      !duplicated(metDataT_colnames_trunc)
    ]

    ## Get get non-duplicated samples out of both the expression and
    ## methylation data:
    expDataT <- expDataT[
      names(expDataT_colnames_trunc_dup_removed)
    ]

    metDataT <- metDataT[
      names(metDataT_colnames_trunc_dup_removed)
    ]

  }

  ## Convert the data frames to matrices:
  expDataN <- as.matrix(expDataN)
  expDataT <- as.matrix(expDataT)

  metDataN <- as.matrix(metDataN)
  metDataT <- as.matrix(metDataT)

  ## Save these objcts, plus clinical data
  ## to an rda file exported to the specified
  ## TENET_directory:
  save(
    clinical,
    expDataN,
    expDataT,
    metDataN,
    metDataT,
    file= paste(
      TENET_directory,
      TCGA_study_abbreviation_all_lower,
      '.expression.methylation.dataset.rda',
      sep=''
    )
  )

}

#' get_diffmeth_regions
#'
#' This is the step2 function of the TENETR package.
#' This function first identifies hm450 DNA methylation probes
#' that fall within enhancer regions based on enhancer and open chromatin
#' datasets from the step1 make_external_datasets function as well as
#' a specified distance from GENCODE v22 transcript TSS.
#' After identifying enhancer DNA methylation probes, the function classifies
#' them as methylated/unmethylated or hypermethylated/hypomethylated based on
#' their differential methylation between the control/normal and
#' experimental/tumor samples, based on input parameters from the user
#'
#'
#' @param TENET_directory Set a path to the directory that contains the step1 directory created by make_external_datasets function as well as the .rda file containing methylation and expression data.
#' @param TSS_dist Set a number to be the buffer in base pairs from GENCODE v22-annotated transcription start sites for DNA methylation probes to not be considered enhancer probes.
#' @param methcutoff Set a number from 0 to 1 to be the beta-value cutoff for methylated probes.
#' @param hypomethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypomethylated probes. Should be set lower than the methcutoff.
#' @param unmethcutoff Set a number from 0 to 1 to be the beta-value cutoff for unmethylated probes.
#' @param hypermethcutoff Set a number from 0 to 1 to be the beta-value cutoff for hypermethylated probes. Should be set higher than the unmethcutoff.
#' @param minExp Sets the minimum number of experimental/tumor samples to be considered for the hypo/hypermethylated groups.
#' @param use_purity_data Set TRUE or FALSE to use purity datasets, as .rda files containing DNA methylation values supplied by the user in a 'purity' subdirectory in the TENET_directory, to select for DNA methylation probes not potentially due to differences in cell type purity.
#' @return Returns several objects including a .rda file with matrices of methylation data for the four quadrants in control/normal and experimental/tumor samples. Also outputs a simple .txt file containing metrics on the number of probes found in different categories.
#' @export

## Write the get_diffmeth_regions function:
get_diffmeth_regions <- function(
  TENET_directory,
  TSS_dist,
  methcutoff,
  hypomethcutoff,
  unmethcutoff,
  hypermethcutoff,
  minExp,
  use_purity_data
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

  ## Create a step2 directory to deposit overlapped probe files:
  dir.create(
    paste(
      TENET_directory,
      'step2/',
      sep=''
    )
  )

  ## Load the hg38 450k annotations:
  hg38_hm450_df <- TENETR.data::hm450_hg38_annotations
  rownames(hg38_hm450_df) <- hg38_hm450_df$probeID

  ## Create a modified dataframe of the hg38 450k annotations
  ## to later convert to granges:
  hg38_450_annotations_granges_df <- data.frame(
    'chr'= hg38_hm450_df$CpG_chrm,
    'start'= hg38_hm450_df$CpG_beg,
    'end'= hg38_hm450_df$CpG_end,
    'strand'= rep(
      '*',
      nrow(hg38_hm450_df)
    ),
    'names' = hg38_hm450_df$probeID,
    stringsAsFactors = FALSE
  )

  ## Remove the big hm450 dataframe:
  rm(hg38_hm450_df)

  ## Remove the probes that have NA values:
  hg38_450_annotations_no_NA_granges_df <- hg38_450_annotations_granges_df[
    !is.na(hg38_450_annotations_granges_df$chr),
  ]

  rownames(hg38_450_annotations_no_NA_granges_df) <- hg38_450_annotations_no_NA_granges_df$names

  ## Remove the dataset with NA probes:
  rm(hg38_450_annotations_granges_df)

  ## Create a granges object from the new
  ## hg38 450k annotations with no NA df:
  hg38_450_annotations_granges <- GenomicRanges::makeGRangesFromDataFrame(
    df= hg38_450_annotations_no_NA_granges_df,
    keep.extra.columns = FALSE,
    starts.in.df.are.0based = TRUE
  )
  names(hg38_450_annotations_granges) <- hg38_450_annotations_no_NA_granges_df$names

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

  # find probes not overlap with TSS windows
  nonTSS_probes <- rownames(
    hg38_450_annotations_no_NA_granges_df[
      setdiff(
        1:nrow(hg38_450_annotations_no_NA_granges_df),
        unique(
          S4Vectors::queryHits(
            GenomicRanges::findOverlaps(
              hg38_450_annotations_granges,
              gencode_v22_TSS_granges
            )
          )
        )
      ),
    ]
  )

  ## Clear the workspace:
  rm(gencode_v22_TSS_df, gencode_v22_TSS_granges, hg38_450_annotations_granges)

  ## Now that non-TSS probes have been identified, let's import the
  ## Enhancer + NDR probes from step1:

  ## List directories in step1 if it is present and produce error message
  ## if not:
  if(
    dir.exists(
      paste(
        TENET_directory,
        "step1/",
        sep=''
      )
    )
  ){

    ## Get the subdirectories in step1 if
    step1_directories <- list.dirs(
      paste(
        TENET_directory,
        "step1",
        sep=''
      ),
      full.names = TRUE
    )

    ## Seperate the directories into ones with ENH and NDR files:
    ENH_directories <- grep(
      'ENH',
      step1_directories,
      value=TRUE
    )

    NDR_directories <- grep(
      'NDR',
      step1_directories,
      value=TRUE
    )

    ## For each directory, go in and list the probelist.txt files
    ## inside and create a list:
    ENH_probelist_files <- character()
    NDR_probelist_files <- character()

    for(i in ENH_directories){

      ## list all files in the directory:
      ENH_filelist_placeholder <- list.files(
        i,
        full.names = TRUE
      )

      ## Make sure only probelist.txt files are grabbed:
      ENH_filelist_probelist <- grep(
        'probelist.txt',
        ENH_filelist_placeholder,
        value=TRUE
      )

      ## Add the files to ENH_probelist_files
      ENH_probelist_files <- c(
        ENH_probelist_files,
        ENH_filelist_probelist
      )
    }

    for(i in NDR_directories){

      ## list all files in the directory:
      NDR_filelist_placeholder <- list.files(
        i,
        full.names = TRUE
      )

      ## Make sure only probelist.txt files are grabbed:
      NDR_filelist_probelist <- grep(
        'probelist.txt',
        NDR_filelist_placeholder,
        value=TRUE
      )

      ## Add the files to NDR_probelist_files
      NDR_probelist_files <- c(
        NDR_probelist_files,
        NDR_filelist_probelist
      )
    }

    ## Let's check step1 to see if the user
    ## has input any custom probelist.txt files:
    custom_placeholder <- list.files(
      paste(
        TENET_directory,
        "step1/",
        sep=''
      ),
      full.names = TRUE
    )

    ## get just the probelist.txt files:
    custom_filelist_probelist <- grep(
      'probelist.txt',
      custom_placeholder,
      value=TRUE
    )

    ## Separate them into NDR and Enhancer files:
    ## Enhancer files allow the terms "H3K27ac", "enhancer", or "ENH"
    ## Open chromatin files allow "open_chromatin" or "NDR"
    ENH_custom_filelist_probelist <- grep(
      'H3K27ac',
      custom_filelist_probelist,
      ignore.case = TRUE,
      value=TRUE
    )

    ENH_custom_filelist_probelist <- c(
      ENH_custom_filelist_probelist,
      ENH_custom_filelist_probelist <- grep(
        'enhancer',
        custom_filelist_probelist,
        ignore.case = TRUE,
        value=TRUE
      )
    )

    ENH_custom_filelist_probelist <- c(
      ENH_custom_filelist_probelist,
      ENH_custom_filelist_probelist <- grep(
        'ENH',
        custom_filelist_probelist,
        ignore.case = FALSE,
        value=TRUE
      )
    )

    NDR_custom_filelist_probelist <- grep(
      'open_chromatin',
      custom_filelist_probelist,
      ignore.case = TRUE,
      value=TRUE
    )

    NDR_custom_filelist_probelist <- c(
      NDR_custom_filelist_probelist,
      grep(
        'NDR',
        custom_filelist_probelist,
        ignore.case = FALSE,
        value=TRUE
      )
    )

    ## Add any of the custom files to the file list:
    ENH_probelist_files <- c(
      ENH_probelist_files,
      ENH_custom_filelist_probelist
    )

    NDR_probelist_files <- c(
      NDR_probelist_files,
      NDR_custom_filelist_probelist
    )

    ## Clear workspace:
    rm(
      custom_filelist_probelist, custom_placeholder, ENH_custom_filelist_probelist,
       ENH_directories, ENH_filelist_placeholder, ENH_filelist_probelist,
       NDR_custom_filelist_probelist, NDR_directories, NDR_filelist_placeholder,
       NDR_filelist_probelist, step1_directories
    )

    ## Now for each of the probelist files assembled,
    ## get the CpGs from it and add them to master list
    ## if probelist files were found
    if(length(ENH_probelist_files)>0){

      ENH_raw_probes <- character()

      for(i in ENH_probelist_files){

        ## Read in the first file:
        probelist_ENH_placeholder <- read.delim(
          i,
          header= FALSE,
          stringsAsFactors = FALSE
        )

        ## Add the probes to the probelist:
        ENH_raw_probes <- c(
          ENH_raw_probes,
          probelist_ENH_placeholder$V1
        )

      }

    } else{

      ## If no files were found, alert the user and quit the function:
      stop(
        'No ENH datasets were found in step1 folder. \nPlease add .probelist.txt files or run make_external_datasets function to generate them.'
      )
    }

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

    } else{

      ## If no files were found, alert the user and quit the function:
      stop(
        'No NDR datasets were found in step1 folder. \nPlease add .probelist.txt files or run make_external_datasets function to generate them.'
      )
    }

    ## Now get the unique probes from each dataset:
    ## If there aren't any probes, return another error:
    if(length(ENH_raw_probes)>0){

      ENH_unique_probes <- unique(ENH_raw_probes)

    } else{

      ## If no files were found, alert the user and quit the function:
      stop(
        'No ENH probes were identified from probelist.txt files. \nPlease check the files and consider re-running make_external_datasets or add your own custom files.'
      )
    }

    if(length(NDR_raw_probes)>0){

      NDR_unique_probes <- unique(NDR_raw_probes)

    } else{

      ## If no files were found, alert the user and quit the function:
      stop(
        'No NDR probes were identified from probelist.txt files. \nPlease check the files and consider re-running make_external_datasets or add your own custom files.'
      )

    }

    ## Now that the probes have all been loaded, let's get the ones
    ## That are in NDR and ENH regions and aren't near TSS:
    ENH_NDR_probes <- intersect(
      ENH_unique_probes,
      NDR_unique_probes
    )

    probes_of_interest <- intersect(
      ENH_NDR_probes,
      nonTSS_probes
    )

    ## Now let's load the combined methylation/expression dataset:
    ## Which should be contained in the TENET directory:
    LS <- list.files(
      TENET_directory,
      pattern=".rda",
      full.names = TRUE
    )

    ## If the file was found, load it. Otherwise return an error message:
    if(length(LS)==1){

      load(
        LS
      )

    } else if(length(LS)<1){

      stop(
        "No .rda file containing methylation and expression data was found. /nPlease place a .rda files with methylation and expression data in the TENET_directory."
      )

    } else if(length(LS)>1){

      stop(
        "Multiple .rda files found in TENET_directory. /nPlease ensure the .rda file containing methylation/expression data is the only one present."
      )
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
      hg38_450_annotations_no_NA_granges_df$names,
      rownames(metDataN)
    )

    matched_probes_experimental <- intersect(
      hg38_450_annotations_no_NA_granges_df$names,
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
        "No HM450 probes were found in the methylation datasets. /nPlease reload the methylation datasets with HM450 probe IDs in the rownames."
      )

    } else if(length(matched_probes_control)==0){

      stop(
        "No HM450 probes were found in the control (normal) methylation dataset. /nPlease reload the methylation dataset with HM450 probe IDs in the rownames."
      )

    } else if(length(matched_probes_experimental)==0){

      stop(
        "No HM450 probes were found in the experimental (tumor) methylation dataset. /nPlease reload the methylation dataset with HM450 probe IDs in the rownames."
      )

    } else if(length(matched_probes)==0){

      stop(
        "No matched HM450 probes were found between the two methylation datasets. /nPlease reload the methylation datasets with HM450 probe IDs in the rownames."
      )

    }

    ## First get quadrant probes if purity data isn't specified:
    if(use_purity_data==FALSE){
      
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
      ## less than the minExp value.
      ## These represent the methylated probes:
      temp_t1 <- temp_t[
        which(
          length_temp_t_cat<minExp
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
      ## greater than the minExp value.
      ## These represent the hypomethylated probes:
      temp_t1 <- temp_t[
        which(
          length_temp_t_cat>minExp
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
      ## greater than the minExp value.
      ## These represent the unmethylated probes:
      temp_t1 <- temp_t[
        which(
          length_temp_t_cat<minExp
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
      ## greater than the minExp value.
      ## These represent the hypermethylated probes:
      temp_t1 <- temp_t[
        which(
          length_temp_t_cat>minExp
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
      
    } else if(use_purity_data==TRUE){
        
      ## Save a variable with the purity subdirectory:
      purity_subdirectory_path <- paste(
        TENET_directory,
        'purity',
        sep=''
      )
      
      ## Check to make sure the user has supplied a 'purity' subdirectory:
      ## If subdirectory isn't found return an error.
      if(!dir.exists(purity_subdirectory_path)){
        
        stop(
          'The purity subdirectory in the provided TENET_directory was not found, although purity information was specified. Please put .rda files containing methylation data from other cell types for purity analyses in a separate subdirectory in the TENET_directory called "purity"'
        )
        
      }
      
      ## Let's list all the files in the purity subdirectory:
      purity_file_list <- list.files(
        purity_subdirectory_path,
        pattern=".rda",
        full.names = TRUE
      )
      
      purity_file_list_trunc <- list.files(
        purity_subdirectory_path,
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
          hg38_450_annotations_no_NA_granges_df$names,
          rownames(get(purity_rda_content_vector[i]))
        )
        
        ## If no probes are found, return an error message
        if(length(matched_probes_purity_dataset)==0){
          
          stop(
            paste(
              'No HM450 probes were found in the',
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
              'No HM450 probes were found overlapping between the',
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
              'No non-NA HM450 probes were found overlapping between the',
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
        
        ## Get the probes which have fewer than minExp samples and save them 
        ## as the methylated probes:
        purity_methylated_probe_vector[[i]] <- rownames(
          temp_t[
            which(
              length_temp_t_cat<minExp
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
        
        ## Get the probes which have more than minExp samples and save them 
        ## as the hypomethylated probes:
        purity_hypomethylated_probe_vector[[i]] <- rownames(
          temp_t[
            which(
              length_temp_t_cat>minExp
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
        
        ## Get the probes which have less than minExp samples and save them 
        ## as the unmethylated probes:
        purity_unmethylated_probe_vector[[i]] <- rownames(
          temp_t[
            which(
              length_temp_t_cat<minExp
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
        
        ## Get the probes which have more than minExp samples and save them 
        ## as the hypermethylated probes:
        purity_hypermethylated_probe_vector[[i]] <- rownames(
          temp_t[
            which(
              length_temp_t_cat>minExp
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

    ## Now let's save the final dataset:
    save(
      metDataT, metDataN, expDataT, expDataN, enhmetDataT, enhmetDataN, methDataT, methDataN, hypomethDataT, hypomethDataN, unmethDataT, unmethDataN, hypermethDataT, hypermethDataN, hypermethcutoff, hypomethcutoff,
      file= paste(
        TENET_directory,
        'step2/',
        "diff.methylated.datasets.rda",
        sep=''
      )
    )

  } else{

    ## Step1 directory does not exist. Return an error:
    stop(
      'Step1 directory was not found. \nPlease rerun make_external_datasets or create a directory called step1 \nand add your own custom .probelist.txt files to it.'
    )
  }
}

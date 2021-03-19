#' step7_top_genes_survival
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates survival plots and information for the expression level of each gene
#' as well as the DNA methylation of each enhancer probe linked to them,
#' using percentile cutoffs as specified by the user.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with a further subdirectories for each of the four analysis types selected, ending with '_survival' containing the results of this function.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create survival plots for the top genes/TFs by most hypermeth probes with G+ links, and these linked DNA methyation probes if specified.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypermeth probes with G- links, and these linked DNA methyation probes if specified.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G+ links, as well as their linked DNA methyation probes if specified.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create survival plots for the top genes/TFs by most hypometh probes with G- links, as well as their linked DNA methyation probes if specified.
#' @param top_gene_number Specify a number to generate survival plots for that many of the top genes/TFs, as well as their linked enhancer probes if specified, based on the most linked enhancer probes.
#' @param visualize_survival_plots_genes Set to TRUE/FALSE depending on if you want to create .pdfs displaying the survival results for the genes/TFs of interest.
#' @param visualize_survival_plots_probes Set to TRUE/FALSE depending on if you want to create .pdfs displaying the survival results for the probes linked to the genes/TFs of interest.
#' @param high_thresh Set a number ranging from 0 to 1, as a threshold for proportion of samples above that number to include in the high expression/methylation group, and should be greater than or equal to low_thresh to prevent samples from appearing in both groups.
#' @param low_thresh Set a number ranging from 0 to 1, as a threshold for proportion of samples below that number to include in the low expression/methylation group, and should be less than or equal to high_thresh to prevent samples from appearing in both groups.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Returns survival information in the form of .tsv files, as well as .pdfs if selected by the user, showing survival information for the expression of the top gene/TFs, as well as the methylation of the enhancer DNA methylation probes linked to them.
#' @export

step7_top_genes_survival <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  visualize_survival_plots_genes,
  visualize_survival_plots_probes,
  high_thresh,
  low_thresh,
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

  ## Check to make sure the clinical dataset was loaded in the .rda file:
  if(!exists('clinical')){

    stop("clinical data was not loaded from diff.methylated.datasets.rda in step2 of TENET directory. Please either rerun step2 get_diffmeth_regions functions referencing a dataset with clinical data, or check the .rda file to make sure the clinical data is included as an object named 'clinical'.")
  }

  ## Process the expression/methylation/clinical data so the functions later don't have to:

  ## combine datasets:
  expData <- cbind(expDataN, expDataT)
  metData <- cbind(metDataN, metDataT)

  ## Get the list of normal and tumor samples:
  exp_normal_samples <- colnames(expDataN)
  exp_tumor_samples <- colnames(expDataT)

  met_normal_samples <- colnames(metDataN)
  met_tumor_samples <- colnames(metDataT)

  ## Get the relevant columns from the clinical data:
  relevant_clinical <- clinical[
    ,c(
      "bcr_patient_barcode",
      "days_to_death",
      "days_to_last_followup",
      "vital_status"
    )
  ]
  rownames(relevant_clinical) <- relevant_clinical$bcr_patient_barcode

  # If there are subjects that are alive in the
  # dataset, set their days to death info to '-Inf'
  relevant_clinical[
    grep(
      "alive",
      relevant_clinical$vital_status,
      ignore.case = TRUE
    ),
    'days_to_death'
  ] <- '-Inf'

  # If there are subjects that are dead in the
  # dataset, set their days to last followup info to '-Inf'
  relevant_clinical[
    grep(
      "dead",
      relevant_clinical$vital_status,
      ignore.case = TRUE
    ),
    'days_to_last_followup'
  ] <- '-Inf'

  ## Remove subjects that remain unaccounted for:
  ## i.e. NA values remain in previous two columns:
  relevant_clinical <- relevant_clinical[
    !is.na(
      relevant_clinical[,"days_to_death"]
    ),
  ]

  relevant_clinical <- relevant_clinical[
    !is.na(
      relevant_clinical[,"days_to_last_followup"]
    ),
  ]

  # Changing the days to last followup info in the cancer clinical data
  # to be numeric values through character values and remove any NAs that are induced
  # as a final check:
  relevant_clinical$days_to_death <- as.numeric(
    as.character(relevant_clinical$days_to_death)
  )

  relevant_clinical <- relevant_clinical[
    !is.na(
      relevant_clinical[,"days_to_last_followup"]
    ),
  ]

  # Changing the days to death info in the cancer clinical data
  # to be numeric values through character values and remove any NAs that are induced
  # as a final check:
  relevant_clinical$days_to_last_followup <- as.numeric(
    as.character(relevant_clinical$days_to_last_followup)
  )

  relevant_clinical <- relevant_clinical[
    !is.na(
      relevant_clinical[,"days_to_death"]
    ),
  ]

  # Add the relevant number as a final "time" variable
  # Equal to days to death for dead individuals, and
  # days to last followup for alive individuals
  relevant_clinical$time <- ifelse(
    relevant_clinical$vital_status=='Alive',
    relevant_clinical$days_to_last_followup,
    ifelse(
      relevant_clinical$vital_status=='Dead',
      relevant_clinical$days_to_death,
      NA
    )
  )

  ## Create a function to get the survival p-value or graph
  ## for each gene of interest:
  expression_survival_function_graph <- function(
    gene_of_interest,
    high_cutoff,
    low_cutoff,
    graph,
    analysis_type,
    gene_or_TF
  ){

    ## Get actual gene name depending on input:
    if(substring(gene_of_interest,1,4)=='ENSG' & nchar(gene_of_interest)==15){

      ## Input is in ENSG, get the gene name:
      gene_ENSG <- gene_of_interest

      gene_name <- gencode_v22_genes[
        gene_ENSG, 'gene_name'
      ]

      ## Get gene ENSG assuming a name is plugged in:
    } else{

      ## Assume what was given was the gene name, get the ENSG:
      gene_name <- gene_of_interest

      gene_ENSG <- rownames(
        gencode_v22_genes[gencode_v22_genes$gene_name==gene_name,]
      )
    }

    ## Get expression values for gene of interest:
    expression_values <- unlist(
      expData[
        gene_ENSG,
      ]
    )
    names(expression_values) <- colnames(expData)

    # Split gene expression values into normal and tumor samples
    tumor_expression_values <- expression_values[exp_tumor_samples]
    normal_expression_values <- expression_values[exp_normal_samples]

    # Changing the sample names to match the subject names:
    names(tumor_expression_values)  <- substr(
      names(tumor_expression_values),
      1,
      12
    )

    names(normal_expression_values)  <- substr(
      names(normal_expression_values),
      1,
      12
    )

    ## Calculate some basic data:
    normal_sample_n <- length(normal_expression_values)
    tumor_sample_n <- length(tumor_expression_values)

    ## Count the number of NA samples:
    NA_normal <- as.numeric(
      unname(
        table(
          is.na(normal_expression_values)
        )['TRUE']
      )
    )

    NA_tumor <- as.numeric(
      unname(
        table(
          is.na(tumor_expression_values)
        )['TRUE']
      )
    )

    if(is.na(NA_normal)){

      NA_normal <- 0

    }

    if(is.na(NA_tumor)){

      NA_tumor <- 0

    }

    ## Calculate mean expression for all samples present:
    mean_tumor_expression <- mean(
      tumor_expression_values,
      na.rm= TRUE
    )

    mean_normal_expression <- mean(
      normal_expression_values,
      na.rm= TRUE
    )

    ## Calculate the number of samples with present clinical data:
    present_clinical_normal_sample_n <- as.numeric(
      unname(
        table(
          names(normal_expression_values) %in% rownames(relevant_clinical)
        )['TRUE']
      )
    )

    present_clinical_tumor_sample_n  <- as.numeric(
      unname(
        table(
          names(tumor_expression_values) %in% rownames(relevant_clinical)
        )['TRUE']
      )
    )

    ## Identify the samples which are present in the clinical data:
    present_normal_samples <- names(normal_expression_values)[
      names(normal_expression_values) %in% rownames(relevant_clinical)
    ]

    present_tumor_samples <- names(tumor_expression_values)[
      names(tumor_expression_values) %in% rownames(relevant_clinical)
    ]


    # Subsetting clinical patient data for individuals with
    # gene expression data in the tumor set:
    function_relevant_clinical <- relevant_clinical[
      present_tumor_samples,
    ]

    ## Calculate quantiles:
    high_cutoff_quantile= quantile(
      tumor_expression_values,
      high_cutoff,
      na.rm= TRUE
    )[1]

    low_cutoff_quantile= quantile(
      tumor_expression_values,
      low_cutoff,
      na.rm= TRUE
    )[1]

    ## Determine if a sample is in the high, low, or intermediate quartiles:
    sample_expression_grouping <- ifelse(
      tumor_expression_values >  high_cutoff_quantile,
      "High",
      ifelse(
        tumor_expression_values <= low_cutoff_quantile,
        "Low",
        ifelse(
          is.na(tumor_expression_values),
          'Missing',
          'Intermediate'
        )
      )
    )

    ## Add the expression and relevant clinical information to the relevant clinical dataframe:
    function_relevant_clinical$tumor_expression <- tumor_expression_values[
      rownames(function_relevant_clinical)
    ]

    ## Add the expression grouping to the clinical data:
    function_relevant_clinical$expression_grouping <- sample_expression_grouping[
      rownames(function_relevant_clinical)
    ]

    ## Remove samples missing expression data:
    function_relevant_clinical_complete <- na.omit(
      function_relevant_clinical,
      cols=c(
        "tumor_expression",
        "expression_grouping"
      )
    )

    function_relevant_clinical_complete <- function_relevant_clinical

    ## Count the number of samples of each group remaining:
    present_tumor_sample_high_n <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='High',
      ]
    )

    present_tumor_sample_intermediate_n  <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='Intermediate',
      ]
    )

    present_tumor_sample_low_n  <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='Low',
      ]
    )

    present_tumor_sample_missing_n <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='Missing',
      ]
    )

    ## Get the mean expression in each group:
    present_tumor_sample_high_mean <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='High',
        'tumor_expression'
      ],
      na.rm= TRUE
    )

    present_tumor_sample_intermediate_mean  <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='Intermediate',
        'tumor_expression'
      ],
      na.rm= TRUE
    )

    present_tumor_sample_low_mean  <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$expression_grouping=='Low',
        'tumor_expression'
      ],
      na.rm= TRUE
    )

    ## Get clinical results for only samples in the high and low group with info:
    function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
      !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
    ]

    ## Calculate proportion of deceased individuals in each group
    proportion_dead_high <- as.numeric(
      nrow(
        function_relevant_clinical_complete_high_low[
          (function_relevant_clinical_complete_high_low$expression_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
        ]
      ) /
        nrow(
          function_relevant_clinical_complete_high_low[
            (function_relevant_clinical_complete_high_low$expression_grouping=='High'),
          ]
        )
    )

    proportion_dead_low <- as.numeric(
      nrow(
        function_relevant_clinical_complete_high_low[
          (function_relevant_clinical_complete_high_low$expression_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
        ]
      ) /
        nrow(
          function_relevant_clinical_complete_high_low[
            (function_relevant_clinical_complete_high_low$expression_grouping=='Low'),
          ]
        )
    )

    ## Determine which group had greater proportion of fatalities:
    highest_dead_proportion_group <- ifelse(
      proportion_dead_high > proportion_dead_low,
      'high_expression_low_survival',
      ifelse(
        proportion_dead_high < proportion_dead_low,
        'low_expression_low_survival',
        'unclear'
      )
    )

    ## Create a high_low survival object:
    high_low_survival_object <- survival::Surv(
      function_relevant_clinical_complete_high_low$time,
      ifelse(
        function_relevant_clinical_complete_high_low$vital_status=='Alive',
        FALSE,
        TRUE
      )
    )

    # Set the rownames of the survival object to be the
    # patient names for the high/low group smaples:
    rownames(high_low_survival_object) <- rownames(
      function_relevant_clinical_complete[
        !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
      ]
    )

    # Combine the top and down vectors into a single vector
    # in the order top-down
    expression_group <- function_relevant_clinical_complete[
      !(function_relevant_clinical_complete$expression_grouping=='Intermediate'  | function_relevant_clinical_complete$expression_grouping=='Missing'),
      'expression_grouping'
    ]

    # Creating legend names for high/low expression groups:
    legend_name_high <- paste(gene_name,"high")

    legend_name_low <- paste(gene_name,"low")

    # Perform the survival analysis.
    # This uses the combined vector as the x variable
    # and the survival data as the y variable to create
    # a table with information about the test
    # including chi-squared p-value
    survival_table <- survival::survdiff(
      high_low_survival_object  ~ expression_group
    )

    # Get the chi-squared test statistic from the analysis above:
    survival_table_chi_squared <- survival_table$chisq

    # Calculating a p value based on the test statistic to get
    # a precise p-value
    survival_pvalue <- as.numeric(
      1 - pchisq(
        abs(survival_table_chi_squared), df = 1
      )
    )

    ## Round the p-value on the graph to 3 digits:
    survival_pvalue_formatted <- formatC(
      survival_pvalue,
      format='e',
      digits=3
    )

    ## Determine where to save the various results:

    ## Set the main directory to save the file to:
    if(analysis_type=='hypermeth_Gplus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypermeth_Gminus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypometh_Gplus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypometh_Gminus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        sep=''
      )
    }

    ## Set the subdirectory to save the files to:
    if(gene_or_TF=='gene'){

      subdirectory <- 'top_genes/'

    } else if(gene_or_TF=='TF'){

      subdirectory <- 'top_TFs/'
    }

    ## Create the graph if specified:
    if(graph==TRUE){

      # Create the plot title
      # with gene name and p-value included:
      survival_title <- paste(
        gene_name,
        "\nKaplan-Meier Survival analysis\nP = ",
        survival_pvalue_formatted,
        sep=''
      )

      ## Create a title for the survival plot pdf:
      survival_plot_pdf_title <- paste(
        main_directory,
        subdirectory,
        gene_name,
        '_survival_plot.pdf',
        sep=''
      )

      ## Open a pdf for saving the plot:
      pdf(survival_plot_pdf_title)

      ## Now actually create the survival plot
      ## Using a similar structure to what was used to generate the p-value:
      plot(
        survival::survfit(
          high_low_survival_object ~ expression_group
        ),

        # Color the lines (high in red first!)
        col = c('red', 'black'),

        ## Add thickness to the lines:
        lwd=3,

        # Use the title that was created earlier as the title
        # of the plot:
        main= survival_title,

        # Set titles of the x and y axis:
        # Note: TCGA measures survival in days as noted:
        xlab="Time (days)",
        ylab="Probability",

        # Change axis size:
        cex.axis=1,
        cex.main=1.5,
        cex.lab=1.25
      )

      ## Add a legend to the plot:
      legend(

        # Set X position of legend in graph:
        x= (
          max(function_relevant_clinical_complete_high_low$days_to_last_followup)*(2/3)
        ),

        # Set Y position of legend in graph
        y= 1,

        # Use the legend titles that were created earlier
        legend= c(legend_name_high, legend_name_low),

        # As above, use black for low and red for high:
        col= c('red', 'black'),

        # Coloring the text in the legend as well
        text.col= c('red', 'black'),

        # Change the shape of the labels in the legend:
        pch= 15
      )

      ## Close the plot:
      dev.off()

    } else if(graph==FALSE){

      ## create vector of relevant info:
      relevant_vector <- c(
        as.numeric(normal_sample_n),
        as.numeric(tumor_sample_n),
        as.numeric(NA_normal),
        as.numeric(NA_tumor),
        as.numeric(mean_normal_expression),
        as.numeric(mean_tumor_expression),
        as.numeric(normal_sample_n),
        as.numeric(tumor_sample_n),
        as.numeric(present_tumor_sample_missing_n),
        as.numeric(present_tumor_sample_low_n),
        as.numeric(present_tumor_sample_intermediate_n),
        as.numeric(present_tumor_sample_high_n),
        as.numeric(present_tumor_sample_low_mean),
        as.numeric(present_tumor_sample_intermediate_mean),
        as.numeric(present_tumor_sample_high_mean),
        as.numeric(proportion_dead_low),
        as.numeric(proportion_dead_high),
        highest_dead_proportion_group,
        as.numeric(survival_pvalue)
      )
      names(relevant_vector) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_expression',
        'mean_tumor_expression',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_expression',
        'mean_tumor_with_clinical_intermediate_expression',
        'mean_tumor_with_clinical_high_expression',
        'proportion_dead_in_low_expression',
        'proportion_dead_in_high_expression',
        'survival_direction_of_effect',
        'survival_p_value'
      )

      ## Return the vector
      return(relevant_vector)
    }

  }

  ## Create a function to get the survival p-value or graph for
  ## each probe of interest:
  methylation_survival_function_graph <- function(
    probe_of_interest,
    high_cutoff,
    low_cutoff,
    graph,
    analysis_type,
    gene_or_TF
  ){

    ## Get methylation values for probe of interest:
    methylation_values <- unlist(
      metData[
        probe_of_interest,
      ]
    )
    names(methylation_values) <- colnames(metData)

    # Split methylation values into normal and tumor samples
    # If both are present in the dataset
    tumor_methylation_values <- methylation_values[met_tumor_samples]
    normal_methylation_values <- methylation_values[met_normal_samples]

    # Changing the sample names to match the subject names:
    names(tumor_methylation_values)  <- substr(
      names(tumor_methylation_values),
      1,
      12
    )

    names(normal_methylation_values)  <- substr(
      names(normal_methylation_values),
      1,
      12
    )

    ## Calculate some basic data:
    normal_sample_n <- length(normal_methylation_values)
    tumor_sample_n <- length(tumor_methylation_values)

    ## Count the number of NA samples:
    NA_normal <- as.numeric(
      unname(
        table(
          is.na(normal_methylation_values)
        )['TRUE']
      )
    )

    NA_tumor <- as.numeric(
      unname(
        table(
          is.na(tumor_methylation_values)
        )['TRUE']
      )
    )

    if(is.na(NA_normal)){

      NA_normal <- 0

    }

    if(is.na(NA_tumor)){

      NA_tumor <- 0

    }

    ## Calculate mean methylation for all samples present:
    mean_tumor_methylation <- mean(
      tumor_methylation_values,
      na.rm= TRUE
    )

    mean_normal_methylation <- mean(
      normal_methylation_values,
      na.rm= TRUE
    )

    ## Calculate the number of samples with present clinical data:
    present_clinical_normal_sample_n <- as.numeric(
      unname(
        table(
          names(normal_methylation_values) %in% rownames(relevant_clinical)
        )['TRUE']
      )
    )

    present_clinical_tumor_sample_n  <- as.numeric(
      unname(
        table(
          names(tumor_methylation_values) %in% rownames(relevant_clinical)
        )['TRUE']
      )
    )

    ## Identify the samples which are present in the clinical data:
    present_normal_samples <- names(normal_methylation_values)[
      names(normal_methylation_values) %in% rownames(relevant_clinical)
    ]

    present_tumor_samples <- names(tumor_methylation_values)[
      names(tumor_methylation_values) %in% rownames(relevant_clinical)
    ]


    # Subsetting clinical patient data for individuals with
    # DNA methylation data in the tumor set:
    function_relevant_clinical <- relevant_clinical[
      present_tumor_samples,
    ]

    ## Calculate quantiles:
    high_cutoff_quantile= quantile(
      tumor_methylation_values,
      high_cutoff,
      na.rm= TRUE
    )[1]

    low_cutoff_quantile= quantile(
      tumor_methylation_values,
      low_cutoff,
      na.rm= TRUE
    )[1]

    ## Determine if a sample is in the high, low, or intermediate quartiles:
    sample_methylation_grouping <- ifelse(
      tumor_methylation_values >  high_cutoff_quantile,
      "High",
      ifelse(
        tumor_methylation_values <= low_cutoff_quantile,
        "Low",
        ifelse(
          is.na(tumor_methylation_values),
          'Missing',
          'Intermediate'
        )
      )
    )

    ## Add the methylation and relevant clinical information to the relevant clinical dataframe:
    function_relevant_clinical$tumor_methylation <- tumor_methylation_values[
      rownames(function_relevant_clinical)
    ]

    ## Add the methylation grouping to the clinical data:
    function_relevant_clinical$methylation_grouping <- sample_methylation_grouping[
      rownames(function_relevant_clinical)
    ]

    ## Remove samples missing methylation data:
    function_relevant_clinical_complete <- na.omit(
      function_relevant_clinical,
      cols=c(
        "tumor_methylation",
        "methylation_grouping"
      )
    )
    function_relevant_clinical_complete <- function_relevant_clinical

    ## Count the number of samples of each group remaining:
    present_tumor_sample_high_n <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='High',
      ]
    )

    present_tumor_sample_intermediate_n  <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='Intermediate',
      ]
    )

    present_tumor_sample_low_n  <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='Low',
      ]
    )

    present_tumor_sample_missing_n <- nrow(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='Missing',
      ]
    )

    ## Get the mean methylation in each group:
    present_tumor_sample_high_mean <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='High',
        'tumor_methylation'
      ],
      na.rm= TRUE
    )

    present_tumor_sample_intermediate_mean  <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='Intermediate',
        'tumor_methylation'
      ],
      na.rm= TRUE
    )

    present_tumor_sample_low_mean  <- mean(
      function_relevant_clinical_complete[
        function_relevant_clinical_complete$methylation_grouping=='Low',
        'tumor_methylation'
      ],
      na.rm= TRUE
    )

    ## Get clinical results for only samples in the high and low group with info:
    function_relevant_clinical_complete_high_low <- function_relevant_clinical_complete[
      !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
    ]

    ## Calculate proportion of deceased individuals in each group
    proportion_dead_high <- as.numeric(
      nrow(
        function_relevant_clinical_complete_high_low[
          (function_relevant_clinical_complete_high_low$methylation_grouping=='High' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
        ]
      ) /
        nrow(
          function_relevant_clinical_complete_high_low[
            (function_relevant_clinical_complete_high_low$methylation_grouping=='High'),
          ]
        )
    )

    proportion_dead_low <- as.numeric(
      nrow(
        function_relevant_clinical_complete_high_low[
          (function_relevant_clinical_complete_high_low$methylation_grouping=='Low' & function_relevant_clinical_complete_high_low$vital_status=='Dead'),
        ]
      ) /
        nrow(
          function_relevant_clinical_complete_high_low[
            (function_relevant_clinical_complete_high_low$methylation_grouping=='Low'),
          ]
        )
    )

    ## Determine which group had greater proportion of fatalities:
    highest_dead_proportion_group <- ifelse(
      proportion_dead_high > proportion_dead_low,
      'high_methylation_low_survival',
      ifelse(
        proportion_dead_high < proportion_dead_low,
        'low_methylation_low_survival',
        'unclear'
      )
    )

    ## Create a high_low survival object:
    high_low_survival_object <- survival::Surv(
      function_relevant_clinical_complete_high_low$time,
      ifelse(
        function_relevant_clinical_complete_high_low$vital_status=='Alive',
        FALSE,
        TRUE
      )
    )

    # Set the rownames of the survival object to be the
    # patient names for the high/low group smaples:
    rownames(high_low_survival_object) <- rownames(
      function_relevant_clinical_complete[
        !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
      ]
    )

    # Combine the top and down vectors into a single vector
    # in the order top-down
    methylation_group <- function_relevant_clinical_complete[
      !(function_relevant_clinical_complete$methylation_grouping=='Intermediate'  | function_relevant_clinical_complete$methylation_grouping=='Missing'),
      'methylation_grouping'
    ]

    # Creating legend names for high/low methylation groups:
    legend_name_high <- paste(probe_of_interest,"high")

    legend_name_low <- paste(probe_of_interest,"low")

    # Perform the survival analysis.
    # This uses the combined vector as the x variable
    # and the survival data as the y variable to create
    # a table with information about the test
    # including chi-squared p-value
    survival_table <- survival::survdiff(
      high_low_survival_object  ~ methylation_group
    )

    # Get the chi-squared test statistic from the analysis above:
    survival_table_chi_squared <- survival_table$chisq

    # Calculating a p value based on the test statistic to get
    # a precise p-value
    survival_pvalue <- as.numeric(
      1 - pchisq(
        abs(survival_table_chi_squared), df = 1
      )
    )

    ## Round the p-value on the graph to 3 digits:
    survival_pvalue_formatted <- formatC(
      survival_pvalue,
      format='e',
      digits=3
    )

    ## Determine where to save the various results:

    ## Set the main directory to save the file to:
    if(analysis_type=='hypermeth_Gplus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypermeth_Gminus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypometh_Gplus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        sep=''
      )

    } else if(analysis_type=='hypometh_Gminus'){

      main_directory <- paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        sep=''
      )
    }

    ## Set the subdirectory to save the files to:
    if(gene_or_TF=='gene'){

      subdirectory <- 'top_genes/'

    } else if(gene_or_TF=='TF'){

      subdirectory <- 'top_TFs/'
    }

    ## Create the graph if specified:
    if(graph==TRUE){

      # Create the plot title
      # with gene name and p-value included:
      survival_title <- paste(
        probe_of_interest,
        "\nKaplan-Meier Survival analysis\nP = ",
        survival_pvalue_formatted,
        sep=''
      )

      ## Create a title for the survival plot pdf:
      ## This is a comination of the probe name with the linked gene:
      survival_plot_pdf_title <- paste(
        main_directory,
        subdirectory,
        probe_of_interest,
        '_survival_plot.pdf',
        sep=''
      )

      ## Open a pdf for saving the plot:
      pdf(survival_plot_pdf_title)

      ## Now actually create the survival plot
      ## Using a similar structure to what was used to generate the p-value:
      plot(
        survival::survfit(
          high_low_survival_object ~ methylation_group
        ),

        # Color the lines (high in red first!)
        col = c('red', 'black'),

        ## Add thickness to the lines:
        lwd=3,

        # Use the title that was created earlier as the title
        # of the plot:
        main= survival_title,

        # Set titles of the x and y axis:
        # Note: TCGA measures survival in days as noted:
        xlab="Time (days)",
        ylab="Probability",

        # Change axis size:
        cex.axis=1,
        cex.main=1.5,
        cex.lab=1.25
      )

      ## Add a legend to the plot:
      legend(

        # Set X position of legend in graph:
        x= (
          max(function_relevant_clinical_complete_high_low$days_to_last_followup)*(2/3)
        ),

        # Set Y position of legend in graph
        y= 1,

        # Use the legend titles that were created earlier
        legend= c(legend_name_high, legend_name_low),

        # As above, use black for low and red for high:
        col= c('red', 'black'),

        # Coloring the text in the legend as well
        text.col= c('red', 'black'),

        # Change the shape of the labels in the legend:
        pch= 15
      )

      ## Close the plot:
      dev.off()

    } else if(graph==FALSE){

      ## create vector of relevant info:
      relevant_vector <- c(
        as.numeric(normal_sample_n),
        as.numeric(tumor_sample_n),
        as.numeric(NA_normal),
        as.numeric(NA_tumor),
        as.numeric(mean_normal_methylation),
        as.numeric(mean_tumor_methylation),
        as.numeric(normal_sample_n),
        as.numeric(tumor_sample_n),
        as.numeric(present_tumor_sample_missing_n),
        as.numeric(present_tumor_sample_low_n),
        as.numeric(present_tumor_sample_intermediate_n),
        as.numeric(present_tumor_sample_high_n),
        as.numeric(present_tumor_sample_low_mean),
        as.numeric(present_tumor_sample_intermediate_mean),
        as.numeric(present_tumor_sample_high_mean),
        as.numeric(proportion_dead_low),
        as.numeric(proportion_dead_high),
        highest_dead_proportion_group,
        as.numeric(survival_pvalue)
      )
      names(relevant_vector) <- c(
        'normal_sample_count',
        'tumor_sample_count',
        'normal_sample_count_missing',
        'tumor_sample_count_missing',
        'mean_normal_methylation',
        'mean_tumor_methylation',
        'normal_sample_with_clinical_count',
        'tumor_sample_with_clinical_count',
        'tumor_sample_with_clinical_NA_count',
        'tumor_sample_with_clinical_low_count',
        'tumor_sample_with_clinical_intermediate_count',
        'tumor_sample_with_clinical_high_count',
        'mean_tumor_with_clinical_low_methylation',
        'mean_tumor_with_clinical_intermediate_methylation',
        'mean_tumor_with_clinical_high_methylation',
        'proportion_dead_in_low_methylation',
        'proportion_dead_in_high_methylation',
        'survival_direction_of_effect',
        'survival_p_value'
      )

      ## Return the vector
      return(relevant_vector)
    }
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
          'hyper_Gplus_survival',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_survival',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hyper_Gplus_survival
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_survival/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_survival/',
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
          'hyper_Gplus_survival/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gplus_survival/',
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

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hyper_Gplus_all_gene_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gplus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hyper_Gplus_all_TF_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gplus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }

    ## Get survival info for each set of genes:
    ## Get survival information for each gene:
    hyper_gene_survival_results_list <- parallel::mclapply(
      X= top_hyper_Gplus_all_gene_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    hyper_TF_survival_results_list <- parallel::mclapply(
      X= top_hyper_Gplus_all_TF_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into a data frames and set the column names
    ## to be the names that should be generated by the function:
    hyper_gene_survival_results_df <- data.frame(
      matrix(
        unlist(hyper_gene_survival_results_list),
        nrow=length(hyper_gene_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    hyper_TF_survival_results_df <- data.frame(
      matrix(
        unlist(hyper_TF_survival_results_list),
        nrow=length(hyper_TF_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(hyper_gene_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(hyper_TF_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the gene ENSG and names to the data frames:
    hyper_gene_survival_results_df$gene_ENSG <- top_hyper_Gplus_all_gene_ENSG

    hyper_gene_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hyper_Gplus_all_gene_ENSG,
        'gene_name'
      ]
    )

    hyper_TF_survival_results_df$gene_ENSG <- top_hyper_Gplus_all_TF_ENSG

    hyper_TF_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hyper_Gplus_all_TF_ENSG,
        'gene_name'
      ]
    )

    ## Write out the files:
    write.table(
      hyper_gene_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        'top_genes/',
        'hyper_Gplus_top_genes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      hyper_TF_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        'top_TFs/',
        'hyper_Gplus_top_TFs_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes/TFs:
    top_genes_CpGs_linked <- hyper_Gplus_sig_link_zscores[
      hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_gene_ENSG,
      'probeID'
    ]

    top_TFs_CpGs_linked <- hyper_Gplus_sig_link_zscores[
      hyper_Gplus_sig_link_zscores$geneID %in% top_hyper_Gplus_all_TF_ENSG,
      'probeID'
    ]

    ## Get the unique probes from each list:
    unique_top_genes_CpGs_linked <- unique(top_genes_CpGs_linked)

    unique_top_TFs_CpGs_linked <- unique(top_TFs_CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    unique_top_genes_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_genes_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    unique_top_TFs_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_TFs_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into data frames and set the columns
    ## according to the function specifications:
    unique_top_genes_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_genes_CpGs_linked_survival_results_list),
        nrow=length(unique_top_genes_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    unique_top_TFs_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_TFs_CpGs_linked_survival_results_list),
        nrow=length(unique_top_TFs_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(unique_top_genes_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(unique_top_TFs_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the names of the probes to the dataframes:
    unique_top_genes_CpGs_linked_survival_results_df$probe_name <- unique_top_genes_CpGs_linked

    unique_top_TFs_CpGs_linked_survival_results_df$probe_name <- unique_top_TFs_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(
      probe_of_interest,
      gene_or_TF,
      ID_or_name
    ){

      all_listed_genes <- hyper_Gplus_sig_link_zscores[
        hyper_Gplus_sig_link_zscores$probeID==probe_of_interest,
        'geneID'
      ]

      ## Identify what the top genes or TFs are depending on the
      ## analysis you want:
      ## Set the subdirectory to save the files to:
      if(gene_or_TF=='gene'){

        top_gene_IDs <- top_hyper_Gplus_all_gene_ENSG

      } else if(gene_or_TF=='TF'){

        top_gene_IDs <- top_hyper_Gplus_all_TF_ENSG
      }

      ## Get the IDs of the top genes that were linked to this probe:
      top_gene_IDs_linked_to_probe <- top_gene_IDs[
        top_gene_IDs %in% all_listed_genes
      ]

      ## Get the names of the probes:
      top_gene_names_linked_to_probe <- gencode_v22_genes[
        top_gene_IDs_linked_to_probe,
        'gene_name'
      ]

      ## Return the listed genes:
      if(ID_or_name=='ID'){

        return(
          paste(
            top_gene_IDs_linked_to_probe,
            collapse=','
          )
        )

      } else if(ID_or_name=='name'){

        return(
          paste(
            top_gene_names_linked_to_probe,
            collapse=','
          )
        )
      }
    }

    ## Add the listed top genes to the survival results dataframes
    ## for the probes:
    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'ID'
      )
    )

    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'name'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'ID'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'name'
      )
    )

    ## Write out the files:
    write.table(
      unique_top_genes_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        'top_genes/',
        'hyper_Gplus_top_genes_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      unique_top_TFs_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gplus_survival/',
        'top_TFs/',
        'hyper_Gplus_top_TFs_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now, generate plots for these probes if the user has specified it:

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gplus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gplus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }
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
          'hyper_Gminus_survival',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_survival',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hyper_Gminus_survival
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_survival/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_survival/',
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
          'hyper_Gminus_survival/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hyper_Gminus_survival/',
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

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hyper_Gminus_all_gene_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gminus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hyper_Gminus_all_TF_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gminus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }

    ## Get survival info for each set of genes:
    ## Get survival information for each gene:
    hyper_gene_survival_results_list <- parallel::mclapply(
      X= top_hyper_Gminus_all_gene_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    hyper_TF_survival_results_list <- parallel::mclapply(
      X= top_hyper_Gminus_all_TF_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into a data frames and set the column names
    ## to be the names that should be generated by the function:
    hyper_gene_survival_results_df <- data.frame(
      matrix(
        unlist(hyper_gene_survival_results_list),
        nrow=length(hyper_gene_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    hyper_TF_survival_results_df <- data.frame(
      matrix(
        unlist(hyper_TF_survival_results_list),
        nrow=length(hyper_TF_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(hyper_gene_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(hyper_TF_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the gene ENSG and names to the data frames:
    hyper_gene_survival_results_df$gene_ENSG <- top_hyper_Gminus_all_gene_ENSG

    hyper_gene_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hyper_Gminus_all_gene_ENSG,
        'gene_name'
      ]
    )

    hyper_TF_survival_results_df$gene_ENSG <- top_hyper_Gminus_all_TF_ENSG

    hyper_TF_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hyper_Gminus_all_TF_ENSG,
        'gene_name'
      ]
    )

    ## Write out the files:
    write.table(
      hyper_gene_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        'top_genes/',
        'hyper_Gminus_top_genes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      hyper_TF_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        'top_TFs/',
        'hyper_Gminus_top_TFs_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes/TFs:
    top_genes_CpGs_linked <- hyper_Gminus_sig_link_zscores[
      hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_gene_ENSG,
      'probeID'
    ]

    top_TFs_CpGs_linked <- hyper_Gminus_sig_link_zscores[
      hyper_Gminus_sig_link_zscores$geneID %in% top_hyper_Gminus_all_TF_ENSG,
      'probeID'
    ]

    ## Get the unique probes from each list:
    unique_top_genes_CpGs_linked <- unique(top_genes_CpGs_linked)

    unique_top_TFs_CpGs_linked <- unique(top_TFs_CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    unique_top_genes_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_genes_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    unique_top_TFs_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_TFs_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypermeth_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into data frames and set the columns
    ## according to the function specifications:
    unique_top_genes_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_genes_CpGs_linked_survival_results_list),
        nrow=length(unique_top_genes_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    unique_top_TFs_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_TFs_CpGs_linked_survival_results_list),
        nrow=length(unique_top_TFs_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(unique_top_genes_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(unique_top_TFs_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the names of the probes to the dataframes:
    unique_top_genes_CpGs_linked_survival_results_df$probe_name <- unique_top_genes_CpGs_linked

    unique_top_TFs_CpGs_linked_survival_results_df$probe_name <- unique_top_TFs_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(
      probe_of_interest,
      gene_or_TF,
      ID_or_name
    ){

      all_listed_genes <- hyper_Gminus_sig_link_zscores[
        hyper_Gminus_sig_link_zscores$probeID==probe_of_interest,
        'geneID'
      ]

      ## Identify what the top genes or TFs are depending on the
      ## analysis you want:
      ## Set the subdirectory to save the files to:
      if(gene_or_TF=='gene'){

        top_gene_IDs <- top_hyper_Gminus_all_gene_ENSG

      } else if(gene_or_TF=='TF'){

        top_gene_IDs <- top_hyper_Gminus_all_TF_ENSG
      }

      ## Get the IDs of the top genes that were linked to this probe:
      top_gene_IDs_linked_to_probe <- top_gene_IDs[
        top_gene_IDs %in% all_listed_genes
      ]

      ## Get the names of the probes:
      top_gene_names_linked_to_probe <- gencode_v22_genes[
        top_gene_IDs_linked_to_probe,
        'gene_name'
      ]

      ## Return the listed genes:
      if(ID_or_name=='ID'){

        return(
          paste(
            top_gene_IDs_linked_to_probe,
            collapse=','
          )
        )

      } else if(ID_or_name=='name'){

        return(
          paste(
            top_gene_names_linked_to_probe,
            collapse=','
          )
        )
      }
    }

    ## Add the listed top genes to the survival results dataframes
    ## for the probes:
    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'ID'
      )
    )

    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'name'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'ID'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'name'
      )
    )

    ## Write out the files:
    write.table(
      unique_top_genes_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        'top_genes/',
        'hyper_Gminus_top_genes_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      unique_top_TFs_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hyper_Gminus_survival/',
        'top_TFs/',
        'hyper_Gminus_top_TFs_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now, generate plots for these probes if the user has specified it:

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gminus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypermeth_Gminus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }
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
          'hypo_Gplus_survival',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_survival',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hypo_Gplus_survival
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_survival/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_survival/',
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
          'hypo_Gplus_survival/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gplus_survival/',
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

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hypo_Gplus_all_gene_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gplus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hypo_Gplus_all_TF_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gplus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }

    ## Get survival info for each set of genes:
    ## Get survival information for each gene:
    hypo_gene_survival_results_list <- parallel::mclapply(
      X= top_hypo_Gplus_all_gene_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    hypo_TF_survival_results_list <- parallel::mclapply(
      X= top_hypo_Gplus_all_TF_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into a data frames and set the column names
    ## to be the names that should be generated by the function:
    hypo_gene_survival_results_df <- data.frame(
      matrix(
        unlist(hypo_gene_survival_results_list),
        nrow=length(hypo_gene_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    hypo_TF_survival_results_df <- data.frame(
      matrix(
        unlist(hypo_TF_survival_results_list),
        nrow=length(hypo_TF_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(hypo_gene_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(hypo_TF_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the gene ENSG and names to the data frames:
    hypo_gene_survival_results_df$gene_ENSG <- top_hypo_Gplus_all_gene_ENSG

    hypo_gene_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hypo_Gplus_all_gene_ENSG,
        'gene_name'
      ]
    )

    hypo_TF_survival_results_df$gene_ENSG <- top_hypo_Gplus_all_TF_ENSG

    hypo_TF_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hypo_Gplus_all_TF_ENSG,
        'gene_name'
      ]
    )

    ## Write out the files:
    write.table(
      hypo_gene_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        'top_genes/',
        'hypo_Gplus_top_genes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      hypo_TF_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        'top_TFs/',
        'hypo_Gplus_top_TFs_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes/TFs:
    top_genes_CpGs_linked <- hypo_Gplus_sig_link_zscores[
      hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_gene_ENSG,
      'probeID'
    ]

    top_TFs_CpGs_linked <- hypo_Gplus_sig_link_zscores[
      hypo_Gplus_sig_link_zscores$geneID %in% top_hypo_Gplus_all_TF_ENSG,
      'probeID'
    ]

    ## Get the unique probes from each list:
    unique_top_genes_CpGs_linked <- unique(top_genes_CpGs_linked)

    unique_top_TFs_CpGs_linked <- unique(top_TFs_CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    unique_top_genes_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_genes_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    unique_top_TFs_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_TFs_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gplus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into data frames and set the columns
    ## according to the function specifications:
    unique_top_genes_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_genes_CpGs_linked_survival_results_list),
        nrow=length(unique_top_genes_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    unique_top_TFs_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_TFs_CpGs_linked_survival_results_list),
        nrow=length(unique_top_TFs_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(unique_top_genes_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(unique_top_TFs_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the names of the probes to the dataframes:
    unique_top_genes_CpGs_linked_survival_results_df$probe_name <- unique_top_genes_CpGs_linked

    unique_top_TFs_CpGs_linked_survival_results_df$probe_name <- unique_top_TFs_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(
      probe_of_interest,
      gene_or_TF,
      ID_or_name
    ){

      all_listed_genes <- hypo_Gplus_sig_link_zscores[
        hypo_Gplus_sig_link_zscores$probeID==probe_of_interest,
        'geneID'
      ]

      ## Identify what the top genes or TFs are depending on the
      ## analysis you want:
      ## Set the subdirectory to save the files to:
      if(gene_or_TF=='gene'){

        top_gene_IDs <- top_hypo_Gplus_all_gene_ENSG

      } else if(gene_or_TF=='TF'){

        top_gene_IDs <- top_hypo_Gplus_all_TF_ENSG
      }

      ## Get the IDs of the top genes that were linked to this probe:
      top_gene_IDs_linked_to_probe <- top_gene_IDs[
        top_gene_IDs %in% all_listed_genes
      ]

      ## Get the names of the probes:
      top_gene_names_linked_to_probe <- gencode_v22_genes[
        top_gene_IDs_linked_to_probe,
        'gene_name'
      ]

      ## Return the listed genes:
      if(ID_or_name=='ID'){

        return(
          paste(
            top_gene_IDs_linked_to_probe,
            collapse=','
          )
        )

      } else if(ID_or_name=='name'){

        return(
          paste(
            top_gene_names_linked_to_probe,
            collapse=','
          )
        )
      }
    }

    ## Add the listed top genes to the survival results dataframes
    ## for the probes:
    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'ID'
      )
    )

    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'name'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'ID'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'name'
      )
    )

    ## Write out the files:
    write.table(
      unique_top_genes_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        'top_genes/',
        'hypo_Gplus_top_genes_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      unique_top_TFs_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gplus_survival/',
        'top_TFs/',
        'hypo_Gplus_top_TFs_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now, generate plots for these probes if the user has specified it:

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gplus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gplus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }
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
          'hypo_Gminus_survival',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_survival',
          sep=''
        )
      )
    }

    ## Create a subdirectories in the new hypo_Gminus_survival
    ## directory to hold results from the top genes and TFs:
    if(
      !dir.exists(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_survival/',
          'top_genes',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_survival/',
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
          'hypo_Gminus_survival/',
          'top_TFs',
          sep=''
        )
      )
    ){

      dir.create(
        paste(
          TENET_directory,
          'step7/',
          'hypo_Gminus_survival/',
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

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_genes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hypo_Gminus_all_gene_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gminus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= top_hypo_Gminus_all_TF_ENSG,
        FUN= expression_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gminus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }

    ## Get survival info for each set of genes:
    ## Get survival information for each gene:
    hypo_gene_survival_results_list <- parallel::mclapply(
      X= top_hypo_Gminus_all_gene_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    hypo_TF_survival_results_list <- parallel::mclapply(
      X= top_hypo_Gminus_all_TF_ENSG,
      FUN= expression_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into a data frames and set the column names
    ## to be the names that should be generated by the function:
    hypo_gene_survival_results_df <- data.frame(
      matrix(
        unlist(hypo_gene_survival_results_list),
        nrow=length(hypo_gene_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    hypo_TF_survival_results_df <- data.frame(
      matrix(
        unlist(hypo_TF_survival_results_list),
        nrow=length(hypo_TF_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(hypo_gene_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(hypo_TF_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_expression',
      'mean_tumor_expression',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_expression',
      'mean_tumor_with_clinical_intermediate_expression',
      'mean_tumor_with_clinical_high_expression',
      'proportion_dead_in_low_expression',
      'proportion_dead_in_high_expression',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the gene ENSG and names to the data frames:
    hypo_gene_survival_results_df$gene_ENSG <- top_hypo_Gminus_all_gene_ENSG

    hypo_gene_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hypo_Gminus_all_gene_ENSG,
        'gene_name'
      ]
    )

    hypo_TF_survival_results_df$gene_ENSG <- top_hypo_Gminus_all_TF_ENSG

    hypo_TF_survival_results_df$gene_name <- c(
      gencode_v22_genes[
        top_hypo_Gminus_all_TF_ENSG,
        'gene_name'
      ]
    )

    ## Write out the files:
    write.table(
      hypo_gene_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        'top_genes/',
        'hypo_Gminus_top_genes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      hypo_TF_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        'top_TFs/',
        'hypo_Gminus_top_TFs_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now for each of the probes, we will need to calculate survival info and graphs:

    ## Get the list of probes linked to the significant genes/TFs:
    top_genes_CpGs_linked <- hypo_Gminus_sig_link_zscores[
      hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_gene_ENSG,
      'probeID'
    ]

    top_TFs_CpGs_linked <- hypo_Gminus_sig_link_zscores[
      hypo_Gminus_sig_link_zscores$geneID %in% top_hypo_Gminus_all_TF_ENSG,
      'probeID'
    ]

    ## Get the unique probes from each list:
    unique_top_genes_CpGs_linked <- unique(top_genes_CpGs_linked)

    unique_top_TFs_CpGs_linked <- unique(top_TFs_CpGs_linked)

    ## Get the survival info for each probe by mclapplying the
    ## methylation survival function to all of the probes
    ## linked to each gene:
    unique_top_genes_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_genes_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    unique_top_TFs_CpGs_linked_survival_results_list <- parallel::mclapply(
      X= unique_top_TFs_CpGs_linked,
      FUN= methylation_survival_function_graph,
      high_cutoff= high_thresh,
      low_cutoff= low_thresh,
      graph= FALSE,
      analysis_type='hypometh_Gminus',
      gene_or_TF='gene',
      mc.cores= core_count
    )

    ## Transform the lists into data frames and set the columns
    ## according to the function specifications:
    unique_top_genes_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_genes_CpGs_linked_survival_results_list),
        nrow=length(unique_top_genes_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    unique_top_TFs_CpGs_linked_survival_results_df <- data.frame(
      matrix(
        unlist(unique_top_TFs_CpGs_linked_survival_results_list),
        nrow=length(unique_top_TFs_CpGs_linked_survival_results_list),
        byrow=T
      ),
      stringsAsFactors = FALSE
    )

    colnames(unique_top_genes_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    colnames(unique_top_TFs_CpGs_linked_survival_results_df) <- c(
      'normal_sample_count',
      'tumor_sample_count',
      'normal_sample_count_missing',
      'tumor_sample_count_missing',
      'mean_normal_methylation',
      'mean_tumor_methylation',
      'normal_sample_with_clinical_count',
      'tumor_sample_with_clinical_count',
      'tumor_sample_with_clinical_NA_count',
      'tumor_sample_with_clinical_low_count',
      'tumor_sample_with_clinical_intermediate_count',
      'tumor_sample_with_clinical_high_count',
      'mean_tumor_with_clinical_low_methylation',
      'mean_tumor_with_clinical_intermediate_methylation',
      'mean_tumor_with_clinical_high_methylation',
      'proportion_dead_in_low_methylation',
      'proportion_dead_in_high_methylation',
      'survival_direction_of_effect',
      'survival_p_value'
    )

    ## Add the names of the probes to the dataframes:
    unique_top_genes_CpGs_linked_survival_results_df$probe_name <- unique_top_genes_CpGs_linked

    unique_top_TFs_CpGs_linked_survival_results_df$probe_name <- unique_top_TFs_CpGs_linked

    ## For each of the CpGs, list which of the top X genes it was listed to:

    ## Wite a function to get the top X gene names a probe is assigned to:
    top_gene_assignment_function <- function(
      probe_of_interest,
      gene_or_TF,
      ID_or_name
    ){

      all_listed_genes <- hypo_Gminus_sig_link_zscores[
        hypo_Gminus_sig_link_zscores$probeID==probe_of_interest,
        'geneID'
      ]

      ## Identify what the top genes or TFs are depending on the
      ## analysis you want:
      ## Set the subdirectory to save the files to:
      if(gene_or_TF=='gene'){

        top_gene_IDs <- top_hypo_Gminus_all_gene_ENSG

      } else if(gene_or_TF=='TF'){

        top_gene_IDs <- top_hypo_Gminus_all_TF_ENSG
      }

      ## Get the IDs of the top genes that were linked to this probe:
      top_gene_IDs_linked_to_probe <- top_gene_IDs[
        top_gene_IDs %in% all_listed_genes
      ]

      ## Get the names of the probes:
      top_gene_names_linked_to_probe <- gencode_v22_genes[
        top_gene_IDs_linked_to_probe,
        'gene_name'
      ]

      ## Return the listed genes:
      if(ID_or_name=='ID'){

        return(
          paste(
            top_gene_IDs_linked_to_probe,
            collapse=','
          )
        )

      } else if(ID_or_name=='name'){

        return(
          paste(
            top_gene_names_linked_to_probe,
            collapse=','
          )
        )
      }
    }

    ## Add the listed top genes to the survival results dataframes
    ## for the probes:
    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'ID'
      )
    )

    unique_top_genes_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'gene',
        ID_or_name= 'name'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_ID <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'ID'
      )
    )

    unique_top_TFs_CpGs_linked_survival_results_df$top_genes_linked_name <- unname(
      sapply(
        unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        top_gene_assignment_function,
        gene_or_TF= 'TF',
        ID_or_name= 'name'
      )
    )

    ## Write out the files:
    write.table(
      unique_top_genes_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        'top_genes/',
        'hypo_Gminus_top_genes_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    write.table(
      unique_top_TFs_CpGs_linked_survival_results_df,
      file= paste(
        TENET_directory,
        'step7/',
        'hypo_Gminus_survival/',
        'top_TFs/',
        'hypo_Gminus_top_TFs_linked_probes_survival_info.tsv',
        sep=''
      ),
      row.names = FALSE,
      quote= FALSE,
      sep='\t'
    )

    ## Now, generate plots for these probes if the user has specified it:

    ## If visualize_survival_plots_genes is TRUE, generate the plots:
    if(visualize_survival_plots_probes==TRUE){

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_genes_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gminus',
        gene_or_TF='gene',
        mc.cores= core_count
      )

      ## Create survival plots for each of the genes designated:
      parallel::mclapply(
        X= unique_top_TFs_CpGs_linked_survival_results_df$probe_name,
        FUN= methylation_survival_function_graph,
        high_cutoff= high_thresh,
        low_cutoff= low_thresh,
        graph= TRUE,
        analysis_type='hypometh_Gminus',
        gene_or_TF='TF',
        mc.cores= core_count
      )
    }
  }
}

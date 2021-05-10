#' step7_top_genes_met_heatmaps
#'
#' This is a step7 function of the TENETR package.
#' This function takes the top genes/TFs by number of linked probes identified from
#' the step6_probe_per_gene_tabulation function up to the number as specified by the user
#' and generates heatmaps showing the methylation level of the enhancer
#' DNA methylation probes linked to those genes of the specified analysis types.
#'
#' @param TENET_directory Set a path to the TENET directory containing the 'step6' subdirectory and results created by the step6_probe_per_gene_tabulation function. This function will also create a new 'step7' subdirectory there, if not already created, with further subdirectories for each of the four analysis types selected, ending with '_met_heatmaps' containing the results of this function.
#' @param hypermeth_Gplus_analysis Set to TRUE/FALSE depending on if you want to create heatmaps showing DNA methylation levels of enhancer probes linked to the top genes/TFs with the most hypermeth probes with G+ links.
#' @param hypermeth_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of enhancer probes linked the top genes/TFs with the most hypermeth probes with G- links.
#' @param hypometh_Gplus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of enhancer probes linked the top genes/TFs with the most hypometh probes with G+ links.
#' @param hypometh_Gminus_analysis Set to TRUE/FALSE depending on if you want to to create heatmaps showing DNA methylation levels of enhancer probes linked the top genes/TFs with the most hypometh probes with G- links.
#' @param top_gene_number Specify a number of the top genes/TFs based on the most linked enhancer probes of a given analysis type to generate heatmaps for with their linked enhancer probes methylation levels.
#' @param core_count Argument passed as mc.cores argument for mclapply. See ?mclapply from the parallel package for more details.
#' @return Currently returns .pdf files with the heatmaps showing the DNA methylation levels of all probes of a specific analysis type linked to those genes, as well as the expression of said genes in the column labels
#' @export

step7_top_genes_met_heatmaps <- function(
  TENET_directory,
  hypermeth_Gplus_analysis,
  hypermeth_Gminus_analysis,
  hypometh_Gplus_analysis,
  hypometh_Gminus_analysis,
  top_gene_number,
  core_count
){

  ## Load latest version of heatmap.3 function
  ## This is from: "https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R"
  ## Acquired 3/17/2021
  heatmap.3 <- function(x,
                        Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                        distfun = dist,
                        hclustfun = hclust,
                        dendrogram = c("both","row", "column", "none"),
                        symm = FALSE,
                        scale = c("none","row", "column"),
                        na.rm = TRUE,
                        revC = identical(Colv,"Rowv"),
                        add.expr,
                        breaks,
                        symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                        col = "heat.colors",
                        colsep,
                        rowsep,
                        sepcolor = "white",
                        sepwidth = c(0.05, 0.05),
                        cellnote,
                        notecex = 1,
                        notecol = "cyan",
                        na.color = par("bg"),
                        trace = c("none", "column","row", "both"),
                        tracecol = "cyan",
                        hline = median(breaks),
                        vline = median(breaks),
                        linecol = tracecol,
                        margins = c(5,5),
                        ColSideColors,
                        RowSideColors,
                        side.height.fraction=0.3,
                        cexRow = 0.2 + 1/log10(nr),
                        cexCol = 0.2 + 1/log10(nc),
                        labRow = NULL,
                        labCol = NULL,
                        key = TRUE,
                        keysize = 1.5,
                        density.info = c("none", "histogram", "density"),
                        denscol = tracecol,
                        symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                        densadj = 0.25,
                        main = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        lmat = NULL,
                        lhei = NULL,
                        lwid = NULL,
                        ColSideColorsSize = 1,
                        RowSideColorsSize = 1,
                        KeyValueName="Value",...){

    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
        return(TRUE)
      if (is.list(x))
        return(all(sapply(x, invalid)))
      else if (is.vector(x))
        return(all(is.na(x)))
      else return(FALSE)
    }

    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
      col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
      warning("Using scale=\"row\" or scale=\"column\" when breaks are",
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
      Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
      Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                   c("both", "row"))) {
        if (is.logical(Colv) && (Colv))
          dendrogram <- "column"
        else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                   c("both", "column"))) {
        if (is.logical(Rowv) && (Rowv))
          dendrogram <- "row"
        else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      ddr <- reorder(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
      rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      ddc <- reorder(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
      colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
      labRow <- if (is.null(rownames(x)))
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
      labCol <- if (is.null(colnames(x)))
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
      if (missing(col) || is.function(col))
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)

      if (!missing(ColSideColors)) {
        #if (!is.matrix(ColSideColors))
        #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
          stop("'ColSideColors' must be a matrix of nrow(x) rows")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
        lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
      }

      if (!missing(RowSideColors)) {
        #if (!is.matrix(RowSideColors))
        #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
          stop("'RowSideColors' must be a matrix of ncol(x) columns")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat))
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors)){
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      } else {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = t(RowSideColors[,rowInd, drop=F])
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
          rsc.colors[rsc.i] = rsc.name
          rsc[rsc == rsc.name] = rsc.i
          rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(rownames(RowSideColors)) > 0) {
          axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
        }
      }
    }

    if (!missing(ColSideColors)) {

      if (!is.matrix(ColSideColors)){
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      } else {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop=F]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
          csc.colors[csc.i] = csc.name
          csc[csc == csc.name] = csc.i
          csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
          axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
        }
      }
    }

    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr"))
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
      retval$rowDendrogram <- ddr
    if (exists("ddc"))
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
      eval(substitute(add.expr))
    if (!missing(colsep))
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in colInd) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol,
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in rowInd) {
        if (!is.null(hline)) {
          abline(h = i + hline, col = linecol, lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote))
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
           col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      par(mar = c(5, 4, 2, 1), cex = 0.75)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min(x, na.rm = TRUE)
        max.raw <- max(x, na.rm = TRUE)
      }

      z <- seq(min.raw, max.raw, length = length(col))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(1, at = xv, labels = lv)
      if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
      else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
      else mtext(side = 1, KeyValueName, line = 2)
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
              lwd = 1)
        axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
        title("Color Key\nand Density Plot")
        par(cex = 0.5)
        mtext(side = 2, "Density", line = 2)
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
              col = denscol)
        axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
        title("Color Key\nand Histogram")
        par(cex = 0.5)
        mtext(side = 2, "Count", line = 2)
      }
      else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
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
      stop('hyper_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hyper_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hyper_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hyper_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hypo_Gplus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hypo_Gplus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hypo_Gminus_links_all_gene_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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
      stop('hypo_Gminus_links_all_TF_freq.txt in step6 of TENET directory was not found. Please check that the file exists and consider rerunning the step6_probe_per_gene_tabulation.')

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

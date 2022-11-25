# Script to be launched from the command line to build the logos and make a logo
# report after running MoDec.

# This script needs 3 command line arguments:
#   1) The name of the sample
#   2) The name of the folder where the results have been saved (and where logos
#     will be made).
#   3) The max y value for the y-axis (i.e. logos will be drawn with the y-axis
#     between 0 and this y-max have common axes limits for the various plots;
#     note however that if the logos go beyond this value a new y-max will
#     be set for the given logo to still see it fully). Can also use a value
#     of 0 to use different values for each logo based on its data.
#   4) The color scheme to use for the logos. (default to "auto").

# Getting input arguments -------------------------------------------------
args <- commandArgs(trailingOnly=T)
sample <- args[1]
resFolder <- args[2]
ymax <- as.numeric(args[3])
if (length(args) >= 4){
  col_scheme <- args[4]
} else {
  col_scheme <- "auto"
}

# -+-+-+-+-+-+
# Indicate below the path where the custom R libraries are installed
.libPaths("PATH/TO/R/LIBRARIES")
# .libPaths("~/RLibraries")

library(ggplot2)

cat("Drawing logos for", sample, "from the folder\n ", resFolder, "\n")

# Defining the function making individual logos based on my visual preferences ----------------------
#' Creating a sequence logo with help of ggseqlogo.
#'
#' This is some intermediate function calling ggseqlogo (the forked version found
#' at https://github.com/GfellerLab/ggseqlogo), and some default parameters
#' have been changed to adapt the visual output.
#'
#' @param data A list of sequences or a position weight matrix.
#' @param smallSampleCorr Include small-sample correction in information content or not.
#' @param col_scheme Color scheme of plot.
#' @param legendText Plot legend or not.
#' @param ylim Range of y-axis.
#' @param titlePlot Title of plot.
#' @param titleSize Size of title.
#' @param titlePos Horizontal position of title.
#' @param axisTextSizeX Size of x tick labels.
#' @param axisTextSizeY Size of y tick labels.
#' @param axisTitleSize Size of axis title.
#' @param seq_type Sequence type of input.
#' @param font Font for plot.
#'
#' @import ggplot2
#'
#' @export
seq_logo <- function( data, smallSampleCorr = T, col_scheme = 'auto',
  legendText = FALSE, ylim = c(0,4.32),
  title = NULL, titleSize = 18, titlePos = 0.5,
  axisTextSizeX = 12, axisTextSizeY = 12, axisTitleSize = 14,
  seq_type = 'aa', font = 'helvetica_bold', ...){

  if (typeof(data) == "character") {
    lengthP = nchar(data[[1]])
  }
  else if (typeof(data) == "double") {
    lengthP = ncol(data)
  }

  #### plot sequence logo
  p = ggplot() +

    #### do the plotting of the sequence logo
    ggseqlogo::geom_logo(data=data, col_scheme=col_scheme,
      smallSampleCorr=smallSampleCorr, font=font, legendText=legendText, ...) +

    #### size of title, size of x and y tick marks, size of y axis description
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5),
      axis.ticks = element_line(colour = 'black', size = 0.5),
      plot.title = element_text(size=titleSize, hjust=titlePos, vjust = 0.1, family="sans"),
      axis.text.x = element_text(size=axisTextSizeX, color = 'black', family="sans"),
      axis.text.y = element_text(size=axisTextSizeY, color = 'black', family="sans"),
      axis.title.y = element_text(size=axisTitleSize, family="sans")) +

    #### remove space between axis and plot space, define x ticks
    scale_y_continuous(expand = c(0, 0)) +

    #### costum y axis range
    coord_cartesian(ylim = ylim) +

    #### add title if given
    ggtitle(title)

  p <- suppressMessages(p +
      scale_x_continuous(breaks = seq(1, lengthP, 1), expand=c(0,0)))

  return(p)
}

# Running the code to make the logo report. -------------------------------
report <- paste0(resFolder, sample, "_report.html")
if (file.exists(report)){
  warning("The html report for ", sample, " already existed, it'll be ",
    "replaced by the new one.")
}

subFolder <- "Logos/"
subFolder_fullPath <- paste0(resFolder, "/Logos/")
if (dir.exists(subFolder_fullPath)){
  prevLogos <- list.files(path=subFolder_fullPath, pattern=".png", full.names=T)
  if (length(prevLogos) > 0){
    file.remove(prevLogos)
    warning("Some logos already existed in `", subFolder_fullPath, "`. Deleted ",
      "all of them and created new logos.")
  }
} else
  dir.create(subFolder_fullPath, recursive = T)
if (is.na(ymax) || (ymax == 0)){
  ylim <- NULL
} else {
  ylim <- c(0, ymax)
}

c_dpi <- 300
if (grepl("Darwin", Sys.info()["sysname"])){
  c_dpi <- 100
  # On mac the resolution is better already with less dpi... so keep it
  # like that to make files of not too big size (it seems this info isn't
  # totally correctly used there...).
}

# First creating all the logos with ggseqlogo
PWMfiles <- gsub(".txt", "", list.files(path=paste0(resFolder, "PWMs"),
  pattern="PWM_"))
n_pep_all <- list()
for (cPWM in PWMfiles){
  cPWM_fullPath <- paste0(paste0(resFolder, "PWMs/"), cPWM, ".txt")
  n_pep <- scan(file=cPWM_fullPath, what=list(character(), character(), numeric()),
    nlines=1, quiet=T)[[3]]
  n_pep_all[[cPWM]] <- n_pep
  figFile <- paste0(subFolder_fullPath, cPWM, "-", round(n_pep), ".png")
  tPWM <- as.matrix(read.delim(cPWM_fullPath, skip=1, row.names=1))
  tplot <- seq_logo(tPWM, ylim=ylim, col_scheme=col_scheme)
  if ((! (is.na(ymax) || ymax==0) )  &&
      (max(ggplot2::ggplot_build(tplot)$data[[1]]$y) > ymax))
    tplot <- suppressMessages(tplot + ggplot2::coord_cartesian(expand=F))
  if (ncol(tPWM) < 30){
  ggplot2::ggsave(figFile, width=6, height=4, units="in", dpi=c_dpi)
  } else {
    # For very long logos (i.e. big L_motif), we use other figure dimensions to
    # make the letters more visible and we rotate the x-axis labels.
    tplot <- tplot + ggplot2::theme(axis.text.x=element_text(angle=90, hjust=1))
    ggplot2::ggsave(figFile, width=ncol(tPWM)/5, height=4, units="in", dpi=c_dpi,
      limitsize=F)
    # limitsize is used so that ggplot doesn't complain for big figure dimensions.
  }
}

# ########'
# Now reading some additional run parameters and making the html report.
tHeadStr <- paste("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"",
  "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">",
  "<html>\n<head>\n<style type='text/css'>",
  "div.float {\n\tfloat: left;\n\tmargin-right: 30px;\n\tmargin-left: 30px;\n}",
  "div.float p {\n\ttext-align: center;\n}",
  "div.spacer {\n\tclear: both;\n}", "</style>\n</head>\n<body>\n", sep="\n")

n_motifs <- sort(unique(as.numeric(gsub("PWM_K(\\d*)_.*", "\\1", PWMfiles))))

if (any(grepl("_run", PWMfiles))){
  report_multiRuns <- paste0(resFolder, sample, "_report_multiRuns.html")
  with_report_multiRuns <- TRUE
} else
  with_report_multiRuns <- FALSE


scoreFile <- paste0(resFolder, "Summary/res_K", n_motifs[1], "_out.txt")
# The same parameters are used in principle for all runs, so we'll read the
# values from the first file (but if --outAdd was used, it's possible that
# different parameter values were used for runs in a same folder; we don't
# check it here).
if (file.exists(scoreFile)){
  tStr_pars <- "<h2>Parameters of the deconvolution</h2>\n<p>"

  temp <- stringr::str_split(readLines(scoreFile), pattern="\t", simplify=T)
  methodName <- temp[1, 2]
  param_remove <- c("n_motifs")
  # This value is different in each summary file and is already visible from
  # the motifs showed so we won't put it in the param summary header.
  temp <- temp[! (temp[,1]  %in% param_remove),]
  params_lastRow <- which(temp[,1] == "dLogL_conv")

  tStr_pars <- paste0(tStr_pars,
    paste(apply(temp[2:params_lastRow,], MARGIN=1, FUN=function(x){
      paste(x[x!=""], collapse=": ")}),  collapse="<br>\n"), "</p>\n")
} else {
  warning("Couldn't find a summary score file to determine the run parameters used.")
  methodName <- "MoDec"
  tStr_pars <- ""
}

tStr <- paste0(tHeadStr, "<h1>", methodName, " - best results for ", sample,
    "</h1>\n", tStr_pars)
cat(tStr, file=report)

if (with_report_multiRuns)
  cat(gsub("best results for",
    "results from multiple runs with highest log-likelihood for",
    tStr), file=report_multiRuns)

nScoreFiles <- list.files(path=paste0(resFolder, "Summary/"),
  pattern="nMotScores", full.names=T)
if (length(nScoreFiles) > 0){
  tStr <- "<h2>Information content table</h2>\n"
  cat(tStr, file=report, append=T)
  tData <- lapply(nScoreFiles, FUN=function(cfile){
    nonEmpty <- length(readLines(cfile, n=2)) > 1
    if (nonEmpty){
      temp <- read.delim(cfile, as.is=T)
    } else
      temp <- NULL
    return(temp)
  })
  tData <- Reduce(rbind, tData)
  colKept <- c("n", "logL", "AIC")
  tData <- tData[,colKept] # Keep only the most useful ones.
  best_n <- which.min(tData[,"AIC"]) # Determine best num of motifs
  # Make a table of the logL and AIC.
  tData_r <- htmlTable::txtRound(as.matrix(tData), digits=2, excl.cols="^n$")
  col.rgroup <- rep_len(c("none", "gray90"), nrow(tData_r))
  col.rgroup[best_n] <- "seagreen1"
  tStr <- paste("Best number of motifs (based on AIC):", best_n)
  tTable <- htmlTable::htmlTable(tData_r, rnames=F, align="r",
    css.cell="padding: 3px 20px", col.rgroup=col.rgroup, tfoot=tStr)
  cat(tTable, file=report, append=T)
}
for (k in n_motifs){
  tStr <- paste0("<h2>", k, " motifs", "</h2>\n")
  cat(tStr, file=report, append=T)
  if (with_report_multiRuns)
    cat(tStr, file=report_multiRuns, append=T)
  clust_ids <- sort(as.numeric(na.omit(stringr::str_match(PWMfiles,
    paste0("PWM_K", k, "_([[:digit:]]*)"))[,2])))
  scoreFile <- paste0(resFolder, "Summary/res_K", k, "_out.txt")
  if (file.exists(scoreFile)){
    temp <- stringr::str_split(readLines(scoreFile), pattern="\t", simplify=T)
    n_pep_tot <- as.numeric(temp[grep("^n_peptides$", temp[,1]), 2])
    n_runs <- as.numeric(temp[grep("^n_runs_real$", temp[,1]), 2])
    best_score <- as.numeric(temp[grep("^logL_bestRun$", temp[,1]), 2])
    run_scores <- as.numeric(temp[grep("^logL_eachRun", temp[,1]), 1+(1:n_runs)])
  } else {
    warning("The output summary file '", scoreFile, "' wasn't found.")
    n_pep_tot <- sum(c(n_pep_all[paste0("res_K", k, "_PWM_",
      clust_ids)], recursive=T))
    best_score <- NA
    run_scores <- NA
  }
  if (with_report_multiRuns){
    tStr_multi <- paste0("<h4>Best run (logL: ", best_score, ")</h4>\n")
    cat(tStr_multi, file=report_multiRuns, append=T)
  }
  for (i in c(clust_ids, 0)){
    if (i == 0){
      cPWM <- paste0("PWM_K", k, "_flat")
    } else
      cPWM <- paste0("PWM_K", k, "_", i)
    if (!(cPWM %in% names(n_pep_all))){
      if (i != 0){
        warning("The logo relating to ", cPWM, "doesn't seem to exist.")
      }
      next
      # In case we're not using a flat motif.
    }
    cLogoFile <- paste0(cPWM, "-", round(n_pep_all[[cPWM]]), ".png")
    tStr <- paste0("<div class=\"float\"><a href=\"", subFolder, cLogoFile,
      "\" ", "target=\"_blank\"><img src=\"", subFolder, cLogoFile, "\"",
      " height=\"150\"></a><p>", ifelse(i==0, "flat motif", i),
      " - ", # width=\"200\"
      n_pep_all[[cPWM]], " (",
      round(n_pep_all[[cPWM]]/n_pep_tot, digits=2), ")</p></div>\n")
    cat(tStr, file=report, append=T)
    if (with_report_multiRuns)
      cat(tStr, file=report_multiRuns, append=T)
  }
  cat("<div class='spacer'></div>\n", file=report, append=T)
  if (with_report_multiRuns)
    cat("<div class='spacer'></div>\n", file=report_multiRuns, append=T)

  # To make subgroups for each run (the results from the best run have
  # already been saved above).
  if (with_report_multiRuns){
    run_ids <- sort(unique(as.numeric(
      stringr::str_match(PWMfiles, paste0("PWM_K", k,
        "_run([[:digit:]]*)_"))[,2])))
    for (cRun in run_ids){
      tStr <- paste0("<h4>Run ", cRun, " (logL: ", run_scores[cRun], ")</h4>\n")
      cat(tStr, file=report_multiRuns, append=T)

      for (i in c(clust_ids, 0)){
        if (i == 0){
          cPWM <- paste0("PWM_K", k, "_run", cRun, "_flat")
        } else
          cPWM <- paste0("PWM_K", k, "_run", cRun, "_", i)
        if (!(cPWM %in% names(n_pep_all))){
          if (i != 0){
            warning("The logo relating to ", cPWM, "doesn't seem to exist.")
          }
          next
          # In case we're not using a flat motif.
        }
        cLogoFile <- paste0(cPWM, "-", round(n_pep_all[[cPWM]]), ".png")
        cat("<div class=\"float\"><a href=\"", subFolder, cLogoFile, "\" ",
          "target=\"_blank\"><img src=\"", subFolder, cLogoFile, "\"",
          " height=\"150\"></a><p>", ifelse(i==0, "flat motif", i),
          " - ", # width=\"200\"
          n_pep_all[[cPWM]], " (",
          round(n_pep_all[[cPWM]]/n_pep_tot, digits=2), ")</p></div>\n",
          sep="", file=report_multiRuns, append=T)
      }
      cat("<div class='spacer'></div>\n", file=report_multiRuns, append=T)
    }
  }
}
tStr <- "<div class='spacer'></div>\n</body>\n</html>"
cat(tStr, file=report, append=T)
if (with_report_multiRuns)
  cat(tStr, file=report_multiRuns, append=T)



# If there were warnings, we output them for the logs ---------------
warn_msg <- unique(names(warnings()))
if (length(warn_msg) > 0){
  cat("Last warning messages were:", paste("\n -(", 1:length(warn_msg), "): ",
    warn_msg, sep="", collapse=""), "\n")
}


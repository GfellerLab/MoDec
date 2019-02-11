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

# Getting input arguments -------------------------------------------------
args <- commandArgs(trailingOnly=T)
sample <- args[1]
resFolder <- args[2]
ymax <- as.numeric(args[3])

# -+-+-+-+-+-+
# Indicate below the path where the custom R libraries are installed
.libPaths("PATH/TO/R/LIBRARIES")
# .libPaths("~/RLibraries")
# .libPaths("C:/Users/jracle/OneDrive/Documents/R/win-library/3.5")

library(ggplot2)

methodName <- "MoDec v1.0"

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
seq_logo <- function( data, smallSampleCorr = T, col_scheme = 'chemistry',
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
  tplot <- seq_logo(tPWM, ylim=ylim)
  if ((! (is.na(ymax) || ymax==0) )  &&
      (max(ggplot2::ggplot_build(tplot)$data[[1]]$y) > ymax))
    tplot <- suppressMessages(tplot + ggplot2::coord_cartesian(expand=F))
  ggplot2::ggsave(figFile, width=6, height=4, units="in", dpi=c_dpi)
}

# ########'
# Now making the html report.
n_motifs <- sort(unique(as.numeric(gsub("PWM_K(\\d*)_.*", "\\1", PWMfiles))))

if (any(grepl("_run", PWMfiles))){
  report_multiRuns <- paste0(resFolder, sample, "_report_multiRuns.html")
  with_report_multiRuns <- TRUE
} else
  with_report_multiRuns <- FALSE
tStr <- paste("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"",
  "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">",
  "<html>\n<head>\n<style type='text/css'>",
  "div.float {\n\tfloat: left;\n\tmargin-right: 30px;\n\tmargin-left: 30px;\n}",
  "div.float p {\n\ttext-align: center;\n}",
  "div.spacer {\n\tclear: both;\n}", "</style>\n</head>\n<body>",
  paste0("<h1>", methodName, " - best results for ", sample,
    "</h1>"),
  sep="\n")

cat(tStr, file=report)
if (with_report_multiRuns)
  cat(gsub("best results for",
    "results from multiple runs with highest log-likelihood for",
    tStr), file=report_multiRuns)

nScoreFiles <- list.files(path=paste0(resFolder, "Summary/"),
  pattern="nMotScores", full.names=T)
if (length(nScoreFiles) > 0){
  tStr <- "<h2>Information content table"
  tStr <- paste0(tStr, "</h2>\n")
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
  tData <- tData[,colKept]
  # Keep only the most useful ones.
  tData_r <- htmlTable::txtRound(tData, digits=2, excl.cols=1)
  tTable <- htmlTable::htmlTable(tData_r, rnames=F, align="r",
    css.cell="padding: 3px 20px", col.rgroup=c("none", "gray90"))
  cat(tTable, file=report, append=T)
  best_n <- which.min(tData[,"AIC"])
  cat("<h4>Best number of motifs (based on AIC): ", best_n, "</h4>\n",
    file=report, append=T)
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
    n_pep_tot_full <- as.numeric(temp[grep("^n_peptides$", temp[,1]), 2])
    n_pep_tot <- as.numeric(temp[grep("^n_weighted_peptides$", temp[,1]), 2])
    n_runs <- as.numeric(temp[grep("^n_runs$", temp[,1]), 2])
    best_score <- as.numeric(temp[grep("^logL_bestRun$", temp[,1]), 2])
    run_scores <- as.numeric(temp[grep("^logL_eachRun", temp[,1]), 1+(1:n_runs)])
  } else {
    warning("The output summary file '", scoreFile, "' wasn't found.")
    n_pep_tot <- sum(c(n_pep_all[paste0("res_K", k, "_PWM_",
      clust_ids)], recursive=T))
    n_pep_tot_full <- NA
    best_score <- NA
    run_scores <- NA
  }
  run_sc_tab <- table(run_scores, useNA="ifany")
  run_sc_tab <- run_sc_tab[order(as.numeric(names(run_sc_tab)), decreasing=T)]
  tStr <- paste0("<h4>logL_bestRun: ", best_score, " (logL from each run: ",
    paste(ifelse(run_sc_tab>1, paste0(run_sc_tab, "*"), ""),
      names(run_sc_tab), sep="", collapse=", "), ")")
  if (with_report_multiRuns)
    tStr_multi <- paste0("<h4>Best run (logL: ", best_score)
  tStr <- paste0(tStr, "<br>N weighted peptides: ", n_pep_tot, " (",
    n_pep_tot_full, " peptides without weighting)")
  if (with_report_multiRuns)
    tStr_multi <- paste0(tStr_multi, "<br>N weighted peptides: ", n_pep_tot,
      " (", n_pep_tot_full, " peptides without weighting)")
  tStr <- paste0(tStr, "</h4>\n")
  cat(tStr, file=report, append=T)
  if (with_report_multiRuns){
    tStr_multi <- paste0(tStr_multi, ")</h4>\n")
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


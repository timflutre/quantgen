#!/usr/bin/env Rscript

## Aim: use 'qvalue' from the command-line
## Author: Timothee Flutre
## Not copyrighted -- provided to the public domain

suppressPackageStartupMessages(require(qvalue))
prog.name <- "qvalue.R"

##-----------------------------------------------------------------------------
## list of functions

help <- function(){
  txt <- paste0("`", prog.name, "' launches 'qvalue'.\n")
  txt <- paste0(txt, "\nUsage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\nOptions:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --in\tinput file (can be gzipped)\n")
  txt <- paste0(txt, "      --col\tcolumn(s) to keep (e.g. '3' or '1-2-4', all if empty)\n")
  txt <- paste0(txt, "\t\t1 column assumed to be p-value\n")
  txt <- paste0(txt, "\t\t2 columns assumed to be gene and p-value (in this order)\n")
  txt <- paste0(txt, "\t\t3 columns assumed to be gene, snp and p-value (in this order)\n")
  txt <- paste0(txt, "      --head\tif there is a header line\n")
  txt <- paste0(txt, "      --fdr\tthreshold on the FDR (default=0.05)\n")
  txt <- paste0(txt, "      --out\toutput file (optional, gzipped if extension is 'gz')\n")
  txt <- paste0(txt, "      --hist\tfile in which to plot the histogram of p-values (optional)\n")
  message(txt)
}

parseArgs <- function(){
  params <- list(verbose=1,
                 in.file=NULL,
                 columns=NULL,
                 header=FALSE,
                 fdr=0.05,
                 out.file=NULL,
                 hist.file=NULL)
  
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)
  
  i <- 0
  while(i < length(args)){ # use "while" loop for options with no argument
    i <- i + 1
    if(args[i] == "-h" || args[i] == "--help"){
      help()
      quit("no", status=0)
    }
    if(args[i] == "-v" || args[i] == "--verbose"){
      params$verbose <- as.numeric(args[i+1])
      i <- i + 1
    }
    if(args[i] == "--in"){
      params$in.file <- args[i+1]
      i <- i + 1
    }
    if(args[i] == "--col"){
      params$columns <- as.numeric(strsplit(args[i+1], "-")[[1]])
      i <- i + 1
    }
    if(args[i] == "--head"){
      params$header <- TRUE
    }
    if(args[i] == "fdr"){
      params$fdr <- as.numeric(args[i+1])
      i <- i + 1
    }
    if(args[i] == "--out"){
      params$out.file <- args[i+1]
      i <- i + 1
    }
    if(args[i] == "--hist"){
      params$hist.file <- args[i+1]
      i <- i + 1
    }
  }
  
  if(params$verbose > 0){
    message("parameters:")
    print(params)
  }
  stopifnot(! is.null(params$in.file),
            file.exists(params$in.file))
  if(! is.null(params$columns) &&
     ((length(params$columns) != 1) &&
      (length(params$columns) != 2) &&
      (length(params$columns) != 3)))
    stop("--col should specify 1, 2 or 3 columns", call.=FALSE)
  
  return(params)
}

loadData <- function(file, header, columns, verbose=0){
  if(verbose > 0)
    message("load input file ...")
  tmp <- read.table(file=file, header=header, nrows=5)
  colClasses <- sapply(tmp, class)
  data <- read.table(file=file, header=header, colClasses=colClasses,
                     stringsAsFactors=FALSE)
  if(! is.null(columns)){
    if(length(columns) > ncol(data))
      stop("--col specified more columns than there are in the file",
           call.=FALSE)
    if(length(columns) == 1){
      data <- list(pvalue=data[,columns])
    } else if(length(columns) == 2){
      data <- list(gene=data[,columns[1]], pvalue=data[,columns[2]])
    } else
      data <- list(gene=data[,columns[1]], snp=data[,columns[2]],
                   pvalue=data[,columns[3]])
  }
  if(verbose > 0){
    message(paste0("nb of p-values: "), length(data$pvalue))
  }
    return(data)
}

callSignifViaStoreyMethod <- function(data, fdr, verbose=0){
  if(verbose > 0)
    message("call significant tests via Storey's method ...")
  set.seed(1859)
  qobj <- qvalue(p=data$pvalue, fdr.level=fdr, robust="TRUE",
                 pi0.method="bootstrap")
  if(verbose > 0){
    message(paste0("pFDR\t", fdr))
    message(paste0("pi0\t", format(qobj$pi0, digits=7)))
    if(! "gene" %in% names(data) || ! "snp" %in% names(data)){
      message(paste0("nbSignifTests\t", sum(qobj$significant)))
    } else{
      message(paste0("nbSignifGeneSnpPairs\t", sum(qobj$significant)))
      message(paste0("nbGenes\t", length(unique(data$gene[qobj$significant]))))
    }
  }
  return(qobj)
}

saveRes <- function(data, qobj, file, verbose=0){
  if(verbose > 0)
    message("save results ...")
  data.saved <- qobj$qvalues
  if("gene" %in% names(data) && ! "snp" %in% names(data)){
    data.saved <- cbind(data$gene, data.saved)
    colnames(data.saved) <- c("gene", "qvalue")
  } else if("gene" %in% names(data) && "snp" %in% names(data)){
    data.saved <- cbind(data$gene, data$snp, data.saved)
    colnames(data.saved) <- c("gene", "snp", "qvalue")
  }
  extension <- (tmp <- strsplit(file, "\\.")[[1]])[length(tmp)]
  if(extension == "gz"){
    write.table(x=data.saved, file=gzfile(file), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=ifelse(is.null(ncol(data.saved)),
                                   FALSE, TRUE))
  } else
    write.table(x=data.saved, file=file, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=ifelse(is.null(ncol(data.saved)),
                                   FALSE, TRUE))
}

getFdr <- function(pvalues=NULL, threshold=0.01){
  stopifnot(! is.null(pvalues), length(pvalues) > 0)
  res <- list()
  set.seed(1859)
  qobj <- qvalue(p=pvalues, robust="TRUE", pi0.method="bootstrap")
  res$prop.true.null.tests <- qobj$pi0
  nb.tests <- length(pvalues)
  res$nb.signif.tests <- length(pvalues[pvalues <= threshold])
  res$fdr <- (res$prop.true.null.tests * nb.tests * threshold) / res$nb.signif.tests
  return(res)
}

plotHistPvalues <- function(data=NULL,
                            maintitle=NULL,
                            subtitle=NULL,
                            xlab=NULL,
                            ylab=NULL,
                            writePi0=TRUE,
                            ylim=NULL,
                            file=NULL,
                            verbose=0){
  stopifnot(! is.null(data), length(data) > 0)
  if(verbose > 0)
    message("plot the histogram of p-values ...")
  
  if(! is.null(file))
    pdf(file)
  
  if(! is.null(maintitle)){
    if(! is.null(subtitle)){
      par(mar=c(4.5, 4.1, 2, 0))
    } else
      par(mar=c(4.5, 4.1, 1, 0))
  } else
    par(mar=c(4.5, 4.1, 0, 0))
  
  tmp <- hist(data, breaks=99, col="darkgrey", border="white", main="",
              xlab="", ylab="", yaxt="n", ylim=ylim)
  if(! is.null(ylim)){
    custom.at <- seq(from=0, to=ylim[2], length.out=5)
    axis(side=2, at=custom.at, labels=TRUE, las=1)
    step <- (custom.at[2] - custom.at[1]) / 5
  } else{
    axis(side=2, labels=TRUE, las=1)
    step <- tmp$counts[1] / 10
  }
  if(is.null(xlab))
    xlab <- expression(paste(italic(p),"-values"))
  if(is.null(ylab))
    ylab <- paste0("Number of tests (out of ", length(data), ")")
  title(main=maintitle, xlab=xlab, ylab=ylab)
  if (! is.null(subtitle)){
    mtext(subtitle, line=-0.9)
  }
  
  fdr <- getFdr(data)
  pi0.perc <- sprintf("%.2f", 100 * fdr$prop.true.null.tests)
  txt <- bquote(paste("proportion of null tests (", pi[0], "): ",
                      .(pi0.perc), "%"))
  if(writePi0){
    if(is.null(ylim)){
      text(x=0.5, y=floor(tmp$counts[1] - 1 * step), labels=txt)
    } else{
      text(x=0.5, y=custom.at[4] - 1 * step, labels=txt)
    }
  }
  
  fdrs <- c(0.05, 0.1, 0.15)
  i <- 2
  cex.fdr <- 1 # 0.8
  tmp2 <- lapply(fdrs, function(fdr){
    qres <- qvalue(p=data, fdr.level=fdr, robust=TRUE, pi0.method="bootstrap")
    txt <- paste("FDR=", fdr, " : ",
                 sum(qres$significant), " tests", sep="")
    if(is.null(ylim)){
      txt <- sprintf("FDR=%.2f", fdr)
      text(x=0.33, y=floor(tmp$counts[1] - i * step), pos=4, labels=txt, cex=cex.fdr)
      arrows(x0=0.51, y0=1.005 * floor(tmp$counts[1] - i * step),
             x1=0.59, y1=1.005 * floor(tmp$counts[1] - i * step),
             length=0.03, lwd=0.5)
      txt <- paste0(sum(qres$significant), " tests")
      text(x=0.6, y=floor(tmp$counts[1] - i * step), pos=4, labels=txt, cex=cex.fdr)
    } else{
      txt <- sprintf("FDR=%.2f", fdr)
      text(x=0.33, y=custom.at[4] - i * step, pos=4, labels=txt, cex=cex.fdr)
      arrows(x0=0.51, y0=custom.at[4] - i * step,
             x1=0.59, y1=custom.at[4] - i * step,
             length=0.03, lwd=0.5)
      txt <- paste0(sum(qres$significant), " tests")
      text(x=0.6, y=custom.at[4] - i * step, pos=4, labels=txt, cex=cex.fdr)
    }
    i <<- i + 1
    qres
  })
  
  if(! is.null(file))
    dev.off()
  
  res <- list(pi0=fdr$prop.true.null.tests)
  res[["fdr0.05"]] <- tmp2[[1]]
  res[["fdr0.1"]] <- tmp2[[2]]
  res[["fdr0.15"]] <- tmp2[[3]]
  invisible(res)
}

run <- function(fdr=0.05, verbose=1){
  params <- parseArgs()
  if(params$verbose > 0)
    message(paste0("START ", prog.name, " (", date(), ")"))
  data <- loadData(params$in.file, params$header, params$columns,
                   params$verbose)
  qobj <- callSignifViaStoreyMethod(data, params$fdr, params$verbose)
  if(! is.null(params$out.file))
    saveRes(data, qobj, params$out.file, params$verbose)
  if(! is.null(params$hist.file))
    plotHistPvalues(data=data$pvalue, file=params$hist.file,
                    verbose=params$verbose)
  if(params$verbose > 0)
    message(paste0("END ", prog.name, " (", date(), ")"))
}

##-----------------------------------------------------------------------------
## Run the program

system.time(run())
print(object.size(x=lapply(ls(), get)), units="Kb")

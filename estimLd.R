#!/usr/bin/env Rscript

# Aim: estimate pairwise linkage disequilibrium correcting for kinship and/or structure
# Copyright (C) 2015-2016 INRA
# License: GPL-3+
# Persons: Timothée Flutre [cre,aut]
# Versioning: none

rm(list=ls())
prog.name <- "estimLd.R"
prog.version <- "1.1.1" # http://semver.org/

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 2 || (R.v.maj == 2 && R.v.min.1 < 15))
    stop("require R >= 2.15 (for paste0)", call.=FALSE)

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
help <- function(){
  txt <- paste0("`", prog.name, "' estimates pairwise linkage disequilibrium correcting for kinship and/or structure.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Usage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --genos\tpath to the genotype file encoded as 0/1/2 (can be gzipped)\n")
  txt <- paste0(txt, "      --snps\tpath to the file with SNP coordinates (can be gzipped)\n")
  txt <- paste0(txt, "\t\tno header, 1 row per SNP, 3 columns: snp<tab>pos<tab>chr (as in GEMMA)\n")
  txt <- paste0(txt, "      --out\tpath to the output gzipped file\n")
  txt <- paste0(txt, "      --mmaf\tminimum for the minor allele frequency (default=0.01)\n")
  txt <- paste0(txt, "      --ck\tcorrect LD estimates with kinship\n")
  txt <- paste0(txt, "      --cs\tcorrect LD estimates with population structure\n")
  txt <- paste0(txt, "      --magr\tmethod to estimate the additive genomic relationships\n")
  txt <- paste0(txt, "\t\tdefault=vanraden1/astle-balding/habier/yang/zhou\n")
  txt <- paste0(txt, "      --chr\tonly chromosome to analyze\n")
  txt <- paste0(txt, "      --pop\tonly population to analyze\n")
  txt <- paste0(txt, "      --threads\tnumber of threads (default=1)\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Remarks:\n")
  txt <- paste0(txt, "  The LDcorSV and rutilstimflutre packages are required.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Examples:\n")
  txt <- paste0(txt, "  ", prog.name, " --genos genos.txt.gz --snps snp_coords.txt.gz --out ld.txt.gz\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <timothee.flutre@supagro.inra.fr>.")
  write(txt, stdout())
}

##' Display version and license information on stdout
##'
##' The person roles comply with R's guidelines (The R Journal Vol. 4/1, June 2012). To comply with help2man (http://www.gnu.org/s/help2man), use --no-discard-stderr.
##' @title Version
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Copyright (C) 2015-2016 INRA.\n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Timothée Flutre [cre,aut].")
  write(txt, stdout())
}

##' Parse the command-line arguments
##'
##' Allow short and long options
##' @title Command-line
##' @param params list of parameters initialized with default values
##' @return List of parameters
parseCmdLine <- function(params){
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)

  i <- 0
  while(i < length(args)){ # use "while" loop for options with no argument
    i <- i + 1
    if(args[i] == "-h" || args[i] == "--help"){
      help()
      quit("no", status=0)
    }
    else if(args[i] == "-V" || args[i] == "--version"){
      version()
      quit("no", status=0)
    }
    else if(args[i] == "-v" || args[i] == "--verbose"){
      params$verbose <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--genos"){
      params$file.genos <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--snps"){
      params$file.snps <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--out"){
      params$file.out <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--mmaf"){
      params$min.maf <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--ck"){
      params$correct.kinship <- TRUE
    }
    else if(args[i] == "--cs"){
      params$correct.structure <- TRUE
    }
    else if(args[i] == "--magr"){
      params$meth.add.gen.rel <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--chr"){
      params$only.chr <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--pop"){
      params$only.pop <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--threads"){
      params$nb.threads <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--utils"){
      params$utils <- args[i+1]
      i <- i + 1
    }
    else{
      write(paste0(prog.name, ": invalid option -- ", args[i], "\n"), stderr())
      help()
      quit("no", status=1)
    }
  }

  return(params)
}

##' Check the values of the command-line parameters
##'
##' @param params list of parameters
checkParams <- function(params){
  if(is.null(params$file.genos)){
    write("ERROR: missing compulsory option --genos\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$file.genos)){
    write(paste0("ERROR: can't find file ", params$file.genos, "\n"), stderr())
    help()
    quit("no", status=1)
  }
  if(is.null(params$file.snps)){
    write("ERROR: missing compulsory option --snps\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$file.genos)){
    write(paste0("ERROR: can't find file ", params$file.snps, "\n"), stderr())
    help()
    quit("no", status=1)
  }
  if(is.null(params$file.out)){
    write("ERROR: missing compulsory option --out\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! params$meth.add.gen.rel %in% c("zhou", "astle-balding", "vanraden1", "habier", "yang")){
    write(paste0("ERROR: unknown option --magr", params$meth.add.gen.rel, "\n"),
          stderr())
    help()
    quit("no", status=1)
  }

  suppressPackageStartupMessages(library(rutilstimflutre))
  if(any(params$correct.kinship, params$correct.structure))
    suppressPackageStartupMessages(library(LDcorSV))
}

run <- function(params){
  if(params$verbose > 0)
    write(paste0("load SNP genotypes (and discard if MAF < ", params$min.maf,
                 ") ..."), stdout())
  X <- as.matrix(read.table(file=params$file.genos, sep="\t"))
  print(dim(X))
  mafs <- estimSnpMaf(X=X)
  if(any(mafs < params$min.maf)){
    X <- X[, - which(mafs < params$min.maf)]
    stopifnot(ncol(X) > 0)
    print(dim(X))
  }

  if(params$verbose > 0)
    write("load SNP coordinates (and order as columns of X) ...", stdout())
  snp.coords <- read.table(file=params$file.snps, sep="\t",
                           stringsAsFactors=FALSE)
  stopifnot(ncol(snp.coords) == 3,
            nrow(snp.coords) >= ncol(X))
  rownames(snp.coords) <- gsub("-", ".", snp.coords[,1])
  snp.coords[,1] <- NULL
  colnames(snp.coords) <- c("pos", "chr")
  snp.coords <- snp.coords[colnames(X), c("chr", "pos")]
  print(dim(snp.coords))

  A <- NULL
  if(params$correct.kinship){
    if(params$verbose > 0)
      write("estimate the additive genomic relationships ...", stdout())
    A <- estimGenRel(X=X, thresh=0, method=params$meth.add.gen.rel)
  }

  pops <- NULL
  if(params$correct.structure){
    if(params$verbose > 0)
      write("extract the population assignments ...", stdout())
    stopifnot(grepl("_", rownames(X)[1]))
    pops <- sapply(strsplit(rownames(X), "_"), function(x){x[1]})
    names(pops) <- rownames(X)
    print(table(pops))
  }

  ld <- estimLd(X=X, K=A, pops=pops, snp.coords=snp.coords,
                only.chr=params$only.chr, only.pop=params$only.pop,
                use.ldcorsv=ifelse(any(! is.null(A), ! is.null(pops)),
                                   TRUE, FALSE),
                verbose=params$verbose)

  write.table(x=ld, file=gzfile(params$file.out), quote=FALSE, row.names=FALSE,
              sep="\t")
}

main <- function(){
  params <- list(verbose=1,
                 file.genos=NULL,
                 file.snps=NULL,
                 file.out=NULL,
                 min.maf=0.01,
                 correct.kinship=FALSE,
                 correct.structure=FALSE,
                 meth.add.gen.rel="vanraden1",
                 only.chr=NULL,
                 only.pop=NULL,
                 nb.threads=1)

  params <- parseCmdLine(params)

  checkParams(params)

  if(params$verbose > 0){
    start.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write(paste0("START ", prog.name, " ", prog.version, " ", start.time),
          stdout())
    args <- commandArgs(trailingOnly=TRUE)
    write(paste("cmd-line:", prog.name, paste(args, collapse=" ")), stdout())
    write(paste0("cwd: ", getwd()), stdout())
  }

  system.time(run(params))

  if(params$verbose > 0){
    end.time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    difft <- as.numeric(
        difftime(as.POSIXct(end.time, format="%Y-%m-%d %H:%M:%S"),
                 as.POSIXct(start.time, format="%Y-%m-%d %H:%M:%S"),
                 units="days"))
    ## difft <- 1 + 2/24 + 3/(24*60) + 3/(24*3600) # 1d 2h 3m 4s in days
    difft.d <- floor(difft)
    difft.h <- floor(((difft - difft.d) * 24) %% 24)
    difft.m <- floor(((difft - difft.d - difft.h/24) * 24*60) %% (24 * 60))
    difft.s <- floor(((difft - difft.d - difft.h/24 - difft.m/(24*60)) *
                      24*60*60) %% (24 * 60 * 60))
    run.length <- sprintf("%02i:%02i:%02i", difft.h, difft.m, difft.s)
    write(paste0("END ", prog.name, " ", prog.version, " ", end.time,
                 " (", run.length, ")"), stdout())
    ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
  }
}

if(! interactive())
    main()

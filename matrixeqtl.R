#!/usr/bin/env Rscript

## Aim: launch MatrixEQTL
## Author: Timothee Flutre
## Not copyrighted -- provided to the public domain

suppressPackageStartupMessages(require(MatrixEQTL))
prog.name <- "matrixeqtl.R"

##-----------------------------------------------------------------------------
## list of functions

## see also ?Matrix_eQTL_main
help <- function(){
  txt <- paste0("`", prog.name, "' launches MatrixEQTL.\n")
  txt <- paste0(txt, "\nUsage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\nOptions:\n")
  txt <- paste0(txt, "  -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, "  -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, "  -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, "      --genos-file\tfile with genotypes\n")
  txt <- paste0(txt, "      --phenos-file\tfile with phenotypes (columns will be re-ordered if necessary)\n")
  txt <- paste0(txt, "      --cvrt-file\tfile with covariates (columns will be re-ordered if necessary)\n")
  txt <- paste0(txt, "      --out-file\t\n")
  txt <- paste0(txt, "      --pv-threshold\t10^{-5}\n")
  txt <- paste0(txt, "      --model\tdefault=linear\n")
  txt <- paste0(txt, "      --err-cvrt-file\t\n")
  txt <- paste0(txt, "      --cis-out-file\t\n")
  txt <- paste0(txt, "      --cis-pv-threshold\tdefault=0\n")
  txt <- paste0(txt, "      --snpspos-file\tBED file\n")
  txt <- paste0(txt, "      --genepos-file\tBED file\n")
  txt <- paste0(txt, "      --cis-dist\tdefault=1000000\n")
  txt <- paste0(txt, "      --pvalue-hist\t\n")
  txt <- paste0(txt, "      --only-tss\t\n")
  txt <- paste0(txt, "      --fdr-storey\trequires 'qvalue'\n")
  txt <- paste0(txt, "      --perm-pheno-samples\t\n")
  txt <- paste0(txt, "      --seed\tdefault=1859\n")
  message(txt)
}

parseArgs <- function(){
  params <- list(verbose=1,
                 genos.file=NULL,
                 phenos.file=NULL,
                 cvrt.file=NULL,
                 out.file=NULL,
                 pv.threshold=10^(-5),
                 model="linear",
                 err.cvrt.file=NULL,
                 cis.out.file=NULL,
                 cis.pv.threshold=0,
                 snpspos.file=NULL,
                 genepos.file=NULL,
                 cis.dist=10^6,
                 pvalue.hist=FALSE,
                 only.tss=FALSE,
                 fdr.storey=FALSE,
                 perm.pheno.samples=FALSE,
                 seed=1859)
  
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
    else if(args[i] == "--genos-file"){
      params$genos.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--phenos-file"){
      params$phenos.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--cvrt-file"){
      params$cvrt.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--out-file"){
      params$out.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--pv-threshold"){
      params$pv.threshold <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--model"){
      params$model <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--err-cvrt-file"){
      params$err.cvrt.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--cis-out-file"){
      params$cis.out.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--cis-pv-threshold"){
      params$cis.pv.threshold <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--snpspos-file"){
      params$snpspos.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--genepos-file"){
      params$genepos.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "--cis-dist"){
      params$cis.dist <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--pvalue-hist"){
      params$pvalue.hist <- as.numeric(args[i+1])
      i <- i + 1
    }
    else if(args[i] == "--only-tss")
      params$only.tss <- TRUE
    else if(args[i] == "--fdr-storey")
      params$fdr.storey <- TRUE
    else if(args[i] == "--perm-pheno-samples")
      params$perm.pheno.samples <- TRUE
    else if(args[i] == "--seed"){
      params$seed <- as.numeric(args[i+1])
      i <- i + 1
    }
  }
  
  if(params$verbose > 0){
    message("parameters:")
    print(params)
  }
  stopifnot(! is.null(params$genos.file),
            file.exists(params$genos.file),
            ! is.null(params$phenos.file),
            file.exists(params$phenos.file),
            (params$pv.threshold >= 0),
            (params$model %in% c("linear","anova")),
            (params$cis.pv.threshold >= 0),
            (params$cis.dist >= 0))
  if(! is.null(params$cvrt.file))
    stopifnot(file.exists(params$cvrt.file))
  if(is.null(params$out.file)){
    stopifnot(! is.null(params$cis.out.file))
    params$pv.threshold <- 0
  } else{
    if(params$pv.threshold == 0)
      params$out.file <- NULL
  }
  if(! is.null(params$cis.out.file)){
    stopifnot(! is.null(params$snpspos.file),
              file.exists(params$snpspos.file),
              ! is.null(params$genepos.file),
              file.exists(params$genepos.file))
  }
  if(! is.null(params$cvrt.file)){
    tmp <- strsplit(params$cvrt.file, "\\.")[[1]]
    if(tmp[length(tmp)] =="gz")
      stop("--cvrt-file should not be gzipped", call.=FALSE)
  }
  if(params$perm.pheno.samples)
    warning("--perm-pheno-samples is not yet implemented", call.=FALSE)
  if(params$fdr.storey)
    suppressPackageStartupMessages(require(qvalue))
  
  return(params)
}

loadGenotypes <- function(genos.file, verbose=0){
  if(verbose > 0)
    message("load genotypes ...")
  genos <- SlicedData$new()
  genos$fileDelimiter <- "\t"
  genos$fileOmitCharacters <- "-1"
  genos$fileSkipRows <- 1 # one row of column labels
  genos$fileSkipColumns <- 1 # one column of row labels
#  genos$fileSliceSize <- 2000 # read file in pieces of 2,000 rows
  genos$LoadFile(genos.file)
  return(genos)
}

loadPhenotypes <- function(phenos.file, genos, perm.pheno.samples, seed,
                           verbose=0){
  if(verbose > 0)
    message("load phenotypes ...")
  ## tmp <- strsplit(basename(params$phenos.file), "\\.")[[1]]
  ## if(tmp[length(tmp)] == "gz"){
  ##   new.phenos.file <- paste(tmp[1:(length(tmp)-1)], collapse=".")
  ##   if(! file.exists(new.phenos.file)){
  ##     cmd <- paste0("zcat ", params$phenos.file, " > ", new.phenos.file)
  ##     system(cmd)
  ##   }
  ##   params$phenos.file <- new.phenos.file
  ## }
  phenos <- SlicedData$new()
  phenos$fileDelimiter <- "\t"
  phenos$fileOmitCharacters <- "NA"
  phenos$fileSkipRows <- 1
  phenos$fileSkipColumns <- 1
#  phenos$fileSliceSize <- 2000
  phenos$LoadFile(phenos.file)
  
  if(phenos$nCols() != genos$nCols())
    stop("different number of columns in genotype and phenotype files",
         call.=FALSE)
  
  ## reorder phenos columns if necessary
  nb.cols <- genos$nCols()
  same.order <- (sum(colnames(phenos) == colnames(genos)) == nb.cols)
  if(! same.order){
    message("re-order the columns of the phenotype object")
    new.col.order <- rep(NA, nb.cols)
    coln.g <- colnames(genos)
    coln.p <- colnames(phenos)
    for(i in 1:nb.cols)
      new.col.order[i] <- which(coln.p == coln.g[i])
    phenos$ColumnSubsample(new.col.order)
  }
  
  ## for permutations
  if(perm.pheno.samples){
    set.seed(seed)
    new.col.order <- sample(x=1:nb.cols, size=nb.cols)
    phenos$ColumnSubsample(new.col.order)
  }
  
  return(phenos)
}

loadCovariates <- function(cvrt.file, genos, verbose=0){
  cvrt <- SlicedData$new()
  if(! is.null(cvrt.file)){
    if(verbose > 0)
      message("load covariates ...")
    cvrt$fileDelimiter <- "\t"
    cvrt$fileOmitCharacters <- "NA"
    cvrt$fileSkipRows <- 1
    cvrt$fileSkipColumns <- 1
#    cvrt$fileSliceSize <- 2000
    cvrt$LoadFile(cvrt.file)
    
    ## reorder phenos columns if necessary
    nb.cols <- genos$nCols()
    same.order <- (sum(colnames(cvrt) == colnames(genos)) == nb.cols)
    if(! same.order){
      message("re-order the columns of the covariate object")
      new.col.order <- rep(NA, nb.cols)
      coln.g <- colnames(genos)
      coln.p <- colnames(cvrt)
      for(i in 1:nb.cols)
        new.col.order[i] <- which(coln.p == coln.g[i])
      cvrt$ColumnSubsample(new.col.order)
    }
  }
  return(cvrt)
}

loadErrorCovariance <- function(err.cvrt.file, verbose=0){
  errorCovariance <- numeric()
  if(! is.null(err.cvrt.file)){
    if(verbose > 0)
      message("load the error covariance ...")
    errorCovariance <- as.matrix(read.table(err.cvrt.file,
                                            stringsAsFactors=FALSE))
  }
  return(errorCovariance)
}

## BED format: chr\tstart\tend\tname ("start" is 0-based)
loadSnpCoords <- function(snpspos.file, verbose=0){
  snpspos <- NULL
  if(! is.null(snpspos.file)){
    if(verbose > 0)
      message("load SNP coordinates (BED file) ...")
    tmp <- read.table(snpspos.file, nrows=5)
    colClasses <- sapply(tmp, class)
    snpspos <- read.table(snpspos.file, colClasses=colClasses,
                          stringsAsFactors=FALSE)[,c(4,1,3)]
    if(verbose > 0)
      message(paste0("dim(snpspos): ", nrow(snpspos), " x ", ncol(snpspos)))
  }
  return(snpspos)
}

## BED format: chr\tstart\tend\tname ("start" is 0-based)
loadGeneCoords <- function(genepos.file, only.tss, verbose=0){
  genepos <- NULL
  if(! is.null(genepos.file)){
    if(verbose > 0)
      message("load gene coordinates (BED file) ...")
    tmp <- read.table(genepos.file, nrows=5)
    colClasses <- sapply(tmp, class)
    genepos <- read.table(genepos.file, colClasses=colClasses,
                          stringsAsFactors=FALSE)
    if(verbose > 0)
      message(paste0("dim(genepos): ", nrow(genepos), " x ", ncol(genepos)))
    genepos[,2] <- genepos[,2] + 1 # to handle the 0-based start
    if(only.tss) genepos[,3] <- genepos[,2]
    genepos <- genepos[, c(4,1,2,3)] # re-order columns for MatrixEQTL
  }
  return(genepos)
}

launchMatrixeqtl <- function(genos, phenos, cvrt, out.file,
                             pv.threshold, model, errorCovariance,
                             cis.out.file, cis.pv.threshold,
                             snpspos, genepos, cis.dist,
                             pvalue.hist, verbose=0){
  if(verbose > 0)
    message("run all association tests ...")
  res <- Matrix_eQTL_main(snps=genos,
                          gene=phenos,
                          cvrt=cvrt,
                          output_file_name=ifelse(is.null(out.file),
                            "", out.file),
                          pvOutputThreshold=pv.threshold,
                          useModel=ifelse(model == "linear",
                            modelLINEAR, modelANOVA),
                          errorCovariance=errorCovariance,
                          verbose=ifelse(verbose > 0, TRUE, FALSE),
                          output_file_name.cis=ifelse(is.null(cis.out.file),
                            "", cis.out.file),
                          pvOutputThreshold.cis=cis.pv.threshold,
                          snpspos=snpspos,
                          genepos=genepos,
                          cisDist=cis.dist,
                          pvalue.hist=pvalue.hist)
  if(verbose > 0)
    message(paste0("total number of tests performed: ", res$cis$ntests))
  return(res)
}

callEqtlsByStoreyMethod <- function(res, fdr, cis.out.file, verbose){
  set.seed(1859)
  qobj <- qvalue(p=res$pvalue, fdr.level=fdr, pi0.method="bootstrap", robust=TRUE)
  if(verbose > 0){
    message(paste0("pFDR=", fdr))
    message(paste0("pi0=", format(qobj$pi0, digits=7)))
    message(paste0("nb of significant gene-SNP pairs: ", sum(qobj$significant)))
    message(paste0("nb of genes with at least one significant gene-SNP pairs: ",
                   length(unique(res$gene[qobj$significant]))))
  }
  res <- cbind(res, qobj$qvalues)
  colnames(res)[ncol(res)] <- "qvalue.storey"
  
  ## adjustment of gene p-values as min ~ Beta(1,N)
  res.g <- cbind(aggregate(pvalue ~ gene, res, length),
                 aggregate(pvalue ~ gene, res, min)$pvalue)
  names(res.g) <- c("gene", "nb.snps", "pvalue")
  res.g$pvalue.adj <- pbeta(q=res.g$pvalue, shape1=1,
                            shape2=res.g$nb.snps)
  set.seed(1859)
  qobj.adj <- qvalue(p=res.g$pvalue.adj, fdr.level=fdr,
                     pi0.method="bootstrap", robust=TRUE)
  if(verbose > 0)
    message(paste0("nb of significant genes (after adjustment): ",
                   sum(qobj.adj$significant)))
  
  extension <- (tmp <- strsplit(cis.out.file, "\\.")[[1]])[length(tmp)]
  if(extension == "gz"){
    write.table(x=res, file=gzfile(cis.out.file), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
  } else
    write.table(x=res, file=cis.out.file, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
}

run <- function(){
  params <- parseArgs()
  if(params$verbose > 0)
    message(paste0("START ", prog.name, " (", date(), ")"))
  genepos <- loadGeneCoords(params$genepos.file, params$only.tss,
                            params$verbose)
  snpspos <- loadSnpCoords(params$snpspos.file, params$verbose)
  errorCovariance <- loadErrorCovariance(params$err.cvrt.file, params$verbose)
  genos <- loadGenotypes(params$genos.file, params$verbose)
  phenos <- loadPhenotypes(params$phenos.file, genos,
                           params$perm.pheno.samples, params$seed,
                           params$verbose)
  cvrt <- loadCovariates(params$cvrt.file, genos, params$verbose)
  res <- launchMatrixeqtl(genos, phenos, cvrt, params$out.file,
                          params$pv.threshold, params$model, errorCovariance,
                          params$cis.out.file, params$cis.pv.threshold,
                          snpspos, genepos, params$cis.dist,
                          params$pvalue.hist, params$verbose)
  if(params$fdr.storey)
    callEqtlsByStoreyMethod(res=res$cis$eqtls, fdr=0.05, params$cis.out.file,
                            params$verbose)
  if(params$verbose > 0)
    message(paste0("END ", prog.name, " (", date(), ")"))
}

##-----------------------------------------------------------------------------
## Run the program

system.time(run())
print(object.size(x=lapply(ls(), get)), units="Kb")

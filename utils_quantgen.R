## `utils_quantgen.R' contains functions for quantitative genetics/genomics
## Copyright (C) 2013-2015 Institut National de la Recherche Agronomique (INRA)
## License: GPL-3+
## Persons: Timothée Flutre [cre,aut]
## Version: see below
## Download: https://github.com/timflutre/quantgen
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

utils_quantgen.version <- "1.14.0" # http://semver.org/

##' Read a large file as fast as possible
##'
##'
##' @param file
##' @param header
##' @param sep
##' @param ... optional arguments
##' @return data.frame
##' @author Timothée Flutre
read.table.fast <- function(file, header, sep="", ...){
  tmp <- read.table(file, header=header, nrows=5)
  colClasses <- sapply(tmp, class)
  if((tmp <- strsplit(p2f, "\\.")[[1]])[length(tmp)] == "gz"){
    nb.lines <- as.numeric(system(paste0("zcat ", file, " | wc -l"), intern=TRUE))
  } else
    nb.lines <- as.numeric(system(paste0("wc -l < ", file), intern=TRUE))
  read.table(file, header=header, nrows=nb.lines, colClasses=colClasses, sep=sep)
}

##' Random generation for the matrix normal distribution
##'
##' https://stat.ethz.ch/pipermail/r-help/2012-February/302442.html
##' http://en.wikipedia.org/wiki/Matrix_normal_distribution
##' @title Matrix normal distribution
##' @param nrow number of rows
##' @param ncol number of columns
##' @param n number of observations
##' @param M matrix of mean
##' @param U covariance matrix among rows
##' @param V covariance matrix among columns
##' @return matrix
##' @author Timothée Flutre
rmatvnorm <- function(nrow, ncol, n, M=NULL, U=NULL, V=NULL){
  if(is.null(M))
    M <- matrix(data=0, nrow=nrow, ncol=ncol)
  if(is.null(U))
    U <- diag(nrow)
  if(is.null(V))
    V <- diag(ncol)
  Z <- matrix(data=rnorm(n=nrow*ncol, mean=0, sd=1),
              nrow=nrow, ncol=ncol)
  return(M + sqrt(U) %*% Z %*% sqrt(V))
}

read.table.gct <- function(file=NULL){
  ## GCT format: http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats
  stopifnot(! is.null(file))
  tmp <- read.table(file, header=TRUE, row.names=1, skip=2)
  return(as.matrix(tmp[,-1], nrow=nrow(tmp), ncol=ncol(tmp)-1))
}

write.table.gct <- function(x=NULL, file=NULL, gzipped=TRUE){
  ## Write with column names ("id" followed by sample names)
  ## and row names (gene names)
  stopifnot(! is.null(x),
            ! is.null(rownames(x)),
            ! is.null(colnames(x)),
            ! is.null(file))

  ## check suffix
  f.split <- strsplit(file, "\\.")[[1]]
  if(length(f.split) == 1)
    stop("file has no suffix")
  suffix <- f.split[length(f.split)]
  if(gzipped && suffix != "gz"){
    stop("option 'gzipped' is set but file suffix is not 'gz'")
  } else if(! gzipped && suffix == "gz"){
    stop("option 'gzipped' is not set but file suffix is 'gz'")
  }

  ## make temporary data.frame
  tmp <- rbind(colnames(x), x)
  tmp <- cbind(c("id", rownames(x)), tmp)

  ## write file
  if(gzipped){
    write.table(x=tmp, file=gzfile(file), quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep="\t")
  } else
    write.table(x=tmp, file=file, quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep="\t")
}

pca.genexp <- function(X=NULL, algo="svd", method="nb", cutoff=10,
                       scree.file=NULL){
  ## Apply PCA on a gene expression matrix and identify
  ## the leading principal components
  ##
  ## Args:
  ##  X: matrix with genes in rows and samples in columns
  ##  algo: can be 'svd' or 'eigen' (same results modulo a rotation)
  ##  method: to choose which PCs to remove, i.e. to interpret cutoff
  ##          'nb': directly the nb of PCs to remove
  ##          'pve': based on prop of variance explained
  ##          'tw': based on p-value from Tracy-Widom test statistic
  ##          'shriner': min the avg squared partial correlation
  ##  cutoff: to choose which PCs to remove, depends on 'method'
  ##  scree.file: file in which to save scree plot if non NULL
  stopifnot(! is.null(X), is.matrix(X),
            (algo %in% c("eigen","svd")),
            ! is.null(method), (method %in% c("nb","pve","tw")),
            ! is.null(cutoff))
  if(nrow(X) < ncol(X))
    warning("input matrix doesn't seem to have genes in rows and samples in columns")
  if(method == "tw")
    require(RMTstat)

  N <- nrow(X) # nb of genes
  P <- ncol(X) # nb of samples -> find PCs as linear combinations of them

  ## center and scale the input matrix
  ## (centering is not essential, but scaling prevents few samples with
  ## big variance to influence the PCs too much)
  X.cs <- scale(X, center=TRUE, scale=TRUE)

  if(algo == "eigen"){
    ## method 1: get the empirical unbiased covariance matrix between samples
    ## and perform its eigendecomposition
    S <- cov(X.cs) # same as 1/(N-1) * t(X.cs) %*% X.cs
    S.evd <- eigen(S)
    PCs <- S.evd$vectors
    e.vals <- S.evd$values
  } else if(algo == "svd"){
    ## method 2: perform the singular value decomposition of X.cs
    X.cs.svd <- svd(X.cs)
    PCs <- X.cs.svd$v
    e.vals <- X.cs.svd$d
  }
  rownames(PCs) <- colnames(X)
  colnames(PCs) <- paste0("PC", 1:P)

  ## choose the nb of PCs to remove
  if(method == "nb"){
    nb.pcs <- cutoff
  } else if(method == "pve"){
    prop.var.exp <- e.vals / sum(e.vals)
    up.to.which.pc <- which(abs(diff(prop.var.exp)) < cutoff)
    stopifnot(length(up.to.which.pc) != 0)
    nb.pcs <- up.to.which.pc[1]
  } else if(method == "tw"){
    pvals <- ptw(q=e.vals, beta=1, lower.tail=FALSE)
    up.to.which.pc <- which(pvals > cutoff)
    stopifnot(length(up.to.which.pc) != 0)
    nb.pcs <- up.to.which.pc[1] # correct for multiple testing?
  } else if(method == "shriner"){
    nb.pcs <- getNbPCsMinimAvgSqPartCor(X)
  }
  message(paste0("nb of PCs to remove: ", nb.pcs))

  ## scree plot: cumulative PVE versus sorted eigenvalues
  if(! is.null(scree.file)){
    pdf(scree.file)
    plot(x=1:length(e.vals),
         y=cumsum(e.vals/sum(e.vals)),
         type="b", ylim=c(0,1),
         main="Scree plot from PCA",
         xlab="Eigenvalues sorted in decreasing order",
         ylab="Cumulative proportion of variance explained")
    abline(v=nb.pcs)
    dev.off()
    embedFonts(scree.file)
  }

  return(list(pcs=PCs, vars=e.vals, nb.pcs=nb.pcs))
}

rm.confound.genexp <- function(X=NULL, confounders=NULL){
  ## Return residuals of linear regressions used to remove a set
  ## of confounders (e.g. PCs or PEER factors) from a matrix of
  ## gene expression levels
  ##
  ## Args:
  ##  X: matrix with samples in rows and genes in columns
  ##     (it will be centered and scaled before PCs are removed)
  ##  confounders: matrix with samples in rows and confounders in columns
  stopifnot(! is.null(X), is.matrix(X),
            ! is.null(confounders), is.matrix(confounders),
            nrow(X) == nrow(confounders))
  if(nrow(X) > ncol(X))
    warning("input matrix doesn't seem to have samples in rows and genes in columns")

  res <- lm.fit(x=confounders, y=scale(X, center=TRUE, scale=TRUE))
  return(t(res$residuals))
}

##' Impute missing expression levels per subgroup and per gene using the mean
##' (only genes expressed in all subgroups are considered)
##'
##'
##' @param list.mat list of matrices, one per subgroup with genes in rows
##' and samples in columns
##' @return
##' @author Timothée Flutre
imp.miss.genexp <- function(list.mat=NULL){
  stopifnot(! is.null(list.mat), is.list(list.mat))
  for(subgroup in names(list.mat)){
    X <- list.mat[[subgroup]]
    stopifnot(is.matrix(X), ! is.null(rownames(X)), ! is.null(colnames(X)))
    if(nrow(X) < ncol(X))
      warning("input matrix doesn't seem to have genes in rows and samples in columns")
  }

  ## identify all individuals and genes expressed in all subgroups
  all.inds <- sort(unique(do.call(c, lapply(list.mat, colnames))))
  message(paste0("total nb of individuals: ", length(all.inds)))
  all.genes <- table(do.call(c, lapply(list.mat, rownames)))
  message(paste0("total nb of genes: ", length(all.genes)))
  com.genes <- sort(names(all.genes[which(all.genes == length(list.mat))]))
  message(paste0("nb of genes expressed in all subgroups: ", length(com.genes)))

  ## impute per subgroup and per gene
  lapply(list.mat, function(X){
    X.impM <- matrix(nrow=length(com.genes), ncol=length(all.inds))
    rownames(X.impM) <- com.genes
    colnames(X.impM) <- all.inds
    X.impM[rownames(X)[which(rownames(X) %in% com.genes)],colnames(X)] <-
      X[rownames(X)[which(rownames(X) %in% com.genes)],]
    X.impM <- t(apply(X.impM, 1, function(x){
      imp.explevels <- x
      imp.explevels[which(is.na(x))] <- mean(x[which(! is.na(x))])
      imp.explevels
    }))
    X.impM
  })
}

##' Simulate a covariance matrix by drawing random numbers from a uniform distribution
##'
##'
##' @param d dimension of the matrix (number of rows and columns)
##' @param u.min minimum for runif()
##' @param u.max maximum for runif()
##' @param names names of rows and columns
##' @return matrix
##' @author Timothée Flutre
simul.covar.mat <- function(d, u.min=0, u.max=0.5, names=NULL){
	suppressPackageStartupMessages(library(Matrix))
  if(! is.null(names))
    stopifnot(length(names) == d)
  mat <- round(nearPD(matrix(runif(n=d*d, min=0, max=0.5),
                             nrow=d))$mat, 2)
  diag(mat) <- diag(mat) / d + 1
  if(! is.null(names))
    rownames(mat) <- colnames(mat) <- names
  return(mat)
}

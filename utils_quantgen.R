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

utils_quantgen.version <- "1.11.0" # http://semver.org/

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

##' Return the Root Mean Squared Error
##'
##'
##' @title Root Mean Squared Error
##' @param error vector \hat{\theta}_i - \theta_i
##' @return numeric
##' @author Timothée Flutre
rmse <- function(error){
  sqrt(mean(error^2))
}

##' Return the Mean Absolute Error
##'
##'
##' @title Mean Absolute Error
##' @param error vector \hat{\theta}_i - \theta_i
##' @return numeric
##' @author Timothée Flutre
mae <- function(error){
  mean(abs(error))
}

##' Return the Mean Signed Difference
##'
##'
##' @title Mean Signed Difference
##' @param error vector \hat{\theta}_i - \theta_i
##' @return numeric
##' @author Timothée Flutre
msd <- function(error){
  mean(error)
}

##' Return the number of true positives, false positives, true negatives,
##' false negatives, true positive proportion (sensitivity), false positive
##' proportion, accuracy, true negative proportion (specificity), false
##' discovery proportion, false negative proportion and positive predictive
##' value (precision)
##'
##' Both input vectors should be sorted beforehand
##' @title Binary classification
##' @param known.nulls vector of booleans (TRUE if the null is true)
##' @param called.nulls vector of booleans (TRUE if the null is accepted)
##' @return vector with names
##' @author Timothée Flutre
binary.classif <- function(known.nulls, called.nulls){
  ## http://en.wikipedia.org/wiki/Sensitivity_and_specificity
  ##
  ##                                  CALLED
  ##                     Accepted null     Rejected null
  ##
  ##       true null         TN (U)            FP (V)          n0
  ## TRUTH
  ##       false null        FN (T)            TP (S)          n1
  ##
  ##                         a                 r               n
  stopifnot(is.vector(known.nulls), is.vector(called.nulls),
            length(known.nulls) == length(called.nulls),
            sum(! is.logical(known.nulls)) == 0,
            sum(! is.logical(called.nulls)) == 0)

  n <- length(known.nulls) # total number of tests
  n0 <- sum(known.nulls)   # nb of true nulls
  n1 <- n - n0             # nb of "false nulls" (i.e. "true alternatives")
  a <- sum(called.nulls)   # nb of accepted nulls ("called not significant")
  r <- n - a               # nb of rejected nulls ("called significant", "discoveries")

  ## true positive = reject a false null
  tp <- sum(which(! called.nulls) %in% which(! known.nulls))

  ## false positive = reject a true null (type I error, "false alarm")
  fp <- sum(which(! called.nulls) %in% which(known.nulls))

  ## true negatives = accept a true null
  tn <- sum(which(called.nulls) %in% which(known.nulls))

  ## false negatives = accept a false null (type II error, "miss")
  fn <- sum(which(called.nulls) %in% which(! known.nulls))

  tpp <- tp / n1        # true positive prop (sensitivity)
  fpp <- fp / n0        # false positive prop
  acc <- (tp + tn) / n  # accuracy
  tnp <- tn / n0        # true negative prop (specificity), = 1 - fpp
  fdp <- fp / r         # false discovery prop
  fnp <- fn / a         # false negative prop
  ppv <- tp / r         # positive predictive value (precision)

  return(c(n=n, n0=n0, n1=n1, a=a, r=r,
           tp=tp, fp=fp, tn=tn, fn=fn,
           tpp=tpp, fpp=fpp, acc=acc, tnp=tnp, fdp=fdp, fnp=fnp, ppv=ppv))
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

##' Stable computation of log_{10}(\sum_i w_i 10^x_i)
##'
##' Use equal weights if not specified
##' @title Log of weighted sum
##' @param x vector
##' @param weights weights
##' @return numeric
##' @author Timothée Flutre
log10.weighted.sum <- function(x, weights=NULL){
  if(is.null(weights))
    weights <- rep(1/length(x), length(x))
  max <- max(x)
  max + log10(sum(weights * 10^(x - max)))
}

##' Return the Moore-Penrose pseudo-inverse of a matrix
##'
##'
##' @title Pseudo-inverse
##' @param mat matrix
##' @return matrix
##' @author Timothée Flutre
mp.inv <- function(mat){
  mat.svd <- svd(mat)
  mat.svd$v %*% diag(1/mat.svd$d) %*% t(mat.svd$u)
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

##' Return the number of PCs that minimizes the average squared partial
##' correlation
##'
##' Shriner (Heredity, 2011)
##' @param X genotype matrix (0,1,2) with SNPs in rows and individuals in columns
##' @return integer
##' @author Shriner
getNbPCsMinimAvgSqPartCor <- function(X){
  if(nrow(X) < ncol(X))
    warning("input matrix doesn't seem to have genes/snps in rows and samples in columns")
  mu <- apply(X, 1, mean, na.rm=TRUE)
  X <- X - mu
  X2 <- cor(X, use="complete.obs")
  a <- eigen(X2)
  a$values[a$values<0] <- 0
  b <- diag(a$values, nrow=length(a$values))
  loadings <- a$vectors %*% sqrt(b)
  partial <- function(x) {
    c <- loadings[,1:x]
    partcov <- X2 - (c %*% t(c))
    d <- diag(partcov)
    if(any(is.element(NaN,d), is.element(0,d), length(d[d<0])!=0)) {
      map <- 1
    } else {
      d <- 1/(sqrt(d))
      e <- diag(d, nrow=length(d))
      pr <- e %*% partcov %*% e
      map <- (sum(pr^2) - ncol(X2)) / (ncol(X2) * (ncol(X2) - 1))
    }
    return(map)
  }
  fm <- sapply(1:(ncol(X2) - 1), partial)
  fm <- c((sum(X2^2) - ncol(X2))/(ncol(X2) * (ncol(X2) - 1)), fm)
  return(max(1, which.min(fm) - 1))
}

##' Quantile-normalize a vector of numbers to a standard normal distribution.
##'
##' TODO: add ref
##' @param x vector of numeric data
##' @param break.ties.rand break ties randomly (default=TRUE)
##' @param seed see for the pseudo-random number generator (default=1859)
##' @return vector
##' @author Timothée Flutre
quant.norm <- function(x, break.ties.rand=TRUE, seed=1859){
  stopifnot(is.vector(x), is.numeric(x), is.logical(break.ties.rand),
            is.numeric(seed))

  out <- setNames(object=rep(NA, length(x)), nm=names(x))

  if(break.ties.rand){
    if(! is.null(seed))
      set.seed(seed)
    idx <- sample.int(n=length(x))
    tmp <- qqnorm(y=x[idx], plot.it=FALSE)$x
    out <- tmp[sort(idx, index.return=TRUE)$ix]
  } else
    out <- qqnorm(y=x, plot.it=FALSE)$x

  return(out)
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

##' Scales a correlation matrix into the corresponding covariance matrix efficiently.
##'
##' Use sweep. https://en.wikipedia.org/wiki/Covariance_matrix#Correlation_matrix
##' @param x correlation matrix
##' @param sd standard deviations
##' @return matrix
##' @author Timothée Flutre
cor2cov <- function(x, sd){
  ## D <- diag(sd); return(D %*% x %*% D)
  return(sweep(sweep(x, 1, sd, "*"), 2, sd, "*"))
}

##' Convert the SFS of independent replicates into a matrix of allele dosage.
##'
##' SFS stands for site frequency spectrum
##' @param seg.sites list returned by scrm()
##' @return matrix with diploid individuals in rows and SNPs in columns
##' @author Timothée Flutre
seg.sites2all.doses <- function(seg.sites){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)

  nb.inds <- nrow(seg.sites[[1]]) / 2 # nb of diploid individuals
  nb.chrs <- sum(sapply(seg.sites, ncol)) # nb of SNPs
  X <- matrix(data=NA, nrow=nb.inds, ncol=nb.chrs)

  j <- 1
  for(x in seq_along(seg.sites)){
    X[,j:(j+ncol(seg.sites[[x]])-1)] <-
      do.call(rbind, lapply(seq(1, 2*nb.inds, by=2), function(i){
        colSums(seg.sites[[x]][c(i,i+1),])
      }))
    j <- j + ncol(seg.sites[[x]])
  }

  return(X)
}

##' Make a data.frame of SNP coordinates from the SFS of independent replicates
##'
##' SFS stands for site frequency spectrum
##' @param seg.sites list returned by scrm()
##' @param snp.ids vector of identifiers (one per SNP)
##' @return data.frame with SNPs in rows and 2 columns (chr, pos)
##' @author Timothée Flutre
seg.sites2snp.coords <- function(seg.sites, snp.ids){
  stopifnot(is.list(seg.sites),
            length(unique(sapply(seg.sites, nrow))) == 1)

  nb.chrs <- length(seg.sites) # nb of chromosomes
  nb.snps.per.chr <- sapply(seg.sites, ncol)
  nb.snps <- sum(nb.snps.per.chr) # nb of SNPs

  snp.coords <- data.frame(chr=rep(NA, nb.snps),
                           pos=-1,
                           row.names=snp.ids)
  snp.coords$chr <- rep(paste0("chr", 1:nb.chrs), nb.snps.per.chr)
  snp.coords$pos <- as.numeric(do.call(c, lapply(seg.sites, colnames)))

  ## convert genomic positions from float to integer
  i <- 0
  while(TRUE){
    tmp <- c(by(floor(snp.coords$pos * 10^i), factor(snp.coords$chr),
                function(x){anyDuplicated(x)}))
    if(! any(tmp))
      break
    i <- i + 1
  }
  snp.coords$pos <- floor(snp.coords$pos * 10^i)

  return(snp.coords)
}

##' Simulate according to an approximation to the coalescent with recombination named the Sequential Coalescent with Recombination Model.
##'
##' Requires the scrm package (Staab et al, 2014).
##' @param nb.inds diploids (thus nb of haplotypes is 2 * nb.inds)
##' @param ind.ids vector of identifiers (one per individual)
##' @param nb.reps number of independent loci that will be produced (could be seen as distinct chromosomes)
##' @param pop.mut.rate theta = 4 N0 mu
##' @param pop.recomb.rate rho = 4 N0 r
##' @param chrom.len in bp
##' @param nb.pops number of populations
##' @param mig.rate migration rate = 4 N0 m (symmetric)
##' @param verbose verbosity level (default=0=nothing, 1=few, 2=more)
##' @return list with haplotypes (list), genotypes as allele doses (matrix) and SNP coordinates (data.frame)
##' @author Timothée Flutre
simul.coalescent <- function(nb.inds=100,
                             ind.ids=NULL,
                             nb.reps=20,
                             pop.mut.rate=50,
                             pop.recomb.rate=5,
                             chrom.len=10^3,
                             nb.pops=1,
                             mig.rate=5,
                             verbose=0){
  suppressPackageStartupMessages(library(scrm))
  stopifnot(nb.inds > nb.pops)

  if(is.null(ind.ids))
    ind.ids <- sprintf(fmt=paste0("ind%0", floor(log10(nb.inds))+1, "i"),
                       1:nb.inds)

  ## simulate according to the SCRM
  nb.samples <- nb.inds * 2 # e.g. 2 chr1 in ind1, 2 chr1 in ind2, etc
  cmd <- paste0(nb.samples, " ", nb.reps)
  cmd <- paste0(cmd, " -t ", pop.mut.rate)
  cmd <- paste0(cmd, " -r ", pop.recomb.rate, " ", chrom.len)
  cmd <- paste0(cmd, " -T") # print genealogies in newick
  cmd <- paste0(cmd, " -L") # print TMRCA and local tree lengths
  cmd <- paste0(cmd, " -SC abs") # absolute seq positions in bp
  cmd <- paste0(cmd, " -oSFS") # print site freq spectrum, requires -t
  if(nb.pops > 1){
    cmd <- paste0(cmd, " -I ", nb.pops)
    nb.inds.per.pop <- rep(0, nb.pops)
    for(p in 1:(nb.pops-1)){
      nb.inds.per.pop[p] <- floor(nb.inds / nb.pops)
      cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[p])
    }
    nb.inds.per.pop[nb.pops] <- nb.inds - sum(nb.inds.per.pop)
    cmd <- paste0(cmd, " ", 2 * nb.inds.per.pop[nb.pops])
    cmd <- paste0(cmd, " ", mig.rate)
  }
  if(verbose > 1)
    message(cmd)
  sum.stats <- scrm(cmd)
  if(verbose > 1)
    print(str(sum.stats))

  ## make a data.frame with SNP coordinates
  nb.snps.per.chr <- sapply(sum.stats$seg_sites, ncol)
  nb.snps <- sum(nb.snps.per.chr)
  snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                     1:nb.snps)
  snp.coords <- seg.sites2snp.coords(sum.stats$seg_sites, snp.ids)
  for(c in 1:nb.reps){
    colnames(sum.stats$seg_sites[[c]]) <-
      snp.ids[(ifelse(c == 1, 1, 1 + cumsum(nb.snps.per.chr)[c-1])):
                (cumsum(nb.snps.per.chr)[c])]
    rownames(sum.stats$seg_sites[[c]]) <-
      paste0(rep(ind.ids, each=2), rep(c("_h1", "_h2"), nb.inds))
  }

  ## make a matrix with genotypes as allele doses
  X <- seg.sites2all.doses(sum.stats$seg_sites)
  rownames(X) <- ind.ids
  colnames(X) <- snp.ids
  if(verbose > 0){
    txt <- paste0("nb of SNPs: ", nb.snps)
    write(txt, stdout())
    print(sapply(sum.stats$seg_sites, ncol))
  }

  return(list(haplos=sum.stats$seg_sites,
              genos=X,
              snp.coords=snp.coords))
}

##' Calculate the distances between SNPs, assuming they are sorted.
##'
##' Useful before estimating pairwise linkage disequilibrium.
##' @param snp.coords data.frame with SNP identifiers as row names, and with two columns "chr" and "pos"
##' @param nb.cores the number of cores to use (default=1)
##' @return list with one component per chromosome
##' @author Timothée Flutre
snp.distances <- function(snp.coords, nb.cores=1){
  stopifnot(is.data.frame(snp.coords),
            colnames(snp.coords) == c("chr", "pos"),
            ! is.null(rownames(snp.coords)))

  chr.names <- unique(snp.coords$chr)
  snp.dists <- mclapply(chr.names, function(chr.name){
    pos <- snp.coords$pos[snp.coords$chr == chr.name]
    names(pos) <- rownames(snp.coords)[snp.coords$chr == chr.name]
    dis <- pos[2:length(pos)] - pos[1:(length(pos)-1)]
    names(dis) <- paste(names(dis), names(pos)[-length(pos)], sep="-")
    dis
  }, mc.cores=nb.cores)
  names(snp.dists) <- chr.names

  return(snp.dists)
}

##' Estimate kinship matrix from SNPs.
##'
##' SNPs with missing data are ignored.
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows
##' @param mafs vector with minor allele frequencies (calculated with `maf.from.dose` if NULL)
##' @param thresh threshold on allele frequencies below which SNPs are ignored (default=0.01, NULL to skip this step)
##' @param method default is "astle-balding"; "animal-model"; "center", "center-std"
##' @param verbose verbosity level (default=1)
##' @return matrix
##' @author Timothée Flutre
estim.kinship <- function(X, mafs=NULL, thresh=0.01,
                          method="astle-balding", verbose=1){
  stopifnot(is.matrix(X),
            method %in% c("astle-balding", "animal-model", "center", "center-std"))
  if(! is.null(thresh))
    stopifnot(thresh >= 0, thresh <= 0.5)
  N <- nrow(X)
  P <- ncol(X)
  if(P < N)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")

  idx.rm <- c()

  ## discard SNPs with missing data
  snps.na <- apply(X, 2, function(x){
    any(is.na(x))
  })
  if(any(snps.na)){
    if(verbose > 0)
      message(paste0("skip ", sum(snps.na), " SNPs with missing data"))
    idx.rm <- which(snps.na)
    X <- X[, -idx.rm]
    P <- ncol(X)
  }

  ## estimate MAFs
  if(is.null(mafs)){
    mafs <- maf.from.dose(X)
    if(verbose > 1)
      message(paste0("allele freqs: ",
                     "min=", format(min(mafs), digits=2),
                     " Q1=", format(quantile(mafs, 0.25), digits=2),
                     " med=", format(median(mafs), digits=2),
                     " mean=", format(mean(mafs), digits=2),
                     " Q3=", format(quantile(mafs, 0.75), digits=2),
                     " max=", format(max(mafs), digits=2)))
  }

  ## discard SNPs with low MAFs
  if(! is.null(thresh)){
    snps.low <- mafs < thresh
    if(any(snps.low)){
      if(verbose > 0)
        message(paste0("skip ", sum(snps.low), " SNPs with freq below ", thresh))
      idx.rm <- which(snps.low)
      X <- X[, -idx.rm]
      P <- ncol(X)
      mafs <- mafs[-idx.rm]
    }
  }

  ## estimate kinship
  if(method == "astle-balding"){
    tmp <- sweep(x=X, MARGIN=2, STATS=2 * mafs, FUN="-")
    tmp <- sweep(x=tmp, MARGIN=2, STATS=sqrt(4 * mafs * (1 - mafs)), FUN="/")
    K <- tcrossprod(tmp, tmp) / P
  } else if(method == "animal-model"){
    K <- tcrossprod(X, X) / (2  * sum(mafs * (1 - mafs)))
  } else if(method == "center"){
    ## tmp <- sweep(x=X, MARGIN=2, STATS=colMeans(X), FUN="-")
    tmp <- scale(x=X, center=TRUE, scale=FALSE)
    K <- tcrossprod(tmp, tmp) / P
  } else if(method == "center-std"){
    ## X.cs <- sweep(x=X, MARGIN=2, STATS=colMeans(X), FUN="-")
    ## tmp <- sweep(x=X.cs, MARGIN=2, STATS=apply(X=X.cs, MARGIN=2, sd), FUN="/")
    tmp <- scale(x=X, center=TRUE, scale=TRUE)
    K <- tcrossprod(tmp, tmp) / P
  }

  return(K)
}

##' Return estimates of linkage disequilibrium between pairs of SNPs.
##'
##' Requires package LDcorSV
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in columns and individuals in rows
##' @param K matrix of kinship
##' @param pops vector of characters indicating the population of each individual
##' @param snp.coords data.frame with SNP identifiers as row names, and two columns, "chr" and "pos"
##' @param only.chr identifier of a given chromosome
##' @param only.pop identifier of a given population
##' @return
##' @author Timothée Flutre
estim.ld <- function(X, K=NULL, pops=NULL, snp.coords,
                     only.chr=NULL, only.pop=NULL, verbose=0){
  suppressPackageStartupMessages(library(LDcorSV))
  stopifnot(is.matrix(X),
            ! is.null(dimnames(X)),
            sum(is.na(X)) == 0,
            is.data.frame(snp.coords),
            colnames(snp.coords) == c("chr", "pos"))
  if(! is.null(K))
    stopifnot(is.matrix(K),
              nrow(K) == ncol(K),
              nrow(K) == nrow(X),
              ! is.null(dimnames(K)),
              all(rownames(K) == colnames(K)),
              all(rownames(K) == rownames(X)))
  W.s <- NA
  if(! is.null(pops)){
    stopifnot(length(pops) == nrow(X),
              ! is.null(names(pops)),
              names(pops) == rownames(X))
    W.s <- model.matrix(~ as.factor(pops))[, -1]
    rownames(W.s) <- names(pops)
  }
  if(! is.null(only.chr))
    if(! only.chr %in% snp.coords$chr)
      stop(paste0("chr '", only.chr, "' absent from snp.coords"))
  if(! is.null(only.pop))
    if(! only.pop %in% pops)
      stop(paste0("pop '", only.pop, "' absent from pops"))

  ld <- NULL

  subset.snps <- 1:ncol(X)
  if(! is.null(only.chr))
    subset.snps <- which(snp.coords$chr == only.chr)
  subset.inds <- 1:nrow(X)
  if(! is.null(only.pop))
    subset.inds <- which(pops == only.pop)

  if(verbose > 0)
    write("estimate pairwise LD ...", stdout())
  if(is.null(K)){
    if(is.null(only.pop)){
      ld <- LD.Measures(donnees=X[subset.inds, subset.snps],
                        V=NA,
                        S=W.s,
                        data="G", supinfo=FALSE, na.presence=FALSE)
    } else
      ld <- LD.Measures(donnees=X[subset.inds, subset.snps],
                        V=NA,
                        S=NA,
                        data="G", supinfo=FALSE, na.presence=FALSE)
  } else{
    if(is.null(only.pop)){
      ld <- LD.Measures(donnees=X[subset.inds, subset.snps],
                        V=K[subset.inds, subset.inds],
                        S=W.s,
                        data="G", supinfo=FALSE, na.presence=FALSE)
    } else
      ld <- LD.Measures(donnees=X[subset.inds, subset.snps],
                        V=K[subset.inds, subset.inds],
                        S=NA,
                        data="G", supinfo=FALSE, na.presence=FALSE)
  }

  return(ld)
}

##' Simulate a data set from a basic animal model.
##'
##' y = mu 1_n + X b + Z u + e = W a + Z u + e
##' y is n x 1; X is n x P; Z is n x Q; W is n x (P+1)
##' u ~ Norm_Q(0, sigma_u^2 A); e ~ Norm_n(0, sigma^2 I_n)
##' @param n number of individuals (default is 300)
##' @param mu global mean (default is 4)
##' @param P number of fixed effects (default is 1)
##' @param b fixed effects (default is 2)
##' @param nb.snps number of SNPs (default is 1000; ignored if A is given)
##' @param maf minor allele frequency (default is 0.3; ignored if A is given)
##' @param A matrix of additive relationships
##' @param sigma2 variance component of the errors (default is 5)
##' @param lambda ratio of variance components as sigma_u^2 /sigma^2 (default is 3)
##' @return list with all input variables and the data set ready to be analyzed
##' @author Timothée Flutre
simul.animal.model <- function(n=300, mu=4, P=1, b=2, nb.snps=1000, maf=0.3,
                               A=NULL, sigma2=5, lambda=3){
  suppressPackageStartupMessages(library(MASS))
  suppressPackageStartupMessages(library(Matrix))

  animal.ids <- sprintf(fmt=paste0("ind%0", floor(log10(n))+1, "i"), 1:n)
  X <- matrix(data=rnorm(n=n), nrow=n, ncol=P)
  b <- matrix(data=rep(b, P), nrow=P, ncol=1)
  W <- cbind(rep(1, n), X)
  a <- matrix(c(mu, b))
  Q <- n
  if(is.null(A)){
    stopifnot(nrow(A) == n, ncol(A) == n)
    snp.ids <- sprintf(fmt=paste0("snp%0", floor(log10(nb.snps))+1, "i"),
                       1:nb.snps)
    M <- matrix(data=rbinom(n=Q*nb.snps, size=2, prob=maf),
                nrow=Q, ncol=nb.snps, dimnames=list(animal.ids, snp.ids))
    A <- (1/nb.snps) * M %*% t(M)
  }
  Z <- diag(Q)
  sigmau2 <- lambda * sigma2
  h2 <- sigmau2 / (sigmau2 + sigma2)
  G <- as.matrix(nearPD(sigmau2 * A)$mat)
  u <- matrix(mvrnorm(n=1, mu=rep(0, Q), Sigma=G))
  R <- sigma2 * diag(n)
  e <- matrix(mvrnorm(n=1, mu=rep(0, n), Sigma=R))
  y <- W %*% a + Z %*% u + e
  dat <- data.frame(fix=W[,2],
                    animal=factor(animal.ids),
                    response=y[,1])
  return(list(X=X, W=W, Z=Z, G=G, a=a, u=u, sigmau2=sigmau2, sigma2=sigma2,
              h2=h2, dat=dat))
}

##' Simulate phenotypes according to the BSLMM model.
##'
##' See Zhou, Carbonetto & Stephens (2013).
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows (SNPs with missing values or low MAF
##' should be discarded beforehand)
##' @param Q number of covariates (including the intercept)
##' @param pi proportion of beta-tilde values that are non-zero
##' @param h approximation to E[PVE] (h and rho should be NULL or not together)
##' @param rho approximation to E[PGE]
##' @param seed seed for the pseudo-random number generator
##' @return list
##' @author Timothée Flutre
simul.bslmm <- function(X, Q=1, pi=NULL, h=NULL, rho=NULL, seed=NULL){
  stopifnot(xor(is.null(h) & is.null(rho), ! (is.null(h) & is.null(rho))),
            sum(is.na(X)) == 0,
            ! is.null(rownames(X)),
            ! is.null(colnames(X)))
  suppressPackageStartupMessages(library(MASS))
  if(! is.null(seed))
    set.seed(seed)

  ## genetic incidence/design matrices and kinship matrix
  N <- nrow(X)
  P <- ncol(X)
  if(P < N)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")
  X.c <- scale(x=X, center=TRUE, scale=FALSE)
  Z <- diag(N)
  K <- tcrossprod(X.c, X.c) / P

  ## non-genetic covariates
  W <- matrix(data=rep(1, N), nrow=N, ncol=1,
              dimnames=list(rownames(X), c("mu")))
  if(Q > 1)
    W <- cbind(W, matrix(data=rnorm(n=N*(Q-1), mean=0, sd=1), nrow=N, ncol=Q-1,
                         dimnames=list(rownames(X), paste0("c", 1:(Q-1)))))
  alpha <- rnorm(n=Q, mean=50, sd=5)

  ## hyper-parameters
  s.a <- (1 / (N*P)) * sum(colSums(X.c^2))
  s.b <- (1/N)  * sum(diag(K))
  if(is.null(pi))
    pi <- exp(runif(n=1, min=log(1/P), max=log(1)))
  if(is.null(h))
    h <- runif(n=1, min=0, max=1)
  if(is.null(rho))
    rho <- runif(n=1, min=0, max=1)
  sigma.a2 <- (h * rho) / ((1 - h) * P * pi * s.a)
  sigma.b2 <- (h * (1 - rho)) / ((1 - h) * s.b)
  tau <- 1

  ## sparse genetic effects
  betat <- setNames(object=rep(0, P), nm=colnames(X))
  gamma <- setNames(object=rbinom(n=P, size=1, prob=pi), nm=colnames(X))
  betat[gamma == 1] <- rnorm(n=sum(gamma == 1), mean=0,
            sd=sqrt(sigma.a2 * tau^(-1)))

  ## polygenic effects
  u <- setNames(object=mvrnorm(n=1, mu=rep(0, N),
                    Sigma=sigma.b2 * tau^(-1) * K),
                nm=rownames(X))

  ## errors
  epsilon <- setNames(object=matrix(rnorm(n=N, mean=0, sd=sqrt(tau^(-1)))),
                      nm=rownames(X))

  ## phenotypes
  y <- setNames(object=W %*% alpha + X.c %*% betat + Z %*% u + epsilon,
                nm=rownames(X))

  return(list(y=y, W=W, alpha=alpha, X.c=X.c, s.a=s.a, s.b=s.b, pi=pi, h=h,
              rho=rho, sigma.a2=sigma.a2, sigma.b2=sigma.b2, tau=tau,
              betat=betat, u=u))
}

##' Make a scatter plot of y as a function of x, along with the regression
##' line from lm() as well as both confidence and both prediction lines from
##' predict().
##'
##'
##' @param x vector of points
##' @param y vector of points
##' @return nothing
##' @author Timothée Flutre
regplot <- function(x, y, ...){
  x <- as.numeric(x)
  y <- as.numeric(y)
  plot(x, y, ...)
  fit <- lm(y ~ x)
  abline(fit, col="red")
  newx <- seq(min(x), max(x), length.out=length(x))
  pred.ci <- predict(fit, newdata=data.frame(x=newx), interval="confidence")
  lines(newx, pred.ci[,"lwr"], lty=2)
  lines(newx, pred.ci[,"upr"], lty=2)
  pred.pi <- predict(fit, newdata=data.frame(x=newx), interval="prediction")
  lines(newx, pred.pi[,"lwr"], lty=3)
  lines(newx, pred.pi[,"upr"], lty=3)
}

##' Produce a quantile-quantile plot for p values and display its confidence
##' interval
##'
##' A quantile is an order statistic, and the j-th order statistic from a
##' Uniform(0,1) sample has a Beta(j,N-j+1) distribution (Casella & Berger,
##' 2001, 2nd edition, p230).
##' Let us assume we have N independent p values, $\{p_1,\ldots,p_N\}$, for
##' instance: pvalues <- c(runif(99000,0,1), rbeta(1000,0.5,1)). Under the
##' null, they are identically uniformly distributed:
##' $\forall i \; p_i \overset{\text{i.i.d.}{\sim}} \mathcal{U}_{[0,1]}$.
##' Therefore, the 95% confidence interval for the j-th quantile of the set
##' of p values can be calculated with: qbeta(0.95, j, N-j+1).
##' TODO: look at this https://github.com/stephenturner/qqman/blob/v0.0.0/qqman.r
##' @param pvalues vector of raw p values
##' @param plot.conf.int show the confidence interval (default=TRUE)
##' @param xlab a title for the x axis (see default)
##' @param ylab a title for the x axis (see default)
##' @param main an overall title for the plot (default: "Q-Q plot (<length(pvalues)> p-values)")
##' @param col plotting color for the points (default is all points in black)
##' @param ... graphical parameters other than xlim, ylim, xlab, ylab, las and col
##' @author Timothée Flutre (inspired from an anonymous comment to http://gettinggeneticsdone.blogspot.fr/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html)
qqplot.pval <- function(pvalues, plot.conf.int=TRUE,
                        xlab=expression(Expected~~-log[10](italic(p)~values)),
                        ylab=expression(Observed~~-log[10](italic(p)~values)),
                        main=NULL, col=NULL){
  N <- length(pvalues)
  expected <- - log10(1:N / N)
  observed <- - log10(pvalues)
  MAX <- max(c(expected, observed))

  if(plot.conf.int){
    c95 <- rep(0, N)
    c05 <- rep(0, N)
    for(j in 1:N){
      c95[j] <- qbeta(0.95, j, N-j+1)
      c05[j] <- qbeta(0.05, j, N-j+1)
    }
    c95 <- - log10(c95)
    c05 <- - log10(c05)
    plot(expected, c95, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    par(new=T)
    plot(expected, c05, ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    par(new=T)
  }

  if(is.null(main))
    main <- paste0("Q-Q plot (", N, " p-values)")

  if(is.null(col)){
    col <- rep(1, N)
  } else if(length(col) != N)
    stop("param 'col' should have the same length as 'pvalues'")

  plot(x=sort(expected), y=sort(observed),
       xlim=c(0,MAX), ylim=c(0,MAX),
       las=1, col=col[order(observed)],
       xlab=xlab, ylab=ylab, main=main)
  abline(0, 1, col="red")
}

##' Plot a Hinton diagram
##'
##' Modified from http://www.cs.princeton.edu/~mimno/factor-analysis.R
##' @title Hinton diagram
##' @param m matrix
##' @param main main title
##' @param max.sqrt.m to play with the scaling
##' @author Timothée Flutre
hinton <- function(m, main="", max.sqrt.m=NULL){
  rows <- dim(m)[1]
  cols <- dim(m)[2]

  left <- rep(0, rows * cols)
  right <- rep(0, rows * cols)
  bottom <- rep(0, rows * cols)
  top <- rep(0, rows * cols)

  box.colors <- rep("white", rows * cols)

  if(is.null(max.sqrt.m))
    max.sqrt.m <- max(sqrt(abs(m)))
  scale <- 0.9 / (2 * max.sqrt.m)

  position <- 1

  for(row in 1:rows){
    for(col in 1:cols){
      if(m[row,col] < 0)
        box.colors[position] <- "black"
      x <- sqrt(abs(m[row,col]))
      left[position] <- col - (x * scale)
      right[position] <- col + (x * scale)
      top[position] <- -(row - (x * scale))
      bottom[position] <- -(row + (x * scale))
      position <- position + 1
    }
  }

  xlab <- ""
  ylab <- ""
  if(! is.null(names(dimnames(m)))){
    xlab <- names(dimnames(m))[2]
    ylab <- names(dimnames(m))[1]
  }

  par(mar=c(ifelse(xlab == "", 1, 3),
        ifelse(ylab == "", 1, 3), 5, 1) + 0.1)

  plot(0, xlim=c(0.25,cols+0.75), ylim=c(-rows-0.75, -0.25),
       type="n", xaxt="n", yaxt="n", xlab="", ylab="")

  rect(left, bottom, right, top, col=box.colors)

  if(main != "")
    title(main=main, line=3)

  if(! is.null(colnames(m))){
    axis(side=3, at=1:cols, labels=FALSE)
    text(x=1:cols, y=par("usr")[4] + 0.25, labels=colnames(m), adj=0, srt=45, xpd=TRUE)
  } else
    axis(side=3, at=1:cols, labels=1:cols)

  if(! is.null(rownames(m))){
    axis(side=2, at=-(1:rows), labels=rownames(m), las=1)
  } else
    axis(side=2, at=-(1:rows), labels=1:rows, las=1)

  if(xlab != "")
    mtext(xlab, side=1, line=1)
  if(ylab != "")
    mtext(ylab, side=2, line=2)
}

##' Plot a scale, e.g. to add on the side of image()
##'
##' Takes some time to draw (there is one polygon per break...)
##' http://menugget.blogspot.de/2011/08/adding-scale-to-image-plot.html
##' @param z
##' @param zlim
##' @param col
##' @param breaks
##' @param horiz
##' @param ylim
##' @param xlim
##' @param ...
##' @author Timothée Flutre
plot.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(! missing(breaks))
    if(length(breaks) != (length(col)+1))
      stop("must have one more break than colour")

  if(missing(breaks) & ! missing(zlim))
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))

  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2] + c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1] - c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }

  poly <- vector(mode="list", length(col))
  for(i in seq(poly))
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])

  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){
    YLIM <- c(0,1)
    XLIM <- range(breaks)
  } else{
    YLIM <- range(breaks)
    XLIM <- c(0,1)
  }
  if(missing(xlim))
    xlim <- XLIM
  if(missing(ylim))
    ylim <- YLIM

  plot(1, 1, t="n", ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt,
       xaxs="i", yaxs="i", bty="n", ...)

  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    } else
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
}

##' Plot a matrix as a heatmap in its natural orientation, with a colored
##' scale on the right side, and optionally using its dimension names for
##' rows and columns
##'
##' To print all row names, choose idx.rownames=1:nrow(z). To print a subset
##' of 10 row names, choose idx.rownames=floor(seq(1, nrow(z), length.out=10)).
##' Similarly for column names.
##' @param z matrix to be plotted
##' @param main title to appear above the heatmap
##' @param idx.rownames vector giving the indices of the row names of z to be added on the left side of the plot
##' @param idx.colnames vector giving the indices of the column names of z to be added on top of the plot
##' @param breaks vector (default=seq(min(z), max(z), length.out=100))
##' @author Timothée Flutre
image.scale <- function(z, main=NULL, idx.rownames=NULL, idx.colnames=NULL,
                        breaks=NULL){
  if(! is.null(idx.rownames) & is.null(rownames(z)))
    stop("non-null idx.rownames requires z to have row names")
  if(! is.null(idx.colnames) & is.null(colnames(z)))
    stop("non-null idx.colnames requires z to have column names")
  if(is.null(breaks))
    breaks <- seq(min(z), max(z), length.out=100)

  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(7,1))
  ## layout.show(2) # for debugging purposes

  col.pal <- colorRampPalette(c("black", "red", "yellow"), space="rgb")

  ## plot the heatmap
  custom.mar <- c(1, 5, 6, 1)
  if(is.null(idx.rownames))
      custom.mar[2] <- 1
  if(is.null(idx.colnames))
      custom.mar[3] <- 3
  par(mar=custom.mar)
  image(t(z)[,nrow(z):1], axes=FALSE, col=col.pal(length(breaks)-1))
  if(! is.null(main))
    mtext(text=main, side=3, line=ifelse(is.null(idx.colnames), 1, 4),
          font=2, cex=1.3)
  if(! is.null(idx.colnames))
      text(x=seq(0,1,length.out=length(idx.colnames)), y=par("usr")[4]+0.02,
           srt=45, adj=0, labels=colnames(z)[idx.colnames], xpd=TRUE)
  if(! is.null(idx.rownames))
      mtext(text=rev(rownames(z)[idx.rownames]), side=2, line=1,
            at=seq(0,1,length.out=length(idx.rownames)),
            las=2)

  ## plot the scale
  par(mar=c(1,0,6,3))
  plot.scale(z, col=col.pal(length(breaks)-1), breaks=breaks, horiz=FALSE,
             yaxt="n")
  axis(4, at=format(breaks[seq.int(from=1,to=100,length.out=5)], digits=2),
       las=2, lwd=0, lwd.ticks=1)
}

##' Returns the genetic map contained in a BioMercator TXT file.
##'
##' http://moulon.inra.fr/index.php/en/tranverse-team/atelier-de-bioinformatique/projects/projets/135
##' @param file the name of the file which the data are to be read from
##' @return list
##' @author Timothée Flutre
read.biomercator <- function(file){
    stopifnot(file.exists(file))

    gmap <- list()

    lines <- readLines(file)

    ## load meta-data
    i <- 1
    while(! grepl("chr=", lines[i])){
        tokens <- strsplit(lines[i], "=")[[1]]
        key <- gsub(" ", "", tokens[1])
        if(key %in% names(gmap))
            key <- paste(key, sum(key == names(gmap)) + 1, sep=".")
        val <- tokens[2]
        gmap[[key]] <- val
        i <- i + 1
    }

    ## load chromosomes (can have several linkage groups)
    chrs <- split(lines[-(1:(i-1))], cumsum(grepl("^chr=", lines[-(1:(i-1))])))
    gmap[["map"]] <- lapply(chrs, function(chr){
        lgs <- split(chr[-1], cumsum(grepl("^lg=", chr[-1])))
        data <- lapply(lgs, function(lg){
            tmp <- as.data.frame(do.call(rbind, strsplit(lg[-1], "\t")))
            tmp <- tmp[,-1] # discard column of identifiers
            tmp <- tmp[! duplicated(tmp[,1]),] # discard redundant marker names
            df <- as.numeric(as.character(tmp[,2]))
            names(df) <- as.character(tmp[,1])
            df[order(df)]
        })
        names(data) <- paste("lg",
                             sapply(lgs, function(lg){
                                 strsplit(lg[1], "=")[[1]][2]
                             }),
                             sep=".")
        data
    })
    names(gmap[["map"]]) <- paste("chr",
                                  sapply(chrs, function(chr){
                                      strsplit(chr[1], "=")[[1]][2]
                                  }),
                                  sep=".")

    ## add useful numbers
    gmap$nbMarkers <- sum(sapply(gmap$map, function(chr){
        sapply(chr, function(lg){length(lg)})
    }))
    gmap$nbLinkageGroups <- sum(sapply(gmap$map, function(chr){length(chr)}))
    gmap$mapSize <- sum(sapply(gmap$map, function(chr){
        sapply(chr, function(lg){lg[length(lg)]})
    }))

    txt <- paste0("map '", gmap$mapName, "':")
    txt <- paste0(txt, "\n\tnb of individuals: ", gmap$popSize)
    txt <- paste0(txt, "\n\tnb of markers: ", gmap$nbMarkers)
    txt <- paste0(txt, "\n\tnb of chromosomes: ", length(gmap$map))
    txt <- paste0(txt, "\n\tnb of linkage groups: ", gmap$nbLinkageGroups)
    txt <- paste0(txt, "\n\tmap size: ", gmap$mapSize, " ", gmap$mapUnit)
    mean.dist <- do.call(c, lapply(gmap$map, function(chr){
        sapply(chr, function(lg){
            rev(lg)[1:(length(lg)-1)] - lg[(length(lg)-1):1]
        })}))
    txt <- paste0(txt, "\n\tmean distances per linkage group:")
    message(paste0(txt))
    print(summary(mean.dist))

    return(gmap)
}

##' Convert SNP genotype data from alleles (say, "AA" and "AT") to minor allele
##' doses (here, 0 and 1 if "T" is the minor allele).
##'
##' Not particularly efficient, but at least it exists.
##' @param x data.frame with SNPs in rows and individuals in columns, the SNP
##' identifiers being in the first column
##' @param na.string a character to be interpreted as NA values
##' @return list of a matrix (allele doses, SNPs in columns and individuals in
##' rows) and a vector (minor alleles)
##' @author Timothée Flutre
genotypes.alleles2dose <- function(x, na.string="--"){
  stopifnot(is.data.frame(x),
            ! is.null(colnames(x)))

  snp.names <- x[,1]
  P <- length(snp.names)
  ind.names <- colnames(x)[-1]
  N <- length(ind.names)
  message(paste0(P, " SNPs and ", N, " individuals"))

  geno.doses <- matrix(data=NA, nrow=N, ncol=P,
                       dimnames=list(ind=ind.names, snp=snp.names))
  alleles <- matrix(data=NA, nrow=P, ncol=2,
                    dimnames=list(snp.names, c("minor", "major")))

  for(p in 1:P){ # for each SNP
    raw.genos <- unlist(x[p, -1])
    raw.genos[raw.genos == na.string] <- NA
    if(all(is.na(raw.genos))){
      next
    }
    tmp <- do.call(c, strsplit(raw.genos[! is.na(raw.genos)], ""))
    distinct.alleles <- sort(unique(tmp))
    allele.counts <- sort(table(tmp))
    if(length(distinct.alleles) > 2){ # SNP with more than 2 alleles
      stop("SNP ", paste0(x[p,1], " has more than 2 alleles"))
    } else if(length(distinct.alleles) == 2){ # SNP with exactly 2 alleles
      alleles[p, "minor"] <- names(allele.counts)[1]
      alleles[p, "major"] <- names(allele.counts)[2]
      raw.genos <- gsub(pattern=paste0(alleles[p, "minor"], alleles[p, "minor"]),
                        replacement="2", x=raw.genos)
      raw.genos <- gsub(pattern=paste(paste0(alleles[p, "minor"], alleles[p, "major"]),
                            paste0(alleles[p, "major"], alleles[p, "minor"]), sep="|"),
                        replacement="1", x=raw.genos)
      raw.genos <- gsub(pattern=paste0(alleles[p, "major"], alleles[p, "major"]),
                        replacement="0", x=raw.genos)
    } else if(length(distinct.alleles) == 1){ # SNP with only 1 allele
      alleles[p, "major"] <- names(allele.counts)[1]
      raw.genos <- gsub(pattern=paste0(alleles[p, "major"], alleles[p, "major"]),
                        replacement="0", x=raw.genos)
    }
    raw.genos <- as.numeric(raw.genos)
    if(! all(names(table(raw.genos, useNA="no")) %in% c("0", "1", "2")))
      stop("check SNP ", paste0(x[p,1], " (row ", p, ")"))
    geno.doses[,p] <- raw.genos
  }

  return(list(geno.doses=geno.doses,
              alleles=alleles))
}

##' Plot missing SNP genotypes as a grid.
##'
##' Data will be represented in black if missing, white otherwise.
##' @param x matrix with SNP genotypes as allele doses (NA if missing) with
##' SNPs in columns and individuals in rows
##' @param main an overall title for the plot (default="Missing genotypes")
##' @param xlab a title for the x axis (default="Individuals")
##' @param ylab a title for the y axis (default="SNPs")
##' @return nothing
##' @author Timothée Flutre
plotGridMissGenos <- function(x, main="Missing genotypes", xlab="Individuals",
                              ylab="SNPs"){
  if(ncol(x) < nrow(x))
    message("did you put SNPs in columns and individuals in rows?")
  image(1:nrow(x), 1:ncol(x), is.na(x), col=c("white","black"),
        main=main, xlab=xlab, ylab=ylab)
}

##' Estimate minor allele frequencies of SNPs.
##'
##' Missing values should be encoded as NA.
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows
##' @return vector
##' @author Timothée Flutre
maf.from.dose <- function(X){
  stopifnot(is.matrix(X))
  N <- nrow(X)
  P <- ncol(X)
  if(P < N)
    warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")

  maf <- apply(X, 2, function(x){
    x <- x[complete.cases(x)]
    tmp <- sum(x) / (2 * length(x))
    ifelse(tmp <= 0.5, tmp, 1 - tmp)
  })

  return(maf)
}

##' Plot the histogram of the minor allele frequency per SNP
##'
##' Missing values (encoded as NA) are discarded.
##' @title Histogram of minor allele frequencies
##' @param X matrix of SNP genotypes encoded as allele doses, with SNPs in
##' columns and individuals in rows (optional if maf is not null)
##' @param maf vector of minor allele frequencies (optional if X is not null)
##' @param main string for the main title (default="")
##' @param xlim default=c(0,0.5)
##' @param col
##' @param border
##' @param las
##' @param breaks
##' @param ...
##' @return nothing
##' @author Timothée Flutre
plotHistMinAllelFreq <- function(X=NULL, maf=NULL, main="", xlim=c(0,0.5),
                                 col="grey", border="white", las=1,
                                 breaks="FD", ...){
  stopifnot(! is.null(X) || ! is.null(maf))

  if(! is.null(X) & is.null(maf)){
    N <- nrow(X)
    P <- ncol(X)
    if(P < N)
      warning("input matrix doesn't seem to have SNPs in columns and individuals in rows")
    message(paste0(P, " SNPs and ", N, " individuals"))
    maf <- maf.from.dose(X)
  }

  tmp <- hist(x=maf, xlab="Minor allele frequency", ylab="Number of SNPs",
              main=main, xlim=xlim, col=col, border=border, las=las,
              breaks=breaks, ...)
}

##' Convert genotype data to the "mean genotype" file format from BimBam
##'
##' The format is specified in BimBam's manual http://www.haplotype.org/download/bimbam-manual.pdf#page=6
##' @param X matrix with individuals in rows and SNPs in columns
##' @param tX matrix with SNPs in rows and individuals in columns
##' @param alleles data.frame with SNPs in rows (names as row names) and
##' alleles in columns (first is "minor", second is "major")
##' @param file prints the genotype data to this file if non NULL (for instance
##' 'genotypes_bimbam.txt' or gzfile('genotypes_bimbam.txt.gz'))
##' @return data.frame
##' @author Timothée Flutre
genotypes.dose2bimbam <- function(X=NULL, tX=NULL, alleles, file=NULL){
    stopifnot(xor(is.null(X), is.null(tX)),
              ! is.null(row.names(alleles)),
              colnames(alleles) == c("minor","major"))
    if(is.null(tX))
        tX <- t(X)
    tmp <- cbind(alleles, tX)
    if(! is.null(file))
        write.table(x=tmp, file=file, quote=FALSE, sep="\t", row.names=TRUE,
                    col.names=FALSE)
    return(tmp)
}

##' Return the additive relationship matrix, A, from the output of the coancestry() function in the "related" package.
##'
##' The "related" package can be found here: http://frasierlab.wordpress.com/software/. The A matrix is also known as the numerator relationship matrix. It is calculated as explained in chapter 2 from Mrode (2005).
##' @param x list returned by coancestry()
##' @param estim.coancestry name of the coancestry estimator (e.g. "dyadml")
##' @param estim.inbreeding name of the inbreeding estimator (e.g. "LR")
##' @param debug boolean (TRUE to check the output matrix is indeed symmetric)
##' @return matrix
##' @author Timothée Flutre
coancestry2addrel <- function(x, estim.coancestry, estim.inbreeding=NULL,
                              debug=FALSE){
  stopifnot(estim.coancestry %in% colnames(x$relatedness))
  if(! is.null(estim.inbreeding)){
    stopifnot("inbreeding" %in% names(x))
    stopifnot(estim.inbreeding %in% colnames(x$inbreeding))
  }

  ind.ids <- unique(c(x$relatedness$ind1.id,
                      x$relatedness$ind2.id))
  nb.inds <- length(ind.ids)
  A <- matrix(NA, nrow=nb.inds, ncol=nb.inds,
              dimnames=list(ind.ids, ind.ids))

  diag(A) <- 1
  if(! is.null(estim.inbreeding)){
    stopifnot(nrow(x$inbreeding) == nb.inds)
    for(i in 1:nrow(x$inbreeding))
      A[x$inbreeding$ind.id[i], x$inbreeding$ind.id[i]] <-
        1 + x$inbreeding[[estim.inbreeding]][i]
  }

  for(i in 1:nrow(x$relatedness))
    A[x$relatedness$ind1.id[i], x$relatedness$ind2.id[i]] <-
      x$relatedness[[estim.coancestry]][i]
  idx.na.upper <- which(upper.tri(A) & is.na(A))
  idx.equiv.lower <- which(lower.tri(t(A)) & is.na(t(A)))
  A[idx.na.upper] <- A[idx.equiv.lower]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]

  if(debug){ # check that the matrix is symmetric
    for(i in 1:(nrow(A)-1))
      for(j in (i+1):ncol(A))
        if(A[i,j] != A[j,i])
          stop(paste0("matrix not symmetric at (", i, ",", j, ")"))
  }

  return(A)
}

##' Calculate the GC content of a set of sequences.
##'
##' Requires the Biostrings package.
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @return vector
##' @author Timothée Flutre
gc.content <- function(x){
  stopifnot(is.character(x))
	suppressPackageStartupMessages(library(Biostrings))

  sapply(x, function(xi){
    sum(alphabetFrequency(DNAString(xi),
                          baseOnly=TRUE, as.prob=TRUE)[c("C","G")])
  })
}

##' Align each sequence against each other (and itself).
##'
##' Requires the Biostrings package.
##' @title All pairwise alignments
##' @param x vector of sequences (e.g. "AGGT"), possibly with names
##' @param type type of alignment (default="global", i.e. Needleman-Wunsch)
##' @return list of instances of class PairwiseAlignments
##' @author Timothée Flutre
all.pair.aligns <- function(x, type="global", ...){
	stopifnot(is.character(x))
	suppressPackageStartupMessages(library(Biostrings))

	aligns <- list()

	for(i in 1:(length(x)-1)){
		for(j in i:length(x)){
			aligns[[paste0(i,"-",j)]] <-
				pairwiseAlignment(pattern=x[j],
													subject=x[i],
													type=type,
													substitutionMatrix=
													nucleotideSubstitutionMatrix(match=1,
																											 mismatch=0,
																											 baseOnly=FALSE,
																											 type="DNA"),
													...)
		}
	}

	return(aligns)
}

##' Extract statistics from all pairwise alignments
##'
##' Requires the Biostrings package.
##' @param aligns list of instances of class PairwiseAlignments (see all.pair.aligns())
##' @param nb.sequences number of sequences
##' @return list of matrices
##' @author Timothée Flutre
stats.all.pair.aligns <- function(aligns, nb.sequences){
	stopifnot(is.list(aligns))
	suppressPackageStartupMessages(library(Biostrings))

	scores <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	dists <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	pids <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	nmatchs <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	nmismatchs <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	ninss <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)
	ndels <- matrix(data=NA, nrow=nb.sequences, ncol=nb.sequences)

	for(i in 1:(nb.sequences-1)){
		for(j in i:nb.sequences){
			## message(paste0(i,"-",j))
			scores[i,j] <- score(aligns[[paste0(i,"-",j)]])
			dists[i,j] <- nedit(aligns[[paste0(i,"-",j)]])
			pids[i,j] <- pid(aligns[[paste0(i,"-",j)]])
			nmatchs[i,j] <- nmatch(aligns[[paste0(i,"-",j)]])
			nmismatchs[i,j] <- nmismatch(aligns[[paste0(i,"-",j)]])
			ninss[i,j] <- sum(nindel(aligns[[paste0(i,"-",j)]])@insertion[,"Length"] != 0) # insertions in pattern wrt subject
			ndels[i,j] <- sum(nindel(aligns[[paste0(i,"-",j)]])@deletion[,"Length"] != 0) # idem
		}
	}

	return(list(scores=scores, dists=dists, pids=pids, nmatchs=nmatchs,
							nmismatchs=nmismatchs, ninss=ninss, ndels=ndels))
}

##' Initialize plates as data.frames with missing data.
##'
##' Useful for molecular biologists.
##' @param n number of plates
##' @param nrow vector of number of rows for each plate
##' @param ncol vector of number of columns for each plate
##' @param names vector of names for each plate
##' @return list of data.frame, one per plate, in the "wide" format
##' @author Timothée Flutre
init.plates <- function(n, nrow, ncol, names){
  plates <- list()
  for(i in 1:n){
    plates[[names[i]]] <-
      as.data.frame(matrix(data=NA, nrow=nrow[i],
                           ncol=ncol[i],
                           dimnames=list(LETTERS[1:nrow[i]], 1:ncol[i])))
  }
  return(plates)
}

##' Lengthen a "wide" plate into 3 columns for easier processing.
##'
##' Useful for molecular biologists.
##' @param plate data.frame of a plate in the "wide" format
##' @return data.frame of a plate in the "long" format (1 well per row)
##' @author Timothée Flutre
lengthen.plate <- function(plate.w){
  stopifnot(is.data.frame(plate.w),
            ! is.null(rownames(plate.w)),
            ! is.null(colnames(plate.w)))
  nb.samples <- nrow(plate.w) * ncol(plate.w)
  plate.l <- data.frame(sample=rep(NA, nb.samples),
                      row=rep(NA, nb.samples),
                      col=rep(NA, nb.samples))
  sample.id <- 1
  for(i in 1:nrow(plate.w)){
    for(j in 1:ncol(plate.w)){
      plate.l$sample[sample.id] <- plate.w[i,j]
      plate.l$row[sample.id] <- rownames(plate.w)[i]
      plate.l$col[sample.id] <- colnames(plate.w)[j]
      sample.id <- sample.id + 1
    }
  }
  return(plate.l)
}

##' Identify empty wells, if any, in a plate.
##'
##' Useful for molecular biologists.
##' @param plates data.frame of a plate in the "wide" format
##' @return 2 column data.frame (row;col) corresponding to empty wells
##' @author Timothée Flutre
empty.wells <- function(plate.w){
  plate.l <- lengthen.plate(plate.w)
  empty.idx <- is.na(plate.l$sample)
  return(plate.l[empty.idx, c("row", "col")])
}

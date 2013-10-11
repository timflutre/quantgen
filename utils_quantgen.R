## Gather useful functions in quantitative genetics/genomics
## Author: Timothee Flutre
## License: GPL-3

read.table.fast <- function(file, header, sep="", ...){
  ## Read a large file as fast as possible
  tmp <- read.table(file, header=header, nrows=5)
  colClasses <- sapply(tmp, class)
  if((tmp <- strsplit(p2f, "\\.")[[1]])[length(tmp)] == "gz"){
    nb.lines <- as.numeric(system(paste0("zcat ", file, " | wc -l"), intern=TRUE))
  } else
    nb.lines <- as.numeric(system(paste0("wc -l < ", file), intern=TRUE))
  read.table(file, header=header, nrows=nb.lines, colClasses=colClasses, sep=sep)
}

rmse <- function(error){
  ## Return the Root Mean Squared Error
  ##
  ## Args:
  ##  error: vector \hat{\theta}_i - \theta_i
  sqrt(mean(error^2))
}

mae <- function(error){
  ## Return the Mean Absolute Error
  ##
  ## Args:
  ##  error: vector \hat{\theta}_i - \theta_i
  mean(abs(error))
}

msd <- function(error){
  ## Return the Mean Signed Difference
  ##
  ## Args:
  ##  error: vector \hat{\theta}_i - \theta_i
  mean(error)
}

binary.classif <- function(known.nulls, called.nulls){
  ## Return the number of true positives, false positives, true negatives,
  ## false negatives, true positive proportion (sensitivity), false positive
  ## proportion, accuracy, true negative proportion (specificity), false
  ## discovery proportion, false negative proportion and positive predictive
  ## value (precision)
  ##
  ## Args:
  ##  known.nulls: vector of booleans (TRUE if the null is true)
  ##  called.nulls: vector of booleans (TRUE if the null is accepted)
  ##
  ## Notes:
  ##  both vectors should be sorted beforehand
  ##
  ## Returns:
  ##  vector with names
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

rmatvnorm <- function(nrow, ncol, n, M=NULL, U=NULL, V=NULL){
  ## Random generation for the matrix-variate normal distribution
  ## with mean equal to M, among-row covariance equal to U and
  ## among-column covariance equal to V
  ## https://stat.ethz.ch/pipermail/r-help/2012-February/302442.html
  ## http://en.wikipedia.org/wiki/Matrix_normal_distribution
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

log10.weighted.sum <- function(x, weights=NULL){
  ## Stable computation of log_{10}(\sum_i w_i 10^x_i)
  ## use equal weights if not specified
  if(is.null(weights))
    weights <- rep(1/length(x), length(x))
  max <- max(x)
  max + log10(sum(weights * 10^(x - max)))
}

mp.inv <- function(mat){
  ## Return the Moore-Penrose pseudo-inverse of a matrix
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

getNbPCsMinimAvgSqPartCor <- function(X){
  ## Return the number of PCs that minimizes the average squared partial
  ##  correlation (Shriner, Heredity, 2011)
  ##
  ## Args:
  ##  dat: genotype matrix (0,1,2) with SNPs in rows and individuals in columns
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

imp.miss.genexp <- function(list.mat=NULL){
  ## Impute missing expression levels per subgroup and per gene using the mean
  ## (only genes expressed in all subgroups are considered)
  ##
  ## Args:
  ##  list.mat: list of matrices, one per subgroup with genes in rows
  ##            and samples in columns
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

hinton <- function(m, main="", max.sqrt.m=NULL){
  ## Plot a Hinton diagram
  ##
  ## Args:
  ##  m: matrix
  ##  main: main title (default="")
  ##  max.sqrt.m: to play with the scaling
  ##
  ## Refs:
  ##  modified from http://www.cs.princeton.edu/~mimno/factor-analysis.R
  
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

## `utils_quantgen.R' contains utility functions for quantitative genetics/omics
## Copyright (C) 2012-2014 Timothée Flutre
## License: GPL-3+
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

##' Read a large file as fast as possible
##'
##' .. content for \details{} ..
##' @title 
##' @param file 
##' @param header 
##' @param sep 
##' @param ... 
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
##' .. content for \details{} ..
##' @title Root Mean Squared Error
##' @param error vector \hat{\theta}_i - \theta_i
##' @return numeric
##' @author Timothée Flutre
rmse <- function(error){
  sqrt(mean(error^2))
}

##' Return the Mean Absolute Error
##'
##' .. content for \details{} ..
##' @title Mean Absolute Error
##' @param error vector \hat{\theta}_i - \theta_i
##' @return numeric
##' @author Timothée Flutre
mae <- function(error){
  mean(abs(error))
}

##' Return the Mean Signed Difference
##'
##' .. content for \details{} ..
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
##' .. content for \details{} ..
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
##' @title 
##' @param X genotype matrix (0,1,2) with SNPs in rows and individuals in columns
##' @return 
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
##' .. content for \details{} ..
##' @title 
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

##' Estimate kinship matrix K via the Astle-Balding model
##'
##' eqn 2.2 of Astle & Balding (Statistical Science, 2009)
##' @title Astle-Balding matrix
##' @param genos.dose matrix with SNPs in rows and individuals in columns,
##' and genotypes coded as allele dose (i.e. in [0,2]),
##' but no missing data is allowed
##' @return Kinship matrix
estimKinshipAstleBalding <- function(genos.dose){
    allele.freqs <- apply(genos.dose, 1, function(x){
        sum(x) / (2 * length(x))
    })
    tmp <- sweep(genos.dose, 1, 2 * allele.freqs, FUN="-")
    tmp <- sweep(tmp, 1, 2 * sqrt(allele.freqs * (1 - allele.freqs)), FUN="/")
    K <- (1/nrow(genos.dose)) * crossprod(tmp, tmp)
    return(K)
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
##' @title 
##' @param pvalues vector of raw p values
##' @param main main title (default="Q-Q plot")
##' @return 
##' @author Timothée Flutre (inspired from an anonymous comment to http://gettinggeneticsdone.blogspot.fr/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html)
qqplot.pval <- function(pvalues, main="Q-Q plot", ...){
    N <- length(pvalues)
    c95 <- rep(0, N)
    c05 <- rep(0, N)
    for(j in 1:N){
        c95[j] <- qbeta(0.95, j, N-j+1)
        c05[j] <- qbeta(0.05, j, N-j+1)
    }
    
    pvalues <- -log10(pvalues)
    null <- -log10(1:N/N)
    MAX <- max(c(pvalues, null))
    
    plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    par(new=T)
    plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
         axes=FALSE, xlab="", ylab="")
    abline(0, 1, col="red") # diagonal
    par(new=T)
    qqplot(x=null, y=pvalues, ylim=c(0,MAX), xlim=c(0,MAX), las=1,
           main=main, xlab=expression(Expected~~-log[10](italic(p)~values)),
           ylab=expression(Observed~~-log[10](italic(p)~values)), ...)
}

##' Plot a Hinton diagram
##'
##' Modified from http://www.cs.princeton.edu/~mimno/factor-analysis.R
##' @title Hinton diagram
##' @param m matrix
##' @param main main title
##' @param max.sqrt.m to play with the scaling
##' @return 
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
##' @title 
##' @param z 
##' @param zlim 
##' @param col 
##' @param breaks 
##' @param horiz 
##' @param ylim 
##' @param xlim 
##' @param ... 
##' @return 
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

##' Plot a matrix as a heatmap in its natural orientation,
##' using its dimension names for rows and columns,
##' and with a colored scale on the right side
##'
##' .. content for \details{} ..
##' @title 
##' @param z matrix
##' @param main title to appear above the heatmap
##' @param breaks vector (default=seq(min(z), max(z), length.out=100))
##' @return 
##' @author Timothée Flutre
image.scale <- function(z, main=NULL, breaks=NULL){
  
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(7,1))
  ## layout.show(2)
  
  col.pal <- colorRampPalette(c("black", "red", "yellow"), space="rgb")
  if(is.null(breaks))
    breaks <- seq(min(z), max(z), length.out=100)
  
  ## plot the heatmap
  par(mar=c(1, 5, 6, 1))
  image(t(z)[,nrow(z):1], axes=FALSE, col=col.pal(length(breaks)-1))
  if(! is.null(main))
    mtext(text=main, side=3, line=4, font=2, cex=1.3)
  text(x=seq(0,1,length.out=nrow(z)), y=par("usr")[4]+0.02, srt=45, adj=0,
       labels=rownames(z), xpd=TRUE)
  mtext(text=rev(colnames(z)), side=2, line=1, at=seq(0,1,length.out=nrow(z)),
        las=2) # left side
  
  ## plot the scale
  par(mar=c(1,0,6,3))
  plot.scale(z, col=col.pal(length(breaks)-1), breaks=breaks, horiz=FALSE,
             yaxt="n")
  axis(4, at=format(breaks[seq.int(from=1,to=100,length.out=5)], digits=2),
       las=2, lwd=0, lwd.ticks=1)
}

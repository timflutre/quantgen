## Gather useful functions in quantitative genetics/genomics
## Author: Timothee Flutre
## License: GPL-3

## Random generation for the matrix-variate normal distribution
## with mean equal to M, among-row covariance equal to U and
## among-column covariance equal to V
rmatvnorm <- function(nrow, ncol, n, M=NULL, U=NULL, V=NULL){
  if(is.null(M))
    M <- matrix(data=0, nrow=nrow, ncol=ncol)
  if(is.null(U))
    U <- diag(nrow)
  if(is.null(V))
    V <- diag(ncol)
  Z <- matrix(data=rnorm(n=nrow*ncol, mean=1, sd=1), nrow=nrow, ncol=ncol)
  return(M + U %*% Z %*% V)
}

## Stable computation of log_{10}(\sum_i w_i 10^x_i)
## use equal weights if not specified
log10.weighted.sum <- function(x, weights=NULL){
  if(is.null(weights))
    weights <- rep(1/length(x), length(x))
  max <- max(x)
  max + log10(sum(weights * 10^(x - max)))
}

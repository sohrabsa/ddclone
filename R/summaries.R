# summary
load.original.data <- function(dataID) {
  dat <- simulated.data.for.id(dataID)
  nMut <- length(dat$mutPrevalence)
  data.frame(id=1:nMut, prev=as.vector(dat$mutPrevalence))
}

#' Computes a point estimate (using MAXPEAR method) for clustering of the mutations from the results of a ddClone analysis run
#'
#' @param expPath path to the directory where the ddClone results are stored
#' @param MCMCOptions a listing containing MCMC options including \code{thinning} and \code{burnIn} to be used.
#' @return A matrix each row of which is a point estimate of clustering assignment for each mutation. See mcclust::maxpear for more details.
#' @export
estimatePointClustering <- function(expPath, MCMCOptions=list(thinning=10, burnIn=100)) {
  filePath <- file.path(expPath, 'clust-trace.csv')
  clust <- read.table(filePath, stringsAsFactors=F)
  step <- MCMCOptions$thinning
  clust <- clust[seq(MCMCOptions$burnIn + 1, nrow(clust), by=step), ]
  # remove na at the end with a warning
  if (any(is.na(clust))) {
    print('warning - there are NA entries in the trace. Maybe an early finish?')
    clust <- na.omit(clust)
  }
  psm <- mcclust::comp.psm(as.matrix(clust))
  tempMat <- diag(ncol(clust))
  if (all(psm == tempMat))  {
    mpear <- mcclust::maxpear(psm, method='avg')
    mpear$cl <- as.matrix(t(data.frame(best=mpear$cl, avg=mpear$cl, comp=mpear$cl, draws=mpear$cl)))
  } else {
    mpear <- mcclust::maxpear(psm, method='all', cls.draw = as.matrix(clust))
  }
  mpear$cl
}

#' Computes a Monte Carlo estimate for the cellular prevalences of each mutation
#'
#' @param expPath path to the directory where the ddClone results are stored
#' @param MCMCOptions a listing containing MCMC options including \code{thinning} and \code{burnIn} to be used.
#' @return a vector containingcontaining the Monte Carlo estimate for each mutation
#' @export
estimatePointPhi <- function(expPath, MCMCOptions=list(thinning=10, burnIn=100)) {
  filePath <- file.path(expPath, 'phi-trace.csv')
  phi <- read.table(filePath)
  phi <- phi[-(1:(MCMCOptions$burnIn)), ]
  phi <- phi[seq(1, nrow(phi), by=MCMCOptions$thinning), ]
  colMeans(phi, na.rm = T)
}

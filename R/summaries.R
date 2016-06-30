# summary
load.original.data <- function(dataID) {
  dat <- simulated.data.for.id(dataID)
  nMut <- length(dat$mutPrevalence)
  data.frame(id=1:nMut, prev=as.vector(dat$mutPrevalence))
}

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
  psm <- comp.psm(as.matrix(clust))
  tempMat <- diag(ncol(clust))
  if (all(psm == tempMat))  {
    mpear <- maxpear(psm, method='avg')
    mpear$cl <- as.matrix(t(data.frame(best=mpear$cl, avg=mpear$cl, comp=mpear$cl, draws=mpear$cl)))
  } else {
    mpear <- maxpear(psm, method='all', cls.draw = as.matrix(clust))
  }
  mpear$cl
}


estimatePointPhi <- function(expPath, MCMCOptions=list(thinning=10, burnIn=100)) {
  filePath <- file.path(expPath, 'phi-trace.csv')
  phi <- read.table(filePath)
  phi <- phi[-(1:(MCMCOptions$burnIn)), ]
  phi <- phi[seq(1, nrow(phi), by=MCMCOptions$thinning), ]
  colMeans(phi, na.rm = T)
}

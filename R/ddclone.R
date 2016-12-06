# ddClone starting point
source('R/env-setup.R')
source('R/driver.R')
source('R/summaries.R')

#' Runs the ddClone MCMC-based inference algorithm over the genotype and allele-counts of a set of muations
#'  and computes point estimates for clustering assignment and cellular prevalences.
#'
#' @param dataObj The data object, containing genotpe matrix and allele counts
#' @param outputPath a directory where to store the results and temp files
#' @return A list containing \code{df} the estimated clustering and cellular prevalences; \code{expPath} where the results are stored
#' @export
ddclone <- function(dataObj, outputPath='.', tumourContent = 1.0,
                    numOfIterations = 100, thinning = 1, burnIn,
                    a = 0.01, alpha = 1, s = 1, seed = 10, useTraditionalCRP = F,
                    dist.fn = jaccardDist,
                    grid.mA = 10, grid.mS = 10, grid.mAlpha = 10) {
  # run the sampler
  if (is.null(burnIn)) burnIn <- floor(niter/10)
  MCMCOptions <- list(thinning = thinning, niter = numOfIterations, burnIn = burnIn)

  hyperParams <- list(hyperParamAlpha=alpha, hyperParamA=a, hyperParamS=s)
  decay.fn.name <- 'exp.decay.s'

  if (useTraditionalCRP) {
    decay.fn.name <- 'identity.decay'
    dist.fn <- identity.s
  }
  set.seed(seed)
  res <- driver(niter = numOfIterations, dist.fn = dist.fn,
                decay.fn.name =  decay.fn.name,
                decay.fn = getFunction(decay.fn.name),
                dataObj = dataObj,
                outputPath = outputPath,
                tumourContent = tumourContent,
                genotype.prior.scheme ='PCN',
                hyperParamAlpha=alpha, hyperParamA=a, hyperParamS=s,
                MCMCOptions = MCMCOptions,
                grid.mA=grid.mA, grid.mS=grid.mS, grid.mAlpha=grid.mAlpha)

  res$hyperParams <- hyperParams
  res$TraditionalCRP <- useTraditionalCRP
  res
  expPath <- res$expPath

  # compute point estiamtes (use MCMC options)
  estClust <- estimatePointClustering(expPath = expPath, MCMCOptions = MCMCOptions)
  estPhi <- estimatePointPhi(expPath = expPath, MCMCOptions = MCMCOptions)

  write.table(estClust, file.path(expPath, 'clust-est.csv'))
  write.table(estPhi, file.path(expPath, 'phi-est.csv'))

  # wrap up into a data.frame
  df <- data.frame(mutID = colnames(dataObj$filteredMutMatrix), phi = unname(estPhi), clusterID = estClust[1, ])
  list(df = df, expPath = expPath)
}


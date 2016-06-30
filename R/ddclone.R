# ddClone starting point
source('R/env-setup.R')
source('R/driver.R')
source('R/summaries.R')

# dataPath = './data/dollo.10.48.5.f0.gl10-u.dat'; outputPath = './output/out1'
# burnIn = 10; thinning = 1; tumourContent = 1.0; numOfIterations = 100; a = 0.01; alpha = 0.01; s = 1000; seed = 1; useTraditionalCRP = F;
# grid.mA = 100; grid.mS = 100; grid.mAlpha = 100
ddclone <- function(dataPath, outputPath, tumourContent = 1.0,
                    numOfIterations = 100, thinning = 1, burnIn,
                    a = 0.01, alpha = 1, s = 1, seed = 10, useTraditionalCRP = F,
                    grid.mA = 10, grid.mS = 10, grid.mAlpha = 10) {
  # run the sampler
  if (is.null(burnIn)) burnIn <- floor(niter/10)
  MCMCOptions <- list(thinning = thinning, niter = numOfIterations, burnIn = burnIn)

  hyperParams <- list(hyperParamAlpha=alpha, hyperParamA=a, hyperParamS=s)
  decay.fn.name <- 'exp.decay.s'
  dist.fn <- jaccardDist
  if (useTraditionalCRP) {
    decay.fn.name <- 'identity.decay'
    dist.fn <- identity.s
  }
  set.seed(seed)
  res <- driver(niter = numOfIterations, dist.fn = dist.fn,
                decay.fn.name =  decay.fn.name,
                decay.fn = getFunction(decay.fn.name),
                dataPath = dataPath,
                outputPath = outputPath,
                tumourContent = tumourContent,
                genotype.prior.scheme ='PCN',
                hyperParamAlpha=alpha, hyperParamA=a, hyperParamS=s,
                MCMCOptions = MCMCOptions,
                grid.mA=grid.mA, grid.mS=grid.mS, grid.mAlpha=grid.mAlpha)

  res$hyperParams <- hyperParams
  res$dataID <- basename(dataPath)
  res$TraditionalCRP <- useTraditionalCRP
  res
  expPath <- res$expPath

  # compute point estiamtes (use MCMC options)
  estClust <- estimatePointClustering(expPath = expPath, MCMCOptions = MCMCOptions)
  estPhi <- estimatePointPhi(expPath = expPath, MCMCOptions = MCMCOptions)

  write.table(estClust, file.path(expPath, 'clust-est.csv'))
  write.table(estPhi, file.path(expPath, 'phi-est.csv'))

  # wrap up into a data.frame
  df <- data.frame(mutID = colnames(readRDS(dataPath)$filteredMutMatrix), phi = unname(estPhi), clusterID = estClust[1, ])
  list(df = df, expPath = expPath, dataPath = dataPath)
}

# ddclone(dataPath = './data/dollo.10.48.4.f0.gl0-u.dat', outputPath = './output/out1', burnIn = 10, seed = 10, numOfIterations = 100)

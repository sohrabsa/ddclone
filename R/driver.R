source('R/helper.R')
source('R/ddcrp-inference.R')
source('R/env-setup.R')

driver <- function(niter=100, decay.fn=window.fn.s, decay.fn.name='window.fn.s',
                   dataPath='', outputPath, genotype.prior.scheme='AB',
                   dist.fn=jaccardDist, hyperParamAlpha=1, hyperParamA=.5, hyperParamS=1000, tumourContent=1,
                   MCMCOptions,
                   grid.mS=10, grid.mA=1000, grid.mAlpha=1000,
                   rndShuffleHyperParams = T, resampleHyperParams = T, permuteCustomers = T) {
  mPhi <- 100
  if (resampleHyperParams) {
    mS <- grid.mS
    mA <- grid.mA
    mAlpha <- grid.mAlpha
    sRange <- list(min=10, max=max(hyperParamS, 1000))
    aRange <- list(min=0.01, max=max(hyperParamA, 1))
    alphaRange <- list(min=0.01, max=max(10, hyperParamAlpha))
  } else {
    mS <- 1
    mA <- 1
    mAlpha <- 1
    sRange <- list(min=hyperParamS, max=hyperParamS)
    aRange <- list(min=hyperParamA, max=hyperParamA)
    alphaRange <- list(min=hyperParamAlpha, max=hyperParamAlpha)
  }

  # create experiment DIR
  if (is.null(outputPath)) {
    basePath <- get.path('.')
    timeTag <- format(Sys.time(), "%Y-%m-%d-%H-%M-%OS6")
    expPath <- file.path(basePath, timeTag)
  } else {
    expPath <- outputPath
  }

  dir.create(expPath, recursive = T, showWarnings = F)
  print(expPath)

  dataID <- basename(dataPath)
  simulatedData <- readRDS(dataPath)
  distMat <- dist.fn(simulatedData)
  dat <- make.pyclone.input(simulatedData)

  # dirichlet process concentration alpha, decay function parameter a, and betabinomial precision s
  hyperParams <- list(alpha=hyperParamAlpha, a=hyperParamA, s=hyperParamS, t=tumourContent)

  # Jaccard distance
  dist.fn <- matrix.dist.fn(distMat)

  datM <- Matrix(as.matrix(data.frame(seq(nrow(dat)))))

  profile <- list(mPhi=mPhi, mS=mS, mA=mA, mAlpha=mAlpha,
                  MCMCOptions, hyperParams, decay.fn.name=decay.fn.name,
                  simulated.data.id=basename(dataPath), genotype.prior.scheme=genotype.prior.scheme,
                  result.path=expPath)
  write.table(t(as.data.frame(profile)), file.path(expPath, 'config.csv'))

  psi.priors <- make.psi.priors(cpDat = data.frame(minor_cn=dat$minor_cn, major_cn=dat$major_cn), scheme = genotype.prior.scheme)

  print('Generating cached matrices')
  LCACHED <<- make.cached.genotype.aware.likelihood(mPhi, dat, mS, dataID, psi.priors, sRange, tumourContent = tumourContent)
  Decay.CACHED <<- make.cached.decay(mA, distMat, decay.fn, decay.fn.name, dataID, aRange)
  AlphaCACHED  <<- make.cached.alpha(mAlpha, distMat, decay.fn, decay.fn.name, dataID, alphaRange)

  # log-likelihood
  lhood.fn <- cached.pyclone.dd.crp.likelihood(dat)

  ddcrp <- ddcrp.gibbs(dat=datM,
                       lhood.fn=lhood.fn,
                       dist.fn=dist.fn,
                       decay.fn=decay.fn,
                       summary.fn=ncomp.summary,
                       log.prior.thresh = -10,
                       hyperParams=hyperParams,
                       MCMCOptions=MCMCOptions,
                       phi.traj=TRUE,
                       clust.traj=TRUE,
                       processor.fn=sample.processor,
                       resampleHyperParams = resampleHyperParams,
                       permuteCustomers = permuteCustomers,
                       expPath = expPath)

  saveRDS(ddcrp, file.path(expPath, 'clust-result-obj.dat'))
  append.traj(ddcrp$phi.traj, ddcrp$clust.traj, ddcrp$hyperParamTraj, expPath)

  print(expPath)
  msg('Done!')
  list(expPath=expPath)
}

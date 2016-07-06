source('R/genotype.R')
source('R/env-setup.R')

library(VGAM)
library(vegan)
library(matrixStats)

EmptyCache <- T

#' Generate a ddClone input object from bulk allele counts and single cell data and saves it in \code{bulkDataPath}.
#'
#' @param bulkDataPath Path to a .tsv file that contains allele counts and copy number data. Expects colnames to be "mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", and "major_cn". "minor_cn" and "major_cn" should be integer values for the minimum and maxmum estimated copy number at that locus respectively.
#' @param genotypeMatrixPath Path to a .tsv file of binary entries where rows correspond to genotypes and coloumns to mutation IDs respectively.
#' @param outputPath What directory should the output be saved into.
#' @param nameTag An optional string to be appended to the name of the output object
#' @return A list with appropriate members that could be given as input for ddClone analysis, for instance \code{ddclone::ddclone}
#' @export
make.ddclone.input <- function(bulkDataPath, genotypeMatrixPath, outputPath, nameTag='') {
#   $ mutCounts        : int [1:2, 1:36]  row1: total_counts', row2: 'var_counts'
#   $ psi              :List of 36
#   $ filteredMutMatrix: num [1:11, 1:36] 0 0 0 0 0 0 0 1 1 1 ...

  require(xlsx)
  # example data set
  # 1. read the genotype-mutation matrix
  inputPath <- '/Users/sohrab/Google Drive/Masters/Thesis/scripts/ddcrppaper/additional_files/additional_file_4_inputs_simulated.xlsx'
  genDat <- read.xlsx(file = inputPath, sheetName = 'seed1_genotypes', row.names=T)
  genDatMutList <- colnames(genDat)

  # 2. read the bulk data
  bulkDat <- read.xlsx(file = inputPath, sheetName = 'seed_1_allele_counts', row.names=T)
  bulkMutList <- as.vector(bulkDat$mutation_id)
  rownames(bulkDat) <- bulkMutList

  # keep and sort by loci shared between genitypes and bulk data
  sharedMutList <- intersect(genDatMutList, bulkMutList)
  genDat <- genDat[, sharedMutList]
  bulkDat <- bulkDat[sharedMutList, ]

  stopifnot(sharedMutList == rownames(bulkDat))

  # mutCount
  nMut <- length(sharedMutList)
  mutCounts <- matrix(nrow = 2, ncol=nMut)
  mutCounts[1, ] <- bulkDat$ref_counts + bulkDat$var_counts
  mutCounts[2, ] <- bulkDat$var_counts

  # set the copy numbers
  psi <- make.psi.priors(bulkDat, scheme = 'PCN')
  psi <- lapply(psi, function(x) x[[1]])

  # wrap the values in a list
  dataObj <- list()
  dataObj$mutCounts <- mutCounts
  dataObj$psi <- psi
  dataObj$filteredMutMatrix <- mutMatrix

  timeTag <- format(Sys.time(), "%Y-%m-%d-%H-%M-%OS6")
  saveRDS(dataObj, file.path(outputPath, paste0(nameTag, timeTag, '.dat')))

  dataObj
}




make.pyclone.input <- function(mutDat) {
  # mutation_id  ref_counts  var_counts	normal_cn	minor_cn	major_cn
  dat <- data.frame(t(mutDat$mutCounts), stringsAsFactors=F)
  # reference counts should be total counts - variant counts (in mutCounts, first row is d, total number of reads)
  colnames(dat) <- c('total_counts', 'var_counts')
  dat$ref_counts <- dat$total_counts - dat$var_counts
  stopifnot(dat$ref_counts >= 0)
  nMut <- length(mutDat$mutPrevalence)
  # mutationId should be read off the mutDat
  dat$mutation_id <- colnames(mutDat$filteredMutMatrix)

  for (i in seq(nMut)) {
    dat$normal_cn[i] <- genotype.c(mutDat$psi[[i]]$gN)
    # should we use min for minor?
    b <- genotype.b(mutDat$psi[[i]]$gV)
    a <- genotype.a(mutDat$psi[[i]]$gV)
    dat$minor_cn[i] <- min(b, a)
    dat$major_cn[i] <- max(b, a)
  }

  dat
}

relabel.clusters <- function(numericVector) {
  t <- as.factor(x = numericVector)
  levels(t) <- 1:length(levels(t))
  as.numeric(t)
}

# -- Clustering evaluation and prevalence accuracy
# returns -sum(pk * log(pk), axis=0).
# ref: http://pydoc.net/Python/scikit-learn/0.14.1/sklearn.metrics.cluster.supervised/
entropy <- function(labels) {
  if (length(labels) == 0)
    return(1.0)
  pi <- table(labels)
  pi_sum <- sum(pi)
  return(-sum((pi / pi_sum) * (log(pi) - log(pi_sum))))
}

# MI(U,V)=\sum_{i=1}^R \sum_{j=1}^C P(i,j)\log\\frac{P(i,j)}{P(i)P'(j)}
mutual.info.score <- function(trueLabels, predLabels) {
  # contingency table
  contingency <- ftable(trueLabels, predLabels)
  contingency_sum <- sum(contingency)
  pi <- rowSums(contingency) # row
  pj <- colSums(contingency) # column
  outer <-  pi %o% pj
  nnz <- contingency != 0.0
  # normalized contingency
  contingency_nm <- contingency[nnz]
  log_contingency_nm <- log(contingency_nm)
  contingency_nm <- contingency_nm / contingency_sum
  log_outer <- -log(outer[nnz]) + log(sum(pi)) + log(sum(pj))
  mi <- (contingency_nm * (log_contingency_nm - log(contingency_sum))
         + contingency_nm * log_outer)
  sum(mi)
}

# V-measure for clustering
evaluate.clustering <- function(labels_true, labels_pred) {
  if (length(labels_true) == 0)
    return(list(homogeneity=1.0, completeness=1.0, v_measure_score=1.0))

  entropy_C <- entropy(labels_true)
  entropy_K <- entropy(labels_pred)

  MI <- mutual.info.score(labels_true, labels_pred)

  homogeneity <- ifelse(entropy_C != 0, MI / (entropy_C), 1.0)
  completeness <- ifelse(entropy_K != 0, MI / (entropy_K), 1.0)

  if (homogeneity + completeness == 0.0)
    v_measure_score <- 0.0
  else
    v_measure_score <- (2.0 * homogeneity * completeness/ (homogeneity + completeness))

  list(homogeneity=homogeneity, completeness=completeness, v_measure_score=v_measure_score)
}

evaluate.prevalence <- function(truePhi, predPhi) {
  mean(abs(truePhi - predPhi))
}


# Absolute difference of mean MCMC-samples and true value
evaluate.prevalence.trace <- function(state, trace, trueState, burn.in) {
  # trace: matrix[values , mutations]
  trace <- trace[-c(1:burn.in), ]
  means <- colMeans(trace)
  abs(means - trueState)
}



## cached methods

# mS: grid size for precision
# mPhi: grid size for monte carlo estiamte of prevalence
# THIS IS THE NATURAL LOGARITHM OF THE FUNCTION
make.cached.genotype.aware.likelihood <- function(mPhi=1000, dat, mS=10, simulated.data.id, psi.priors, range, tumourContent) {
  fileName <- paste0('l', mPhi, '-', mS, '-', simulated.data.id)
  baseDir <- get.path('likelihoods')
  filePath <- file.path(baseDir, fileName)
  if (file.exists(filePath) && !EmptyCache) {
    readRDS(filePath)
  } else {
    phi <- seq(0, 1, length=mPhi)
    psi <- psi.priors

    # make grid for s
    s <- seq(range$min, range$max, length=mS)
    # for each s, phi in rows, and F in columns over observations [F(phi, obs)] sum(prod(rows[obs]))
    ncol <- nrow(dat)
    dims <- c(mS, mPhi, ncol)

    distMat <- array(0, dim = dims)

    sArray <- array(rep(s, mS), dim=c(mS, mPhi))
    for (i in 1:ncol) { # per mutation
      bArray <- array(dat$var_counts[i], dim=c(mS, mPhi))
      dArray <- array(dat$total_counts[i], dim=c(mS, mPhi))
      for (psi_i in psi[[i]]) {
        eArray <- array(rep(pyclone.eta(phi = phi, psi = psi_i, t = tumourContent), each=mS), dim=c(mS, mPhi))
        distMat[, , i] <- distMat[, , i] + beta.binomial(bArray, dArray, eArray, sArray)
      }
    }

    distMat <- safelog(distMat)

    res <- list(mat=distMat, phi=phi, s=s)
    saveRDS(res, filePath)
    res
  }
}

# returns logged likelihood
cached.pyclone.dd.crp.likelihood <- function(dat, LCACHED) {
  # for now return a monte carlo estimate
  function(obs.indices, precision) {
    discMat <- LCACHED$mat[which.min(abs(precision-LCACHED$s)), , ]

    if (length(obs.indices) == 1) {
      val <- discMat[, as.vector(obs.indices)]
      if (any(val == -Inf)) return(-Inf)
      matrixStats::logSumExp(val) - log(nrow(discMat))
    } else {
      val <- rowSums(discMat[, as.vector(obs.indices)])
      if (any(val == -Inf)) return(-Inf)
      matrixStats::logSumExp(val) - log(nrow(discMat)) # makes NaN
    }
  }
}


# should be done for current value of a, the decay function parameter
cached.resample.alpha <- function(state, a, Decay.CACHED, AlphaCACHED) {
  alpha <- AlphaCACHED$alpha
  l.p.c.Grid <- AlphaCACHED$grid[[which.min(abs(a - Decay.CACHED$a))]]
  K <- sum(state$idx == state$customer)
  prob <- exp(K * safelog(alpha) + l.p.c.Grid)
  sample(alpha, 1, prob=prob)
}


# TODO: Do I need to include the f(d_ij) in the cluster probabilities?
# for now I'm assuming no, because f is not influenced by s
resample.neal.precision <- function(state, cluster.fn, curPrecision, LCACHED) {
  alphaS <- 1.0
  betaS <- 0.0001
  # use the full conditional S according to Neal's paper (prior X F)
  s <- LCACHED$s
  lprior <- dgamma(s, alphaS, betaS, log =T)
  llhood <- rep(0, length(s))

  for (sVal in s) {
    tmpLHood <- 0
    for (i in unique(state$cluster)) {
      obs.indexes <- which(state$cluster == i)
      tmpLHood <- tmpLHood + cluster.fn(obs.indexes, precision = sVal)
    }
    llhood[length(llhood + 1)] <- tmpLHood
  }

  lprob <- lprior + llhood
  sample(x = s, size=1, replace = F, prob = exp(lprob))
}

resample.pyclone.precision <- function(state, cluster.fn, curPrecision) {
  alphaS <- 1.0
  betaS <- 0.0001
  # use the full conditional S according to Neal's paper (prior X F)
  oldS <- curPrecision
  proposalS <- .01
  newS <- sample.gamma.proposal(oldS, proposalS)
  if (is.nan(newS) || is.na(newS)) return(oldS)

  l.oldTarget <- dgamma(oldS, alphaS, betaS, log =T)
  l.newTarget <- dgamma(newS, alphaS, betaS, log =T)

  for (i in unique(state$cluster)) {
    obs.indexes <- which(state$cluster == i)
    # cluster.fn which is the likelihood.fn is already in log space
    l.oldTarget <- l.oldTarget + cluster.fn(obs.indexes, precision = oldS)
    l.newTarget <- l.newTarget + cluster.fn(obs.indexes, precision = newS)
  }

  # x -> x'
  params <- pyclone.gamma.params(oldS, proposalS)
  l.forwardProposal <- dgamma(newS, params$a, params$b, log = T)

  # x' -> x
  params <- pyclone.gamma.params(newS, proposalS)
  l.reverseProposal <- dgamma(oldS, params$a, params$b, log = T)

  # use a simple MH ratio, just newJoint/oldJoint
  MHRatio <- exp(l.newTarget-l.oldTarget)
  if (is.nan(MHRatio) || is.na(MHRatio)) {
    return(oldS)
  } else {
    return(ifelse(rbinom(1, 1, min(MHRatio, 1) ) == 1, newS, oldS))
  }
}

# resamples a
cached.resmple.decay.fn.param <- function(state, alpha, Decay.CACHED) {
  a <- Decay.CACHED$a

  noSelfLoopIndices <- which(state$idx != state$customer)

  theSum <- rep(0, length(a))
  for (n_a in seq(length(a))) {
    for (i in noSelfLoopIndices) {
      theSum[n_a] <- theSum[n_a] + safelog(Decay.CACHED$decays[[n_a]][i, state$customer[i]])
    }
  }

  N <- nrow(Decay.CACHED$decays)
  for (n_a in seq(length(a))) {
    for (i in seq(N)) {
      temp <- ifelse(i > 1, safelog(alpha + sum(Decay.CACHED$decays[[n_a]][i, 1:(i-1)])), safelog(alpha))
      theSum[n_a] <- theSum[n_a] - temp
    }
    # using an exponential prior
    lambda <- 10
    theSum[n_a] <- theSum[n_a] + log(lambda) - lambda * a[n_a]
  }

  # just added - log.sum(theSum), check if it's any better
  sample(a, 1, replace = F, prob = exp(theSum - log.sum(theSum)))
}


# FOR EACH VALUE OF a, the decay function parameter
# Dirichlet process concentration parameter
# returns log.p(alpha) - sum(log(alpha + sum(f(dij))))
# lacks the term for k * log(alpha)
make.cached.alpha <- function(M=1000, distMat, decay.fn, decay.fn.name, simulated.data.id, range, Decay.CACHED) {
  fileName <- paste0('alpha.',  M, '-', simulated.data.id)
  baseDir <- get.path('likelihoods')
  filePath <- file.path(baseDir, fileName)
  if (file.exists(filePath) && !EmptyCache) {
    readRDS(filePath)
  } else {
    # pick an equally distanced grid in (, 10), how about sampling from gamma?
    alpha <- seq(range$min, range$max, length=M)
    #cachedAlpha <- make.cached.decay(M, distMat, decay.fn, decay.fn.name)
    l.p.c <- list()
    # over different values of a (decay.param)
    for (i in seq(length(Decay.CACHED$a))) {
      # set f(i,i)=alpha, then there's one operation less at lesat
      cachedDecay <- Decay.CACHED$decays[[i]] #f(d_ij)
      totalDist <- rowSums(cachedDecay) - diag(cachedDecay) # vector of size N  (sum(f(dij)) over j!=i)
      # should later add k.log.alpha
      xt1 <- safelog(t(sapply(seq(alpha), function(j) {alpha[j] + totalDist})))
      xt <- rowSums(xt1)
      l.p.c[[i]] <- dgamma(alpha, shape=1, rate=0.4, log=T) - xt
    }

    res <- list(grid=l.p.c, alpha=alpha)
    saveRDS(res, filePath)

    res
  }
}


inspect.alpha.prior <- function() {
  y=dgamma(alpha, shape=1, rate=0.4, log=T)
  plot(alpha, exp(y), xlim=c(0,20))
  lines(alpha, y=rep(0, length(alpha)))
}

make.cached.decay <- function(M=1000, distMat, decay.fn, decay.fn.name, simulated.data.id, range) {
  fileName <- paste0('decay.fn.', decay.fn.name,  M, '-', simulated.data.id)
  baseDir <- get.path('likelihoods')
  filePath <- file.path(baseDir, fileName)
  if (file.exists(filePath) && !EmptyCache) {
    readRDS(filePath)
  } else {
    # make the grid
    a <- seq(range$min, range$max, length=M)
    decays <- lapply(a, FUN = function(a_i) { apply(distMat, c(1,2), function(d_ij) (decay.fn(d_ij, a_i)))} )
    res <- list(a=a, decays=decays)
    saveRDS(res, filePath)
    res
  }
}


# use http://cran.r-project.org/web/packages/Rmpfr/vignettes/Rmpfr-pkg.pdf
# to fix underflow
# given current state (table assignments) sample parameters for each cluster
# # ref: Sudderth PhD Thesis: Alg2.1. step 3:  θ^t_k ∼ p(θ_k |{x_i | z^t_i = k},λ)  {page 88}
# why isn't the prior being used? IT"S all based on the likelihood -> because p(theta | lambda) is uniform (uniform prior)
cached.sample.pyclone.cluster.parameters <- function(state, precision, LCACHED) {
  state$phi <- seq(nrow(state))
  probs <- list()
  discMat <- LCACHED$mat[which.min(abs(precision-LCACHED$s)), , ]
  for (clust in unique(state$cluster)) {
    obs.indices <- which(state$cluster == clust)
    if (length(obs.indices) == 1) {
      t <- discMat[, obs.indices]
      prob <- exp(t - matrixStats::logSumExp(t))
      if (any(is.nan(prob))) prob <- 0
      if (all(prob == 0)) {
        print('warning - phi prob NaN or Zero (single)')
      }
    }
    else {
      t <- rowSums(discMat[, obs.indices])
      prob <- exp(t - matrixStats::logSumExp(t))
      if (any(is.nan(prob))) prob <- 0
      if (all(prob == 0)) {
        print(obs.indices)
      }
    }

        if (all(prob == 0))
      state$phi[obs.indices] <- NA
    else
      state$phi[obs.indices] <- sample(LCACHED$phi, 1, replace = F, prob = prob)
  }

    state
}


# mean precision
log.beta.binomial <- function(b, d, m, s) {
  # alpha = sm
  # beta = s(1-m)
  # (d b) B(b + sm, d - b + s(1-m)) / B(sm, s(1-m))
  r <- VGAM::dbetabinom.ab(x=b, size=d, shape1=s*m, shape2=s*(1-m), log=T)
  r
}

beta.binomial <- function(b, d, m, s) {
  # alpha = sm
  # beta = s(1-m)
  # (d b) B(b + sm, d - b + s(1-m)) / B(sm, s(1-m))
  VGAM::dbetabinom.ab(x=b, size=d, shape1=s*m, shape2=s*(1-m), log=F)
}


# --- Other ---

pyclone.gamma.params <- function(data, proposalS) {
  b = data * proposalS
  a = b * data

  list(a=a, b=b)
}


sample.gamma.proposal <- function(mean, proposalS=.01) {
  params <- pyclone.gamma.params(mean, proposalS)
  # shape and rate have to be positive
  rgamma(n=1, params$a, params$b)
}



r.beta.binomail <- function(d, m, s) {
  VGAM::rbetabinom.ab(1, size=d, shape1=s*m, shape2=s*(1-m))
}

# example ddcrp summary functions
ncomp.summary <- function(dat, iter, state, lhood, alpha)
{
  c(iter = iter, ncomp = length(unique(state$cluster)))
}

crp.score.summary <- function(dat, iter, state, lhood, alpha)
{
  c(iter = iter,
    ncomp = length(unique(state$cluster)),
    crp.score = (crp.log.prior(alpha, state$cluster) + sum(lhood)))
}


# returns customers that are connected to i, directly or indirectly, inclusing i
# NOT THOSE THAT i is connected to!!! (3->i<-4, i->5) returns c(3,4) and NOT 5
connections <- function(i, links)
{
  visited <- c()
  to.visit <- c(i)
  while (length(to.visit) > 0) {
    curr <- to.visit[1]
    visited <- c(visited, curr)
    to.visit <- to.visit[-1]
    pointers <- which(links == curr)
    for (p in pointers) {
      if (!(p %in% visited))
        to.visit <- c(to.visit, p)
    }
  }
  visited
}

sample.processor <- function(state, hyperParams, iter) {
  state <- cached.sample.pyclone.cluster.parameters(state, hyperParams$s)
  fileName <- paste0(iter, '-clust.csv')
  write.table(state[, c(2, 6)], file.path(expDir, fileName), sep='\t', col.names=F, row.names=F, append = T)
}


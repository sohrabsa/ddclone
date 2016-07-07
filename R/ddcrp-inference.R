source("R/decay.R")
source("R/helper.R")
source("R/helper-math.R")
source('R/env-setup.R')

library(plyr)
library(Matrix)
library(mcclust)

##
append.traj <- function(phi.traj, clust.traj, hyper.param.traj, expPath) {
  write.table(phi.traj, file.path(expPath, 'phi-trace.csv'))
  write.table(clust.traj, file.path(expPath, 'clust-trace.csv'))
  write.table(hyper.param.traj, file.path(expPath, 'param-trace.tsv'), sep='\t')
}

# ddcrp sampler
ddcrp.gibbs <- function(dat, dist.fn, decay.fn, lhood.fn, summary.fn = ncomp.summary,
                        log.prior.thresh=-10, clust.traj=FALSE, phi.traj=FALSE,
                        hyperParams, MCMCOptions=NULL,
                        permuteCustomers = T, resampleHyperParams = T, expPath,
                        LCACHED, AlphaCACHED, Decay.CACHED)
{
  ### set up summary statistics and trajectories
  if (is.null(MCMCOptions)) {
    MCMCOptions <- list(thinning=1, niter=100, burnIn=10)
  }

  niter <- MCMCOptions$niter
  ndata <- dim(dat)[1]
  msg.inc <- 10^(floor(log10(dim(dat)[1]))-1)
  pb <- txtProgressBar(style=3)

  if (clust.traj) {
    clust.traj <- matrix(NA, nrow=niter, ncol=ndata)
  }
  if (phi.traj) {
    phi.traj <- matrix(NA, nrow=niter, ncol=ndata)
  }

  hyperParamTraj <- matrix(NA, nrow=niter, ncol=3)
  colnames(hyperParamTraj) <- c('a', 'alpha', 's')

  score <- numeric(niter)
  map.score <- 0

  ### set up initial state, summaries, and cluster likelihoods
  msg("setting up the initial state")
  st <- data.frame(idx=1:ndata, cluster=1:ndata, customer=1:ndata)
  lhood <- plyr::daply(st, plyr::.(cluster), function(df) lhood.fn(dat[df$idx,], hyperParams$s))
  summary <- summary.fn(dat, 0, st, lhood, hyperParams$alpha)
  # precompute log.prior
  ## log.prior.lst is indexed first by alpha, then by a
  N <- nrow(Decay.CACHED$decays[[1]])
  n_a <- length(Decay.CACHED$a)
  n_alpha <- length(AlphaCACHED$alpha)
  ldecay <- array(.1, dim=c(n_a, N, N))
  l.prior.mat <- array(.1, dim=c(n_alpha, n_a, N, N))
  l.alpha <- safelog(AlphaCACHED$alpha)
  for (aIndex in 1:n_a) {
    ldecay[aIndex, ,] <- safelog(Decay.CACHED$decays[[aIndex]])
    for (alphaIndex in 1:n_alpha) {
      diag(ldecay[aIndex, ,]) <- l.alpha[alphaIndex]
      diag(Decay.CACHED$decays[[aIndex]]) <- AlphaCACHED$alpha[alphaIndex]
      l.prior.mat[alphaIndex, aIndex, ,] <- ldecay[aIndex, ,] - array(safelog(rowSums(Decay.CACHED$decays[[aIndex]])), dim = c(N, N))
    }
  }
  l.prior.mat[, , ,][which(l.prior.mat[, , ,] < log.prior.thresh)] <- NA

  msg("Done setting up the initial state")


  ### run for niter iterations
  for (iter in 1:niter) {
    iter.score <- 0

    ## set the customer order for this iteration
    permuted <- seq(1, ndata)
    if (permuteCustomers) {
      permuted <- sample(permuted)
    }

    ## run the sitting sampler for each customer
    tableRes <- ddcrp.resample.tables.assignments(permuted, st, lhood, lhood.fn, hyperParams, l.prior.mat, dat, Decay.CACHED, AlphaCACHED)
    iter.score <- tableRes$iter.score
    lhood <- tableRes$lhood
    st <- tableRes$st

    ### update the summary
    iter.score <- iter.score + sum(lhood)
    score[iter] <- iter.score
    if ((score[iter] > map.score) || (iter==1)) {
      map.score <- score[iter]
      map.state <- st
    }
    summary <- rbind(summary, summary.fn(dat, iter, st, lhood, hyperParams$alpha))


    ### sample table parameters
    if (!is.null(dim(clust.traj))) clust.traj[iter,] <- st$cluster
    if (!is.null(dim(phi.traj))) phi.traj[iter,] <- cached.sample.pyclone.cluster.parameters(st, hyperParams$s, LCACHED)$phi

    ## update all hyperparameters
    hyperParamTraj[iter, 'a'] <- hyperParams$a
    hyperParamTraj[iter, 'alpha'] <- hyperParams$alpha
    hyperParamTraj[iter, 's'] <- hyperParams$s

    if (iter %% 20 == 0)
      append.traj(phi.traj, clust.traj, hyperParamTraj, expPath)

    if (resampleHyperParams) {
      shuffle <- c(1,2,3) # 1. a.decay.f, 2. alpha.dirichlet , 3. s.betaBinom
      shuffleList <- sample(shuffle)
      if (length(shuffle) == 1) shuffleList <- shuffle

      for (i in shuffleList) {
        if (i == 1) {
          hyperParams$a <- cached.resmple.decay.fn.param(state=map.state, alpha=hyperParams$alpha, Decay.CACHED)
        } else if (i == 2) {
          hyperParams$alpha <- cached.resample.alpha(state=map.state, a=hyperParams$a, Decay.CACHED, AlphaCACHED)
        } else if (i == 3) {
          hyperParams$s <- resample.neal.precision(state=map.state, cluster.fn = lhood.fn, hyperParams$s, LCACHED)
        }
      }
    }

    setTxtProgressBar(pb, iter/niter)
  }

  close(pb)
  ### return everything
  list(summary=summary, clust.traj=clust.traj, phi.traj=phi.traj, score=score,
       map.score = map.score, map.state = map.state, hyperParams = hyperParams, hyperParamTraj=hyperParamTraj)
}


ddcrp.resample.tables.assignments <- function(customerOrder, st, lhood, lhood.fn, hyperParams, l.prior.mat, dat, Decay.CACHED, AlphaCACHED) {
  iter.score <- 0
  for (i in customerOrder) { # note: index i = 1 is correct at the outset
    ### "remove" the i-th data point from the state
    ### to do this, set its cluster to i, and set its connected data to i
    old.cluster <- st$cluster[i]
    old.customer <- st$customer[i]
    conn.i <- connections(i, st$customer) # all who are connected to i from back and i, NOT THOSE that i is connected to
    st$cluster[conn.i] <- i   # why? splitting the table for customers that remain connected to this customer
    # why set to i? what if there is other cluster called i?
    st$customer[i] <- i # self loop

    if (old.customer == i) {
      lhood[char(old.cluster)] <- 0
    } else if (old.customer != i) {
      old.idx <- st$idx[st$cluster==old.cluster]
      # if old connections are still connected to i, then old.idx will be null and there'll be no need to update likelihoods
      if (!is.null(old.idx) & length(old.idx) != 0) {
        lhood[char(old.cluster)] <- lhood.fn(dat[old.idx,], hyperParams$s)
      }
    }

    ### get log prior
    aIndex <- which.min(abs(hyperParams$a - Decay.CACHED$a))
    alphaIndex <- which.min(abs(hyperParams$alpha - AlphaCACHED$alpha))

    cand.links <- which(!is.na(l.prior.mat[alphaIndex, aIndex, i,]))
    log.prior <- l.prior.mat[alphaIndex, aIndex, i, cand.links]

    ### compute the likelihood of data point i (and its connectors)
    ### with all other tables (!!! do we need to use "idx"?)
    cand.clusts <- unique(st$cluster[cand.links]) # different values for c.j

    # new likelihood for clusters if you were to add link i
    new.lhood <- plyr::daply(subset(st, cluster %in% cand.clusts), plyr::.(cluster),
                       function (df)
                         lhood.fn(dat[unique(c(df$idx,st[conn.i,"idx"])),], hyperParams$s))

    if (length(new.lhood)==1) names(new.lhood) <- cand.clusts

    ### set up the old likelihoods
    old.lhood <- lhood[char(cand.clusts)] # p(k)
    sum.old.lhood <- sum(old.lhood)

    ### compute the conditional distribution
    log.prob <-
      log.prior +
      plyr::laply(cand.links,
            .fun=function (j) {
              # This is just reusing the same likelihood
              sum.old.lhood - old.lhood[char(st$cluster[j])] + new.lhood[char(st$cluster[j])] }) # for each i at a table, update the likelihood

    prob <- exp(log.prob - log.sum(log.prob))
    if (length(prob)==1)
      st$customer[i] <- cand.links[1]
    else
      st$customer[i] <- sample(cand.links, 1, prob=prob)

    ### update the score with the prior and update the clusters
    iter.score <- iter.score + log.prior[which(cand.links == st$customer[i])]
    st$cluster[conn.i] <- st$cluster[st$customer[i]]
    clust.i.idx <- subset(st, cluster == st$cluster[i])$idx
    lhood[char(st$cluster[i])] <- lhood.fn(dat[clust.i.idx,], hyperParams$s)
  } # for i in 1:nCustomers

  list(iter.score=iter.score, lhood=lhood, st=st)
}


# Genotype methods
# psi: genptype string
# t: tumour content
# phi: mutation prevalence (i.e., percentage of cells that have it)
pyclone.eta <- function(phi, psi, t=1) {
  gN <- psi$gN
  gR <- psi$gR
  gV <- psi$gV

  eta <- (1-t) * genotype.c(gN) * genotype.mu(gN) +
    t*(1-phi)*genotype.c(gR)*genotype.mu(gR) +
    t*phi*genotype.c(gV) * genotype.mu(gV)
  Z <- (1-t)*genotype.c(gN) + t*(1-phi)*genotype.c(gR) + t * phi * genotype.c(gV)
  eta/Z
}

noushdarou.eta <- function(Phi, Psi, t=1) {
  gN <- genotypeFromCN(2,0)
  N <- length(Phi)
  Z <- 0
  eta <- 0
  for (i in 1:N) {
    gC <- genotypeFromCN(Psi$ref_counts[i], Psi$alt_counts[i])
    eta <- eta + Phi[i] * genotype.c(gC) * genotype.mu(gC)
    Z <- Z + Phi[i] * genotype.c(gC)
  }
  eta <- (1-t) * genotype.c(gN) * genotype.mu(gN)  + t * eta
  Z <- (1-t) * genotype.c(gN) + t * Z

  ifelse(Z == 0, 0, eta/Z)
}

pyclone.generic.eta <- function(phi) {
  0.001 + 0.499*phi
}

test.genereic.eta <- function() {
  s <- seq(0, 1, length=1000)
  s1 <- pyclone.generic.eta(s)
  s2 <- pyclone.eta(s, generic.psi(), 1)
  s3 <- which(s1 != s2)
  (s1[s3][1]) == (s2[s3][1])

  if (any(s1 != s2))
    stop('Error!')
}

genotype.c <- function(genotype) {
  nchar(genotype)
}

genotype.b <- function(genotype) {
  length(grep('B', strsplit(genotype, '')[[1]]))
}

genotype.mu <- function(genotype) {
  epsilon = .001
  C <- genotype.c(genotype)
  if (C == 0) return(0)
  m <- genotype.b(genotype)/C
  if (m == 0) m <- epsilon
  else if (m == 1) m <- 1 - epsilon
  m
}

genotype.a <- function(genotype) {
  length(grep('A', strsplit(genotype, '')[[1]]))
}

# AB prior
generic.psi <- function() {
  gN <- genotypeFromCN(2, 0)
  gR <- genotypeFromCN(2, 0)
  gV <- genotypeFromCN(1, 1)

  psi <- list()
  psi$gV <- gV
  psi$gN <- gN
  psi$gR <- gR

  psi
}

sample.psi<- function(normalMajorCN=2, normalMinorCN=0) {
  # gV
  maxCN <- 5
  totalCN <- sample(1:maxCN, 1)
  starCN <- sample(0:totalCN, 1)
  c1 <- max(starCN, totalCN - starCN)
  c2 <- totalCN - c1

  gN <- genotypeFromCN(normalMajorCN, normalMinorCN)
  starG <- genotypeFromCN(totalCN, 0)
  gR <- sample(c(gN, starG), 1)

  # gV
  if (gN == gR) {
    gV.c2 <- sample(c(c1, c2), 1)
    gV <- genotypeFromCN(totalCN - gV.c2, gV.c2)
  } else {
    gV <- genotypeFromCN(totalCN - 1, 1)
  }

  psi <- list()
  psi$gV <- gV
  psi$gN <- gN
  psi$gR <- gR

  psi
}

sample.psi.with.gV <- function(varientCNs) {
  variantMajorCN <- varientCNs[1]
  variantMinorCN <- varientCNs[2]
  # gV
  maxCN <- 5
  totalCN <- variantMajorCN + variantMinorCN
  if (totalCN > maxCN) print(paste0('WARNING! totalCN (', totalCN, ') > maxCN (', maxCN,')'))

  if (totalCN == 0) {
    psi <- list()
    psi$gV <- genotypeFromCN(0, 0)
    psi$gN <- genotypeFromCN(2, 0)
    psi$gR <- psi$gN

    return(psi)
  }

  starCN <- sample(0:totalCN, 1)
  c1 <- max(starCN, totalCN - starCN)
  c2 <- totalCN - c1

  gN <- genotypeFromCN(2, 0)

  if (variantMinorCN == 1)
    gR <- genotypeFromCN(totalCN, 0)
  else
    gR <- gN

  # gV
  gV <- genotypeFromCN(variantMajorCN, variantMinorCN)

  psi <- list()
  psi$gV <- gV
  psi$gN <- gN
  psi$gR <- gR

  psi
}

genotypeFromCN <- function(majorCN, minorCN) {
  paste0(c(rep('A', majorCN), rep('B', minorCN)), collapse = '')
}


pick.variant.genotype <- function(cnAtMut, Phi, mode='max') {
  nonZeros <- cnAtMut$ref_counts + cnAtMut$alt_counts != 0
  cnAtMut <- cnAtMut[nonZeros, ]
  if (nrow(cnAtMut) == 0)
    c(0, 0)
  else {
    theRow <- cnAtMut[which.max(Phi[nonZeros]), ]
    c(theRow$ref_counts, theRow$alt_counts)
  }
}

# cpDat is a data.frame with two columns majorCN, and minorCN
# scheme the way priors are ought to be made
# 5 different versions from PyClone + a very vague one from tutorial at
# https://bitbucket.org/aroth85/pyclone/wiki/Priors
make.psi.priors <- function(cpDat, scheme='PCN') {
  res <- list()
  if (scheme == 'PCN') { # Parental Copy Number Prior
    for (i in 1:nrow(cpDat)) {
      if (cpDat$major_cn[i] + cpDat$minor_cn[i] > 1) {
        c_mn <- cpDat$minor_cn[i]
        c_mj <- cpDat$major_cn[i]
        c_t <- c_mn + c_mj
        #if (mn == 0) mn <- 1

        res[[i]] <- list(list(gV=genotypeFromCN(c_mn, c_mj),  gN='AA', gR='AA'),
                         list(gV=genotypeFromCN(c_t - 1, 1),  gN='AA', gR='AA'),
                         list(gV=genotypeFromCN(c_t - 1, 1),  gN='AA', gR=genotypeFromCN(c_t, 0))
                        )

        # Don't add the case where c_minor = 0 as b(gV) = c_minor
        if (c_mn != 0) res[[i]][[4]] <- list(gV=genotypeFromCN(c_mj, c_mn),  gN='AA', gR='AA')

        res[[i]] <- unique(res[[i]])
      } else if (cpDat$major_cn[i] + cpDat$minor_cn[i] == 1) {
        ## If genotype of the variant has only one B, reference would have
        ## sum minor and major letters but all A
        # we're not allowing deletion (gV has to have at least one B)
        res[[i]] <- list(list(gV='B', gN='AA', gR='AA'),
                         list(gV='B', gN='AA', gR='A')) # more likely that first you get mutated and survive  the copy number change
      } else {
        res[[i]] <- list(list(gV=genotypeFromCN(cpDat$major_cn[i], cpDat$minor_cn[i]), gN='AA', gR='AA'))
      }
    }
  } else if (scheme == 'Dollo') { # Dollo?
    for (i in 1:nrow(cpDat)) {
      if (cpDat$major_cn[i] + cpDat$minor_cn[i] > 1) {
        mn <- cpDat$minor_cn[i]
        if (mn == 0) mn <- 1
        res[[i]] <- list(list(gV=genotypeFromCN(cpDat$major_cn[i] + cpDat$minor_cn[i] - 1, mn), gN='AA', gR='AA'),
                         list(gV=genotypeFromCN(cpDat$minor_cn[i], cpDat$major_cn[i]), gN='AA', gR='AA'),
                         list(gV=genotypeFromCN(cpDat$major_cn[i] + cpDat$minor_cn[i] - 1, 1), gN='AA',
                              gR=genotypeFromCN(cpDat$major_cn[i] + cpDat$minor_cn[i], 0)))

        if (cpDat$major_cn[i] + cpDat$minor_cn[i] == 2) {
          res[[i]] <- list(list(gV=genotypeFromCN(cpDat$major_cn[i] + cpDat$minor_cn[i] - 1, mn), gN='AA', gR='AA'),
                           list(gV=genotypeFromCN(cpDat$minor_cn[i], cpDat$major_cn[i]), gN='AA', gR='AA'))
        }
      } else if (cpDat$major_cn[i] + cpDat$minor_cn[i] == 1) {
        res[[i]] <- list(list(gV='B', gN='AA', gR='A')) # more likely that first you get mutated and survive  the copy number change
      } else {
        res[[i]] <- list(list(gV=genotypeFromCN(cpDat$major_cn[i], cpDat$minor_cn[i]), gN='AA', gR='AA'))
      }
    }
  } else if (scheme == 'TCN') {
    for (i in 1:nrow(cpDat)) {
      totalCN <- cpDat$major_cn[i] + cpDat$minor_cn[i]
      gVs <- lapply(1:totalCN, function(i) genotypeFromCN(totalCN - i, i))
      tempList <- list()
      for (j in 1:length(gVs)) {
        tempList[[length(tempList) + 1]] <- list(gV=gVs[[j]], gN='AA', gR='AA')
        tempList[[length(tempList) + 1]] <- list(gV=gVs[[j]], gN='AA', gR=genotypeFromCN(totalCN,0))
      }
      res[[i]] <- tempList

    }
  } else if (scheme == 'AB') {
    for (i in 1:nrow(cpDat)) {
      res[[i]] <- list(list(gV='AB', gN='AA', gR='AA'))
    }
  } else if (scheme == 'BB') {
    for (i in 1:nrow(cpDat)) {
       res[[i]] <-list(list(gV='BB', gN='AA', gR='AA'))
    }
  } else if (scheme == 'NZ') { # no_zygosity
    for (i in 1:nrow(cpDat)) {
      res[[i]] <- list(list(gV=genotypeFromCN(cpDat$major_cn[i] + cpDat$minor_cn[i] - 1, 1), gN='AA', gR='AA'))
    }
  }

  res
}

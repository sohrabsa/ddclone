library(lsa)
## decay functions and distance functions

indexed.window.decay <- function(dat, threshold) {
  function(i) {
    as.integer(dat[i] <= threshold)
  }
}

window.decay <- function(w) {
  function (x) {
    as.integer(x <= w)
  }
}

window.decay.s <- function(x, w) {
  as.integer(x <= w)# || x >=(1-w))
}

logistic.decay.s <- function(x, a, b=.1) {
 logistic((-x+a)/b)
}

logistic.decay <- function(a, b=1)
  function (x) logistic((-x+a)/b)

exp.decay <- function(a, b = 1.0)
  function (x) exp(-x/a) * b + 1e-6

exp.decay.s <- function(x, a, b = 1) {
  exp(-x/a) * b + 1e-6
}


# TODO: change back to 1
identity.decay <- function (x, dummy)
{
  ifelse(x == Inf, 0, 1)
}

euclidean.dist <- function(p1, p2)
{
  v <- p1 - p2
  sqrt(v %*% v)
}

manhattan.dist <- function(p1, p2)
{
  v <- abs((p1[1] - p2[1])) + abs((p1[2] - p2[2]))
}

subtract.dist <- function(p1, p2)
{
  p1 - p1
}

seq.dist <- function(i,j)
{
  if (j <= i)
    i - j
  else
    Inf
}

input.based.dist.fn <- function(input, dist.fn)
{
  function (i,j) { dist.fn(input[i], input[j]) }
}

matrix.dist.fn <- function(dist.matrix)
{
  function (i,j) {
    dist.matrix[i,j]
  }
}

link.dist.fn <- function(adj)
{
  function (i, j) if (adj[i,j]==0) { Inf } else { 1 }
}

jaccardDist <- function(sDat) {
  d <- as.matrix(vegdist(t(sDat$filteredMutMatrix), method='jaccard', na.rm = T))
  if (nrow(sDat$filteredMutMatrix) == 1) {
    d <- as.matrix(vegdist(t(sDat$filteredMutMatrix), method='euclidean', na.rm = F))
  }
  d[which(is.nan(d))] <- 1
  d
}

identity.s <- function(simulatedData) {
  d <- matrix(NA, nrow=ncol(simulatedData$filteredMutMatrix), ncol=ncol(simulatedData$filteredMutMatrix))

  for (i in seq(nrow(d))) {
    for (j in seq(ncol(d))) {
      d[i, j] <- Inf
      if (i == j) d[i, j] <- 2
      else if (i > j) d[i, j] <- i - j
    }
  }
  d
}

cosineDist <- function(simulatedData) {
  as.matrix(cosine(simulatedData$filteredMutMatrix))
}


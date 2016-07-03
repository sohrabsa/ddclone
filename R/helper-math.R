#
# --- mathematics ---
#
#
# safe log function: smoothly handles zero and infinity
exp300 <- exp(-500)

safelog <- function(x) {
  x[x < exp300] <- exp300
  lx <- log(x)
  lx
}

# given log(v), returns log(sum(v))
log.sum <- function(v) {
  log.sum.pair <- function(x,y)
  {
    if ((y == -Inf) && (x == -Inf)) return(-Inf);
    if (y < x) return(x+log(1 + exp(y-x)))
    else return(y+log(1 + exp(x-y)));
  }

  if (length(v) == 1) return(v)
  r <- v[1];
  for (i in 2:length(v))
    r <- log.sum.pair(r, v[i])
  return(r)
}


# harmonic mean
# !!! what's with all these negatives?

harmonic.mean <- function(x)
{
  -(-log(length(x)) + log.sum(-x))
}


# exchangeable dirichlet likelihood

exch.dirichlet.lhood <- function(counts, hyper)
{
  k <- length(counts)
  idx <- counts > 0
  v <- (lgamma(k*hyper) - sum(idx)*lgamma(hyper) +
          sum(lgamma(hyper+counts[idx])) - lgamma(sum(counts[idx])+k*hyper))


  v
}


dirichlet.lhood <- function(counts, alpha)
{
  k <- length(counts)
  idx <- counts > 0
  v <- (lgamma(sum(alpha)) -
          sum(lgamma(alpha[idx])) +
          sum(lgamma(alpha[idx] + counts[idx])) -
          lgamma(sum(counts)+sum(alpha)))

  v
}

# draw from a dirichlet
rdirichlet <- function(n, alpha)
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}


# the logistic function
logistic <- function(x) exp(x) / (1 + exp(x))




# liu's alpha identity

liu <- function(a, n)
{
  sum(sapply(1:n, function (i) a/(a + i - 1)))
}


# --- other ---

msg <- function(s, ...)
{
  time <- format(Sys.time(), "%X")
  cat(sprintf("%s | %s\n", time, s))
}

char <- as.character

colsum <- function (x)
{
  if (is.null(dim(x)) || (dim(x)[1] == 1))
    x
  else
    colSums(x)
}



safelog.old <- function (x) {d
  safelog.f <- function (x) {
    if (x == Inf)
      Inf
    else if (x == 0)
      -100
    else
      log(x)
  }

  if (length(x) == 1)
    safelog.f(x)
  else{
    #apply(xx, seq(dim(xx)), log)
    apply(x, 1:ifelse(is.null(dim(x)), 1, length(dim(x))), log)
  }
}


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



# the logistic function
logistic <- function(x) exp(x) / (1 + exp(x))




# --- other ---

msg <- function(s, ...)
{
  time <- format(Sys.time(), "%X")
  cat(sprintf("%s | %s\n", time, s))
}

char <- as.character




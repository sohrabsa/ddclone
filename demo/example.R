## show case

# library(ddclone)
source('R/ddclone.R')
source('R/helper.R')

###########################
# 1. Simulated Data
###########################
# Run ddClone over simulated dataone")

datObj <- readRDS(system.file("extdata", "dollo.10.48.4.f0.gl0-u.dat", package = "ddclone"))
  ddCloneRes <- ddclone(dataObj = datObj,
                outputPath = './output', tumourContent = 1.0,
                numOfIterations = 100, thinning = 1, burnIn = 1,
                seed = 1)

# Display the result
df <- ddCloneRes$df
expPath <- ddCloneRes$expPath


# Evaluate against the gold standard
dat <- readRDS(dataPath)
nMut <- length(dat$mutPrevalence)
goldStandard <- data.frame(mutID = 1:nMut,
                           clusterID = relabel.clusters(as.vector(dat$mutPrevalence)),
                           phi = as.vector(dat$mutPrevalence))

# Evaluate clustering
(clustScore <- evaluate.clustering(goldStandard$clusterID, df$clusterID))

# Evaluate prevalence estimates
(phiScore <- mean(abs(goldStandard$phi - df$phi)))

# Save the result
score <- data.frame(clustScore, phiMeanError = phiScore)
write.table(score, file.path(expPath, 'result-scores.csv'))




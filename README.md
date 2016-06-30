```{r}
library(ddclone)
```


# 1. Simulated Data

Run ddClone over simulated data
```{r}
dataPath <- './data/dollo.10.48.4.f0.gl0-u.dat'
ddCloneRes <- ddclone(dataPath = dataPath,
              outputPath = './output', tumourContent = 1.0,
              numOfIterations = 100, thinning = 1, burnIn = 1,
              seed = 1)
```

Display the result
```{r}
df <- ddCloneRes$df
expPath <- ddCloneRes$expPath
```

Evaluate against the gold standard
```{r}
dat <- readRDS(dataPath)
nMut <- length(dat$mutPrevalence)
goldStandard <- data.frame(mutID = 1:nMut,
                           clusterID = relabel.clusters(as.vector(dat$mutPrevalence)),
                           phi = as.vector(dat$mutPrevalence))
```

Evaluate clustering
```{r}
(clustScore <- evaluate.clustering(goldStandard$clusterID, df$clusterID))
```

Evaluate prevalence estimates
```{r}
(phiScore <- mean(abs(goldStandard$phi - df$phi)))
```

Save the result
score <- data.frame(clustScore, phiMeanError = phiScore)
write.table(scores, file.path(expPath, 'result-scores.csv'))

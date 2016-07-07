# ddClone: Joint statistical inference of clonal populations from single cell and bulk tumour sequencing data

A statistical framework leveraging data obtained from both single cell and bulk
sequencing strategies. The ddClone approach
is predicated on the notion that single cell sequencing
data will inform and improve clustering of allele fractions
derived from bulk sequencing data in a joint statistical model.
ddClone combines a Bayesian non-parametric prior informed by single cell
data with a likelihood model based on bulk sequencing data to infer
clonal population architecture. Intuitively, the prior encourages genomic
loci with co-occurring mutations in single cells to cluster together.
Using a cell-locus binary matrix from single cell sequencing,
ddClone computes a distance matrix between mutations using the Jaccard distance with exponential decay.
This matrix is then used as a prior for inference over mutation clusters and their prevalences from deeply
sequenced bulk data in a distance-dependent Chinese restaurant process framework.
The output of the model is the most probable set of clonal genotypes present and the
prevalence of each genotype in the population.


# Install the package


# A simple example



# 1. Simulated Data

Load the library:
```{r}
library(ddclone)
```


Run ddClone over simulated data:
```{r}
data()
ddCloneRes <- ddclone(dataPath = dataPath,
              outputPath = './output', tumourContent = 1.0,
              numOfIterations = 100, thinning = 1, burnIn = 1,
              seed = 1)
```

Display the result:
```{r}
df <- ddCloneRes$df
expPath <- ddCloneRes$expPath
```

Evaluate against the gold standard:
```{r}
dat <- readRDS(dataPath)
nMut <- length(dat$mutPrevalence)
goldStandard <- data.frame(mutID = 1:nMut,
                           clusterID = relabel.clusters(as.vector(dat$mutPrevalence)),
                           phi = as.vector(dat$mutPrevalence))
```
Note that in this example the data was packaged in such a way that it contained the gold standard. 


Evaluate clustering
```{r}
(clustScore <- evaluate.clustering(goldStandard$clusterID, df$clusterID))
```

Evaluate prevalence estimates
```{r}
(phiScore <- mean(abs(goldStandard$phi - df$phi)))
```

Save the result
```{r}
score <- data.frame(clustScore, phiMeanError = phiScore)
write.table(score, file.path(expPath, 'result-scores.csv'))
```

# 2. Create a ddclone input object



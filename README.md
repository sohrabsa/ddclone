---
references:
- id: ddcrp
  title: Distance dependent Chinese restaurant processes
  author:
  - family: Frazier
    given: PI
  - family: Blei
    given: DM
  volume: 12
  publisher: Journal of Machine Learning Research
  page: 2461--2488
  type: article-journal
  issued:
    year: 2012
    month: Aug
- id: ddclone
  title: Joint statistical inference of clonal populations from single cell and bulk tumour sequencing data
  author:
  - family: Salehi
    given: Sohrab
  - family: Steif
    given: Adi
  - family: Andrew
    given: Roth
  - family: Aparicio
    given: Samuel
  - family: Bouchard
    given: Alexandre
  - family: Shah
    given: Sohrab P.
  publisher: (submitted)
  type: article-journal
---



=======
# ddClone: Joint statistical inference of clonal populations from single cell and bulk tumour sequencing data

A statistical framework leveraging data obtained from both single cell and bulk sequencing strategies. 
The ddClone (Salehi et al.) approach is predicated on the notion that single cell sequencing
data will inform and improve clustering of allele fractions
derived from bulk sequencing data in a joint statistical model.  
ddClone combines a Bayesian non-parametric prior informed by single cell
data with a likelihood model based on bulk sequencing data to infer
clonal population architecture. Intuitively, the prior encourages genomic
loci with co-occurring mutations in single cells to cluster together.
Using a cell-locus binary matrix from single cell sequencing,
ddClone computes a distance matrix between mutations using the Jaccard distance with exponential decay.
This matrix is then used as a prior for inference over mutation clusters and their prevalences from deeply
sequenced bulk data in a distance-dependent Chinese restaurant process (Frazier and Blei 2012) framework.
The output of the model is the most probable set of mutational clusters present and the
prevalence of each mutation in the population.
The code is based on the ddCRP model, as introduced and implemented in (Frazier and Blei 2012).


## Install the package

An easy way to install ddclone is as follows:

```{r}
library('devtools')
install_github('sohrabsa/ddClone')
```

## A simple example

### 1. Simulated Data

Load the library:
```{r}
library(ddclone)
```

Run ddClone over simulated data:
```{r}
data(gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da)
ddCloneRes <- ddclone(dataObj = gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da,
              outputPath = './output/dollo.0/', tumourContent = 1.0,
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
data(gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da)
nMut <- length(gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da$mutPrevalence)
goldStandard <- data.frame(mutID = 1:nMut,
                           clusterID = relabel.clusters(as.vector(gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da$mutPrevalence)),
                           phi = as.vector(gd.10.48.0.f0.gl0.u.tips.ADO.00.doublet.00.m.50.nus.0.da$mutPrevalence))
```
Note that in this example the data was packaged in such a way that it contained the gold standard. 


Evaluate clustering:
```{r}
(clustScore <- evaluate.clustering(goldStandard$clusterID, df$clusterID))
```

Evaluate prevalence estimates:
```{r}
(phiScore <- mean(abs(goldStandard$phi - df$phi)))
```

Save the result:
```{r}
score <- data.frame(clustScore, phiMeanError = phiScore)
write.table(score, file.path(expPath, 'result-scores.csv'))
```

### 2. Create a ddclone input object
ddClone's input object is a list of 3 elements, `mutCounts`, `psi`, and `filteredMutMatrix`.
We use the simulated data from the Generalized Dollo model:
```{r}
require(xlsx)
intputFilePath <- system.file("extdata", "inputs_simulated.xlsx", package = "ddclone")
```

Read the genotype-mutation matrix:
```{r}
genDat <- read.xlsx(file = intputFilePath, sheetName = 'seed1_genotypes', row.names = T)
genDatMutList <- colnames(genDat)
```

Read the bulk data:
```{r}
bulkDat <- read.xlsx(file = intputFilePath, sheetName = 'seed_1_allele_counts', row.names = T)
bulkMutList <- as.vector(bulkDat$mutation_id)
rownames(bulkDat) <- bulkMutList
```

Generate the ddClone compatible data object:
```{r}
ddCloneInputObj <- make.ddclone.input(bulkDat = bulkDat, genDat = genDat, outputPath = './output/dollo.0/', nameTag = '')
```

Inspect the data object:
```{r}
str(ddCloneInputObj, max.level = 1)
```
Now we can run the analysis similar to sample 1 above.
```{r}
ddCloneRes <- ddclone(dataObj = ddCloneInputObj,
              outputPath = './output/dollo.0/', tumourContent = 1.0,
              numOfIterations = 100, thinning = 1, burnIn = 1,
              seed = 1)
```

# References
Frazier, PI, and DM Blei. 2012. “Distance Dependent Chinese Restaurant Processes” 12 (Aug). Journal of Machine Learning Research: 2461–88.

Salehi, Sohrab, Adi Steif, Roth Andrew, Samuel Aparicio, Alexandre Bouchard, and Sohrab P. Shah. “Joint Statistical Inference of Clonal Populations ￼from Single Cell and Bulk Tumour Sequencing Data.” (submitted).

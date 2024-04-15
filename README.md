# MCANOVA Package

The MC-ANOVA R package provides:
  
  - **MC_ANOVA**: A function to estimate the Relative Accuracy (RA) of cross-ancestry prediction for short chromosome segments.
  - **Maps of the Relative Accuracy** of European-derived local genomic scores for African, Caribbean, East Asian, and South Asian ancestry groups.
  - **A Shiny App** providing a graphical interface to the Relative Accuracy maps. (Note, a slower but URL-acccessible version of the app can be found [here](https://lupia.github.io/Cross-Ancestry-Portability/))
  - **Tools** that, together with `MC_ANOVA()` can be used develop Relative Accuracy maps.

## Installation

To install the development version from GitHub, first install `remotes`: 
```r 
 install.packages("remotes")
 library(remotes)
```
   
Then install the package:
```r
 # detach("package:MCANOVA", unload = TRUE)
 install_github("lupiA/MCANOVA")
 library(MCANOVA)
```
This should take about one minute for a typical computer when including installation of all non-base dependencies.

This package should be compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Monterey, and Linux CentOS 7.

<div id="MENUE" />

 
## Examples

 - [Loading Relative Accuracy maps in an R session](#DATA).
 - [Shiny App](#APP): Launches a Shiny App for the RA Maps we developed using UK Biobank data.
 - [Segments](#SEGMENTS): Finds disjoint chromosome segments.
 - [MCANOVA](#MCANOVA): Estimate Within- and Cross-ancestry R-squared.



<div id="DATA" />
  
### Loading the Relative Accuracy maps into an R session

```r
 library(MCANOVA)

# Relative Accuracy Map (UK Biobank arrays)
 data(MAP_UKB)
```

[Back](#MENUE)

   

<div id="APP" />


### Launching the Shiny App

```r
 library(MCANOVA)
 PGS_portability_app()
```

[Back](#MENUE)



<div id="SEGMENTS" />


### Creating chromosome segments of a minimum base pair length and size (# of SNPs).


```r
 library(MCANOVA)

# Genotype map
 data(geno_map_example)

# Initialize MAP and define segments
 minSNPs <- 10
 minBP <- 10e3
 MAP_example <- geno_map_example
 MAP_example$segments <- getSegments(MAP_example$base_pair_position, chr = MAP_example$chromosome, minBPSize = minBP, minSize = minSNPs, verbose = TRUE)

```


[Back](#MENUE)



<div id="MCANOVA" />


### Running MC-ANOVA 


This example requires the R package [BGData](https://github.com/QuantGen/BGData/tree/master) which is installed along with the MCANOVA package:

```r
# Load necessary packages
# install.packages("BGData")
 library(MCANOVA)
 library(BGData)

# Set seed
 set.seed(12345)

# Generate genotypes (100 subjects and 500 SNPs)
 n <- 100
 p <- 500
 X <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
 data(geno_map_example)
 colnames(X) <- geno_map_example$SNPs
 minSNPs <- 10
 minBP <- 10e3
 MAP_example <- geno_map_example
 MAP_example$segments <- getSegments(MAP_example$base_pair_position, chr = MAP_example$chromosome, minBPSize = minBP, minSize = minSNPs, verbose = TRUE)

# Assign ancestry IDs (80% to ancestry 1, 20% to ancestry 2)
 n_1 <- round(0.8 * n)
 n_2 <- round(0.2 * n)
 ancestry <- rep(c("Group_1", "Group_2"), times = c(n_1, n_2))
 rownames(X) <- ancestry

# Initialize portability estimates
 MAP_example$correlation_within <- NA
 MAP_example$correlation_within_SE <- NA
 MAP_example$correlation_across <- NA
 MAP_example$correlation_across_SE <- NA
 MAP_example$R_squared_within <- NA
 MAP_example$R_squared_across <- NA

# Set parameters for MC-ANOVA

# a small constant added to the diagonals of X'X to avoid numerical errors when some SNPs are in perfect LD
 lambda <- 1e-8
# number of Monte Carlo simulations
 nRep <- 300
# number of causal variants
 nQTL <- 3

# Loop over segments and run MC-ANOVA
 for (i in min(MAP_example$segments):max(MAP_example$segments)) {
   core <- which(MAP_example$segments == i)
   flank_size <- 10
   chunk_start <- max(min(core) - flank_size, 1)
   chunk_end <- min(max(core) + flank_size, nrow(MAP_example))
   chunk <- chunk_start:chunk_end
   isCore <- chunk %in% core
  
   X_1 <- X[rownames(X) == "Group_1", chunk]
   X_2 <- X[rownames(X) == "Group_2", chunk]
  
   # Run MC-ANOVA
   out <- MC_ANOVA(X = X_1, X2 = X_2, core = which(isCore), lambda = lambda, nQTL = nQTL, nRep = nRep)
  
   # Extract portability estimates
   MAP_example$correlation_within[chunk[isCore]] <- out[1, 1]
   MAP_example$correlation_within_SE[chunk[isCore]] <- out[1, 2]
   MAP_example$correlation_across[chunk[isCore]] <- out[2, 1]
   MAP_example$correlation_across_SE[chunk[isCore]] <- out[2, 2]
   MAP_example$R_squared_within[chunk[isCore]] <- out[1, 1]^2
   MAP_example$R_squared_across[chunk[isCore]] <- out[2, 1]^2
 }

 RA <- MAP_example$R_squared_across/MAP_example$R_squared_within
```

[Back](#MENUE)

### Expected outputs

  - **MAP_example**: a 500 x 10 data frame with columns 1-3 containing variant information (chromosome, RS ID, and base pair position), column 4 containing the numeric segment the SNP belongs to estimated from `getSegments()`, columns 5 and 6 containing the within-ancestry group `MC_ANOVA()` correlation estimates and standard errors, columns 7 and 8 containing the across-ancestry group `MC_ANOVA()` correlation estimates and standard errors, and columns 9 and 10 containing the within-ancestry and across-ancestry R-squared `MC_ANOVA()` estimates.
  - **RA**: a length 500 vector of Relative Accuracy estimates.
  - **Interactive shiny app interface.**

All demos should take only a few seconds to run.

### System Requirements

Depends:
    R (>= 3.5.0)
Imports:
  shiny,
  ggplot2,
  BGData

### References

- [Wang et al.(Nat. Comm., 2020)](https://www.nature.com/articles/s41467-020-17719-y) Theoretical and empirical quantification of the accuracy of polygenic scores in ancestry divergent populations. 
- Lupi A., Vazquez A.I., and G. de los Campos (*submitted*). Mapping the Relative Accuracy of Cross-Ancestry Prediction.


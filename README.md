# MCANOVA Package

The MC-ANOVA R package provides:
  
  - **MCANOVA**: A function = to estimate the Relative Accuracy (RA) of cross-ancestry prediction for short chromosome segments.
  - **Maps of the Relative Accuracy** of European-derived local genomic scores for poulations of African, Caribbean, South Asian, and East Asian ancestry,
  - **A Shiny App** providing a graphical interface to the Relative Accuracy maps.
  - **Tools** that, together with `MCANOVA()` can be used develop Relative Accuracy maps.

## Installation

To install the development version from Github:

```r
 # install.packages("remotes")
  remotes::install_github("lupiA/MCANOVA")
  library(MCANOVA)
```


<div id="MENUE" />

 
## Examples

 - [Loading Relative Accuracy maps in an R session](#DATA).
 - [Shiny App](#APP): Launches a Shiny App for the RA Maps we developed using UK-Biobank data.
 - [Segments](#SEGMENTS): Finds disjoint chromosome segments.
 - [MCANOVA](#MCANOVA): Estimate Within- and Cross-ancestry R-squared.



<div id="DATA" />
  
### Loading the Relative Accuracy maps into an R session

```r
 library(MCANOVA)
 data(AF) # African ancestry

```

[Back](#MENUE)

   

<div id="APP" />


### Launcing the Shiny App

```r
 library(MCANOVA)
 PGS_portability_app()
```

[Back](#MENUE)



<div id="SEGMENTS" />


### Creating chromosome segments of a minimum basepair length and size (# of SNPs).


```r
# Genotype map
 path <- system.file("dat", package = "MCANOVA")
 genotype_map <- read.csv(paste0(path, "/geno_map_example.csv"), header = TRUE)


# Initialize MAP and define segments
 minSNPs <- 10
 minBP <- 10e3
 MAP <- genotype_map
 MAP$segments <- getSegments(MAP$base_pair_position, chr = MAP$chromosome, minBPSize = minBP, minSize = minSNPs, verbose = TRUE)

```


[Back](#MENUE)



<div id="S" />


### Running MC-ANOVA 


This example requires the R package [BGData](https://github.com/QuantGen/BGData/tree/master) which is installed along with the MCANOVA package:

```r
# Load necessary packages
library(MCANOVA)
library(BGData)

# Set seed
set.seed(12345)

# Generate genotypes (100 subjects and 500 SNPs)
n <- 100
p <- 500
X <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
colnames(X) <- genotype_map$SNPs

# Assign ancestry IDs (80% to ancestry 1, 20% to ancestry 2)
n_1 <- round(0.8 * n)
n_2 <- round(0.2 * n)
ancestry <- rep(c("Pop_1", "Pop_2"), times = c(n_1, n_2))
rownames(X) <- ancestry

# Initialize portability estimates
MAP$correlation_within <- NA
MAP$correlation_within_SE <- NA
MAP$correlation_across <- NA
MAP$correlation_across_SE <- NA
MAP$R_squared_within <- NA
MAP$R_squared_across <- NA

# Set parameters for MC-ANOVA
lambda <- 1e-8 # a small constant added to the diagonals of X'X to avoid numerical errors when some SNPs are in perfect LD
nRep <- 300 # number of Monte Carlo simulations
nQTL <- 3 # numbre of causal variants

# Loop over segments and run MC-ANOVA
for (i in min(MAP$segments):max(MAP$segments)) {
  core <- which(MAP$segments == i)
  flank_size <- 10
  chunk_start <- max(min(core) - flank_size, 1)
  chunk_end <- min(max(core) + flank_size, nrow(MAP))
  chunk <- chunk_start:chunk_end
  isCore <- chunk %in% core
  
  X_1 <- X[rownames(X) == "Pop_1", chunk]
  X_2 <- X[rownames(X) == "Pop_2", chunk]
  
  # Run MC-ANOVA
  out <- MC_ANOVA(X1=X_1, X2 = X_2, core = which(isCore), lambda = lambda, nQTL = nQTL, nRep = nRep)
  
  # Extract portability estimates
  MAP$correlation_within[chunk[isCore]] <- out[1, 1]
  MAP$correlation_within_SE[chunk[isCore]] <- out[1, 2]
  MAP$correlation_across[chunk[isCore]] <- out[2, 1]
  MAP$correlation_across_SE[chunk[isCore]] <- out[2, 2]
  MAP$R_squared_within[chunk[isCore]] <- out[1, 1]^2
  MAP$R_squared_across[chunk[isCore]] <- out[2, 1]^2
}
```

[Back](#MENUE)


#### References

- [Wang et al.(Nat. Comm., 2020)](https://www.nature.com/articles/s41467-020-17719-y) Theoretical and empirical quantification of the accuracy of polygenic scores in ancestry divergent populations. 
- Lupi A., Vazquez A.I., and G. de los Campos (*submitted*). Mapping the Relative Accuracy of Cross-Ancestry Prediction.


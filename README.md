# MCANOVA Package

# Installation
To install the development version from Github:

```
# install.packages("remotes")
remotes::install_github("lupiA/MCANOVA")
```

# Cross-Ancestry-Portability

MC-ANOVA is an extension of HD-ANOVA [de los Campos et al., 2020](https://pubmed.ncbi.nlm.nih.gov/33315963/) that predicts the portablity of SNP segments in the context of cross-ancestry Polygenic Risk Scores (PGS). The goal is to estimate the extent of genome differentiation with within and across ancestry R-squared. The [MC_ANOVA.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/R/MC-ANOVA.R) function draws genetic values for QTL from a local core of SNPs and then predicts those values using SNPs not in the core or randomly chosen to be QTL. The R-squared is the squared correlation between the generated genetic values and the predicted values.
\
\
We also provide a function, [getSegments.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/R/getSegments.R), to group SNPs into local segments based on a provided Kbp size (e.g., 10 Kbp) and minimum number of SNPs (e.g., 10 SNPs).
\
\
Finally, we have provided an interactive tool, an [R Shiny App](https://github.com/lupiA/Cross-Population-Portability/blob/main/R-shiny-app), in which users can input a single SNP (base pair [BP] position), range of SNPs (BP positions), or a comma-separated list of SNPs, and the App will output portability and marker information. Users are able to download the main output from the App to a .csv file.


### Example
#### Running MC-ANOVA and obtaining portability map
\
After downloading the [MC_ANOVA.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/R/MC-ANOVA.R) and [getSegments.R](https://github.com/lupiA/Cross-Population-Portability/blob/main/R/getSegments.R) functions and the [geno_map_example.csv](https://github.com/lupiA/Cross-Population-Portability/blob/main/geno_map_example.csv) file (requires the R package [BGData](https://github.com/QuantGen/BGData/tree/master)):

```
# Load necessary packages
library(BGData)

# Set seed
set.seed(12345)

# Genotype map
genotype_map <- read.csv("~/geno_map_example.csv", header = TRUE)

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

# Initialize MAP and define segments
minSNPs <- 10
minBP <- 10e3
MAP <- genotype_map
MAP$segments <- getSegments(MAP$base_pair_position, chr = MAP$chromosome, minBPSize = minBP, minSize = minSNPs, verbose = TRUE)

# Initialize portability estimates
MAP$correlation_within <- NA
MAP$correlation_within_SE <- NA
MAP$correlation_across <- NA
MAP$correlation_across_SE <- NA
MAP$R_squared_within <- NA
MAP$R_squared_across <- NA

# Set parameters for MC-ANOVA
lambda <- 1e-8
nRep <- 300
nQTL <- 3

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
  out <- MC_ANOVA(X_1, X2 = X_2, core = which(isCore), lambda = lambda, nQTL = nQTL, nRep = nRep)
  
  # Extract portability estimates
  MAP$correlation_within[chunk[isCore]] <- out[1, 1]
  MAP$correlation_within_SE[chunk[isCore]] <- out[1, 2]
  MAP$correlation_across[chunk[isCore]] <- out[2, 1]
  MAP$correlation_across_SE[chunk[isCore]] <- out[2, 2]
  MAP$R_squared_within[chunk[isCore]] <- out[1, 1]^2
  MAP$R_squared_across[chunk[isCore]] <- out[2, 1]^2
}
```

# Single Cell Informatics Workshop
## Oz Single Cells 2018
By Joseph Powell, Quan Nguyen and Anne Senabouth

Please follow the instructions below to prepare your R environment, before beginning [this vignette](OzSingleCellsWorkshop.md).

## 1. Preparing the R Environment
Feel free to skip some steps if you have already done those steps.

### 1.1 R installation
Please follow the R installation instructions [here](https://mirror.aarnet.edu.au/pub/CRAN/).
If you are a Windows user, make sure you install Rtools. Please note the `ascend`
package requires R version >= 3.4.3. The latest version of R version 3.5.1 is 
best.

### 1.2 R programming environment
The workshop will be done in [RStudio](https://www.rstudio.com/products/rstudio/download/),
but feel free to set up your own R workspace.

### 1.3 Package Installations
You will need to install the following packages to run the development version
of `ascend`. Feel free to skip these steps if you already have these packages.

#### 1.3.1 Packages from CRAN
You can use the install.packages() to install the packages described in this 
section. The pcakages you require from this repository are as follows:

1. [devtools](https://cran.r-project.org/web/packages/devtools/index.html): This 
package will be used to load the development version of `ascend`.
2. [tidyverse](https://www.tidyverse.org/): This is a series of R packages 
for data science and visualisation. This will install packages such as dplyr,
ggplot2 and tidyr.
3. [data.table](https://github.com/Rdatatable/data.table/wiki/Installation):
Please follow the instructions on this page for your operating system.

Remaining packages can be installed as follows:

```{r install_cran, eval = FALSE}
# List of packages to install
cran_packages <- c("reshape2", "fields", "ggbeeswarm", "gridExtra", 
                   "dynamicTreeCut", "dendextend", "RColorBrewer",
                   "locfit", "KernSmooth")

# Easy command to install all at once
install.packages(cran_packages)
```

#### 1.3.2 Packages from Bioconductor
Bioconductor is a repository for R packages  related to the analysis and 
comprehension of high-throughput genomic data. It uses a separate set of 
commands for the installation of packages.

##### 1.3.2.1 Setting up Bioconductor
Use the following code to retrieve the latest installer from Bioconductor.

```{r setup_bioconductor, eval = FALSE}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

You can then install the Bioconductor packages using `biocLite`.

```{r bioconductor_packages, eval = FALSE}
bioconductor_packages <- c("Biobase", "BiocGenerics", "BiocParallel",
                           "SingleCellExperiment", "GenomeInfoDb", "GenomeInfoDbData")

biocLite(bioconductor_packages)
```

##### 1.3.2.2 Differential expression packages
`ascend` provides wrappers for [DESeq](https://bioconductor.org/packages/release/bioc/html/DESeq.html) 
and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), 
so you may choose to add them to your installation. However, we will only be 
using DESeq for the workshop as DESeq2 will require more time than allocated 
for the workshop.

### 1.4 Installing 'ascend' via devtools
As `ascend` is still under development, we will use devtools to install the
package.

```{r install_ascend}
# Load devtools package
library(devtools)

# Use devtools to install the package
install_github("IMB-Computational-Genomics-Lab/ascend", ref = "devel")

# Load the package in R
library(ascend)
```

#### 1.4 Configuring BiocParallel
This package makes extensive use of [BiocParallel](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html), enabling `ascend` to make the most of your computer's hardware. As each system is different, BiocParallel needs to be configured by the user. Here are some example configurations.

#### Unix/Linux/MacOS (Single Machine)
```{r SetupNix, eval = FALSE}
library(BiocParallel)
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores, progressbar=TRUE), default = TRUE)
```

#### Windows (Single Machine - Quad-core system)
The following commands allows Windows to parallelise functions via BiocParallel.
Unlike multicore processing in *nix systems, Snow creates additional R sessions 
to export tasks to. This requires additional computational resources to run and 
manage the tasks.

We recomend you bypass this step if your machine has lower specs.

```{r SetupWin, eval = FALSE}
library(BiocParallel)
workers <- 3 # Number of cores on your machine - 1
register(SnowParam(workers = workers, 
                   type = "SOCK", 
                   progressbar = TRUE), default = TRUE)
```

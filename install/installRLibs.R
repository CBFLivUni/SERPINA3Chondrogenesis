#function to install libraries from Bioconductor/CRAN if not installed already
installPackage <- function(libName) {
  if (libName %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(libName, ask = FALSE)
  }
}

#need to install BiocManager to install all the other libraries
if (!require(BiocManager)){
  install.packages("BiocManager")
}

#read the libraries needed
packagesToInstall <- read.delim("install/rlibs.txt", header = FALSE)
 
#install all the libraries
sapply(packagesToInstall[,1],installPackage)

# install the packages that will be needed to install other packages

BiocManager::install("devtools")

library("BiocManager")
library("devtools")
library("vegan")
library("locfit")

# POTENTIAL ERROR: ARE YOU MISSING "make"? 
# install Rtools for R 4.0 as in https://cran.r-project.org/bin/windows/Rtools/rtools40.html,
# then, if using windows, run this line to add rtools to the system path
# write('PATH="${RTOOLS40_HOME}\\user\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
# Sys.which("make")


                            

# install some of the basic packages
renv::restore(packages = c( "ggplot2",
                            "dplyr",
                            "tibble",
                            "tidyr"))


# install some of the basic packages
renv::restore(packages = c( "metagenomeSeq",
                            "phyloseq",
                            "vegan"))

# install packages separated
BiocManager::install("phyloseq")
BiocManager::install("vegan")
BiocManager::install("EcolUtils") # failed
BiocManager::install("microbiome")
BiocManager::install("metagMisc") # failed
devtools::install("metagMisc") # failed
devtools::install_github("vmikk/metagMisc")
devtools::install_github("GuillemSalazar/EcolUtils")
BiocManager::install("metagenomeSeq")
BiocManager::install("DESeq2")
BiocManager::install("locfit")
BiocManager::install("DelayedArray")




# update all packages
renv::update()

#save renv state
renv::snapshot()

#check status
renv::status()




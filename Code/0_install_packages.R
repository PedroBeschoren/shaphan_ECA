
# install the packages that will be needed to install other packages

BiocManager::install("devtools")

library("BiocManager")
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
BiocManager::install("metagenomeSeq")
BiocManager::install("phyloseq")
BiocManager::install("vegan")
BiocManager::install("EcolUtils") # failed
BiocManager::install("microbiome")
remotes::install_github("vmikk/metagMisc@v.0.0.4")
BiocManager::install("DESeq2")
BiocManager::install("locfit")
BiocManager::install("DelayedArray")




# update all packages
renv::update()

#save renv state
renv::snapshot()

#check status
renv::status()


#pallete
frass_source_y<-"#3b76c5"
frass_source_x<-"#e68613"

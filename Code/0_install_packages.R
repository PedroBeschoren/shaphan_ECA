# install the packages that will be needed to isntall other packages
renv::restore(packages = c("devtools",
                           "BiocManager",
                           "remotes"))

library("BiocManager")
BiocManager::install("metagenomeSeq")

# install Rtools for R 4.0 as in https://cran.r-project.org/bin/windows/Rtools/rtools40.html,
# then, if using windows, run this line to add rtools to the system path
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
Sys.which("make")

renv::restore(packages = c( "rlang"))
                            

install.packages("rlang")
# install some of the basic packages
renv::restore(packages = c( "ggplot2",
                            "dplyr",
                            "tibble",
                            "tidyr"))

BiocManager::install("phyloseq")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tibble")
install.packages("tidyr")


# install some of the basic packages
renv::restore(packages = c( "metagenomeSeq",
                           
                            # "phyloseq",
                            "vegan"))




# install some of the problematic packages
renv::restore(packages = c( "microbiome",
                            "agricolae",
                            "car"))



# install all other packages
renv::restore()



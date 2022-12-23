
#ignore these:
#testing connection to github


# install the packages that will be needed to isntall other packages
renv::restore(packages = c("devtools",
                           "BiocManager",
                           "remotes"))

library("BiocManager")


# POTENTIAL ERROR: ARE YOU MISSING "mkae"? 
# install Rtools for R 4.0 as in https://cran.r-project.org/bin/windows/Rtools/rtools40.html,
# then, if using windows, run this line to add rtools to the system path
# write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
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


# install all other packages
renv::restore()

# update all packages
renv::update()

#save renv state
renv::snapshot()

#check status
renv::status()




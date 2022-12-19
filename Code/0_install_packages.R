# install the packages that will be needed to isntall other packages
renv::restore(packages = c("devtools",
                           "BiocManager",
                           "remotes"))

# install Rtools for R 4.0 as in https://cran.r-project.org/bin/windows/Rtools/rtools40.html, 
# then, if using windowsn, run this line to add rtools to the system path
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)

# install some of the basic packages
renv::restore(packages = c( "ggplot2",
                            "dplyr",
                            "tibble",
                            "tidyr"))


# install some of the basic packages
renv::restore(packages = c( "metagenomeSeq",
                            "phyloseq",
                            "vegan"))


# install some of the problematic packages
renv::restore(packages = c( "metagMisc",
                            "EcolUtils",
                            "MicEco",
                            "R.rsp",
                            "VGAM",
                            "RCurl"))

# install all other packages
renv::restore(packages != "RCurl")

renv::snapshot(packages = "RCurl")

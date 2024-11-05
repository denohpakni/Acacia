install.packages(c("devtools", "ggplot2", "gridExtra", "gtable", "label.switching",
                   "tidyr", "dplyr"), dependencies = T)

install.packages("adegenet")
install.packages("ade4")
install.packages("poppr")
install.packages("ape")
install.packages("BiocManager")
install.packages("hierfstat")
install.packages("iterpc")
install.packages("expm")
install.packages("leaflet")
install.packages("directlabels")

# From BiocManager

BiocManager::install(c("SNPRelate", "qvalue"))

BiocManager::install("ggtree")

install.packages("dartR")



library(directlabels)
library(devtools)
library(dartR)

# test installation
gl.smearplot(testset.gl)


# make sure to get the STRUCTURE software for your computer from https://web.stanford.edu/group/pritchardlab/home.html


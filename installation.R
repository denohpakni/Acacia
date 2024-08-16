install.packages("devtools")


install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "qvalue"))

install.packages("dartR")

library(devtools)
library(dartR)

# test installation
gl.smearplot(testset.gl)


# make sure to get the STRUCTURE software for your computer from https://web.stanford.edu/group/pritchardlab/home.html


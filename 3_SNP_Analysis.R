#####--------SNPs Input--------####

snp1 <- gl.read.dart(filename ="SNPs.csv",
                       ind.metafile ="acacia_pop.csv")

#To keep the original genlight object untouched.

snp2 <- snp1 


###----- 1. Reports for data quality-----####
# Report for Monomorphs

gl.report.monomorphs(snp2)


# Report for Reproducibility

gl.report.reproducibility(snp2)
nLoc(snp1) #31823

# Report for CallRate

gl.report.callrate(snp2, method="loc")

gl.report.callrate(snp2, method="ind")

####---------------2 .Filtering markers for SNPs ------------------####
#Filtering Monomorphs for snps

snp3 <- gl.filter.monomorphs(snp2, v=0)

nLoc(snp3) #31823


#Filtering Reproducibility for snps

snp4 <- gl.filter.reproducibility(snp3, t=0.95)
nLoc(snp4) #31823

#Filtering Callrate per individual for snps

snp5 <- gl.filter.callrate(snp4, 
                              method = "ind", 
                              threshold = 0.10, #cant afford to loose all my samples
                              recalc = T)

nLoc(snp5)
nInd(snp5)


#Filtering Callrate per marker for snps

snp6 <- gl.filter.callrate(snp5,
                              method = "loc",
                              threshold = 0.70,
                              recalc = T)

nLoc(snp6) # 1816
nInd(snp6)

## Filtering by Minor Allele Frequency

gl.report.maf(snp6)


snp7 <- gl.filter.maf(snp6, threshold = 0.05)

# Check No. of loci, sample and population changes
nLoc(snp7)
nInd(snp7)
nPop(snp6)
nPop(snp7)

# Save filtered data
save(snp7, file="snp7_filtered.rdata")


####----- 3. The Next Day ----####

load("snp7_filtered.rdata") 

#Load the libraries "poppr" and "ape"

library(poppr)

library(ape)

####------------------- 4. From genlight to a genind object ------------------####
""" 
  GLOSSARY:
  
  A genind object

  The package adegenet (dartR is based on this package) uses two types of object to storage genetic information: 
  genind is one of them and it is used for individual genotypes.
  Unfortunately the function aboot only allows genind objects to use “Nei distance”, 
  which is what we need to work with SNP markers in general."""

# Converting genlight into genind.

snp7gi <- gl2gi(snp7)

####---------------------5. Distance and NJ tree for SNPs ----------------####
nInd(snp7gi)

snp8gi <- aboot(snp7gi,
                   tree="nj",
                   dist="nei.dist",
                   sample = 100, #number of iteration for bootstrap
                   missing = "mean",
                   showtree=FALSE,
                   root=FALSE)


library(ggtree)
library(ape)


snpinfo <- read.csv("acacia_pop.csv")

ggtree(snp8gi, layout="circular")%<+% snpinfo %>% +
  geom_point(size=1, aes(color=pop))+
  geom_tree(aes(color=pop))+
  geom_tiplab(aes(label=label), color="black",size=1.5)+
  scale_color_manual(values = c("#fc5623",
                                "#B931FC", 
                                "#EE9322", 
                                "#00ff00"),
                     labels = c("kitui","tsavo","naivasha","rongai"))

# save the Rectangular plot 
ggsave("dendrogram_snp_ind.png", width = 20, height = 20)

# save the Circular plot 
ggsave("dendrogram_Circular_snp_ind.png", width = 20, height = 20)


# Comparing with the normal method of NJ

gl.tree.nj(snp7)


####------------------6. Principal Coordinates Analysis - PCoA ----------------####

""" 
  Using nomral genlight data.
  By using PCoA we can visualize individual and population differences.
  PCoA is a method to analyze and visualize similarities or dissimilarities of data
"""

####----Redefine the population information---####

pop(snp7) <- snp7@other$ind.metrics$pop

####----Pcoa---####

install.packages("directlabels")
library(directlabels)

pcoa_snps <- gl.pcoa(snp7, nfactors=10)


gl.pcoa.plot(pcoa_snps,
             snp7,
             xaxis=1,
             yaxis=2)


#### STRUCTURE #####
rm(str1)

sr <- gl.run.structure(snp7, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
ev <- gl.evanno(sr)

qmat <- gl.plot.structure(sr, K=5)

# Draw a Life like Map of population distribution

head(qmat)
gl.map.structure(qmat, K=5, snp7, scalex=1, scaley=0.5)

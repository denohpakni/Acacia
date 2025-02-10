install.packages(c("dartRverse","devtools","hierfstat","adegenet"))

install.packages("ggplot2")

library(dartRverse)
library(hierfstat)


snpgl1 <- gl.read.dart(filename ="SNPs.csv",
                     ind.metafile ="acacia_pop.csv")
plot(snpgl1)

gl.map.interactive(snp5gl)

gl.report.maf(snp5gl)

gl.report.callrate(snpgl1)


#### 2. Filtering ####
snp2gl <- gl.filter.callrate(snpgl1, method = "loc", threshold = 0.8)
nLoc(snpgl1)
nLoc(snp2gl)


#Always run this after removing individuals – removes loci that are no longer variable
snp3gl <- gl.filter.monomorphs(snp2gl)
nLoc(snp2gl)

gl.report.reproducibility(snp3gl)
snp4gl <- gl.filter.reproducibility(snp3gl)

gl.report.rdepth(snp4gl)
snp5gl <- gl.filter.rdepth(snp4gl, lower = 8, upper = 90)
nLoc(snp5gl)

#look at the data to see if you see any obvious issues and redo if you do.
plot(snp5gl)

####---------------------------- 3. Analysis ---------------------------####


  #####----- Relatedness -------#####
# Run EMIBD9 to detect related individuals. Output saved to save time
ibd9 <- gl.run.EMIBD9(snp5gl, Inbreed = FALSE, 
                    emibd9.path = "C:/Program Files (x86)/IoZ_ZSL/EMIBD/")

# save the results
saveRDS(ibd9, file="./ibd9.rds")

#####-- PCOA --#####
pcoa <- gl.pcoa(snp5gl)
gl.pcoa.plot(pcoa, snp5gl)

#####----------- Allelic Richness -----------------#####

# convert the genlight object to genind format
gi <- gl2gi(snp5gl)

#convert genind object to hierfstat format
hfstat <- genind2hierfstat(gi)

#calculate allelic richness
ar <- allelic.richness(hfstat)
names(ar)
summary(ar$Ar) #gives mean AR for each population

ar <- as.data.frame(ar$Ar)
mean.ar <- colMeans(ar)

#first, extend the margins of the graphing window to fit long axis labels
par(mar=c(8,3,3,3))
boxplot(ar, ylab="Allelic richness", las=2)


#####------ heterozygosity for each population #####------ 
""" observed, expected, and unbiased heterozygosities and Fis (inbreeding coefficient) by population """
hets <- gl.report.heterozygosity(snp5gl, method="pop")
ind.hets <- gl.report.heterozygosity(snp5gl, method="ind")


#####------ Genetic differentiation by Fst #####------ 
#Calculate pairwise Fsts between populations.

fsts <- gl.fst.pop(snp5gl, nboots=100, percent=95)
knitr::kable((round(fsts$Fsts,3)))
summary(fsts$Fsts)



# Save  data and dependencies for next day work
save(snp5gl, file="snp5gl_filtered.rdata")
load("snp5gl_filtered.rdata")



#####------------------ Distance and NJ tree for SNPs ----------------#####
library(poppr)

library(ape)

gl.tree.nj(snp5gl)

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


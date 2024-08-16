library(devtools)
library(dartR)

ScoresAcacia <- gl.read.silicodart(
  filename="SilicoDArT.csv",
  ind.metafile="acacia_pop.csv")


snpAcacia <- gl.read.dart( 
  filename="SNPs.csv", 
  ind.metafile="acacia_pop.csv")

####### Boot strap #####
gl.sample(ScoresAcacia)
gl.sample(snpAcacia)


bssnp <- gl.sample(
  snpAcacia,
  nsample = min(table(pop(snpAcacia))),
  replace = TRUE,
  onepop = FALSE,
  verbose = NULL
)

######################## basic statistics ##################
# for each loci (Hs, Ho, Fis etc.)


basicStatSNP <- gl.basic.stats(snpAcacia)
basicStatScore <- gl.basic.stats(ScoresAcacia)


#Calculates the expected heterozygosities for each population in a genlight object
expectedPopHetz <- gl.test.heterozygosity(snpAcacia)

############### observed Heterozygosity per loci #############

observedHetz <- gl.Ho(snpAcacia)
write.csv(observedHetz,file = "./Results/observedHetz.csv")

gl.report.heterozygosity(snpAcacia)


############################# Diversity Analysis ###########################

pca <- gl.pcoa(snpAcacia)
gl.pcoa.plot(pca,snpAcacia)
gl.select.shapes(x = snpAcacia, select = NULL, verbose = NULL)



gl.report.maf(snpAcacia) # MAF for each locus for SNP data

gl.report.diversity(snpAcacia) #diversity indexes for SNPs
gl.report.heterozygosity(snpAcacia) #Estimates Expected Heterozygosity


gl.diagnostics.hwe(snpAcacia)
gl.report.hwe(snpAcacia) # populations with less than 5 individuals are skipped


# The fixation index ( FST) is a measure of population differentiation due to genetic structure. 
gl.fst.pop(snpAcacia)

#### AMOVA ####
# AMOVA is used to detect whether or not there is significant population structure
AMOVA_snp <- gl.amova(snpkig)
AMOVA_snp

AMOVA_Scores <- gl.amova(Scoreskig)
AMOVA_Scores

# PCA Analysis
pcaSnp <- gl.pcoa(
  snpAcacia,
  nfactors = 5,
  correction = NULL,
  mono.rm = TRUE,
  parallel = FALSE,
  n.cores = 16,
  plot.out = TRUE,
  save2tmp = FALSE,
  verbose = NULL
)

snpgl <- snpAcacia
pca <- gl.pcoa(snpgl[1:6,],verbose=2)
gl.pcoa.plot(pca, gl)

gl.grm(snpAcacia) # mean probability of identity by state (IBS)


# Generate a geographical map
gl.map.interactive(snpAcacia)

###### Population STRUCTURE ##########################
gl.tree.nj(ScoresAcacia)
gl.tree.nj(snpAcacia)# nj tree to summarize genetic similarity among populations

## Run for SNP 
sr <- gl.run.structure(snpAcacia, k.range = 2:5, num.k.rep = 3, exec = './structure.exe')
ev <- gl.evanno(sr)
ev
qmat <- gl.plot.structure(sr, K=3)
qmat
head(qmat)
gl.map.structure(qmat, K=3, snpAcacia, scalex=1, scaley=0.5)

################### Hardy-Weinberg tests over loci and populations ############
gl.hwe.pop(snpAcacia)

################### Mantel test ############
mantelTest <- gl.ibd(snpAcacia)

print(mantelTest)


############## Report private alleles in one population compared with a second population ##############
gl.report.pa(snpAcacia)

gl.grm(snpAcacia)

########## Distances #########

gl.dist.pop(snpkig)
gl.dist.ind(snpkig)

########### Distance Matrices #############

gl.report.ld.map(snpkig)

snpdm <- gl.propShared(snpkig) #Calculates a similarity (distance) matrix
write.csv(snpdm,file = "./Results/snpdm.csv")

#### Reports ####


gl.report.rdepth(snpAcacia,save2tmp=TRUE)
gl.report.maf(snpAcacia)
gl.report.parent.offspring(snpAcacia)
gl.report.diversity(snpAcacia)
gl.report.hwe(snpAcacia)
gl.report.ld.map(snpAcacia)
gl.report.callrate(snpAcacia)


utils.n.var.invariant(gl, verbose = NULL)

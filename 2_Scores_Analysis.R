library(devtools)
library(dartR)

sil1 <- gl.read.silicodart(
  filename="SilicoDArT.csv",
  ind.metafile="acacia_pop.csv")


snp1 <- gl.read.dart( 
  filename="SNPs.csv", 
  ind.metafile="acacia_pop.csv")

####### Boot strap #####
gl.sample(sil1)
gl.sample(snp1)

######################## basic statistics ##################
# for each loci (Hs, Ho, Fis etc.)

basicStatSNP <- gl.basic.stats(snp1)
basicStatScore <- gl.basic.stats(sil1)


#save genelight(gl) object when you are done for the day and want to go home
save(sil1, file="sil1.rdata")

# next day 
load("sil1.rdata")

#Display names of for the gl object
names(sil1@other$loc.metrics)
indNames(sil1)

####-------- Monomorphs Report----------####
gl.report.monomorphs(sil1)

####---------CallRate Report----------####
gl.report.callrate(sil1,method="loc")


# for individual samples
gl.report.callrate(sil1,method="ind")

####---------Reproducibility Report---------####

gl.report.reproducibility(sil1)



####---------2. Filtering ---------####

sil2 <- gl.filter.reproducibility(
  sil1,
  threshold = 0.96,
  plot.out = TRUE,
)
#check no. of Loci
nLoc(sil2) #188499

sil3 <- gl.filter.callrate(
  sil2,
  method = "loc",
  threshold = 0.95,
  mono.rm = FALSE,
  recalc = TRUE,
  plot.out = TRUE,
)

nLoc(sil3) #187549
nInd(sil3)


# save and load gl filtered data. For using later
save(sil3,file="sil3filtered.rdata")

load("sil3filtered.rdata")

####------------- Geographic Origin ----------####
# you need lat and long co-ordinates
library(maps)
map_sil3 <- borders(database = "world", 
                     regions = c("kenya","kitui","tsavo","naivasha","rongai"),
                     fill = "transparent", colour = "black")

points_sil3 <- read.csv("acacia_pop.csv")


#Map using ggplot and automatic colors
ggplot() + map_sil3 +
  geom_point(data = points_sil3, 
             aes(x = lon, y = lat, colour = pop)) + # Plot points of location of the samples 
  scale_color_discrete("populations") +  # Change legend title
  xlab("Longitude") +  # X-axis label
  ylab("Latitude")+ # Y-axis label
  theme_light()+
  coord_fixed(ratio = 1) # wont work without actual co-ordinates!

####------------- 3. Genetic distance ----------####
gl.dist.pop(sil3)

Jac_sil3_ind <- gl.dist.ind(sil3, method = "jaccard")

####------------- Clustering ----------####
library(ape)

nj_sil3<- nj(Jac_sil3_ind)


# plot the tree
library(ggtree)

ggtree(nj_sil3)+
  geom_tiplab(aes(label=label), color="red", size=3)

# save the plot 
ggsave("dendrogram_ind_Scores.png", width = 10, height = 20)


#### ---------- 4. Imputing missing data ---------####

sil3_imputed <- gl.impute(sil3)

####------Repeat the process to plot a tree of imputed data --------####

JacSil3Imp_ind <- gl.dist.ind(sil3_imputed, method = "jaccard")

####----Solve clustering---####

njImpSil3 <- nj(JacSil3Imp_ind)

####----Plot imputed data tree---####
ggtree(njImpSil3)+
  geom_tiplab(aes(label=label), color="red", size=3)

####------5. Annotating the tree ------####
library(ggtree)

library(ape)

library(dplyr)
#To annotate the tree, we need to load again your ind.meta file. 
#You can add any additional information you want. 
#We will storage this information as "info".

sil3_info <- read.csv("acacia_pop.csv")


# Plotting the NJ tree with ggtree. Switch between rectangular and circular shapes


ggtree(njImpSil3, layout="circular")%<+% sil3_info %>% +
  geom_point(size=1, aes(color=pop))+
  geom_tree(aes(color=pop))+
  geom_tiplab(aes(label=label), color="black",size=1.5)+
  scale_color_manual(values = c("#fc5623",
                                "#B931FC", 
                                "#EE9322", 
                                "#00ff00"),
                     labels = c("kitui","tsavo","naivasha","rongai"))

# save the plot 
ggsave("dendrogram_imputed_ind_Scores.png", width = 15, height = 20)

# save the Plot
ggsave("dendrogram_imputed_ind_Scores_Circular.png", width = 20, height = 20)


####----------6. Calculate the distance between populations---------------####

##----Euclidean distances between populations---##

sil4Imp_pop <- gl.dist.pop(sil3_imputed, method = "euclidean")


# save imputed  data 
save(sil3_imputed,file="sil3imputed.rdata")

load("sil3imputed.rdata")



####---- 7. Principal Coordinates Analysis - PCoA ---####

# Redefine the population information

sil4_imputed <- sil3_imputed

pop(sil4_imputed) <- sil4_imputed@other$ind.metrics$pop

####----PCoA---####

pcoa_silico <- gl.pcoa(sil4_imputed, nfactors=10)

gl.pcoa.plot(pcoa_silico,
             sil4_imputed,
             scale = FALSE,
             ellipse = FALSE,
             plevel = 0.95,
             pop.labels = "pop")


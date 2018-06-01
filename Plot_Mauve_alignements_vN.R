########################################################################################
###############
############### Visualize Mauve alignment
########## 21 juin 2017
# IM.
########################################################################################

#Script to import Mauve Alignements into R. and produce Synteny plots.



########################################################################################
# Install and opend libraries

#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)

##############################################################################
##############################################################################
# Synteny N16961 different assemblies
##############################################################################

###################
# Large chromosome
###################
setwd('~/Documents/Noemie/NoemiePaper_VF/ComparisonsN16961s/1_MAUVE/')

###############################
# import mauve backbone file
align<-read_mauve_backbone('N16961_Blokesch_chr1.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("This study","Heidelberg et al.","to erase")
# plot alignment
pdf("N16961_comparisons_LargeChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Large chromosome",dna_seg_labels = names)
dev.off()

####################
# Small Chromosome
####################

align<-read_mauve_backbone('N16961_Blokesch_chr2.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
  align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("This study","Heidelberg et al.","to erase")
# plot alignment
pdf("N16961_comparisons_SmallChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Small chromosome",dna_seg_labels = names)
dev.off()
#####################################################################################################################

#######################################
# Synteny N16961 A1552
#######################################

####################
# Large chromosome
####################

setwd('~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/1b_Mauve/')
###############################
# import mauve backbone file
align<-read_mauve_backbone('N16961_A1552_comparison_chr1.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
  align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("N16961","A1552","N16961toerase")
# plot alignment
pdf("N16961A1552_comparison_LargeChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Large chromosome",dna_seg_labels = names)
dev.off()

####################
# Small Chromosome
####################

align<-read_mauve_backbone('N16961_A1552_comparison_chr2.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
  align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("N16961","A1552","N16961Toerase")
# plot alignment
pdf("N16961A1552_comparison_SmallChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Small chromosome",dna_seg_labels = names)
dev.off()
#####################################################################################################################

#######################################
# Synteny A1552 Sa5Y
#######################################

#######################################
# Large chromosome
setwd('~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/1b_Mauve/')

###############################
# import mauve backbone file
align<-read_mauve_backbone('A1552_Sa5Y_comparison_chr1.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
  align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("Sa5Y","A1552","Sa5Y","A1552")
# plot alignment
pdf("1A1552Sa5Y_comparison_LargeChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Large chromosome",dna_seg_labels = names)
dev.off()

####################
# Small Chromosome
####################
align<-read_mauve_backbone('A1552_Sa5Y_comparison_chr2.backbone')
# change colors of blocks
colors<-rainbow(length(align[[1]][[1]]$col))

for (i in align[[1]][[1]]$col) {
  for (ca in 1:length(align[[1]])){
    align[[1]][[ca]]$col[align[[1]][[ca]]$col==i]<- colors[match(i,align[[1]][[ca]]$col)]
  }
}

# change colors of correspondences
for (i in 1:length(align[[2]])){
  align[[2]][[i]]$col<-rep("gray33", length(align[[2]][[i]]$col))
}

names<-c("A1552","Sa5Y","A1552","Sa5Y")
# plot alignment
pdf("A1552Sa5Y_comparison_SmallChromosome.pdf",height = 6,width = 8)
plot_gene_map(align$dna_segs, align$comparisons,main = "Small chromosome",dna_seg_labels = names)
dev.off()



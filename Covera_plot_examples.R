
########################################################################################
###############
############### PLot coverage. 
########## 21 juin 2017

########################################################################################

# Libraries
#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)
library(ggplot2)

#######################################
# Synteny N16961 versions
#######################################

#######################################
# Large chromosome
setwd('~/Documents/Noemie/NoemiePaper_VF/Coverage_plot/')

sizeChr1_A1552<- 3015094
sizeChr2_A1552 <- 1070374
sizeChr1_Sa5Y<-2955400
sizeChr2_Sa5Y<-1096954

filesToProcess <- dir(pattern = "*\\_detail.txt.gz$")  #files to process.
p4<-list()
for (i in filesToProcess) {
  
  cov<-read.delim(gzfile(i), h=F,sep = "\t")


#sampling data to reduce memory. but do differentially for A1552 a lot of sampling and no sampling for Sa5Y
covA1552<-cov[1:(sizeChr1_A1552+sizeChr2_A1552),]
covA1552<-covA1552[sample(nrow(covA1552), 1000000),]
covSa5Y<-cov[(sizeChr1_A1552+sizeChr2_A1552):dim(cov)[1],]
covSa5Y<-covSa5Y[sample(nrow(covSa5Y), 1000000),]

cov<-rbind.data.frame(covA1552,covSa5Y)
summary(cov)
tail(cov)
p4[[i]] <- ggplot(aes(y = V3, x = V2), data = cov) + geom_area(stat="identity", fill="black")+ labs(title = unlist(strsplit(i,split = "_"))[1], x="Position in bp", y="Coverage")+ 
  scale_x_continuous(breaks = seq(0,8000000,by = 250000)) + scale_y_continuous(limits = c(-10, 150))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_rect(data=NULL, aes(xmin=0, xmax=sizeChr1_A1552, ymin=-10, ymax=0), fill='red', alpha=0.5)+
geom_rect(data=NULL, aes(xmin=sizeChr1_A1552, xmax=sizeChr1_A1552+sizeChr2_A1552, ymin=-10, ymax=0), fill='pink', alpha=0.25)+
geom_rect(data=NULL, aes(xmin=sizeChr1_A1552+sizeChr2_A1552, xmax=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y, ymin=-10, ymax=0), fill='green', alpha=0.5)+
geom_rect(data=NULL, aes(xmin=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y, xmax=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y+sizeChr2_Sa5Y, ymin=-10, ymax=0), fill='blue', alpha=0.25)

}

pdf("4703_vN.pdf",useDingbats = F,width = 6.692913*3,height = 6.692913/8*3)
p4[[1]]
dev.off()

pdf("4776_vN.pdf",useDingbats = F,width = 6.692913*3,height = 6.692913/8*3)
p4[[2]]
dev.off()

## Zoom chromosome 2 4703


p4<-list()

filesToProcess<-filesToProcess[grep("4703",filesToProcess)]
for (i in filesToProcess) {
  
  cov<-read.delim(gzfile(i), h=F,sep = "\t")
  
  
  #sampling data to reduce memory. but do differentially for A1552 a lot of sampling and no sampling for Sa5Y
  covA1552<-cov[1:(sizeChr1_A1552+sizeChr2_A1552),]
  covA1552<-covA1552[sample(nrow(covA1552), 1000),]
  covSa5Y<-cov[(sizeChr1_A1552+sizeChr2_A1552):dim(cov)[1],]
  covSa5Y<-covSa5Y[sample(nrow(covSa5Y), 100000),]
  
  cov<-rbind.data.frame(covA1552,covSa5Y)
  summary(cov)
  tail(cov)
  p4[[i]] <- ggplot(aes(y = V3, x = V2), data = cov) + geom_area(stat="identity", fill="black")+ labs(title = unlist(strsplit(i,split = "_"))[1], x="Position in bp", y="Coverage")+ 
    scale_x_continuous(limits=c((sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y),(sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y+sizeChr2_Sa5Y)),breaks = seq(0,8000000,by = 250000)) + scale_y_continuous(limits = c(-10, 150))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_rect(data=NULL, aes(xmin=0, xmax=sizeChr1_A1552, ymin=-10, ymax=0), fill='red', alpha=0.5)+
    geom_rect(data=NULL, aes(xmin=sizeChr1_A1552, xmax=sizeChr1_A1552+sizeChr2_A1552, ymin=-10, ymax=0), fill='pink', alpha=0.25)+
    geom_rect(data=NULL, aes(xmin=sizeChr1_A1552+sizeChr2_A1552, xmax=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y, ymin=-10, ymax=0), fill='green', alpha=0.5)+
    geom_rect(data=NULL, aes(xmin=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y, xmax=sizeChr1_A1552+sizeChr2_A1552+sizeChr1_Sa5Y+sizeChr2_Sa5Y, ymin=-10, ymax=0), fill='blue', alpha=0.25)
  
}

pdf("4703_chr2_detail_vN.pdf",useDingbats = F,width = 6.692913,height = 6.692913*1.3)
p4[[1]]
dev.off()




################################################################################################################################################################################
################################################################################################################################################################################
################################################################################################################################################################################


# Gene homology between strains
################################################################################################################################################################################


setwd('~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/3_GeneHomology/')


filesToProcess <- dir(pattern = "*\\.tsv$")  #files to process.

genes<-list()

for (i in filesToProcess) {
  
  genes[[i]]<-read.table(i, h=T,sep = "\t")
}

genes[[1]]

names(genes)  

mer1<-merge(genes[[3]],genes[[1]],by = c("note","Name"),all.y = T)

?sort
mer1[mer1$Minimum.y==NA]

tail(mer1)



mer1$locus_tag.y
  
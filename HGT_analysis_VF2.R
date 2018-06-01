############################################################################################### 
############################################################################################### 
######## Coverage analysis of HGT event between two Vibrio strains
# IM. v1
#16.mai.2017


####################################################################################
## LIBRARIES

library(pheatmap)
library(GenomicRanges)
#library(GenomeGraphs)
library(ggbio)
library(ggplot2)
library("ggtree")
library(circlize)
library(plyr )
require(scales)
####################################################################################




# Template ref LOAD FIRST
########################
#####################CHANGE AT EACH TREATMENT

template<-read.delim("/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/BamAddRG/Coverage_data/Templates_INFO/Template_C_vf.csv",h=F,sep="\t")

Inserts<-template[grep('blokesch',template$V9),]
Inserts<-Inserts[grep('insert',Inserts$V3),]


Coordinates<-template[grep('source',template$V3),]
Coordinates<-Coordinates[grep('Geneious',Coordinates$V2),]

Coor_LengthA1552<-max(Coordinates[grep('A1552',Coordinates$V1),5])

Coor_LengthA1552_CHR1<-Coordinates[grep('A1552',Coordinates$V1),]
Coor_LengthA1552_CHR1<-Coor_LengthA1552_CHR1[grep('ch1',Coor_LengthA1552_CHR1$V1),5]     

Coor_Length_Chimera<-max(Coordinates$V5)

Coor_LengthA1552_Sa5Y_CHR1<-max(Coordinates[grep('ch1',Coordinates$V1),5])

template2<-template[grep("CDS",template$V3),c(1,4,5,7,9)]


template3<-strsplit(as.character(template2$V9),';')
template3<-lapply(template3,function (x) x[grep("locus_tag",x)])
template3<-lapply(template3,function (x) gsub("locus_tag=","",x))
template3<-unlist(lapply(template3,function (x) ifelse (length(x) == 0, "NA", x)))


template4<-cbind.data.frame(template2[,c(1:4)],template3)
colnames(template4)<-c("chr","start","end","strand","id")

head(template4)
summary(template4)
# plot all results in only one track object only one track
template5<-GRanges(template4)
autoplot(template5,aes(color = strand, fill = strand))
########################################################################################################################################################################
########################################################################################################################################################################





########################################################################################################################################################################
####################################################################################
# Coverage CALCULATIONS VF
####################################################################################

setwd("/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/BamAddRG/Coverage_data/TemplateC1/")
filesToProcess <- dir(pattern = "*\\_detail.txt.gz$")  #files to process.

ALLconfirmedV2<-list()
A1552_in_Sa5Y_FF_VF<-list()

for (i in filesToProcess) {

cov_4703<-read.delim(gzfile(i), h=F,sep = "\t")

# gaps in A1552
###############
temp<-cov_4703[c(1:Coor_LengthA1552),]
## plot raw data
#pdf(paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"A1552_cov.pdf",sep = "_"),width = 12, height = 4)
#barplot(temp$V3)
#abline(h = 10,col="red")
#dev.off()
x <- temp$V3


# all coverage bellow threshold transform it to Zero
x[x<10]<-0
#############

gap_A115<-rle(x)[[1]][rle(x)[[2]]==0]# size of gaps in A155
pos_A115<-cumsum(rle(x)[[1]])
end_A115<-pos_A115[rle(x)[[2]]==0] # end position of gap
start_A115<-end_A115-(gap_A115-1) # start position of gap

# making the dataframe
A115_vf<-cbind.data.frame(paste(rep("Gap",length(gap_A115)),seq(1:length(gap_A115)),sep=""),
                          start_A115,end_A115,gap_A115)
colnames(A115_vf)<-c("name","start","end","length")
strain<-ifelse(A115_vf$end < Coor_LengthA1552, "A1552", "SA5Y")
chrom<-ifelse(A115_vf$end < Coor_LengthA1552_CHR1, "chr1", "chr2")
A115_vf<-cbind.data.frame(A115_vf,strain,chrom)


write.table(A115_vf[A115_vf$length>100,],paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Gaps100.txt",sep = "_"),sep="\t",quote =F)
write.table(A115_vf[A115_vf$length>1000,],paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Gaps1000.txt",sep = "_"),sep="\t",quote =F)

#sum(A115_vf[A115_vf$length>1000,4])
#plot gap measurement
#plot(x)
#points(x=start_A115,y=rep(10,length(start_A115)),col="red",pch=16)
#points(x=end_A115,y=rep(10,length(end_A115)),col="pink",pch=16)


##################################################################################
# Islands in Sa5Y

# gaps in A1552
temp<-cov_4703[c(Coor_LengthA1552:Coor_Length_Chimera),]

# plot raw data
#pdf(paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Sa5Y_cov.pdf",sep = "_"),width = 12, height = 4)
#barplot(temp$V3)
#abline(h = 10,col="red")
#dev.off()
x <- temp$V3

# all coverage bellow threshold transform it to Zero
x[x>9]<-1000
#barplot(x)
x
rle(x)
rle(x)[[1]][rle(x)[[2]]==1000]
isle_SA5Y<-rle(x)[[1]][rle(x)[[2]]==1000]# size of gaps in A155
pos_SA5Y<-cumsum(rle(x)[[1]])
end_SA5Y<-pos_SA5Y[rle(x)[[2]]==1000] # end position of gap
start_SA5Y<-end_SA5Y-(isle_SA5Y-1) # start position of gap

# making the dataframe
SA5Y_vf<-cbind.data.frame(paste(rep("Island",length(isle_SA5Y)),seq(1:length(isle_SA5Y)),sep=""),
                          start_SA5Y+Coor_LengthA1552,end_SA5Y+Coor_LengthA1552,isle_SA5Y)
colnames(SA5Y_vf)<-c("name","start","end","length")
strain<-ifelse(SA5Y_vf$end < Coor_LengthA1552, "A1552", "SA5Y")
chrom<-ifelse(SA5Y_vf$end < Coor_LengthA1552_Sa5Y_CHR1, "chr1", "chr2")
SA5Y_vf<-cbind.data.frame(SA5Y_vf,strain,chrom)


#SA5Y_vf

#SA5Y_vf[SA5Y_vf$length>100,]
write.table(SA5Y_vf[SA5Y_vf$length>100,],paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Islands100.txt",sep = "_"),sep="\t",quote =F)
write.table(SA5Y_vf[SA5Y_vf$length>1000,],paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Islands1000.txt",sep = "_"),sep="\t",quote =F)

#################
# CDS info table#
#################

# data islands
interval_1k<-SA5Y_vf[SA5Y_vf$length>1000,]

# gene info inside interval
genes_island<-list()
for (gene in 1:length(interval_1k$name)) {
  template_inside<-template4[template4$start>interval_1k$start[gene] & template4$end<interval_1k$end[gene] ,] # gene complete inside interval
  template_right<-template4[template4$start<interval_1k$end[gene] & template4$end>interval_1k$end[gene] ,] # end inside gene
  template_left<-template4[template4$start<interval_1k$start[gene]& template4$end>interval_1k$start[gene] ,] # start inside gene
  genes_island[[gene]]<-rbind.data.frame(template_inside,template_right,template_left)
  genes_island[[gene]]<-genes_island[[gene]][!duplicated(genes_island[[gene]]$id),]
  names(genes_island)[gene]<-as.character(interval_1k$name[gene])
}

genes_island<-lapply(genes_island, function(df){  df[order(df$id),]})
gene_is_nam<-lapply(genes_island, function(df){ as.character (df[,5])})
gene_is_nam<-lapply(gene_is_nam,function (df) {paste(df,collapse=";")})


SA5Y_1k_CDS<-cbind.data.frame(SA5Y_vf[SA5Y_vf$length>1000,],unlist(gene_is_nam))
colnames(SA5Y_1k_CDS)[length(colnames(SA5Y_1k_CDS))]<-"CDS"
write.table(SA5Y_1k_CDS,paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Islands1kb_CDS.txt",sep = "_"),sep="\t",quote =F)


### in GAPS
interval_1k<-A115_vf[A115_vf$length>1000,]

genes_gap<-list()
for (gene in 1:length(interval_1k$name)) {
  template_inside<-template4[template4$start>interval_1k$start[gene] & template4$end<interval_1k$end[gene] ,] # gene complete inside interval
  template_right<-template4[template4$start<interval_1k$end[gene] & template4$end>interval_1k$end[gene] ,] # end inside gene
  template_left<-template4[template4$start<interval_1k$start[gene]& template4$end>interval_1k$start[gene] ,] # start inside gene
  genes_gap[[gene]]<-rbind.data.frame(template_inside,template_right,template_left)
  genes_gap[[gene]]<-genes_gap[[gene]][!duplicated(genes_gap[[gene]]$id),]
  names(genes_gap)[gene]<-as.character(interval_1k$name[gene])
}

genes_gap<-lapply(genes_gap, function(df){  df[order(df$id),]})
gene_is_nam<-lapply(genes_gap, function(df){ as.character (df[,5])})
gene_is_nam<-lapply(gene_is_nam,function (df) {paste(df,collapse=";")})

A115_1k_CDS<-cbind.data.frame(A115_vf[A115_vf$length>1000,],unlist(gene_is_nam))
colnames(A115_1k_CDS)[length(colnames(A115_1k_CDS))]<-"CDS"
A115_1k_CDS
write.table(A115_1k_CDS,paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Gaps1kb_CDS.txt",sep = "_"),sep="\t",quote =F)


########################################
# Confirmation of HGT by reciprocity CDS 
# in A1552
CDS_Sa5Y<-unlist(strsplit(gsub(";"," ",
     paste(SA5Y_1k_CDS$CDS,collapse= " ")),
                                            " "))
CDS_Sa5Y<-gsub("Sa5Y-","",CDS_Sa5Y)
CDS_Sa5Y<-CDS_Sa5Y[grep("VC",CDS_Sa5Y)]

list_gaps<-strsplit(gsub("A1552-","",as.character(A115_1k_CDS$CDS)),";")
names(list_gaps)<-A115_1k_CDS$name
list_gaps<-lapply(list_gaps, function(x) x[x %in% CDS_Sa5Y])

A1552_in_Sa5Y<-names(Filter(length, list_gaps))

A1552_in_Sa5Y<-A115_1k_CDS[A115_1k_CDS$name %in% A1552_in_Sa5Y,]

# CONDITIONING FOR NON RECIPROCAL BETWEEN TWO TEMPLATES
if (dim (A1552_in_Sa5Y)[1]==0) {
  A1552_in_Sa5Y<-cbind.data.frame(chrom="NOT_RECIPROCAL",start=9999, end=9999,strand="NA")
} else {
  A1552_in_Sa5Y<-cbind.data.frame(chrom=rep("template",length(A1552_in_Sa5Y$start)),start=A1552_in_Sa5Y$start, end=A1552_in_Sa5Y$end,strand=rep("+",length(A1552_in_Sa5Y$start)))
  
}


## in Sa5Y
CDS_A115<-unlist(strsplit(gsub(";"," ",
                               paste(A115_1k_CDS$CDS,collapse= " ")),
                          " "))
CDS_A115<-gsub("A1552-","",CDS_A115)
CDS_A115<-CDS_A115[grep("VC",CDS_A115)]

list_gaps<-strsplit(gsub("Sa5Y-","",as.character(SA5Y_1k_CDS$CDS)),";")
names(list_gaps)<-SA5Y_1k_CDS$name
list_gaps<-lapply(list_gaps, function(x) x[x %in% CDS_A115])

Sa5Y_in_A1552<-names(Filter(length, list_gaps))

Sa5Y_in_A1552<-SA5Y_1k_CDS[SA5Y_1k_CDS$name %in% Sa5Y_in_A1552,]
# CONDITIONING FOR NON RECIPROCAL BETWEEN TWO TEMPLATES
if ( dim (Sa5Y_in_A1552)[1]==0) {
  Sa5Y_in_A1552<-cbind.data.frame(chrom="NOT_RECIPROCAL",start=9999, end=9999,strand="NA")
} else {
  Sa5Y_in_A1552<-cbind.data.frame(chrom=rep("template",length(Sa5Y_in_A1552$start)),start=Sa5Y_in_A1552$start,end=Sa5Y_in_A1552$end,strand=rep("+",length(Sa5Y_in_A1552$start)))
  
}


# summary
ALLconfirmed<-rbind.data.frame(A1552_in_Sa5Y,Sa5Y_in_A1552)
write.table(ALLconfirmed,paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Confirmed.txt",sep = "_"),sep="\t",quote =F)





##############################
# merge when  inferior to 9999. And produce final table of merges, go until 3 independant events
#A1552
A1552_in_Sa5Y_FF_tmp<-list()
A1552_in_Sa5Y_FF<-list()
A1552_in_Sa5Y_FF_VF<-list()
bool<-list()
for ( ii in 1:length(A1552_in_Sa5Y[,1]))  {
  bool[[ii]]<-A1552_in_Sa5Y[ii+1,2]-A1552_in_Sa5Y[ii,3] < 9999}
bool
# if only 1 event 
if (length(A1552_in_Sa5Y[,1]) ==1 ){ 
  A1552_in_Sa5Y_FF<-list()
  for ( lama in 1) {
    A1552_in_Sa5Y_FF[[lama]]<-cbind.data.frame(chrom=paste(unlist(strsplit(i,"_"))[1],"deletion",sep="_"),start=A1552_in_Sa5Y$start, end=A1552_in_Sa5Y$end, strand=A1552_in_Sa5Y$end-A1552_in_Sa5Y$start)
  }
} else {   
  
  inin<-c(1,(cumsum(rle(unlist(bool))$lengths)[-length(cumsum(rle(unlist(bool))$lengths))]  +1 )  )[-length(cumsum(rle(unlist(bool))$lengths))]
  
  finin<-cumsum(rle(unlist(bool))$lengths) +1[-length(cumsum(rle(unlist(bool))$lengths))]
  
  VALAS<-rle(unlist(bool))$values[-length(cumsum(rle(unlist(bool))$lengths))]
  
  for ( numerito in 1:length(VALAS)){
    tryCatch(    if (VALAS[numerito] == TRUE ){
      A1552_in_Sa5Y_FF_tmp[[numerito]]<-tryCatch(A1552_in_Sa5Y[inin[numerito]:finin[numerito],],error= function(e) e[1])
      A1552_in_Sa5Y_FF[[numerito]]<-cbind.data.frame(chrom=paste(unlist(strsplit(i,"_"))[1],"deletion",sep="_"),start=min(A1552_in_Sa5Y_FF_tmp[[numerito]]$start), end=max(A1552_in_Sa5Y_FF_tmp[[numerito]]$end), 
                                                     strand=max(A1552_in_Sa5Y_FF_tmp[[numerito]]$end)-min(A1552_in_Sa5Y_FF_tmp[[numerito]]$start))
    }else{
      A1552_in_Sa5Y_FF_tmp[[numerito]]<-tryCatch(A1552_in_Sa5Y[inin[numerito]:finin[numerito],],error= function(e) e[1])
      A1552_in_Sa5Y_FF[[numerito]]<-cbind.data.frame(chrom=rep(paste(unlist(strsplit(i,"_"))[1],"deletion",sep="_"),length(A1552_in_Sa5Y_FF_tmp[[numerito]]$start)),start=A1552_in_Sa5Y_FF_tmp[[numerito]]$start, end=A1552_in_Sa5Y_FF_tmp[[numerito]]$end, 
                                                     strand=A1552_in_Sa5Y_FF_tmp[[numerito]]$end-A1552_in_Sa5Y_FF_tmp[[numerito]]$start)
    },error= function(e) e[1])
    
  }
  
}
A1552_in_Sa5Y_FF_VF[[i]]<- ldply(A1552_in_Sa5Y_FF, data.frame)
A1552_in_Sa5Y_FF_VF[[i]]<- A1552_in_Sa5Y_FF_VF[[i]][!duplicated(A1552_in_Sa5Y_FF_VF[[i]]$end),]
A1552_in_Sa5Y_FF_VF[[i]]<- A1552_in_Sa5Y_FF_VF[[i]][rev(rownames(A1552_in_Sa5Y_FF_VF[[i]])),]
A1552_in_Sa5Y_FF_VF[[i]]<-   A1552_in_Sa5Y_FF_VF[[i]][!duplicated(A1552_in_Sa5Y_FF_VF[[i]]$start),]
A1552_in_Sa5Y_FF_VF[[i]]<- A1552_in_Sa5Y_FF_VF[[i]][rev(rownames(A1552_in_Sa5Y_FF_VF[[i]])),]


# Sa5Y
Sa5Y_in_A1552_FF_tmp<-list()
Sa5Y_in_A1552_FF<-list()
Sa5Y_in_A1552_FF_VF<-list()
bool<-list()
for ( ii in 1:length(Sa5Y_in_A1552[,1]))  {
  bool[[ii]]<-Sa5Y_in_A1552[ii+1,2]-Sa5Y_in_A1552[ii,3] < 9999}
bool
# if only 1 event 
if (length(Sa5Y_in_A1552[,1]) ==1 ){ 
  Sa5Y_in_A1552_FF<-list()
  for ( lama in 1) {
    Sa5Y_in_A1552_FF[[lama]]<-cbind.data.frame(chrom=paste(unlist(strsplit(i,"_"))[1],"insertion",sep="_"),start=Sa5Y_in_A1552$start, end=Sa5Y_in_A1552$end, strand=Sa5Y_in_A1552$end-Sa5Y_in_A1552$start)
  }
} else {   
  
  inin<-c(1,(cumsum(rle(unlist(bool))$lengths)[-length(cumsum(rle(unlist(bool))$lengths))]  +1 )  )[-length(cumsum(rle(unlist(bool))$lengths))]
  
  finin<-cumsum(rle(unlist(bool))$lengths) +1[-length(cumsum(rle(unlist(bool))$lengths))]
  
  VALAS<-rle(unlist(bool))$values[-length(cumsum(rle(unlist(bool))$lengths))]
  
  for ( numerito in 1:length(VALAS)){
    tryCatch(    if (VALAS[numerito] == TRUE ){
      Sa5Y_in_A1552_FF_tmp[[numerito]]<-tryCatch(Sa5Y_in_A1552[inin[numerito]:finin[numerito],],error= function(e) e[1])
      Sa5Y_in_A1552_FF[[numerito]]<-cbind.data.frame(chrom=paste(unlist(strsplit(i,"_"))[1],"insertion",sep="_"),start=min(Sa5Y_in_A1552_FF_tmp[[numerito]]$start), end=max(Sa5Y_in_A1552_FF_tmp[[numerito]]$end), 
                                                     strand=max(Sa5Y_in_A1552_FF_tmp[[numerito]]$end)-min(Sa5Y_in_A1552_FF_tmp[[numerito]]$start))
    }else{
      Sa5Y_in_A1552_FF_tmp[[numerito]]<-tryCatch(Sa5Y_in_A1552[inin[numerito]:finin[numerito],],error= function(e) e[1])
      Sa5Y_in_A1552_FF[[numerito]]<-cbind.data.frame(chrom=rep(paste(unlist(strsplit(i,"_"))[1],"insertion",sep="_"),length(Sa5Y_in_A1552_FF_tmp[[numerito]]$start)),start=Sa5Y_in_A1552_FF_tmp[[numerito]]$start, end=Sa5Y_in_A1552_FF_tmp[[numerito]]$end, 
                                                     strand=Sa5Y_in_A1552_FF_tmp[[numerito]]$end-Sa5Y_in_A1552_FF_tmp[[numerito]]$start)
    },error= function(e) e[1])
    
  }
  
}
Sa5Y_in_A1552_FF_VF[[i]]<- ldply(Sa5Y_in_A1552_FF, data.frame)
Sa5Y_in_A1552_FF_VF[[i]]<- Sa5Y_in_A1552_FF_VF[[i]][!duplicated(Sa5Y_in_A1552_FF_VF[[i]]$end),]
Sa5Y_in_A1552_FF_VF[[i]]<- Sa5Y_in_A1552_FF_VF[[i]][rev(rownames(Sa5Y_in_A1552_FF_VF[[i]])),]
Sa5Y_in_A1552_FF_VF[[i]]<-   Sa5Y_in_A1552_FF_VF[[i]][!duplicated(Sa5Y_in_A1552_FF_VF[[i]]$start),]
Sa5Y_in_A1552_FF_VF[[i]]<- Sa5Y_in_A1552_FF_VF[[i]][rev(rownames(Sa5Y_in_A1552_FF_VF[[i]])),]


final_Confirmed<-rbind.data.frame(A1552_in_Sa5Y_FF_VF[[i]],Sa5Y_in_A1552_FF_VF[[i]])
write.table(final_Confirmed,paste(unlist(strsplit(i,"_"))[1],"Confirmed_10kbthreshold_VF.txt",sep = "_"),sep="\t",quote =F, row.names = F)
ALLconfirmedV2[[i]]<-final_Confirmed
}
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################


########################################################################################################################################################################
# Linear PLOT by treatment
########################################################################################################################################################################
# 1 do first template begiingin script only

# ATTENTION: Issue with delimitation of Chr2 and A1552 Sa5Y. Some big framents acquired by HGT make that the total chr size became bigger. 
# these could be discarded at separating the data on A1552, Sa5y, Chr1 Chr2 first  step..
# ONLY figure that do not have this issue is when plotting everything together.
# Not allways accurate, Cassette vertical line, division of chromosomes and strains
# Issues with next to end template of chr2 Sa5Y
# Issues with next to start template chr1 Sa5y
# manually curated figures Sa5y chrI and Sa5Y chr2

wd<-"/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/BamAddRG/Coverage_data/TemplateC1/"
setwd(wd)
cassette<-unlist(strsplit(wd,split="/"))[grep("Template",unlist(strsplit(wd,split="/")) )]

# import files
filesToProcess <- dir(pattern = "*\\_VF.txt$")  #files to pro# if event 3 merged
table(unlist(lapply(strsplit(filesToProcess,"_"), function (x) x[4])))
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = T),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA")))
names(listOfFiles)<-gsub("_Confirmed_10kbthreshold_VF.txt","",gsub(".xml","",gsub("Blast_","",filesToProcess))   )


listOfFiles[[2]]$Gene
names(listOfFiles)

HGT<-list()
for (i in 1:length(listOfFiles)){

HGT[[i]]<- cbind.data.frame(listOfFiles[[i]],ymin=rep(i,dim(listOfFiles[[i]])[1])-0.5,ymax=rep(i,dim(listOfFiles[[i]])[1])+0.5  )  
listOfFiles[[i]]$chrom<-as.character(listOfFiles[[i]]$chrom)
listOfFiles[[i]][grep("deletion",listOfFiles[[i]]$chrom),1]<-"darkred"
listOfFiles[[i]][grep("insertion",listOfFiles[[i]]$chrom),1]<-"darkblue"
HGT[[i]]<- cbind.data.frame(HGT[[i]],Colorito=listOfFiles[[i]]$chrom)
names(HGT)[i]<-names(listOfFiles)[i]
}
HGT_vf<-rbind.data.frame(HGT[[1]],HGT[[2]])
for (i in 3:length(HGT)){
HGT_vf<-rbind.data.frame(HGT_vf,HGT[[i]])}
#
#
#
#
#

# PLOT ALL CHIMERA
plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(Coor_Length_Chimera,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -100000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse()

plotlines + geom_rect(aes(xmin=HGT_vf$start, xmax=HGT_vf$end, ymin=HGT_vf$ymin, ymax=HGT_vf$ymax),data=NULL,
                      fill=HGT_vf$Colorito, alpha=2/4)  + geom_vline(xintercept=Coor_LengthA1552)  +
  annotate("text", x = 400000, y = 1-  2, label = "A1552 template",size=2)   +
  annotate("text", x = 5000000, y = 1 - 2, label = "Sa5Y template",size=2)   +
  annotate("text", x = Inserts$V4[1], y = length(names(listOfFiles)) +1, label = "Insert", size =2) +
  theme(axis.line=element_blank(),axis.text.y=element_blank(), axis.ticks=element_blank(),axis.title.x=element_blank(),
         axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma) 

ggsave(paste("Results",cassette,"VF2.pdf",sep = "_"), plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))




# PLOT A1552 template
HGT_vf_A1552<-HGT_vf[HGT_vf$end<(Coor_LengthA1552),]
HGT_vf_A1552$Colorito<-as.character(HGT_vf_A1552$Colorito)
HGT_vf_A1552$Colorito<-"black"

plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(Coor_LengthA1552,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -100000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse() 

plotlines + geom_rect(aes(xmin=HGT_vf_A1552$start, xmax=HGT_vf_A1552$end, ymin=HGT_vf_A1552$ymin, ymax=HGT_vf_A1552$ymax),data=NULL,
                      fill=HGT_vf_A1552$Colorito, alpha=2/4)  + geom_vline(xintercept=Coor_LengthA1552_CHR1) +
  annotate("text", x = 400000, y = 1 - 2, label = "Chromosome I", size=2)   +
  annotate("text", x = 3300000, y = 1 - 2, label = "Chromosome II", size=2)   +
  annotate("text", x = Inserts$V4[1]-Coor_LengthA1552, y = length(names(listOfFiles))+1, label = "Insert", size=2) +
    theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma)

ggsave(paste("Results",cassette,"A1552","VF2.pdf",sep = "_"), plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))

# PLOT Sa5Y template
HGT_vf_Sa5Y<-HGT_vf[HGT_vf$start>Coor_LengthA1552,]
HGT_vf_Sa5Y$Colorito<-as.character(HGT_vf_Sa5Y$Colorito)
HGT_vf_Sa5Y$Colorito<-"black"
plotlines<-  ggplot(aes(xmin=rep(Coor_LengthA1552,length(listOfFiles)), xmax=rep(Coor_Length_Chimera,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = 3850000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse()

plotlines + geom_rect(aes(xmin=HGT_vf_Sa5Y$start, xmax=HGT_vf_Sa5Y$end, ymin=HGT_vf_Sa5Y$ymin, ymax=HGT_vf_Sa5Y$ymax),data=NULL,
                      fill=HGT_vf_Sa5Y$Colorito, alpha=2/4)  + geom_vline(xintercept=Coor_LengthA1552_Sa5Y_CHR1) +
  annotate("text", x = 4400000, y = 1 - 2, label = "Chromosome I", size=2)   +
  annotate("text", x = 7300000, y = 1 - 2, label = "Chromosome II", size=2)   +
  annotate("text", x = Inserts$V4[1], y = length(names(listOfFiles)) +1, label = "Insert", size=2) +
  theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma)

ggsave(paste("Results",cassette,"SA5Y","VF2.pdf",sep = "_"), plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))

###
#Coordinates selection antibio_A1552
VCA0747_A1552<-691610
VC2338_A1552 <-2263614
VCA0107_A1552<-117330

cassette2<-c(VCA0107_A1552,VCA0747_A1552,VCA0747_A1552,VCA0747_A1552,VC2338_A1552,VC2338_A1552,VC2338_A1552,VCA0107_A1552)
names(cassette2)<- c("Template1","Template2","TemplateA","TemplateB","TemplateC1","TemplateC2","TemplateD","TemplateE")

###
# Ploting separately CHR1 CHR2 for A1552
sizeChr1_A1552 <- 3015094
sizeChr2_A1552 <- 1070374

HGT_vf_A1552<-HGT_vf[HGT_vf$end<Coor_LengthA1552,]
HGT_vf_A1552$Colorito<-as.character(HGT_vf_A1552$Colorito)
HGT_vf_A1552$Colorito<-"black"
HGT_vf_A1552_LargeChr<-HGT_vf_A1552[HGT_vf_A1552$start< ( sizeChr1_A1552 + 1),]
HGT_vf_A1552_SmallChr<-HGT_vf_A1552[HGT_vf_A1552$start> ( sizeChr1_A1552 ),]
HGT_vf_A1552_SmallChr$start<-HGT_vf_A1552_SmallChr$start-sizeChr1_A1552
HGT_vf_A1552_SmallChr$end<-HGT_vf_A1552_SmallChr$end-sizeChr1_A1552

#large chromosome
plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(sizeChr1_A1552,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -100000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse()

plotlines + geom_rect(aes(xmin=HGT_vf_A1552_LargeChr$start, xmax=HGT_vf_A1552_LargeChr$end, ymin=HGT_vf_A1552_LargeChr$ymin, ymax=HGT_vf_A1552_LargeChr$ymax),data=NULL,
                      fill=HGT_vf_A1552_LargeChr$Colorito, alpha=2/4) +
  annotate("text", x = 400000, y = 1-2, label = "Large chromosome", size=2)   +
  theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma, limits =c(-100000,sizeChr1_A1552+20000))+ 
  geom_vline(xintercept=cassette2[cassette[cassette %in% names(cassette2)]])
ggsave(paste("Results",cassette,"A1552","LargeChromosome","VF2.pdf",sep = "_") , plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))


#Small chromosome
plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(sizeChr2_A1552,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -33000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse()

plotlines + geom_rect(aes(xmin=HGT_vf_A1552_SmallChr$start, xmax=HGT_vf_A1552_SmallChr$end, ymin=HGT_vf_A1552_SmallChr$ymin, ymax=HGT_vf_A1552_SmallChr$ymax),data=NULL,
                      fill=HGT_vf_A1552_SmallChr$Colorito, alpha=2/4) +
  annotate("text", x = 120000, y = 1- 2, label = "Small chromosome", size=2)   +
  theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma, limits =c(-33000,sizeChr2_A1552))+ 
  geom_vline(xintercept=cassette2[cassette[cassette %in% names(cassette2)]])

ggsave(paste("Results",cassette,"A1552","SmallChromosome","VF2.pdf",sep = "_"), plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))

###

###
#Coordinates selection antibio_Sa5Y
VCA0747_Sa5Y<-721129
VC2338_Sa5Y <-591308
VCA0107_Sa5Y<-119142

cassette2<-c(VCA0107_Sa5Y,VCA0747_Sa5Y,VCA0747_Sa5Y,VCA0747_Sa5Y,VC2338_Sa5Y,VC2338_Sa5Y,VC2338_Sa5Y,VCA0107_Sa5Y)
names(cassette2)<- c("Template1","Template2","TemplateA","TemplateB","TemplateC1","TemplateC2","TemplateD","TemplateE")

# Ploting separately CHR1 CHR2 for Sa5Y
sizeChr1_A1552<- 3015094
sizeChr2_A1552 <- 1070374
sizeChr1_Sa5Y <- 2955400
sizeChr2_Sa5Y <- 1095478


HGT_vf_Sa5Y<-HGT_vf[HGT_vf$end>Coor_LengthA1552,]
HGT_vf_Sa5Y$Colorito<-as.character(HGT_vf_Sa5Y$Colorito)
HGT_vf_Sa5Y$Colorito<-"black"

HGT_vf_Sa5Y_LargeChr<-HGT_vf_Sa5Y[HGT_vf_Sa5Y$start< (sizeChr1_A1552 +sizeChr2_A1552+sizeChr1_Sa5Y + 1),]
HGT_vf_Sa5Y_LargeChr$start<-HGT_vf_Sa5Y_LargeChr$start- sizeChr1_A1552 - sizeChr2_A1552
HGT_vf_Sa5Y_LargeChr$end<-HGT_vf_Sa5Y_LargeChr$end-sizeChr1_A1552 - sizeChr2_A1552 

HGT_vf_Sa5Y_SmallChr<-HGT_vf_Sa5Y[HGT_vf_Sa5Y$start> ( sizeChr1_A1552 +sizeChr2_A1552+sizeChr1_Sa5Y ),]
HGT_vf_Sa5Y_SmallChr$start<-HGT_vf_Sa5Y_SmallChr$start-sizeChr1_Sa5Y - sizeChr1_A1552 - sizeChr2_A1552
HGT_vf_Sa5Y_SmallChr$end<-HGT_vf_Sa5Y_SmallChr$end-sizeChr1_Sa5Y -sizeChr1_A1552 - sizeChr2_A1552 

#large chromosome
plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(sizeChr1_Sa5Y,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -100000, y = (1:length(names(listOfFiles))), label = names(listOfFiles),size =2)+ scale_y_reverse()

plotlines + geom_rect(aes(xmin=HGT_vf_Sa5Y_LargeChr$start, xmax=HGT_vf_Sa5Y_LargeChr$end, ymin=HGT_vf_Sa5Y_LargeChr$ymin, ymax=HGT_vf_Sa5Y_LargeChr$ymax),data=NULL,
                      fill=HGT_vf_Sa5Y_LargeChr$Colorito, alpha=2/4) +
  annotate("text", x = 400000, y = 1-2, label = "Large chromosome", size=2)   +
  theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma, limits =c(-100000,sizeChr1_Sa5Y))+ 
  geom_vline(xintercept=cassette2[cassette[cassette %in% names(cassette2)]])

ggsave( paste("Results",cassette,"Sa5Y","LargeChromosome","_sVF2.pdf",sep = "_") , plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
        width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))


#Small chromosome
plotlines<-  ggplot(aes(xmin=rep(1,length(listOfFiles)), xmax=rep(sizeChr2_Sa5Y,length(listOfFiles)), ymin=1:length(listOfFiles)-0.05, ymax= 1:length(listOfFiles)+0.05),data=NULL)+ geom_rect() + 
  annotate("text", x = -33000, y = (1:length(names(listOfFiles))),  label = names(listOfFiles),size =2)+ scale_y_reverse() 

plotlines + geom_rect(aes(xmin=HGT_vf_Sa5Y_SmallChr$start, xmax=HGT_vf_Sa5Y_SmallChr$end, ymin=HGT_vf_Sa5Y_SmallChr$ymin, ymax=HGT_vf_Sa5Y_SmallChr$ymax),data=NULL,
                      fill=HGT_vf_Sa5Y_SmallChr$Colorito, alpha=2/4) +
  annotate("text", x = 120000, y = 1- 2, label = "Small chromosome", size=2)   +
  theme(axis.line=element_blank(), axis.text.y=element_blank(),      axis.ticks=element_blank(),axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",  panel.border=element_blank() ) + scale_x_continuous(labels = comma, limits =c(-33000,sizeChr2_Sa5Y+500))+ 
  geom_vline(xintercept=cassette2[cassette[cassette %in% names(cassette2)]])

ggsave(paste("Results",cassette,"Sa5Y","SmallChromosome","_sVF2.pdf",sep = "_"), plot = last_plot(), path = "/media/imateus/IvanHD2_Ubuntu/Noemie_HGT_Analysis/Figures_VF2/",
       width = 106, height = 72.5, units = c("mm"),
       dpi = 300, limitsize = TRUE) + theme_set(theme_gray(base_size = 6))
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################




 # scripts not used
########################################################################################################################################################################
################## CONTROL STATISTICS NOEMIE, HGT events

setwd("~/Documents/V_cholerae/all_data/refVCA0107/")

AHGT<-read.table("ALL_HGT_vf.txt",h=F)

plot(AHGT[AHGT[,2]<4019277,4],AHGT[AHGT[,2]>4019277,4], xlab="HGT event in A1552", ylab="HGT event in Sa5Y", pch=16)
hist(AHGT[AHGT[,2]<4019277,4],breaks=30, col="black", xlab="Total length HGT", main="")

NOEMIE<-c("75.234","95.786","76.492","129.153","38.295","90.916","73.366","104.938","63.66","96.463","23.138","65.836","74.207","57.691",
          "89.668","35.662","103.776","142.511","34.03","81.253","12.308","116.171","73.68","60.228","38.361","26.306","99.865","85.804",
          "21.768","124.038","42.773","74.343","71.487","39.684","28.528","10.82","49.103","17.083","47.719","44.915","168.241","23.559",
          "100.092","47.549","95.89","56.84","140.235","33.55","68.975","161.271","65.791","81.765","73.121","24.078","58.881","30.162",
          "24.342","145.152","62.483","71.829")
AHGT[AHGT[,2]<4019277,]
plot(NOEMIE,tapply(AHGT[AHGT[,2]<4019277,4], AHGT[AHGT[,2]<4019277,1], sum),  ylab="HGT event calculation", pch=16)

hist(tapply(AHGT[AHGT[,2]<4019277,4], AHGT[AHGT[,2]<4019277,1], length), xlab="Independent events",col="black", main="")

AHGT_5kb_a<-read.table("ALL_EXP1_VC0107_5kblim_vf.txt",h=F)
AHGT_5kb<-tapply(AHGT_5kb_a[,4], AHGT_5kb_a[,1], sum)

AHGT_10kb_a<-read.table("ALL_EXP1_VC0107_10kblim_vf.txt",h=F)
AHGT_10kb<-tapply(AHGT_10kb_a[,4], AHGT_10kb_a[,1], sum)
tapply(AHGT_5kb_a[,4], AHGT_5kb_a[,1], length)
pdf("~/Documents/Noemie/Summary_results/Controls_HGTcompa_v1.pdf", height =6 ,width=8)
par(mfrow=c(2,2),mar=c(5,7,4,4))
plot(tapply(AHGT_10kb_a[,4], AHGT_10kb_a[,1], length),tapply(AHGT_5kb_a[,4], AHGT_5kb_a[,1], length), pch=16, cex=1.5, 
     ylab="HGT min. sep. distance 5kb", xlab="HGT min. sep. distance 10kb", las=1,cex.axis=1.5, cex.lab=1.5,
     main="Nb. HGT event by strain")
plot(AHGT_5kb,AHGT_10kb, pch=16, cex=1.5, ylab="HGT min. sep. distance 10kb", xlab="HGT min. sep. distance 5kb",las=1,cex.axis=1.5, cex.lab=1.5,
     main="Total HGT event by strain")
plot(NOEMIE,AHGT_5kb[grep("ins",names(AHGT_5kb))], pch=16, cex=1.5, ylab="HGT min. sep. distance 5kb", las=1,xlab="Noemie manual estimation (kb)",
     main="Total HGT event by strain",cex.axis=1.5, cex.lab=1.5)
plot(NOEMIE,AHGT_10kb[grep("ins",names(AHGT_10kb))], pch=16, cex=1.5, ylab="HGT min. sep. distance 10kb",las=1, xlab="Noemie manual estimation (kb)",
     main="Total HGT event by strain",cex.axis=1.5, cex.lab=1.5)
dev.off()

pdf("~/Documents/Noemie/Summary_results/HGT_hist_v1.pdf", height =3 ,width=8)
par(mfrow=c(1,3),mar=c(5,7,4,4))
hist(tapply(AHGT_10kb_a[grep("ins",AHGT_10kb_a[,1]),4], AHGT_10kb_a[grep("ins",AHGT_10kb_a[,1]),1], length), col="black",
     xlab="Nb. HGT events by strain", las=1,cex.axis=1.5, cex.lab=1.5,
     main="HGT min. sep. distance 10kb")
hist(tapply(AHGT_10kb_a[grep("ins",AHGT_10kb_a[,1]),4], AHGT_10kb_a[grep("ins",AHGT_10kb_a[,1]),1], sum), col="black",breaks=30,
                                             xlab="Total length HGT events by strain", las=1,cex.axis=1.5, cex.lab=1.5,
                                             main="HGT min. sep. distance 10kb")
hist(AHGT_10kb_a[grep("ins",AHGT_10kb_a[,1]),4], col="black",breaks=30,
     xlab="Length HGT events by strain", las=1,cex.axis=1.5, cex.lab=1.5,
     main="HGT min. sep. distance 10kb")
dev.off()

# control Reciprocity HGT 

df<-read.table("~/Documents/Noemie/Experiment2/Extra_samples_template_E/Extra_samples_templateE_v1.txt",h=F)
df
plot(df[grep("ins",df$V1),4],df[grep("del",df$V1),4], pch=16, cex=1.5, ylab="Deletions in Sa5y",las=1, xlab="Insertions in A1552",
     main="Total HGT event by strain",cex.axis=1.5, cex.lab=1.5)
df<-read.table("~/Documents/Noemie/Experiment2/template_A/ALL_EXP2_templateA_v1.txt",h=F)
df
plot(df[grep("ins",df$V1),4],df[grep("del",df$V1),4], pch=16, cex=1.5, ylab="Deletions in Sa5y",las=1, xlab="Insertions in A1552",
     main="Total HGT event by strain",cex.axis=1.5, cex.lab=1.5)


########################################################################################################################################################################





########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################




### import VCF info

tmp.vcf<-readLines("4703_S6_R1_Trimm_Pair_sortedQ30RG.vcf")
tmp.vcf.data<-read.table("4703_S6_R1_Trimm_Pair_sortedQ30RG.vcf")
# filter for the columns names
tmp.vcf<-tmp.vcf[-(grep("#CHROM",tmp.vcf)+1):-(length(tmp.vcf))]
vcf.names<-unlist(strsplit(tmp.vcf[length(tmp.vcf)],"\t"))
names(tmp.vcf.data)<-vcf.names

tmp.vcf.data<-tmp.vcf.data[grep("snp",tmp.vcf.data$INFO),1:5]
tmp.vcf.data[150:500,]

# 
tmp.vcf<-readLines("~/Documents/Microsynth/References_single_strains/db_strains/Sa5Y_in_A1552.vcf")
tmp.vcf.data2<-read.table("~/Documents/Microsynth/References_single_strains/db_strains/Sa5Y_in_A1552.vcf")
# filter for the columns names
tmp.vcf<-tmp.vcf[-(grep("#CHROM",tmp.vcf)+1):-(length(tmp.vcf))]
vcf.names<-unlist(strsplit(tmp.vcf[length(tmp.vcf)],"\t"))
names(tmp.vcf.data2)<-vcf.names

tmp.vcf.data2<-tmp.vcf.data2[grep("snp",tmp.vcf.data2$INFO),1:5]
tmp.vcf.data2[150:500,]


####################################################################################
####################################################################################

####################################################################################
### merge fragment when between fragments is less than 1kb.
####################################################################################

test<-head(A115_1k_CDS,20)
test<-A115_1k_CDS

#1. merge consecutive rows and produce dataframe
new_tab<-list()
for (st in 2:dim(test)[1])
ifelse( (test[st,2] - test[st-1,3]  ) < 1000, 
        new_tab[[st]]<-c(paste(test[st-1,1],test[st,1],sep = "_"),c(test[st,2] - test[st-1,3]  ) ,
                   test[st-1,2],test[st-1,3],test[st,2],test[st,3],
                   (test[st,3]-test[st-1,2]) ,
                   as.character(test[st,5]),as.character(test[st,6])), 
        "sup")
new_tab

new_tab2<-data.frame(matrix(unlist(new_tab),ncol = 9,byrow = T))
colnames(new_tab2)<-c("id","length_espace","start_fg1","end_fg1","start_fg2","end_fg2","length","strain","chrom")
new_tab2
test
#2. rows that espace higher than 1000 kb, put in separate. first fragment
CACA<-list()
for (i in 1:length(new_tab2$length_espace)) {
if(as.numeric(levels(new_tab2$length_espace))[new_tab2$length_espace][i] >1000) {
   CACA[[i]]<-c(unlist(strsplit(as.character(new_tab2$id),split="_")[i])[1],
     "NA",
     as.numeric(levels(new_tab2$start_fg1))[new_tab2$start_fg1][i],
     as.numeric(levels(new_tab2$end_fg1))[new_tab2$end_fg1][i],
     "NA",
     "NA",
     as.numeric(levels(new_tab2$end_fg1))[new_tab2$end_fg1][i] - as.numeric(levels(new_tab2$start_fg1))[new_tab2$start_fg1][i],
     as.vector(new_tab2$strain[i]),
     as.vector(new_tab2$chrom[i]))} else
 {   "NA"}
}

CACA

#3. rows that espace higher than 1000 kb, put in separate. Second fragment
CACA2<-list()
for (i in 1:length(new_tab2$length_espace)) {
  if(as.numeric(levels(new_tab2$length_espace))[new_tab2$length_espace][i] >1000) {
    CACA2[[i]]<-c(unlist(strsplit(as.character(new_tab2$id),split="_")[i])[2],
                 "NA",
                 as.numeric(levels(new_tab2$start_fg2))[new_tab2$start_fg2][i],
                 as.numeric(levels(new_tab2$end_fg2))[new_tab2$end_fg2][i],
                 "NA",
                 "NA",
                 as.numeric(levels(new_tab2$end_fg2))[new_tab2$end_fg2][i] - as.numeric(levels(new_tab2$start_fg2))[new_tab2$start_fg2][i],
                 as.vector(new_tab2$strain[i]),
                 as.vector(new_tab2$chrom[i]))} else
                 {   "NA"}
}

CACA2

#4 Reasemble, first merge without badmerge and individual not merged.

CACA3<-rbind.data.frame(data.frame(matrix(unlist(CACA),ncol = 9,byrow = T)),
data.frame(matrix(unlist(CACA2),ncol = 9,byrow = T)))
colnames(CACA3)<-c("id","length_espace","start_fg1","end_fg1","start_fg2","end_fg2","length","strain","chrom")


new_tab3<-new_tab2[as.numeric(levels(new_tab2$length_espace))[new_tab2$length_espace]<1000,]
new_tab4<-rbind.data.frame(new_tab3,CACA3[!as.character(CACA3$id) %in% unlist(strsplit(as.character(new_tab3$id),split = "_")),])
# put into numeric
new_tab4$end_fg1<-as.numeric(levels(new_tab4$end_fg1))[new_tab4$end_fg1]
new_tab4$start_fg1<-as.numeric(levels(new_tab4$start_fg1))[new_tab4$start_fg1]
new_tab4$end_fg2<-as.numeric(levels(new_tab4$end_fg2))[new_tab4$end_fg2]
new_tab4$start_fg2<-as.numeric(levels(new_tab4$start_fg2))[new_tab4$start_fg2]

new_tab5<-new_tab4[order(new_tab4$end_fg1),]
new_tab5<-new_tab5[!duplicated(new_tab5$id),]

#5 Reduce coordinates at fg1 start end , only 1 pair
new_tab5
red_tb<-list()
for (st in 1:dim(new_tab5)[1]) {
if(!is.na(new_tab5$end_fg2[st])) {
       red_tb[[st]]<-c(id=as.character(new_tab5$id[st]),
                       start_bd1=new_tab5$start_fg1[st],end_bd1=new_tab5$end_fg2[st],
                       length_bound=(new_tab5$end_fg2[st]- new_tab5$start_fg1[st]),
                       part= length(new_tab5[st,3:6][!is.na(new_tab5[st,3:6])])/2, 
                       strain=new_tab5$strain[st],chrom=new_tab5$chrom[st])} 
  else {
       red_tb[[st]]<-c(id=as.character(new_tab5$id[st]),
                       start_bd1=new_tab5$start_fg1[st],end_bd1=new_tab5$end_fg1[st],
                       length_bound=(new_tab5$end_fg1[st]- new_tab5$start_fg1[st]),
                       part= length(new_tab5[st,3:6][!is.na(new_tab5[st,3:6])])/2, 
                       strain=new_tab5$strain[st],chrom=new_tab5$chrom[st])}
}
red_tb<-data.frame(matrix(unlist(red_tb),ncol = 7,byrow = T))
colnames(red_tb)<-c("id_merged","start_bd","end_bd","lenght_boundary","nb_parts","strain","chrm")

# re merge gaps

second_round<-list()
for (st in 2:dim(red_tb)[1])
  ifelse( (red_tb[st,2] - red_tb[st-1,3]  ) < 1000, 
          second_round[[st]]<-c(id=paste(red_tb[st-1,1],red_tb[st,1],sep = "_"),length_boundaries=c(red_tb[st,2] - red_tb[st-1,3]  ) ,
                           start_bd1_fg1=red_tb[st-1,2],end_bd1_fg1=red_tb[st-1,3],
                           start_bd2_fg2=red_tb[st,2],end_bd2_fg2=red_tb[st,3],
                           nb_parts=red_tb$nb_parts[st],
                           as.character(red_tb[st,6]),as.character(red_tb[st,7])), 
          "sup")
second_round

second_round2<-data.frame(matrix(unlist(second_round),ncol = 9,byrow = T))
colnames(second_round2)<-c("id","length_espace","start_fg1","end_fg1","start_fg2","end_fg2","length","strain","chrom")
second_round2
new_tab5






#############################################################################################
#############################################################################################3
# CDS in interval.

### only for each file, need to incorporate into loop for all files

template<-read.delim("~/Documents/V_cholerae/test1/ref_VCA0107/reference_template/All_chromosomes-Sa5Y-VCA0107-frt-kan-frt.gff",h=F,sep="\t")

template2<-template[grep("cds",template$V3),c(1,4,5,7,9)]


template3<-strsplit(as.character(template2$V9),';')
template3<-lapply(template3,function (x) x[grep("locus_tag",x)])
template3<-lapply(template3,function (x) gsub("locus_tag=","",x))
template3<-unlist(lapply(template3,function (x) ifelse (length(x) == 0, "NA", x)))


template4<-cbind.data.frame(template2[,c(1:4)],template3)
colnames(template4)<-c("chr","start","end","strand","id")

tail(template4)

# data islands
interval_1k<-SA5Y_vf[SA5Y_vf$length>1000,]

# gene info inside interval
genes_island<-list()
for (gene in 1:length(interval_1k$name)) {
template_inside<-template4[template4$start>interval_1k$start[gene] & template4$end<interval_1k$end[gene] ,] # gene complete inside interval
template_right<-template4[template4$start<interval_1k$end[gene] & template4$end>interval_1k$end[gene] ,] # end inside gene
template_left<-template4[template4$start<interval_1k$start[gene]& template4$end>interval_1k$start[gene] ,] # start inside gene
genes_island[[gene]]<-rbind.data.frame(template_inside,template_right,template_left)
genes_island[[gene]]<-genes_island[[gene]][!duplicated(genes_island[[gene]][all_genes_inside$id,]),]
names(genes_island)[gene]<-as.character(interval_1k$name[gene])
}
genes_island<-lapply(genes_island, function(df){  df[order(df$id),]})
gene_is_nam<-lapply(genes_island, function(df){ as.character (df[,5])})
gene_is_nam<-lapply(gene_is_nam,function (df) {paste(df,collapse=";")})


SA5Y_1k_CDS<-cbind.data.frame(SA5Y_vf[SA5Y_vf$length>1000,],unlist(gene_is_nam))
colnames(SA5Y_1k_CDS)[length(colnames(SA5Y_1k_CDS))]<-"CDS"
write.table(SA5Y_1k_CDS,paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Islands1kb_CDS.txt",sep = "_"),sep="\t",quote =F)


### in GAPS

interval_1k<-A115_vf[A115_vf$length>1000,]

genes_gap<-list()
for (gene in 1:length(interval_1k$name)) {
  template_inside<-template4[template4$start>interval_1k$start[gene] & template4$end<interval_1k$end[gene] ,] # gene complete inside interval
  template_right<-template4[template4$start<interval_1k$end[gene] & template4$end>interval_1k$end[gene] ,] # end inside gene
  template_left<-template4[template4$start<interval_1k$start[gene]& template4$end>interval_1k$start[gene] ,] # start inside gene
  genes_gap[[gene]]<-rbind.data.frame(template_inside,template_right,template_left)
  genes_gap[[gene]]<-genes_gap[[gene]][!duplicated(genes_gap[[gene]][all_genes_inside$id,]),]
  names(genes_gap)[gene]<-as.character(interval_1k$name[gene])
}
genes_gap<-lapply(genes_gap, function(df){  df[order(df$id),]})
gene_is_nam<-lapply(genes_gap, function(df){ as.character (df[,5])})
gene_is_nam<-lapply(gene_is_nam,function (df) {paste(df,collapse=";")})

A115_1k_CDS<-cbind.data.frame(A115_vf[A115_vf$length>1000,],unlist(gene_is_nam))
colnames(A115_1k_CDS)[length(colnames(A115_1k_CDS))]<-"CDS"
A115_1k_CDS
write.table(A115_1k_CDS,paste(gsub("sorted","",gsub("\\RG.bam_detail\\.txt","",i)),"Gaps1kb_CDS.txt",sep = "_"),sep="\t",quote =F)
########################################################################################################################################################################
########################################################################################################################################################################



####################################################################################
############
############ multiple samples circos
####################################################################################

###### 1.put template coordinates and outside ring
pdf("~/Documents/V_cholerae/all_data/refVCA0107/4703_circos.pdf",height =90,width=90)

# template
templat = data.frame(name = c("template"),
                     start = c(1),
                     end= c(8057838))
circos.genomicInitialize(templat)

### add rectangle of chromosomes and genomes
chrom= data.frame(name = c("template","template","template","template","template","template"),
                  start = c(1,2995339,4019277,6960891,7075966,7077442),
                  end= c(2995338,4019276,6960890,7075965,7077441,8057838),
                  value=c(1,0,1,0,1,0))

circos.par("track.height" = 0.05)

circos.genomicTrackPlotRegion(chrom, ylim =c(0, 1),
                              panel.fun =function(region, value, ...) {
                                i =getI(...)
                                circos.genomicRect(region, value, ytop = i + 0.1  , ybottom = i - 1.1 ,
                                                   col =c("black", "darkgrey","coral4","coral2","Blue","coral2"), ...)
                              })
######### FINAL HGT
setwd("~/Documents/Noemie/Experiment2/template_A")
filesToHGT <- dir(pattern = "*\\_Confirmed_VF.txt$")  #files to pro# if event 3 merged
if (bool[2] == TRUE ){
  A1552_in_Sa5Y_FF1<-tryCatch(ACT_A1552_in_Sa5Y[(2-1):2,],error= function(e) e[1])
  A1552_in_Sa5Y_FF1<-cbind.data.frame(chrom=paste(unlist(strsplit(i,"_"))[1],"tSa5Y",sep="_"),start=min(A1552_in_Sa5Y_FF1$start), end=max(A1552_in_Sa5Y_FF1$end), 
                                      strand=max(A1552_in_Sa5Y_FF1$end)-min(A1552_in_Sa5Y_FF1$start))
  ACT_A1552_in_Sa5Y<- rbind.data.frame(A1552_in_Sa5Y_FF1,ACT_A1552_in_Sa5Y[c(-1,-2),])
} else {
  A1552_in_Sa5Y_FF1<-tryCatch(ACT_A1552_in_Sa5Y[(2-1):2,],error= function(e) e[1])
  A1552_in_Sa5Y_FF1<-cbind.data.frame(chrom=rep(paste(unlist(strsplit(i,"_"))[1],"tSa5Y",sep="_"),length(A1552_in_Sa5Y_FF1[,1])),start=A1552_in_Sa5Y_FF1$start, end=A1552_in_Sa5Y_FF1$end, 
                                      strand=A1552_in_Sa5Y_FF1$end-A1552_in_Sa5Y_FF1$start)
  ACT_A1552_in_Sa5Y<- rbind.data.frame(A1552_in_Sa5Y_FF1,ACT_A1552_in_Sa5Y[c(-1,-2),])
}cess.

list_HGT<-list()
for (n in filesToHGT) {
  HGT<-read.delim(n, h=T,sep = "\t")
  HGT<-cbind.data.frame(rep("template",length(HGT[,2])) , HGT[,2:3] , rep("+",length(HGT[,2])))
  colnames(HGT)<- c("chrom","start", "end","strand")
  list_HGT[[n]]<-HGT
}


posTransform.fun = function(region) {
  return(region)
}

### Iniitalize with template
circos.genomicInitialize(templat)

### configure tracks
circos.par("track.height" = 0.05)

### chromosomes tracks
circos.genomicTrackPlotRegion(chrom, ylim =c(0, 1),
                              panel.fun =function(region, value, ...) {
                                i =getI(...)
                                circos.genomicRect(region, value, ytop = i + 0.1  , ybottom = i - 1.1 ,
                                                   col =c("gray35", "indianred","grey","hotpink","Blue","hotpink"), ...)
                              })


circos.genomicTrackPlotRegion(list_HGT,ylim =c(0, 1), stack = T, panel.fun = function(region, value, ...) {
  i =getI(...)
  circos.genomicRect(region, value, col = "black", ytop = i + 0.3  , ybottom = i - 1.3 ,  border = "black", ...)
}, track.height = 0.8, bg.border= 0.2)


circos.clear()

###################################################################################33333

for (i in filesToCov) {
  cov1<-read.delim(i, h=F,sep = "\t")
#cov1<-read.delim("~/Documents/V_cholerae/test1/ref_VCA0107/Bowtie2_bef/4703_detail.txt", h=F,sep = "\t")
head(cov1)
colnames(cov1)<-c("chr","start","cov")
cov1<-cov1[cov1$cov>10,]
covo<-cbind.data.frame(chr=rep("template",length=dim(cov1[cov1$cov>10,])[1]),cov1[,c(2,3)])
#covo$cov<-sqrt(covo$cov)
head(covo)

circos.trackPlotRegion(factors = covo$chr, y = covo$cov,panel.fun =function(x, y) {
})
circos.trackPoints(covo$chr, covo$start, covo$cov, col = "black", pch = 16, cex = 0.4)
}

dev.off()

circos.clear()

##########################################################################################################################################################################
##########################################################################################################################################################################


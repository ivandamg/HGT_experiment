########################################################################################
###############
############### Modify locus Tag with homolog to N16961 using blast .xml (output6)
########## 12 oct 2017

########################################################################################

# Libraries
#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)
library(pheatmap)
#######################################
# Set working directory

setwd('~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/')

###############################
# Preléiminary 
#IN BASH

# 1. Create PROKKA ANNOTATION OF NEW GENOMES

#for i in $(ls Vibrio_*.fasta); do echo $i ; ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus Vibrio --species cholerae --strain $(echo $i | cut -d'_' -f3) --locustag VC-$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f3)_Prokka --rfam --usegenus $i ; done

# 2. Make protein db with the new anotations
#a=0;for i in $(ls Annotation*/*.faa); do echo $(echo $i | cut -d'/' -f2|cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_Prokka_protein/$(echo $i | cut -d'/' -f2|cut -d'_' -f1)_db ; done

# 3. Blastp REFERENCE proteins TO new assemblies proteins.
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1,2,3) ;blastp -db db_Prokka_protein/$(echo $i | cut -d'_' -f1,2,3)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_Prokka_protein/blast_$(echo $i | cut -d'_' -f1,2,3).xml -query REF_*/*.faa ; done

# Import xml files (honmology info)
# xml best blast hit per gene.
# 4 CREATE GENENAMES FILES PER STRAIN
# for i in $(ls *_Prokka.faa); do cat $i | grep '>' | sed 's/ /\t/' | sed 's/>//' | cut -f1 > db_Prokka_protein/$(echo $i | cut -d'_' -f3).GeneNames ; done

###############################################################3333
filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
Homolog_info <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))

names(Homolog_info)<-gsub(".xml","",gsub("blast_","",filesToProcess))

# Create column header

colnam<-strsplit(names(Homolog_info), "_" )
colnam<-lapply(colnam, function (x) x[c(2,1)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))
   
for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(Homolog_info[[i]]) <- NAME }

Homolog_info[[1]]
names(Homolog_info)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
Homolog_info<-lapply(Homolog_info, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

Homolog_info<-lapply(Homolog_info, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

Homolog_info<-lapply(Homolog_info, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )

# select only names of homology
Homolog_info<-lapply(Homolog_info, function(x) x[,c(1,2)])

Homolog_info<-lapply(Homolog_info, function(x) x[!duplicated(x),])

##############################################################################################

###############################

# Import Gene names for file
# xml best blast hit per gene.

filesToProcess <- dir(pattern = "*\\.GeneNames$")  #files to pro# if event 3 merged
Gene_Names <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA")))

names(Gene_Names)<-gsub(".xml","",gsub("blast","",filesToProcess))

# Create column header
head(Gene_Names[[1]])
colnam<-strsplit(names(Gene_Names), "\\." )
colnam<-lapply(colnam, function (x) x[c(1)])

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]))
  
  colnames(Gene_Names[[i]]) <- NAME }

head(Gene_Names[[1]])





##############################################################################################
names(Gene_Names)<-gsub("\\.GeneNames","",names(Gene_Names))
names(Homolog_info)<-gsub("_N16961VCheidelbergNames","", names(Homolog_info) )

head(Homolog_info[[i]])
head(Gene_Names[[i]])

Info_homo<-list()
Merged_name_homolog<- list()
for (i in names(Gene_Names) ){
  
  Info_homo[[i]]<-merge (Gene_Names[[i]],Homolog_info[[i]],by=gsub(".GeneNames","",i), all=T)
  paste("temp",seq(1,dim(Info_homo[[i]])[1]),sep="")
      Info_homo[[i]][,2]<-gsub("^","inASM1542v1_",Info_homo[[i]][,2])
      Info_homo[[i]][,3]<-paste("temp",seq(1,dim(Info_homo[[i]])[1]),sep="")
      c("file",paste("temp",seq(1,dim(Info_homo[[i]])[1]  -1),sep=""))
      Info_homo[[i]]<-cbind.data.frame(c("file",paste("temp",seq(1,dim(Info_homo[[i]])[1]  -1),sep="")),Info_homo[[i]])
      
  write.table(Info_homo[[i]],paste("~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/NewLocusTag/LocusHomology_",gsub(".GeneNames","",i),".txt",sep="" ),sep = "\t",quote = F,row.names = F, col.names = F)
  
  Merged_name_homolog[[i]]<-merge (Gene_Names[[i]],Homolog_info[[i]],by=gsub(".GeneNames","",i), all=T)
  
  Merged_name_homolog[[i]]<-gsub("^","inASM1542v1_",Merged_name_homolog[[i]][,2])
    Merged_name_homolog[[i]]<- as.data.frame(Merged_name_homolog[[i]])
  
  
  write.table(Merged_name_homolog[[i]],paste("~/Documents/Noemie/NoemiePaper_VF/Strains_paper/db_Prokka_protein/NewLocusTag/LocusTag_",gsub(".GeneNames","",i),".txt",sep="" ),sep = "\t",quote = F,row.names = F, col.names = F)
  
  
}

names(Merged_name_homolog)
Merged_name_homolog[[3]]

length(Merged_name_homolog[[1]])

dim(as.data.frame(Merged_name_homolog[[1]]))

##########################################################################
# modify gbk files

# in sublime add a new space to replace info

#1. replace: /locus_tag="Ab_(\w+)"\n  for /locus_tag="Ab_\1"\n\t\t\t\t\t\t\t/note="Ab_\1"\n
# 2. save with ANNO IM

##########################################################################
# FInal modifiaction
# 1. OPEN LOCUSomology in sublime LocusHomology_VC-Sa5Y.txt

# 2. replace file by the respective annotation
# 2. replace :  
        #         ^             with    cat 
        #         \tVC          WITH     | sed 's/\\/note="VC
        #         \tin          with    "/\\/note="in
        #         \t            with    "/' > 
        #         ' > NA"/'     WITH    \\/note="NA"/'

# 3. Add end of file : cat the last temp file  
        # cat temp3610 > VC-A1552_Annotated.gbk 
        #  rm temp*

#4 Copy paste on bash

################################################################################################
# Genome-wide gene homology

# use gene names   and homolog info. 

# use N16961 as reference.


head(Homolog_info[[i]])
head(Gene_Names[[i]])

# Reference

Reference<-Homolog_info[[2]]
Reference

dim(Homolog_info[[2]])
dim(Gene_Names[[2]])

# which of Gene_Names doesnt have correspondance in N16961 heidelberg?

CompN16961<-merge(Gene_Names[[2]],Homolog_info[[2]],by = "N16961Blokesch", all = T)

summary(CompN16961)

ALLN16961Heidelberg_genes<- read.table("~/Documents/Noemie/NoemiePaper_VF/Strains_paper/GeneNames_N16961Heidelberg.txt",h=F)

names(ALLN16961Heidelberg_genes)<-"N16961VCheidelbergNames"

N16961Heid_Blo<-merge( CompN16961, ALLN16961Heidelberg_genes,by = "N16961VCheidelbergNames", all = T)


summary(N16961Heid_Blo)
N16961Heid_Blo$N16961VCheidelbergNames<- as.character(N16961Heid_Blo$N16961VCheidelbergNames)
N16961Heid_Blo$N16961Blokesch<- as.character(N16961Heid_Blo$N16961Blokesch)
N16961Heid_Blo$N16961VCheidelbergNames[N16961Heid_Blo$N16961VCheidelbergNames==NA]<-0
N16961Heid_Blo$N16961Blokesch[N16961Heid_Blo$N16961Blokesch==NA]<-0
N16961Heid_Blo$N16961VCheidelbergNames[N16961Heid_Blo$N16961VCheidelbergNames!=0]<-1
N16961Heid_Blo$N16961Blokesch[N16961Heid_Blo$N16961Blokesch!=0]<-1

tail(N16961Heid_Blo,100)

# Compare A1552 to Reference

Homolog_info[[1]]

################################################################################################
################################################################################################
################################################################################################
#Make homology plot between strains pheatmap

###############################
# USe protein .faa file 
# 2. Make protein db with the new anotations
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_protein/$(echo $i |cut -d'_' -f1)_db ; done

# 3. Blastp REFERENCE proteins TO new assemblies proteins. Change reference
#a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f1) ;blastp -db db_protein/$(echo $i | cut -d'_' -f1)_db -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_protein/Blast_$(echo $i | cut -d'_' -f1)_IN_N16961Blokesch_.xml -query N16961Blokesch_Prokka.faa ; done

# 4. Extract Gene names 
# for i in $(ls *.faa); do cat $i | grep ">" | cut -d' ' -f1 | sed 's/>//' > $(echo $i | cut -d'_' -f1).GeneNames ; done

# 5. Extracting Gene Name of 1st gene on 2nd chromosome
# look at gbk file and extract the name
Chr2A1552_1stGene<-"VC-A1552Blokesch_02796"
Chr2N16961_1stGene<-"VC-N16961Blokesch_02736"
Chr2Sa5Y_1stGene<-"VC-Sa5YBlokech_02729"

setwd('~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/db_protein/')

# Import xml files
# xml best blast hit per gene.

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged



####################################################
############  N16961-REF A1552        ##############
####################################################

filesToProcess<-filesToProcess[grep("IN_N16961Blo",filesToProcess)]
filesToProcess<-filesToProcess[grep("Blast_A1552Blo",filesToProcess)]

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))


# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(4,2)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )


###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles

# Merge names and blast
labREFN16961<-read.table("../N16961Blokesch.GeneNames",h=F)
colnames(labREFN16961)<-"N16961Blokesch"

merged<-merge(Comp1[[1]],labREFN16961,by="N16961Blokesch",all=T)



Comp1<-cbind.data.frame(N16961=as.character(merged$N16961Blokesch),  A1552=as.character(merged$A1552Blokesch) ,stringsAsFactors = FALSE)

summary(Comp1)
legenda<-Comp1$N16961
Comp1[!is.na(Comp1)] <- 1
Comp1[is.na(Comp1)] <- 0


rownames(Comp1)<-legenda

# error if duplicate rownames 
# correct by
legenda[grep("VC-N16961Blokesch_01244",legenda)]<-c("VC-N16961Blokesch_01244","VC-N16961Blokesch_01244b")
rownames(Comp1)<-legenda

Comp1$N16961<-as.numeric(Comp1$N16961)
Comp1$A1552<-as.numeric(Comp1$A1552)

rownames(Comp1[Comp1$A1552==0,])

# Order Comp1 
Comp1  <- Comp1[order(row.names(Comp1 )),] 
# Order Comp1 vector

#

LargeChromosomeN16961<-Comp1[1:match(Chr2N16961_1stGene,rownames(Comp1)),]
SmallChromosomeN16961<-Comp1[match(Chr2N16961_1stGene,rownames(Comp1)):length(rownames(Comp1)),]

rownames(SmallChromosomeN16961)

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/N16961REF_A1552_LargeChromosome.pdf", height=9, width = 4)
pheatmap(LargeChromosomeN16961,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/N16961REF_A1552_SmallChromosome.pdf", height=9, width = 4)
pheatmap(SmallChromosomeN16961,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()
#####################################################

####################################################
############  A1552-REF N16961        ##############
####################################################
filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged

filesToProcess<-filesToProcess[grep("IN_A1552Blo",filesToProcess)]
filesToProcess<-filesToProcess[grep("Blast_N16961Blo",filesToProcess)]

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))


# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(4,2)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )


###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles

# Merge names and blast
labREFA1552<-read.table("../A1552Blokesch.GeneNames",h=F)
colnames(labREFA1552)<-"A1552Blokesch"

merged<-merge(Comp1[[1]],labREFA1552,by="A1552Blokesch",all=T)



Comp1<-cbind.data.frame(A1552=as.character(merged$A1552Blokesch),  N16961=as.character(merged$N16961Blokesch) ,stringsAsFactors = FALSE)

summary(Comp1)
legenda<-Comp1$A1552
Comp1[!is.na(Comp1)] <- 1
Comp1[is.na(Comp1)] <- 0


rownames(Comp1)<-legenda

# error if duplicate rownames 
# correct by
legenda[grep("VC-A1552Blokesch_01412",legenda)]<-c("VC-A1552Blokesch_01412","VC-A1552Blokesch_01412b")
rownames(Comp1)<-legenda

Comp1$N16961<-as.numeric(Comp1$N16961)
Comp1$A1552<-as.numeric(Comp1$A1552)



# Order Comp1 
Comp1  <- Comp1[order(row.names(Comp1 )),] 
# Order Comp1 vector
rownames(Comp1[Comp1$N16961==0,])
#

LargeChromosomeA1552<-Comp1[1:match(Chr2A1552_1stGene,rownames(Comp1)),]
SmallChromosomeA1552<-Comp1[match(Chr2A1552_1stGene,rownames(Comp1)):length(rownames(Comp1)),]

rownames(SmallChromosomeA1552)

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/REFA1552_N16961_LargeChromosome.pdf", height=9, width = 4)
pheatmap(LargeChromosomeA1552,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/REFA1552_N16961_SmallChromosome.pdf", height=9, width = 4)
pheatmap(SmallChromosomeA1552,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F,
         show_colnames = T, show_rownames = F)
dev.off()
# if doesnt work because breaks are equal because all data is equal
pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/REFA1552_N16961_SmallChromosome.pdf", height=9, width = 4)
summary(SmallChromosomeA1552)
SmallChromosomeA1552<-rbind.data.frame(SmallChromosomeA1552,c(0,0))
pheatmap(SmallChromosomeA1552,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F,
         show_colnames = T, show_rownames = F)
dev.off()


####################################################################################################################################################




####################################################
############  REF A1552 - SA5Y        ##############
####################################################
filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged

filesToProcess<-filesToProcess[grep("IN_A1552Blo",filesToProcess)]
filesToProcess<-filesToProcess[grep("Blast_Sa5YBlo",filesToProcess)]

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))


# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(4,2)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )


###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles

# Merge names and blast
labREFA1552<-read.table("../A1552Blokesch.GeneNames",h=F)
colnames(labREFA1552)<-"A1552Blokesch"

merged<-merge(Comp1[[1]],labREFA1552,by="A1552Blokesch",all=T)



Comp1<-cbind.data.frame(A1552=as.character(merged$A1552Blokesch),  Sa5Y=as.character(merged$Sa5YBlokech) ,stringsAsFactors = FALSE)

summary(Comp1)
legenda<-Comp1$A1552
Comp1[!is.na(Comp1)] <- 1
Comp1[is.na(Comp1)] <- 0


rownames(Comp1)<-legenda


Comp1$A1552<-as.numeric(Comp1$A1552)
Comp1$Sa5Y<-as.numeric(Comp1$Sa5Y)

names<-rownames(Comp1[Comp1$Sa5Y==0,])
names[order(names)]
hist(unlist(lapply(strsplit(names,split="_"), function (x) as.numeric(x[2]))),breaks = 200)
# Order Comp1 
Comp1  <- Comp1[order(row.names(Comp1 )),] 
# Order Comp1 vector

#

LargeChromosomeA1552<-Comp1[1:match(Chr2A1552_1stGene,rownames(Comp1)),]
SmallChromosomeA1552<-Comp1[match(Chr2A1552_1stGene,rownames(Comp1)):length(rownames(Comp1)),]

rownames(SmallChromosomeA1552)

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/A1552REF_Sa5Y_LargeChromosome.pdf", height=9, width = 4)
pheatmap(LargeChromosomeA1552,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/A1552REF_Sa5Y_SmallChromosome.pdf", height=9, width = 4)
pheatmap(SmallChromosomeA1552,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = F)
dev.off()
#####################################################

####################################################
############  REF Sa5Y A1552        ##############
####################################################
filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged

filesToProcess<-filesToProcess[grep("IN_Sa5YBlo",filesToProcess)]
filesToProcess<-filesToProcess[grep("Blast_A1552Blo",filesToProcess)]

listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))
names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))


# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(4,2)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))

for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )


###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles

# Merge names and blast
labREFSa5Y<-read.table("../Sa5YBlokech.GeneNames",h=F)
colnames(labREFSa5Y)<-"Sa5YBlokesch"

merged<-merge(Comp1[[1]],labREFSa5Y,by="Sa5YBlokesch",all=T)



Comp1<-cbind.data.frame(Sa5Y=as.character(merged$Sa5YBlokesch),  A1552=as.character(merged$A1552Blokesch) ,stringsAsFactors = FALSE)

summary(Comp1)
legenda<-Comp1$Sa5Y
Comp1[!is.na(Comp1)] <- 1
Comp1[is.na(Comp1)] <- 0


rownames(Comp1)<-legenda



Comp1$A1552<-as.numeric(Comp1$A1552)
Comp1$Sa5Y<-as.numeric(Comp1$Sa5Y)



# Order Comp1 
Comp1  <- Comp1[order(row.names(Comp1 )),] 
# Order Comp1 vector
rownames(Comp1[Comp1$A1552==0,])
names<-rownames(Comp1[Comp1$A1552==0,])
names[order(names)]
#

LargeChromosomeSa5Y<-Comp1[1:match(Chr2Sa5Y_1stGene,rownames(Comp1)),]
SmallChromosomeSa5Y<-Comp1[match(Chr2Sa5Y_1stGene,rownames(Comp1)):length(rownames(Comp1)),]

rownames(SmallChromosomeSa5Y)

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/REFSa5Y_A1552_LargeChromosome.pdf", height=9, width = 4)
pheatmap(LargeChromosomeSa5Y,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F, show_colnames = T, show_rownames = T)
dev.off()

pdf("~/Documents/Noemie/NoemiePaper_VF/comparison_N16961_A1552_Sa5Y/4_BlastHomologyStrains/REFSa5Y_A1552_SmallChromosome.pdf", height=9, width = 4)
pheatmap(SmallChromosomeSa5Y,cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =90 ,cellheight = 0.2,fontsize_col = 20,legend = F,
         show_colnames = T, show_rownames = F)
dev.off()



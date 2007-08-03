### Script that modifies the primary data file generated after the
### segmentation/detection according to Maria's annotations

## File names
file.nucleus <- "nuclei_auto.txt"
file.gene <- "genes_auto.txt"
file.gene.annot <- "genes_annot.csv"
file.nucleus.annot <- "nuclei_annot.csv"
file.extra <- "missing_genes.txt"
file.nucleus.corr <- "nuclei_gold.txt"
file.gene.corr <- "genes_gold.txt"
## Directories
annot.dir <- "manual_corrections"
#data.dir <- "data"

## Some functions
source("R/toolbox.R")

## Compute full paths
file.nucleus <- paste(annot.dir,file.nucleus,sep="/")
file.gene <- paste(annot.dir,file.gene,sep="/")
file.gene.annot <- paste(annot.dir,file.gene.annot,sep="/")
file.nucleus.annot <- paste(annot.dir,file.nucleus.annot,sep="/")
file.extra <- paste(annot.dir,file.extra,sep="/")
file.nucleus.corr <- paste(annot.dir,file.nucleus.corr,sep="/")
file.gene.corr <- paste(annot.dir,file.gene.corr,sep="/")

## Read files
gene <- read.table(file.gene,header=TRUE,as.is="img")
nucleus <- read.table(file.nucleus,header=TRUE,as.is="img")
nucleus.to.del <- read.table(file.nucleus.annot,header=TRUE,sep=",",as.is=c("img","Comment"))
gene.to.corr <- read.table(file.gene.annot,header=TRUE,sep="\t",as.is=c("img","comment"))
gene.extra <- read.table(file.extra,header=TRUE)

## Identify nuclei to ignore
nucleus.to.del <- subset(nucleus.to.del,exclude==1)
row.to.del <- match.from.dist(nucleus,nucleus.to.del,c("x","y","z"),"img")
if(!all.equal(row.to.del$dist,rep(0,dim(row.to.del)[1])))
  warning("Check nucleus matches")
cat(dim(row.to.del)[1]," nuclei ignored\n")
nucleus <- remove.lines(nucleus,row.to.del$row)
## And write corrected file for nucleus data
write.table(nucleus,file.nucleus.corr,row.names=FALSE)

## Identify spots to remove
gene.to.del <- subset(gene.to.corr,multiplicity%in%c(0,2))
gene <- merge(gene,nucleus[,c("img","nucleus")],by=c("img","nucleus"))
row.to.del <- match.from.dist(gene,gene.to.del,c("px","py","pz"),c("img","nucleus","gene"))
row.to.del <- row.to.del[order(row.to.del$dist),]
if(any(row.to.del$dist>0.05))
    warning("Check gene matches")
gene.corr <- remove.lines(gene,row.to.del$row)
gene.extra <- remove.columns(gene.extra,c("x.nucleus","y.nucleus","z.nucleus","threshold"))
gene.extra <- rename.columns(gene.extra,c("x.gene","y.gene","z.gene"),c("x","y","z"))
gene.extra[,names(gene.corr)[!(names(gene.corr)%in%names(gene.extra))]] <- NA
gene.corr <- rbind(gene.corr,gene.extra)
write.table(gene.corr,file.gene.corr,row.names=FALSE)

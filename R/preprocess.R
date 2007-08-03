### Pre-processing of data (saved in a RData file)

file.nucleus <- "data/nuclei_corr.txt"
file.gene <- "data/genes_corr.txt"
truncate.evf <- 0
file.save <- "gene_nucleus.RData"

## Some functions
source("R/toolbox.R")

## Read data files
nucleus <- read.table(file.nucleus,header=TRUE)
gene <- read.table(file.gene,header=TRUE)
## Redundant field
gene <- remove.columns(gene,"stimulation")

## Compute slide number
nucleus$slide <- NA
## First slide: before May 2007, second slide after May 2007
nucleus[get.dates.from.strings(nucleus$img)<=as.Date("2007-05-01"),"slide"] <- 1
nucleus[get.dates.from.strings(nucleus$img)>as.Date("2007-05-01"),"slide"] <- 2
nucleus$slide <- as.factor(nucleus$slide)

## Merge data frames
gene <- merge(nucleus,gene,by.x=c("img","nucleus"),by.y=c("img","nucleus"),
                suffixes=c(".nucleus",".gene"))
## Corrections and minor modifications
levels(gene$gene) <- c("csn","wap")
levels(gene$stimulation) <- c("non-stimulated","stimulated")
## Error on two column labels in a data file
gene <- rename.columns(gene,c("elongation","size"),c("volume","flatness"))
## Use EVF instead of CI
gene <- rename.columns(gene,"ci","evf")
gene$diameter <- 2*(gene$volume/(4/3*pi))^(1/3)
# Add a factor identifying the nuclei
gene$nucleus.id <- paste(gene$img,gene$nucleus)
levels(gene$nucleus.id) <- as.character(1:length(levels(gene$nucleus.id)))

## Compute number of csn and wap spots for each nucleus
gene.nb <- by(gene[,1],list(gene$nucleus.id,gene$gene),length)
gene <- merge(gene,unclass(gene.nb),by.x="nucleus.id",by.y="row.names")
## Merging has modified the order of rows, reset it.
gene <- gene[sort.list(gene$nucleus.id),]
## Joint distribution of numbers of csn spots and wap spots
nucleus.before.filtering <- gene[!duplicated(gene$nucleus.id),c("stimulation","img","nucleus.id","csn","wap")]

## Filtering
gene <- subset(gene,csn%in%3:4 & wap%in%3:4)
## Drop unused levels for img factor in order to avoid problems
## when using function by with a grouping on factor img
levels(gene$img)[!(levels(gene$img)%in%gene$img)] <- NA
levels(gene$nucleus.id)[!(levels(gene$nucleus.id)%in%gene$nucleus.id)] <- NA
nucleus <- gene[!duplicated(gene$nucleus.id),c("nucleus.id","csn","wap","stimulation","diameter",
                                               "flatness","mean.nucleus","sigma","threshold")]

## Truncation set negative EVF to zero
neg.evf.number <- sum(gene$evf<0)
evf.min <- min(gene$evf)
gene[gene$evf<truncate.evf,"evf"] <- truncate.evf

## Distances computations
gene44 <- subset(gene,csn==4 & wap==4,
                 select=c("nucleus.id","stimulation","gene",
                   "px","py","pz","diameter","flatness","slide"))
gene44 <- remove.unused.levels(gene44,"nucleus.id")
## Reorder gene lines in order to group spots within the same nucleus
gene44 <- gene44[order(gene44$nucleus.id,gene44$gene),]
gene44$idx.spot <- 1:4
## Distances entre spots dans un même noyau
inter.spot.dist <- by(gene44,gene44$nucleus.id,edist,c("px","py","pz"),"gene")
## Extract minimal distances
csn2wap.dist <- unlist(lapply(inter.spot.dist,function(x){apply(x,1,min)}))
wap2csn.dist <- unlist(lapply(inter.spot.dist,function(x){apply(x,2,min)}))
lut <- row.names(gene44)
names(lut) <- paste(gene44$nucleus.id,lut,sep=".")
names(csn2wap.dist) <- lut[names(csn2wap.dist)]
names(wap2csn.dist) <- lut[names(wap2csn.dist)]
csn.gene <- subset(gene44,gene=="csn")
csn2csn.dist <- by(csn.gene,csn.gene$nucleus.id,edist,c("px","py","pz"))
csn2csn.dist <- lapply(csn2csn.dist,function(x){diag(x)<-Inf;x})
csn2csn.dist <- unlist(lapply(csn2csn.dist,function(x){apply(x,1,min)}))
names(csn2csn.dist) <- lut[names(csn2csn.dist)]
wap.gene <- subset(gene44,gene=="wap")
wap2wap.dist <- by(wap.gene,wap.gene$nucleus.id,edist,c("px","py","pz"))
wap2wap.dist <- lapply(wap2wap.dist,function(x){diag(x)<-Inf;x})
wap2wap.dist <- unlist(lapply(wap2wap.dist,function(x){apply(x,1,min)}))
names(wap2wap.dist) <- lut[names(wap2wap.dist)]
gene44$homo.dist <- NA
gene44$hetero.dist <- NA
gene44[names(csn2csn.dist),"homo.dist"] <- csn2csn.dist
gene44[names(wap2wap.dist),"homo.dist"] <- wap2wap.dist
gene44[names(csn2wap.dist),"hetero.dist"] <- csn2wap.dist
gene44[names(wap2csn.dist),"hetero.dist"] <- wap2csn.dist
gene44 <- remove.columns(gene44,c("px","py","pz"))
spot.dist <- reshape(gene44,direction="wide",idvar=c("nucleus.id","gene"),
                     v.names=c("homo.dist","hetero.dist"),timevar="idx.spot")
spot.dist$min.homo.dist <- apply(spot.dist[,paste("homo.dist",1:4,sep=".")],1,min)

## Save pre-processed data
save(gene,nucleus.before.filtering,nucleus,neg.evf.number,evf.min,gene44,spot.dist,file=file.save)

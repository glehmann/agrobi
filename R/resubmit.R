## File paths
annotation.file <- "manual_corrections_detections.csv"
nucleus.file <- "caswap_data_nucleus.txt"
output.file <- "missing_genes_slide2.csv"
## Some functions
source("toolbox.R")
## Read data
missing.gene <- read.table(annotation.file,header=TRUE,sep="\t")
nucleus <- read.table(nucleus.file,header=TRUE)
## Select slide. First slide: before May 2007, second slide after May 2007
missing.gene <- subset(missing.gene,
                       get.dates.from.strings(img) > as.Date("2007-05-01") & multiplicity==1)
## For nucleus data, keep only required fields
nucleus <- subset(nucleus,select=c("stimulation","img","nucleus","x","y","z","threshold"))
## Merge gene data and nucleus data into a single frame
missing.gene <- merge(nucleus,missing.gene,by.x=c("img","nucleus"),by.y=c("img","nucleus"),
                      suffixes=c(".nucleus",".gene"))
## Keep only relevant fields
missing.gene <- subset(missing.gene,select=c("stimulation.nucleus","img","nucleus",
                                      "x.nucleus","y.nucleus","z.nucleus","threshold",
                                      "gene","x.gene","y.gene","z.gene"))
## Rename stimulation field which is redundant in gene and nucleus files
cnames <- names(missing.gene)
cnames[cnames=="stimulation.nucleus"] <- "stimulation"
names(missing.gene) <- cnames
## Write file
write.table(missing.gene,output.file,row.names=FALSE)


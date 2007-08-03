## Exploratory analysis of component tree

## Type of output for graphics (2D projections)
#output <- "pdf" # PDF file
output <- "screen" # current graphic device (e.g. x11)

## Read data file
#image.id <- "20070320_19" # image with 4 WAP spots and 2 dubious spots
#image.id <- "20070410_20" # 5 WAP spots
image.id <- "20070410_24" # 5 WAP spots
data.dir <- "data"
file.tree <- paste(data.dir,"/tree_",image.id,".txt",sep="")
tree.var <- read.table(file.tree,header=TRUE)

## Correct first row (root)
tree.var[1,"local_intensity"] <- max(tree.var$intensity)
tree.var[1,"spot_base"] <- 0
tree.var[1,"parent_children"] <- NA

## Principal component analysis
acp <- princomp(tree.var[,c("integrated_intensity","physical_size",
                            "volume_levelling","local_intensity","intensity")],
                cor=TRUE)
summary(acp)
tree.comp <- predict(acp)[,1:3] # Keep only first three components

## Precompute some tree features before plotting
edges <- merge(data.frame(i=1:(dim(tree.comp)[1]),
                          node=tree.var$node,parent=tree.var$parent),
              data.frame(i=1:(dim(tree.comp)[1]),node=tree.var$node),
              by.x="parent",by.y="node",suffixes=c(".child",".parent"))
leaves <- !(tree.var$node%in%unique(tree.var$parent))
spots <- tree.var$spot_base

## function that adds a 2D tree to an existing plot
tree.2d.plot <- function(x,y,edges,leaves,spots) {
  ## edges
  matlines(rbind(x[edges$i.child],x[edges$i.parent]),
           rbind(y[edges$i.child],y[edges$i.parent]),
           col="grey",lty=1)
  ## leaves
  points(x[leaves],y[leaves],col="blue")
  ## spots
  points(x[spots>0],y[spots>0],
         pch=spots[spots>0]+48,col="red")
  ## other nodes
  points(x[!leaves & spots==0],y[!leaves & spots==0],pch=".",cex=2)
}

## Plot tree projections of the 3 first main component planes. Plots
## share common scales.
## Uncomment 3 lines below for full view of the tree
#xlim <- range(tree.comp[,1])
#ylim <- range(tree.comp[,2])
#zlim <- range(tree.comp[,3])
## Zoom on the leaves
#xlim <- c(-12,0) # 20070320_19
#ylim <- c(-4.5,3)
#zlim <- c(-6,5)
#xlim <- c(-8,0.5) # 20070410_20
#ylim <- c(-3,3)
#zlim <- c(-4.5,5)
xlim <- c(-6,0.5) # 20070410_24
ylim <- c(-3,5)
zlim <- c(-4.5,5)
if(output=="pdf")
  pdf(file=paste("tree2D_",image.id,".pdf",sep=""),paper="a4",width=7,height=7.3)
par(mfrow=c(2,2),mar=rep(.1,4),oma=c(4,4,6,4))
plot(tree.comp[,1],tree.comp[,3],type="n",xlim=xlim,ylim=zlim,axes=FALSE,xlab="",ylab="")
axis(3); mtext("Comp.1",3,line=3)
axis(2); mtext("Comp.3",2,line=3)
box()
tree.2d.plot(tree.comp[,1],tree.comp[,3],edges,leaves,spots)
plot.new()
plot(tree.comp[,1],tree.comp[,2],type="n",xlim=xlim,ylim=ylim,,axes=FALSE,xlab="",ylab="")
axis(1); mtext("Comp.1",1,line=3)
axis(2); mtext("Comp.2",2,line=3)
box()
tree.2d.plot(tree.comp[,1],tree.comp[,2],edges,leaves,spots)
plot(tree.comp[,3],tree.comp[,2],type="n",xlim=zlim,ylim=ylim,axes=FALSE,xlab="",ylab="")
axis(1); mtext("Comp.3",1,line=3)
axis(4); mtext("Comp.2",4,line=3)
box()
tree.2d.plot(tree.comp[,3],tree.comp[,2],edges,leaves,spots)
title(main=image.id,outer=TRUE,line=5)
if(output=="pdf")
  dev.off()

## Function that helps to identify a node on a tree projection on the
## plane spanned by components 2 and 3. Left click on the x11 device.
get.pos <- function() {
  plan.pos <- locator(1)
  dist.2.nodes <- (tree.comp[,3]-plan.pos$x)^2+(tree.comp[,2]-plan.pos$y)^2
  tree.var[which(dist.2.nodes==min(dist.2.nodes)),c("x","y","z")]
}

## 3D plot of the tree
require(scatterplot3d)
tree3d <- scatterplot3d(tree.comp[,1],tree.comp[,2],tree.comp[,3],type="n")
for(i in 1:dim(edges)[1])
  tree3d$points(tree.comp[c(edges$i.child[i],edges$i.parent[i]),"Comp.1"],
                tree.comp[c(edges$i.child[i],edges$i.parent[i]),"Comp.2"],
                tree.comp[c(edges$i.child[i],edges$i.parent[i]),"Comp.3"],
                col="grey",type="l")
tree3d$points3d(tree.comp[!(tree.var$node%in%unique(tree.var$parent)),1:3],col="blue")
tree3d$points3d(tree.comp[tree.var$spot_base>0,1:3],
       pch=tree.var$spot_base[tree.var$spot_base>0]+48,col="red")
tree3d$points3d(tree.comp[(tree.var$node%in%unique(tree.var$parent))&tree.var$spot_base==0,1:3],pch=".")

## Display the coefficients of the linear combinations defining the main components
loadings(acp)

## Suivi de la décimation
follow <- read.table("granulometry-volume-levelling.txt",header=TRUE)
branches <- follow[-1,][diff(follow$queue)!=0,]
par(mfrow=c(1,2))
hist(branches$attrib[branches$attrib>1000])
boxplot(branches$attrib[branches$attrib>1000])
i4 <- range(which(follow$queue==5))
plot(follow$total,type="l")
abline(v=i4[1])
abline(v=i4[2])
plot(follow$total+follow$queue,type="l")
plot(follow$attrib[1:i4[2]],type="l")

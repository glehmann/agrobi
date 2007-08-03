### MISCELLEANEOUS FUNCTIONS USED FOR PROCESSING NUCLEUS AND
### GENE DATA

require(xtable)

### Extract dates from strings
get.dates.from.strings <- function(x) {
  ## x: vector of strings containing dates formatted as "xxx20070527xx"
  ## Value: vector of dates
  filenames <- x
  starts <- regexpr("[[:digit:]]{8}",filenames)
  as.Date(substr(filenames,starts,starts+7),"%Y%m%d")
}

### Compute a distance matrix from coordinates stored in a dataframe
edist <- function(x,coord,type=NULL) {
  ## x: dataframe containing coordinates
  ## coord: vector of names of the columns containing the coordinates
  ## type: to be used when distances between points of two types are
  ##       to be computed. Then type should be the name of the column
  ##       of x which identifies the point type. This column must be a factor.
  ## Value: a distance matrix. If argument type is NULL, square matrix
  ##        with as many rows as lines of dataframe x
  mat.dist <- as.matrix(dist(subset(x,select=coord))) # compute distances
  ## If only distances between points of different types are required, extract
  ## a sub-matrix containing only distances between points of different types.
  if(!is.null(type)) {
    ## Unvaluated expression for extracting rows (points of type 1)
    cond.row <- parse(text=paste(type,"==","\"",levels(x[,type])[1],"\"",sep=""))
    ## Unvaluated expression for extracting columns (points of type 2)
    cond.col <- parse(text=paste(type,"==","\"",levels(x[,type])[2],"\"",sep=""))
    ## Use row and column names for matrix extraction
    mat.dist <- mat.dist[row.names(subset(x,eval(cond.row))),
                         row.names(subset(x,eval(cond.col))),drop=FALSE]
    ## Give names to the dimensions
    mat.dist.names <- dimnames(mat.dist)
    names(mat.dist.names) <- levels(x[,type][1:2])
    dimnames(mat.dist) <- mat.dist.names
  }
  mat.dist
}

### Dataframe with fewer columns
remove.columns <- function(x,colnames) {
  ## x: dataframe
  ## colnames: vector of column names to be removed
  ## Value: dataframe with fewer columns
  x[,names(x)[!(names(x)%in%colnames)]]
}

### Dataframe with fewer lines
remove.lines <- function(x,rownames) {
  ## x: dataframe
  ## rownames: vector of row names to be removed
  ## Value: dataframe with fewer rows
  x[!(row.names(x)%in%rownames),]
}

### Rename columns of a dataframe
rename.columns <- function(x,old.colnames,new.colnames) {
  ## x: dataframe
  ## old.colnames: vector of current column names
  ## new.colnmaes: vector of new column names
  ## Value: a dataframe with renamed columns
  cnames <- names(x)
  cnames[match(old.colnames,cnames)] <- new.colnames
  names(x) <- cnames
  x
}

### Munkres implementation of the Hungarian algorithm (John Weaver)
if(!is.loaded("munkres_Rwrap")) {
  dyn.load(paste("R/munkres_Rwrap",.Platform$dynlib.ext, sep=""))
}
munkres <- function(md) {
  ## md: distance matrix. Row and columns may be named.
  ## Value: a dataframe if both row and columns of md are named.
  ##        A three-column matrix otherwise defining pairs. First
  ##        column contains row identifiers (names or indices),
  ##        second column contains column identifiers. Third column
  ##        contains distances.
  nr <- dim(md)[1]
  nc <- dim(md)[2]
  asg <- .C("munkres_Rwrap",as.integer(nr),as.integer(nc),as.double(md))[[3]]
  asg <- matrix(asg,nr,nc,dimnames=dimnames(md))
  asg <- which(asg==0,arr.ind=TRUE)
  asg <- cbind(asg,md[asg])
  if(all(!is.null(dimnames(md)))) {
    data.frame(row=I(rownames(md)[asg[,1]]),col=I(colnames(md)[asg[,2]]),dist=asg[,3])
  } else {
    asg
  }
}

### Remove unused levels of a factor in a dataframe
remove.unused.levels <- function(x,colname) {
  ## x: a dataframe
  ## colname: name of a factor column in dataframe x
  ## Value: dataframe with cleaned factor column
  if(!colname%in%names(x))
    stop(paste(colname,"is not a column name"))
  levels(x[,colname])[!(levels(x[,colname])%in%x[,colname])] <- NA
  x
}

### Match two lists of objects based on distances
match.from.dist <- function(x,y,coord,group,value="paired") {
  ## x: dataframe representing the first list of objects
  ## y: dataframe reprensenting the second list of objects
  ## coord: vector of names of the columns which are used for distance
  ##        computations.
  ## group: vector of names of columns defining groups. Matches are
  ##        performed within groups.
  ## value: either "paired" or "unpaired". If "paired", only rows of
  ##        dataframe x which have been paired are returned. If "unpaired",
  ##        only rows of dataframe which have not been paired are returned.
  ## Value: a subset of dataframe x
  xx <- x
  xx$id <- row.names(x)
  xx$from <- "x"
  xx$by <- if(length(group)==1) xx[,group] else apply(xx[,group],1,paste,collapse="_")
  yy <- y
  yy$id <- row.names(y)
  yy$from <- "y"
  yy$by <- if(length(group)==1) yy[,group] else apply(yy[,group],1,paste,collapse="_")
  ## Reduce xx and yy
  yy <- subset(yy,by%in%unique(xx$by))
  xx <- subset(xx,by%in%unique(yy$by))
  xy <- rbind(xx[,c("id","by","from",coord)],yy[,c("id","by","from",coord)])
  xy$from <- as.factor(xy$from)
  xy <- remove.unused.levels(xy,"by")
  distances <- by(xy,xy$by,edist,coord=coord,type="from")
  pair.dist <- lapply(distances,munkres)
  res <- NULL
  for(row in pair.dist) res <- rbind(res,row)
  res
}

## Require package xtable
## Distribution of gene numbers per nucleus as a LaTeX table
latex.gene.numbers.per.nucleus <- function(data,row.vars,...) {
  ## data: dataframe containing nucleus data.
  ## row.vars: vector of strings containing modalities to put as rows
  tab <-ftable(table(data),row.vars=row.vars)
  cnames <- attributes(tab)$col.vars[[1]]
  rnames <- rep("",nrow(tab)+1)
  rnames[1] <- paste("number of",names(attributes(tab)$col.vars),"genes")
  ng <- length(attributes(tab)$row.vars[[2]])
  for(i in 1:length(attributes(tab)$row.vars[[1]])) {
    k <- 2+(i-1)*ng
    rnames[k] <- paste(names(attributes(tab)$row.vars)[1],"=")
    rnames[k] <- paste(rnames[k],attributes(tab)$row.vars[[1]][i])
    rnames[k] <- paste(rnames[k],",",sep="")
    rnames[k] <- paste(rnames[k],attributes(tab)$row.vars[[2]][1],
                     names(attributes(tab)$row.vars)[2],"genes")
    if(ng>=2) {
      idx <- 2:ng 
      rnames[idx+k-1] <- paste(attributes(tab)$row.vars[[2]][idx],
                               names(attributes(tab)$row.vars)[2],"genes")
    }
  }
  ## Convertir en dataframe
  tab <- cbind(rnames,as.data.frame(rbind(as.numeric(cnames),unclass(tab))))
  dsp <- c("s","s",rep("d",ncol(tab)-1))
  print(xtable(tab,display=dsp,align=c("c","r|",rep("r",ncol(tab)-1)),...),
        include.rownames=FALSE,include.colnames=FALSE,hline.after=1)
}

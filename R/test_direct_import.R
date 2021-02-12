################################################################################
library(tictoc)
library(data.table)
library(dplyr)
library(sparsify)
library(mltools)
library(rtracklayer)
library(hdf5r)
library(datastructures)

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\test.data\\10gen")


#setwd(dirname(getActiveDocumentContext()$path))

#tabix.input.file <- "..\\test.data\\s50cell\\bedgraph"


files <- list.files (getwd())

#files <- files[!grepl("^.*AGOG_merged.xlsx$", files)]

files <- files[!grepl("^(_|~).*$", files)]

df1 <- NULL

tic("Overall time")

batch <- 5

for (i in 0:(length(files)/batch)) {
  
  df2 <- NULL
  
  for (j in 1:batch) {
    
    f <- files[i*batch+j]
    
    time2s <- as.numeric(Sys.time())*1000
    
    data <- read_tsv(f, col_names = FALSE, progress=FALSE, col_types = cols())
    colnames(data) <- c("chr","start","end",strsplit(f, "[.]")[[1]][1])
    
    if (is.null(df2)) {
      df2 <- data
    } else {
      df2 <- full_join(df2, data,by=c("chr","start","end"))
    }
    
    time2e <- as.numeric(Sys.time())*1000
    
    print(paste("Parsing [",i*batch+j,"]:", f,"in",time2e-time2s,"ms"))
      
  }
  
  time1s <- as.numeric(Sys.time())*1000
  
  if (is.null(df1)) {
    df1 <- df2
  } else {
    df1 <- full_join(df1, df2,by=c("chr","start","end"))
  }
  
  time1e <- as.numeric(Sys.time())*1000
  
  print(paste("Joining",i*batch,"to",(i+1)*batch-1,"in",time1e-time1s)
  )
  
}
toc()

### Parse n files, then merge as matrix, then convert to sparse and merge in batch



files <- list.files (getwd())

#files <- files[!grepl("^.*AGOG_merged.xlsx$", files)]

files <- files[!grepl("^(_|~).*$", files)]

df1 <- NULL

tic("Overall time")

files <- files[1:10]

batch <- 5

for (i in 0:(ceiling(length(files)/batch))) {
  
  df2 <- df2.sp <- NULL
  
  for (j in 1:batch) {
    
    f <- files[i*batch+j]
    
    time2s <- as.numeric(Sys.time())*1000
    
    #data <- read_tsv(f, col_names = FALSE, progress=FALSE, col_types = cols())
    data <- data.table::fread(f, header=FALSE,nThread=8)
    
    
    colnames(data) <- c("chr","start","end",strsplit(f, "[.]")[[1]][1])

    if (is.null(df2)) {
      df2 <- data
    } else {
      df2 <- merge(df2, data, by = c("chr","start","end"), all = TRUE)#full_join(df2, data,by=c("chr","start","end"))
    } 

    rm(data)
        
    time2e <- as.numeric(Sys.time())*1000
    
    print(paste("Parsing [",i*batch+j,"]:", f,"in",time2e-time2s,"ms"))
    
    if (i*batch+j >= length(files)) {
     print("hit last index")
       break
    }
    
  }
  
  time1s <- as.numeric(Sys.time())*1000
  
  df2$row <- with(df2, paste(chr,start,end, sep="."))
  df2 <- select(df2,row,everything(),-chr,-start,-end)

  df2[is.na(df2)] <- 0 #replace with more efficient

  df2.sp <- Matrix(as.matrix(df2[, -1]), sparse = TRUE) #sparsify(df2[, -1], sparsifyNAs = TRUE,)
  rownames(df2.sp) <- df2$row
  
  if (is.null(df1)) {
    df1 <- df2.sp
  } else {
    df1 <- merge.sparse(df1,df2.sp)
  }
  
  time1e <- as.numeric(Sys.time())*1000
  
  print(paste("Joining",i*batch+1,"to",(i+1)*batch,"in",time1e-time1s)
  )
  
}
toc()

########################

fread("./foo.csv")


fread("cat ./foo.csv | awk -F ',' 'BEGIN { s = 5 } { for (i=1; i<=NF; i++) printf(\"%s%s\", $(i), i<s ? OFS : i<NF ? \"\" : ORS) }'")
################################################################################

files <- list.files (getwd())

#files <- files[!grepl("^.*AGOG_merged.xlsx$", files)]

files <- files[!grepl("^(_|~).*$", files)]

df1 <- NULL

tic("Overall time")

files <- files[1:10]

batch <- 5

for (i in 0:(ceiling(length(files)/batch))) {
  
  df2 <- df2.sp <- NULL
  
  for (j in 1:batch) {
    
    f <- files[i*batch+j]
    
    time2s <- as.numeric(Sys.time())*1000
    
    #data <- read_tsv(f, col_names = FALSE, progress=FALSE, col_types = cols())
    data <- data.table::fread(f, header=FALSE,nThread=8)
    cell <- strsplit(f, "[.]")[[1]][1]
    data$cell <- cell
    
    
    colnames(data) <- c("chr","start","end","val","cell")
    object.size(data)/1000/1000

    
    if (is.null(df2)) {
      df2 <- data
    } else {
      df2 <- merge(df2, data, by = c("chr","start","end"), all = TRUE)#full_join(df2, data,by=c("chr","start","end"))
    } 
    
    rm(data)
    
    time2e <- as.numeric(Sys.time())*1000
    
    print(paste("Parsing [",i*batch+j,"]:", f,"in",time2e-time2s,"ms"))
    
    if (i*batch+j >= length(files)) {
      print("hit last index")
      break
    }
    
  }
  
  time1s <- as.numeric(Sys.time())*1000
  
  df2$row <- with(df2, paste(chr,start,end, sep="."))
  df2 <- select(df2,row,everything(),-chr,-start,-end)
  
  df2[is.na(df2)] <- 0 #replace with more efficient
  
  df2.sp <- Matrix(as.matrix(df2[, -1]), sparse = TRUE) #sparsify(df2[, -1], sparsifyNAs = TRUE,)
  rownames(df2.sp) <- df2$row
  
  if (is.null(df1)) {
    df1 <- df2.sp
  } else {
    df1 <- merge.sparse(df1,df2.sp)
  }
  
  time1e <- as.numeric(Sys.time())*1000
  
  print(paste("Joining",i*batch+1,"to",(i+1)*batch,"in",time1e-time1s)
  )
  
}
toc()


################ Test converting every file to sparse matrix first

data <- data.table::fread(f, header=FALSE,nThread=8)
colnames(data) <- c("chr","start","end",strsplit(f, "[.]")[[1]][1])
data$row <- with(data, paste(chr,start,end, sep="."))
data <- select(data,-chr,-start,-end)

data.sp <- data.table(data)
data.sp <- sparsify(data.table(data.sp[, -"row"]), sparsifyNAs = TRUE,)
rownames(data.sp) <- data$row


##### Combine using genomic ranges #############################

file.bed.all <- NULL

for (file in files) {

  name <- strsplit(file, "[.]")[[1]][1]
  
  cols <- c("file"="integer")
  names(cols) <- strsplit(file, "[.]")[[1]][1]
  
  file.bed <- rtracklayer::import(file,format="bedGraph",extraCols=cols)
  
  if (is.null(file.bed.all)) {
    file.bed.all <- file.bed
  } else {
    file.bed.all <- combine(file.bed.all, file.bed)
  }
  
}




  
 x <- c("name"="5") 
 
 
  
range1 <- rtracklayer::import(files[1],format="bedGraph",extraCols=c("bed1"="integer"))
head(range1)

range2 <- rtracklayer::import(files[3],format="bedGraph",extraCols=c("bed2"="integer"))
head(range2)

range <- combine(range1,range2)
range.meta <- as.data.frame(mcols(range))
mcols(range) <- NULL

df[df==""]<-NA


range.meta[range.meta==0] <- -1

range.meta[is.na(range.meta)] <- 0

range.sparse <- Matrix(data.matrix(range.meta), sparse=TRUE)






mcols(range) <- NULL

################## Merge sparse matrixes together #########################

merge.sparse <- function(...) {
  
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()
  
  for (M in list(...)) {
    
    cnold <- colnames(M)
    rnold <- rownames(M)
    
    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)
    
    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)
    ind <- unname(which(M != 0,arr.ind=T))
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,M@x)
  }
  
  sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
}



########## Find the consensus regions, then fill matrix######################

hdf5 <- paste(getwd(),"/singlecellexp.h5",sep="")

h5createFile(hdf5)

files <- list.files (getwd())
files <- files[!grepl("^(_|~).*$", files)]

files <- files[1:5]

df <- NULL

for (file in files) {
  
  tic()
  
  df.temp <- data.table::fread(file, header=FALSE, nThread=8, select=c(1:3))
  
  df <- if(is.null(df)) df.temp else merge(df,df.temp,all.x=TRUE,all.y=TRUE)

  cat("\nDone:",file,"  Time:", capture.output(toc()))
  
}

colnames(df) <- c("chr","start","end")

df$index <- 1:nrow(df)


#### hdf5r
file.h5 <- H5File$new(hdf5, mode = "w")

assays.grp <- file.h5$create_group("assays")

assays.grp[["assay1"]] <- df

file.h5$ls()

assays.grp$ls()

file.h5

assay1 <- assays.grp[["assay1"]]

assay1.type <- assay1$get_type()
assay1.type$get_class()


file.h5$close_all()


# rhdf5

h5createFile(hdf5)

h5ls(hdf5)

h5f = H5Fopen(hdf5)

#h5write(df, hdf5,"assay")

h5writeDataset.data.frame(df, h5f,"assay",DataFrameAsCompound = FALSE)

head(h5read(hdf5,"assay"))

assay <- HDF5Array(hdf5, "assay")

h5closeAll()

### hashmap
keys <- apply(df[,1:3],1,paste,sep=".",collapse=".")
values <- df$index

bimap <- bimap("character", "integer")
bimap <- insert(bimap, keys, values)

at(bimap, keys[100], "values")
at(bimap, values[100], "keys")

h5write(bimap, hdf5,"bimap")








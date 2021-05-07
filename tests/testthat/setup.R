message("Test setup starting...")

df1 <- data.table(chr=rep("chr1",5),start=(1:5)*2,end=(1:5)*2+1,value=0,cov=1)
df2 <- data.table(chr=rep("chr1",5),start=(3:7)*2,end=(3:7)*2+1,value=25,cov=2)
df3 <- data.table(chr=rep("chr1",5),start=(6:10)*2,end=(6:10)*2+1,value=50,cov=3)
df4 <- data.table(chr=rep("chr2",8),start=(3:10)*2,end=(3:10)*2+1,value=100,cov=4)

files <- c("df1.bedgraph","df2.bedgraph","df3.bedgraph","df4.bedgraph")
files <- file.path(tempdir(),files)

write.table(df1, file = files[1], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df2, file = files[2], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df3, file = files[3], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)
write.table(df4, file = files[4], row.names=FALSE, sep="\t",col.names=FALSE, quote = FALSE)

h5_dir <- paste0(tempdir(),"/sse")

scm.h5 <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,cov_idx=5)
scm.mem <- read_beds(files,h5=FALSE,cov_idx=5)

# scMethrix_data <- list(h5=scm.h5,mem=scm.mem)
# use_data(scMethrix_data)

message("Test setup completed")

#files <- list.files(system.file('extdata', package = 'scMethrix'))
# files <- c("D:/Git/scMethrix/inst/extdata/C0.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C1.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C2.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C3.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C4.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C5.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C6.bedgraph",
#            "D:/Git/scMethrix/inst/extdata/C7.bedgraph")
# 
# h5_dir <- paste0(tempdir(),"/sse")
# 
# scm.h5.cov <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,cov_idx=c(5,6))
# scm.mem.cov <- read_beds(files,h5=FALSE,cov_idx=c(5,6))
# 
# scm.h5 <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE)
# scm.h5 <- read_beds(files,h5=FALSE)
# 



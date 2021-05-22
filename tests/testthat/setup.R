message("Test setup starting...")

files <- c(system.file("extdata", "C1.bedgraph", package="scMethrix"),
           system.file("extdata", "C2.bedgraph", package="scMethrix"),
           system.file("extdata", "C3.bedgraph", package="scMethrix"),
           system.file("extdata", "C4.bedgraph", package="scMethrix"))

files <- c("D:/Git/scMethrix/inst/extdata/C1.bedgraph","D:/Git/scMethrix/inst/extdata/C2.bedgraph",
           "D:/Git/scMethrix/inst/extdata/C3.bedgraph","D:/Git/scMethrix/inst/extdata/C4.bedgraph")

h5_dir <- paste0(tempdir(),"/sse")

scm.h5 <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,cov_idx=5)
scm.mem <- read_beds(files,h5=FALSE,cov_idx=5)

# scMethrix_data <- list(h5=scm.h5,mem=scm.mem)
# use_data(scMethrix_data)

message("Test setup completed")
### Initialization stuff #################################
list.of.packages <- c("dplyr","magrittr","pryr","recommenderlab","GenomicRanges","SummarizedExperiment",
                      "DelayedArray","HDF5Array","tictoc","readr","sos","rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

tabix.input.file <- ".\\sample.data\\big50cell\\Parsed\\big50cell.tsv"
hdf5.output.folder <- ".\\sample.data\\big50cell\\hdf5"

### Test Rust TSV file input #############################


tic("Total")
tic("Rust TSV file import")

scMethrix.data <- read_tsv(tabix.input.file, col_names = TRUE)

toc()
tic("Generate sparse matrix and genomic ranges")
gc()
scMethrix.mat <- data.matrix(dplyr::select(scMethrix.data, -chr,-start,-end), rownames.force = NA)
gc()
scMethrix.rowRanges <- GenomicRanges::makeGRangesFromDataFrame(scMethrix.data)
rm(scMethrix.data)
gc()
scMethrix.mat.sparse <- Matrix(scMethrix.mat, sparse=TRUE)
rm(scMethrix.mat)                                                   

scMethrix.assays <- SimpleList()
scMethrix.assays[[basename(tabix.input.file)]] <- scMethrix.mat.sparse

toc()
tic("Generate summarized experiment")

scMethrix.exp <- SummarizedExperiment::SummarizedExperiment(assays = scMethrix.assays, rowRanges = scMethrix.rowRanges)
rm(scMethrix.rowRanges,scMethrix.assays)

saveHDF5SummarizedExperiment(scMethrix.exp, hdf5.output.folder, as.sparse=TRUE, replace=TRUE)
rm(scMethrix.exp)

toc()
tic("Load the HDF5 file")

scMethrix.hdf5.exp <- loadHDF5SummarizedExperiment(dir=hdf5.output.folder, prefix="")
hdf5.mat.delayed <- assays(scMethrix.hdf5.exp)[[1]]

toc()
tic("Do ops on hdf5")

c <- colSums(hdf5.mat.delayed)
r <- rowSums(hdf5.mat.delayed)

toc()
tic("Do ops on dgCMatrix")

c <- colSums(scMethrix.mat.sparse)
r <- rowSums(scMethrix.mat.sparse)

toc()
toc()

### Generate semi-sparse object (includes ranges)############

scMethrix.sites <- dplyr::pull(scMethrix.data, chr)
scMethrix.cells <- names(scMethrix.data)[2:length(scMethrix.data)]
scMethrix.mat <- data.matrix(dplyr::select(scMethrix.data, -chr), rownames.force = NA)

scMethrix.data.sparse <- Matrix(scMethrix.mat, dimnames = list(scMethrix.sites,scMethrix.cells), sparse=TRUE)

object.size(scMethrix.data)
object.size(scMethrix.data.sparse)


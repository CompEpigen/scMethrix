list.of.packages <- c("dplyr","magrittr","pryr","recommenderlab","GenomicRanges","SummarizedExperiment",
                      "DelayedArray","HDF5Array","tictoc","readr","sos","rstudioapi","sparseMatrixStats",
                      "RangedSummarizedExperiment")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

input.file <- "..\\test.data\\s50cell\\Parsed\\s50cell.tsv"

scMethrix.data <- read_tsv(input.file, col_names = TRUE)
scMethrix.mat <- data.matrix(dplyr::select(scMethrix.data, -chr,-start,-end), rownames.force = NA)
scMethrix.rowRanges <- GenomicRanges::makeGRangesFromDataFrame(scMethrix.data)
scMethrix.mat.sparse <- Matrix(scMethrix.mat, sparse=TRUE)
scMethrix.assays <- SimpleList()
scMethrix.assays[[basename(input.file)]] <- scMethrix.mat.sparse

se.exp <- SummarizedExperiment::SummarizedExperiment(assays = scMethrix.assays, colData = dimnames(scMethrix.mat)[[2]], rowRanges = scMethrix.rowRanges)

scMethrix.exp <- create_scMethrix(methyl_mat=scMethrix.mat.sparse,
                            colData = dimnames(scMethrix.mat.sparse)[[2]], 
                            rowRanges = scMethrix.rowRanges)

scMethrix.exp
dim(scMethrix.exp)
dimnames(scMethrix.exp)
assayNames(scMethrix.exp)
head(assay(scMethrix.exp))
#assays(scMethrix.exp) <- endoapply(assays(rse), asinh)
head(assay(scMethrix.exp))

rowRanges(scMethrix.exp)
rowData(scMethrix.exp)  # same as 'mcols(rowRanges(rse))'
colData(scMethrix.exp)




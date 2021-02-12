library("HelloRanges")
library("tictoc")
library("beepr")


setwd("D:\\Documents\\School\\Thesis\\scMethrix\\test.data\\100gen")

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\test.data\\10gen")

files <- list.files (getwd())

files <- files[1:3]


samples <- strsplit(files, "[.]")
samples <- lapply(samples, `[[`, 1)
samples <- unlist(samples)

R_bedtools_unionbedg(files,names=samples) #need to set names(i) <- samples

tic()
genome <- Seqinfo(genome = NA_character_)
i <- c("bed1.bedgraph", "bed10.bedgraph", "bed2.bedgraph")
names(i) <- samples
bl <- List(lapply(i, import, genome = genome))
gr_b <- stack(bl, "i")
dj <- disjoin(gr_b, ignore.strand = TRUE, with.revmap = TRUE)
mcols(gr_b)$i <- decode(mcols(gr_b)$i)
dfl <- extractList(mcols(gr_b), mcols(dj)$revmap)
assay <- as.matrix(dfl[, "score"], col.names = dfl[, 
                                                   "i"])
rowData <- granges(dj)
ans <- SummarizedExperiment(list(score = assay), rowData)

head(assay)

toc()
beep()

sse <- as(ans, "SingleCellExperiment")
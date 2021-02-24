library(BRGenomics)
library(SingleCellExperiment)

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\10gen")

files <- list.files (getwd(),full.names = TRUE)
files <- files[1:5]

files <- files[grepl("*.gz$", files)] #Select only gz
files <- files[grepl("*.bedgraph$", files)] #Select only gz

### Make a GRanges object from all ranges in files ##############

grange <- NULL

for (file in files) grange <- rbind(grange, data.table::fread(file, header=FALSE,nThread=8, select = c(1:3)))

grange <- unique(grange)
colnames(grange) <- c("chr","start","end")

data.gr <- makeGRangesFromDataFrame(grange)

#### Make a GRangeList object from all files ##############33

grl <- lapply(files,fread, select = c(1:3))
grl <- lapply(grl, setNames, c("chr","start","end"))
grl <- GRangesList(grl)
grl <- makeGRangesBRG(grl,ncores=1)

dl <- lapply(files,fread, select = c(4))
dl <- lapply(dl, setNames, c("score"))
names(dl) <- files

se <- SummarizedExperiment(assays=list(score=dl[[1]]), rowRanges=grl)

### Combine all Granges objects ####################

tic()

score_idx <- 4

data <- NULL

for (file in files) {
  data.new <- data.table::fread(file, header=FALSE,nThread=8,select=c(1:3,score_idx))
  colnames(data.new) <- c("chr","start","end",strsplit(file, "[.]")[[1]][1])
  data.new.gr <- makeGRangesFromDataFrame(data.new,keep.extra.columns=TRUE)
  data <- c(data,data.new.gr)
  rm(data.new)
}


### Use BRGenomics package

data <- BRGenomics::import_bedGraph(files,ncores=1)
data.gr <- GRangesList(data)
data <- makeGRangesBRG(data,ncores=1)
data.mg <- BRGenomics::mergeGRangesData(data,ncores = 1,multiplex=TRUE)

beep()

assignInNamespace(".multiplex_gr", ns = "BRGenomics",
                  function(data_in, field, ncores) {
                    # data must be *sorted*, base-pair resolution coverage data
                    if (all(vapply(data_in, function(x)
                      all(width(x) == 1L), logical(1L)))) {
                      data_in <- mclapply(data_in, sort, mc.cores = ncores)
                    } else {
                      warning(
                        .nicemsg(
                          "One or more inputs are not 'basepair resolution
                         GRanges' objects. Coercing them using
                         makeGRangesBRG()..."
                        ),
                        immediate. = TRUE
                      )
                      data_in <-
                        mclapply(data_in, makeGRangesBRG, mc.cores = ncores)
                    }
                    
                    # merge ranges
                    gr <- do.call(c, c(data_in, use.names = FALSE))
                    mcols(gr) <- NULL
                    gr <- unique(sort(gr))
                    
                    # (Fastest to keep these next steps separated, esp. for large datasets)
                    
                    # get dataframe of signal counts for each dataset in the output GRanges
                    idx <-
                      mclapply(data_in, function(x)
                        which(gr %in% x), mc.cores = ncores)
                    counts <- mcmapply(function(dat, idx, field) {
                      out <- rep.int(NA, length(gr))
                      out[idx] <- mcols(dat)[[field]]
                      out
                    },
                    data_in,
                    idx,
                    field,
                    mc.cores = ncores,
                    SIMPLIFY = TRUE)
                    
                    mcols(gr)[names(data_in)] <- counts
                    gr
                  })








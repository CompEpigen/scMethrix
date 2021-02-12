library(GenomicRanges)

setwd("D:\\Documents\\School\\Thesis\\scMethrix\\sample.data\\100gen")



files <- list.files (getwd())

files <- files[1:5]

data <- NULL

for (file in files) {

  data.new <- data.table::fread(file, header=FALSE,nThread=8, select = c(1:3))

  if (is.null(data)) {
    data <- data.new
  } else {
    data <- rbind(data, data.new)
  }

}

data <- unique(data)
colnames(data) <- c("chr","start","end")
data.gr <- makeGRangesFromDataFrame(data)

for (file in files) {

  data.new <- data.table::fread(file, header=FALSE,nThread=8)
  colnames(data.new) <- c("chr","start","end",strsplit(file, "[.]")[[1]][1])
  data.new.gr <- makeGRangesFromDataFrame(data.new)
  data.gr.overlap <- findOverlaps(data.gr,data.new.gr)

}







disjoin(data.gr)


rowRanges(data.gr)



width(data.gr)



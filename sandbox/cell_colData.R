colData <- readRDS("D:/Git/sampleData/metadata/colData.rds")
colData$Sample <- substr(colData$Sample,1,nchar(colData$Sample)-2)
colData <- unique(colData)

-----------------------------------------------------------

setwd("D:/Git/sampleData/metadata")
files <- list.files (getwd(),full.names = TRUE)
files <- files[grepl(".*txt$", files,ignore.case = TRUE)]

colData <- data.table()

for (file in files) {
  
  data <- suppressWarnings(data.table::fread(file,header = F,col.names = c("Sample")))
  data$Cell <- sub("\\_.*", "", get_sample_name(file))
  
  colData <- rbind(colData,data)

}

colData2 <- colData
colData$Sample <- paste0(colData$Sample,"_1")
colData2$Sample <- paste0(colData2$Sample,"_2")
colData <- rbind(colData,colData2)
rm(colData2)
colData <- colData[order(Sample),]

table(colData$Cell)

#---------------------------
#Derive from exp

scm = scm.impute

samp = row.names(colData(scm))
samp = sub("final_","",samp)
samp = sub(".sra","",samp)
samp = sub("_CpG","",samp)
samp = data.table(Sample=samp)

colData <- merge(samp,colData,all.x=TRUE)

row.names(colData(scm)) <- colData$Sample
colData(scm)$Cell = colData$Cell

table(colData$Cell)

#--- Generate 3 sample list -------------------------
lst <- colData
lst$Sample <- substr(lst$Sample,1,nchar(lst$Sample)-2)
lst <- unique(lst)
lst <- mPv <- lst[Cell=="mPv"]

lst <- mL23 <- lst[Cell=="mL2-3"][1:200,]
lst <- mL62 <- lst[Cell=="mL6-2"][1:200,]
lst <- rbindlist(list(mL62,mL23,mPv))
lst$Sample <- paste0("final_",lst$Sample,".sra_1*")

write.table(lst, file='3_cell_types.tsv', quote=FALSE, sep='\t', col.names = NA)
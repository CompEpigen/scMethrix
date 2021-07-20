library("GenomicRanges")
library("AnnotationHub")

# Get Grange of promoters
ah = AnnotationHub()
qhs = query(ah, c("RefSeq", "Mus musculus", "mm10"))
genes = qhs[[1]]
proms = promoters(genes)
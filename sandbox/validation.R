# Validation ------------------------------------------------------------------------------------------------------
library(scMethrix)

if(!requireNamespace("GEOquery")) {
  BiocManager::install("GEOquery")
}
library(GEOquery) 

data_dir = paste0("D:\\Git\\sampleData\\Vignette.GSE56879\\raw")
dir.create(path = data_dir, showWarnings = FALSE, recursive = TRUE)

# Download the GSM files from GEO ---------------------------------------
gsm_list <- c(paste0("GSM",1370535+0:11),  paste0("GSM",1370555+0:19))
#gsm_list <- c("GSM1370535", "GSM1370536", "GSM1370537", "GSM1370538", "GSM1370539", "GSM1370540", "GSM1370541", "GSM1370555", "GSM1370556", "GSM1370557", "GSM1370558", "GSM1370559", "GSM1370560", "GSM1370561")

bed_files <- sapply(gsm_list, function(gsm){
  GEOquery::getGEOSuppFiles(GEO = gsm, baseDir = data_dir, makeDirectory = FALSE, filter_regex = ".*.cov.txt.gz")
})

bed_files <- list.files(path = data_dir,  full.names = TRUE, pattern = ".*.cov.txt.gz")
file.rename(bed_files,gsub("(.*GSM.*?)_.*","\\1.bedgraph.gz",bed_files))
bed_files <- list.files(path = data_dir,  full.names = TRUE)

print(basename(bed_files))

#Genome of your preference to work with
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
if(!requireNamespace("BSgenome.Mmusculus.UCSC.mm10")) {
  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
}
library(BSgenome.Mmusculus.UCSC.mm10) 

soft <- GEOquery::getGEOfile("GSE56879")
soft <- GEOquery::getGEO(filename=soft)

GSE_colData <- sapply(names(soft@gsms), function(gsm) {
  soft@gsms[[gsm]]@header$title
})

GSE_colData <- gsub("_.*","",GSE_colData)
GSE_colData <- GSE_colData[which(names(soft@gsms) %in% gsm_list)]
GSE_colData <- data.frame(row.names = names(GSE_colData), Medium = GSE_colData)

colData(scm)$Medium[which(colData(scm)$Medium == "Ser")] <- "Serum" 
colData(scm)$Sample <- row.names(colData(scm))

res <- data.table::fread(bed_files[1])
colnames(res) <- c("chr","start","end","beta","M","U")
head(res)

# Import the BedGraph files into scMethrix ------------------------------
scm <- scMethrix::read_beds(
  files = bed_files, ref_cpgs = mm10_cpgs, h5 = T, h5_dir = "D:\\Git\\sampleData\\Vignette.GSE56879\\obj", 
  chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4, M_idx = 5, U_idx = 6, stranded = FALSE, 
  zero_based = FALSE, colData = GSE_colData, batch_size = 5
)

save_scMethrix(scm,dest="D:/Git/sampleData/Vignette.GSE56879/obj/scm.rds")
scm.raw <- load_scMethrix(dest="D:/Git/sampleData/Vignette.GSE56879/exp/scm.rds")
scm <- scm.raw

colData(scm)$Sample <- row.names(colData(scm))

setwd("D:/Documents/School/Thesis/Report/Validation")
# Mean score ------------------------------------------------------------------------------------------------------

#stats <- plot_stats(scm,base_size = 20)

row.names(colData(scm)) <- paste0("C",1:ncol(scm))
colData(scm)$Sample <- row.names(colData(scm))

# Plots for initial quality checks (Figure 4.) --------------------------
density <- plot_density(scm,base_size = 20,palette="Blues")
coverage <- plot_coverage(scm,base_size=20,type="density",max_cov=10)
sparsity <- plot_sparsity(scm,base_size=20,type="scatter",pheno="Sample")
stat <- plot_stats(scm, per_chr = T, stat = "count", base_size=20)

Cairo(file="stats.png",
      type="png",
      units="px", 
      width=1025, 
      height=300, 
      pointsize=24, 
      dpi="auto")

egg::ggarrange(density, sparsity, ncol = 2, widths = c(300,700),labels = c("  A", "  B"),label.args = list(gp = grid::gpar(fontface = "bold", fontsize =20)))  

dev.off()

Cairo(file="stats2.png",
      type="png",
      units="px", 
      width=1025, 
      height=300, 
      pointsize=24, 
      dpi="auto")

egg::ggarrange(stat, coverage, ncol = 2, widths = c(700,300),labels = c("  C", "  D"),label.args = list(gp = grid::gpar(fontface = "bold", fontsize =20)))  

dev.off()
# -----------------------------------------------------------------------------------------------------------------

get_stats(scm,per_chr=FALSE)

mask_scMethrix(scm, assay="score", threshold=1, by = "row", stat="count", op="<=")
mask_scMethrix(scm, assay="score", threshold=0.05, by = "row", stat="sd", op="<")
mask_scMethrix(scm, assay="counts", threshold=3, by = "row", stat="mean", op=">")

scm <- mask_scMethrix(scm, assay="score", threshold=1, by = "row", stat="count", op="<=")
scm <- remove_uncovered(scm)

scm <- mask_scMethrix(scm, assay="score", threshold=0.05, by = "row", stat="sd", op="<")
scm <- remove_uncovered(scm)

scm <- mask_scMethrix(scm, assay="counts", threshold=2, by = "row", stat="mean", op=">")
scm <- remove_uncovered(scm)

scm.qc <- scm

proms.mm10 <- list(hg38 = fread("D:\\Git\\sampleData\\Vignette.GSE56879\\mouse_epdnew_ktWrF.bed", header = F, col.names = c("chr","start","end","GeneSymbol"), select  = c(1:4)))
proms.mm10 <- makeGRangesFromDataFrame(proms.mm10[[1]],keep.extra.columns = T)
proms.mm10 <- promoters(proms.mm10, upstream=2000, downstream=200)
proms.mm10$GeneSymbol <- gsub("_1","",proms.mm10$GeneSymbol)

scm <- bin_scMethrix(scm.qc,regions = proms.mm10)
scm <- convert_HDF5_scMethrix(scm)

scm.bin <- scm

density <- plot_density(scm.b,base_size = 20,palette="Blues",show_legend = T)
coverage <- plot_coverage(scm,base_size=20,type="density",max_cov=10)
sparsity <- plot_sparsity(scm,base_size=20,type="scatter",pheno="Sample")

Cairo(file="stats.png",
      type="png",
      units="px", 
      width=1025, 
      height=300, 
      pointsize=24, 
      dpi="auto")

egg::ggarrange(density, coverage, sparsity, ncol = 3, widths = c(300,300,350),labels = c("  A", "  B", "  C"),label.args = list(gp = grid::gpar(fontface = "bold", fontsize =20)))  

dev.off()


scm <- impute_regions(scm,k=6)

scm.bin <- scm


# Get high and low fold features ----------------------------------------------------------------------------------




# Reduce to variable features -----------------------------------------------------------------------------------

scm <- reduce_scMethrix(scm,n_cpg = 1000)

scm1 <- scm[,which(colData(scm)$Medium == "2i")]
scm1 <- impute_regions(scm1,k=5)

scm2 <- scm[,which(colData(scm)$Medium != "2i")]
scm2 <- impute_regions(scm2,k=5)

scm <- merge_scMethrix(scm1,scm2,by="col")

anno = "Medium"

colData(scm)$Medium[which(colData(scm)$Medium == "Ser")] <- "Serum" 

scm <- dim_red_scMethrix(scm,assay="impute",type="tSNE",perplexity=50)
dimred <- plot_dim_red(scm,"tSNE",color_anno = anno,axis_labels=list(X="UMAP1",Y="UMAP2"),base_size=20,legend.title = element_blank(),legend.position="bottom",
                       #legend.justification="bottom",
                       legend.margin=margin(5,5,5,5),
                       legend.box.margin=margin(-18,0,0,0),
                       legend.box.background = element_blank(),
                       legend.background = element_blank(),
                       legend.text=element_text(size=20))
dimred

density.bin <- plot_density(scm.bin,base_size = 20,palette="Blues")
density.top <- plot_density(scm,base_size = 20,palette="Blues")

Cairo(file="stats3.png",
      type="png",
      units="px", 
      width=1050, 
      height=300, 
      pointsize=24, 
      dpi="auto")

egg::ggarrange(density.bin, density.top, dimred, ncol = 3, widths = c(300,300,300),labels = c("  A", "  B", "  C"),label.args = list(gp = grid::gpar(fontface = "bold", fontsize =20)))  

dev.off()

scm <- cluster_scMethrix(scm,assay="impute",n_clusters = 2,method="ward.D")
dimred <- plot_dim_red(scm,"tSNE",color_anno = "Medium",axis_labels=list(X="UMAP1",Y="UMAP2"),base_size=20,legend.title = element_blank(),legend.position="bottom",
                       #legend.justification="bottom",
                       legend.margin=margin(5,5,5,5),
                       legend.box.margin=margin(-18,0,0,0),
                       legend.box.background = element_blank(),
                       legend.background = element_blank(),
                       legend.text=element_text(size=20))
dimred

scm <- scm[,paste0("C",1:32)]

scm.red <- scm

generate_heatmap <- function(scm,assay = "score", type_anno="Sample", grouping = NULL, ...) {

  mat <- as.matrix(get_matrix(scm,assay))
  type <- colData(scm)[[type_anno]]
  ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=c("#FF0000","#0000FF")),
                         labels = c("Serum", "2i"),labels_gp = gpar(col = "white", fontsize = 20))#,
   # df = data.frame(type = type)
    #annotation_height = unit(8, "mm")
   )
  
  Heatmap(mat, name = "β-value", km = 5, top_annotation = ha,# cluster_columns = agnes(t(mat)), 
          show_row_names = FALSE, column_dend_side = "top", column_names_side = "top", show_column_names = FALSE,
          column_title = NULL, row_title = NULL, column_split = colData(scm)$Medium) 
}


hmap <- generate_heatmap(scm,type_anno="Medium",assay="impute") 

Cairo::Cairo(file="heatmap.png",
      type="png",
      units="px", 
      width=300, 
      height=300, 
      pointsize=24, 
      dpi="auto")

hmap

dev.off()





library(org.Mm.eg.db)

genes <- as.data.frame(findOverlaps(rowRanges(scm),proms.mm10))
genes <- genes[!duplicated(genes$queryHits),]
rowRanges(scm)$SYMBOL <- proms.mm10$GeneSymbol[genes$subjectHits]
genes <- data.table(select(org.Mm.eg.db,keys=rowRanges(scm)$SYMBOL,columns=c("SYMBOL","GO","ENTREZID","PATH"),keytype="SYMBOL"))
genes <- unique(genes, by = c("SYMBOL","GO"))



get_GO <- function(scm,n=10) {
  
  mtx <- data.table(get_matrix(scm,"impute"))
  col_idx <- which(colData(scm)$Medium == "Serum")
  rowData(scm)$B.serum <- rowMeans(mtx[,..col_idx])
  
  col_idx <- which(colData(scm)$Medium != "Serum")
  rowData(scm)$B.2i <- rowMeans(mtx[,..col_idx])
  rowData(scm)$fold <- foldchange(rowData(scm)$B.serum,rowData(scm)$B.2i)
  
  scm <- scm[(rowData(scm)$B.serum > 0.5 & rowData(scm)$B.2i < 0.5) | (rowData(scm)$B.serum < 0.5 & rowData(scm)$B.2i > 0.5),]
  
  fold <- rowData(scm)$fold
  fold <- fold[which(!is.infinite(fold))]
  low <- sort(fold)[1:10]
  low <- match(low,rowData(scm)$fold)
  
  high <- sort(fold,decreasing = T)[1:10]
  high <- match(high,rowData(scm)$fold)
  
  scm <- scm[c(high,rev(low)),]
  
  genes <- as.data.frame(findOverlaps(rowRanges(scm),proms.mm10))
  genes <- genes[!duplicated(genes$queryHits),]
  rowRanges(scm)$SYMBOL <- proms.mm10$GeneSymbol[genes$subjectHits]
  genes <- data.table(select(org.Mm.eg.db,keys=rowRanges(scm)$SYMBOL,columns=c("SYMBOL","GO","ENTREZID","PATH"),keytype="SYMBOL"))
  genes <- unique(genes, by = c("SYMBOL","GO"))
  
  mm.map <- org.Mm.egPATH
  mm.genes <- mappedkeys(mm.map)
  kegg <- as.list(mm.map[mm.genes])
  
  kegg[genes$ENTREZID]
  
  
  library(GO.db)
  select(GO.db, keys=genes$GO, columns="DEFINITION", keytype="GOID")
  
  
  
  
}









xx <- merge(genes,rowData(scm),by="SYMBOL")


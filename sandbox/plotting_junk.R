
cells <- c("Bcell", "CD4Tcell", "CD8Tcell", "Treg", "NKcell","Eosinophil", "Neutrophil",
           "Granulocyte", "Microglia", "Dendritic", "Monocyte", "Neuron", "Glia", "Glioma",
           "WholeBlood","Endothelial")

col_idx <- sapply(colData(scm)$Cell, `%in%`, cells)
scm <- scm[,col_idx]

group <- c(1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,6)
group <- c("Lymphocyte","Granulocytes","Monocytes","Brain","Blood","Endothelial")[group]
shapes <- c(21,22,23,24,25,21)[group] 
colors <- (as.numeric(ave(group, group, FUN = seq_along)))# - 1)# + shapes
colors <-  c('#ffe119', '#4363d8', '#3cb44b', '#e6194B', '#a9a9a9', '#000000')[colors]

names(shapes) <- cells
names(colors) <- cells

ord <- 1:length(cells)
names(ord) <- cells

colData(feat)$Shape <- shapes[colData(feat)$Cell]
colData(feat)$Color <- colors[colData(feat)$Cell]
colData(feat)$Order <- ord[colData(feat)$Cell]


#feat <- dim_red_scMethrix(feat,type="tSNE",assay="impute",perplexity=85,max_iter=5000)
plot_dim_red(feat,"tSNE",color_anno="Color",shape_anno = "Shape", axis_labels=list(Y="",X=paste0("Perplexity=",i)))
      
      
      
      + geom_mark_hull(expand=0.01,aes(fill=Cell)))
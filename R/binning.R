# trans <- c(score = function(x) mean(x,na.rm=TRUE),coverage = function(x) sum(x,na.rm=TRUE))

bin_scMethrix <- function(m, bin_size = 10000000, trans = NULL) {

  bins <- bin_granges(rowRanges(m),bin_size = bin_size)
  bins <- subsetByOverlaps(bins,rowRanges(m))
  
  sites <- sapply(1:length(bins),function (i) {
    (findOverlaps(rowRanges(m),bins[i]))@from
  })
  
  assys <- List()
  
  for (n in 1:length(assays(m))) {

    name <- assayNames(m)[n]
    
    tryCatch(
      expr = {op <- trans[[name]]},
      error = function(e){op <- function(x) mean(x,na.rm=TRUE)}
    )
    
    vals <- lapply(sites,function (i) {
      apply(get_matrix(m[i,],type=name),2,op)
    })
    
    setDT(vals, key=names(vals[[1]]))
    vals <- data.table::transpose(vals)
    colnames(vals) <- rownames(colData(m))
    
    assys[[name]] <- as(vals,typeof(get_matrix(m,type=name)))
  }
  
  if (is_h5(m)) {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = TRUE, 
                            h5_dir = h5_dir, genome_name = m@metadata$genome,desc = m@metadata$desc,colData = colData(m),
                            replace = replace)
  
  } else {
    
    m_obj <- create_scMethrix(assays = assys, rowRanges=bins, is_hdf5 = FALSE, 
                              genome_name = m@metadata$genome,desc = m@metadata$desc,colData = colData(m),)
    
  }
}
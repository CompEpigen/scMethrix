#####
##### OLD CODE FOR BINNING: ATTEMPTS AT BATCHING AND PARALLELIZATION
####


# mtx <- get_matrix(scm,assay=name)[overlap_indices$xid,] #TODO: Somehow missing rows if not subset, not sure why

# cl <- parallel::makeCluster(n_threads)
# doParallel::registerDoParallel(cl)
# parallel::clusterEvalQ(cl, c(library(DelayedMatrixStats)))
# parallel::clusterExport(cl,list('overlap_indices','worker','op'), envir = environment())
# on.exit(parallel::stopCluster(cl))

# 
# 
# idxs <- split(overlap_indices$xid, ceiling(seq_along(overlap_indices$xid)/ceiling(length(overlap_indices$xid)/n_chunks)))
# 
# dat <- do.call("rbind",lapply(idxs, function(idx) get_matrix(scm[idx,])))
# dat = cbind(overlap_indices, dat)
# 
# dat[, lapply(.SD, mean, na.rm = T), by = yid, .SDcols = colnames(scm)]
# 
# avg <- t(sapply(unique(overlap_indices[,yid]), function (rid) {
#   message("Parsing ",rid)
#   idx <- overlap_indices[yid == rid,]$xid
#   DelayedMatrixStats::colMeans2(mtx,row=idx,na.rm=TRUE)
# }))
# colnames(avg) <- colnames(mtx)




#------  This splits the matrix up into seperate data.tables and calculates, but ends up slower

# if (n_threads > 1) {
# 
#   mtx <- get_matrix(scm, assay=name, n_chunks = floor(ncol(scm)/batch_size), by="col")
# 
#   out = NULL
#   
#   if (verbose) message("Generated ", length(mtx), " chunks...")
#   split_time()
#   
#   for (i in 1:length(mtx)) {
#     sub.mtx <- mtx[[i]]
#     col_idx <- split_vector(1:ncol(sub.mtx),chunks=n_threads)
#     sub.mtx <- lapply(col_idx,function(idx) as.data.table(sub.mtx[,idx,drop=FALSE]))
#     sub.mtx <- parLapply(cl,sub.mtx,worker,op = op)
# 
#     out <- cbind(out,cbindlist(sub.mtx))
#     
#     rm(sub.mtx)
#     gc()
#     
#     if (verbose) message("   Processed chunk ",i," (",split_time(),")")
#   }


# mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] 
# cols <- split_vector(1:ncol(mtx),n_chunks)
# 
# cl <- parallel::makeCluster(n_threads)
# doParallel::registerDoParallel(cl)
# parallel::clusterEvalQ(cl, c(library(data.table)))
# parallel::clusterExport(cl,list('overlap_indices','worker','op'), envir = environment())
# on.exit(parallel::stopCluster(cl))
# 
# 
# mtx <- lapply(cols,function(col) mtx[,..col])
# #mtx <- lapply(mtx,worker,overlap_indices = overlap_indices,op = op)
# mtx <- parLapply(cl,mtx,worker,overlap_indices = overlap_indices,op=op)
# 
# 
# mtx <- setDT(unlist(mtx, recursive = FALSE))

### Split by rids
# worker <- function (mtx, yid, op) {
#   mtx <- mtx[,lapply(.SD,op),by=(yid)]
#   mtx <- mtx[,yid:=NULL]
#   return(mtx)
# }
# mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] 
# 
# #Make two chunked lists for xid and yids
# yid_list <- rle(overlap_indices$yid)
# yid_list <- split(overlap_indices$yid,ceiling(rep(1:length(yid_list$lengths), yid_list$lengths) /
#                                                 (length(yid_list$lengths)/n_chunks)))
# lens <- sapply(yid_list,function(yids) length(yids))      
# xid_list <- split(overlap_indices$xid,rep(1:length(lens),lens))
# 
# cl <- parallel::makeCluster(n_threads)  
# doParallel::registerDoParallel(cl) 
# parallel::clusterEvalQ(cl, c(library(data.table)))
# parallel::clusterExport(cl,list('worker'), envir = environment())
# on.exit(parallel::stopCluster(cl))
# 
# mtx <- clusterMap(cl,worker,mtx=lapply(xid_list,function(xid) mtx[unlist(xid),]), yid = yid_list, MoreArgs = list(op = op))
# mtx <- rbindlist(lapply(mtx, as.data.frame.list))

### Mapply algorithm 
# mtx <- mapply(worker,mtx=lapply(xid_list,function(xid) mtx[unlist(xid),]), yid = yid_list, MoreArgs = list(op = op),simplify=TRUE)
# cols <- row.names(mtx)
# mtx <- lapply(1:nrow(mtx),function(x) unlist(mtx[x,]))
# mtx <- t(rbindlist(lapply(mtx, as.data.frame.list)))
# colnames(mtx) <- cols

### Basic algorithm

if (batch_size != ncol(scm)) {
  
  # } else {
  cols <- split_vector(1:ncol(scm),size=batch_size)
  
  if (verbose) message("Generated ", length(cols), " chunks...")
  
  out <- NULL
  
  for (i in 1:length(cols)) {
    
    col <- cols[[i]]
    mtx <- as.data.table(get_matrix(scm,assay=name)[overlap_indices$xid,col])
    mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
    mtx <- mtx[,overlap_indices:=NULL]
    
    out <- cbind(out,mtx)
    
    rm(mtx)
    gc()
    
    if (verbose) message("   Processed chunk ",i," (",split_time(),")")
  }
} else {
  

#------------------------------------------------------------------------------------------------------------
#' Transforms an assay in an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object. The function is
#' applied column-wise as to optimize how HDF5 files access sample data. 
#' 
#' If HDF5 objects are used, transform functions need to accept 'DelayedMatrix' (e.g., from \pkg{DelayedMatrixStats}).
#' Otherwise, 
#' @inheritParams generic_scMethrix_function
#' @param h5_temp string; temporary directory to store the temporary HDF5 files
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' transform_assay(scMethrix_data,assay="score",new_assay="plus1",trans=function(x){x+1})
#' @export
transform_assay <- function(scm, assay = "score", new_assay = "new_assay", trans = NULL, h5_temp = NULL) {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  .validateType(new_assay,"string")
  .validateType(trans, "function")
  .validateType(h5_temp, c("string","null"))
  
  if (!.validateAssay(scm,new_assay,check.absent=T))
    #new_assay %in% SummarizedExperiment::assayNames(scm)) 
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)

  #- Function code -----------------------------------------------------------------------------
  if (is_h5(scm)) {
    
    if (is.null(h5_temp)) {h5_temp <- tempdir()}
    
    grid <- DelayedArray::RegularArrayGrid(refdim = dim(scm),
                                           spacings = c(length(scm), 1L)) 
    
    trans_sink <- HDF5Array::HDF5RealizationSink(dim = dim(scm),
                                                 dimnames = list(NULL,sampleNames(scm)), type = "double",
                                                 filepath = tempfile(pattern="trans_sink_",tmpdir=h5_temp),
                                                 name = new_assay, level = 6)
    
    blocs <- DelayedArray::blockApply(get_matrix(scm,assay=assay), grid = grid, FUN = trans)
    
    for(i in 1:length(blocs)) {
      DelayedArray::write_block(block = as.matrix(blocs[[i]]), viewport = grid[[as.integer(i)]], sink = trans_sink)
    }
    
    rm(blocs)
    mtx <- as(trans_sink, "HDF5Matrix")
    
  } else {
    mtx <- get_matrix(scm,assay=assay)
    dims <- dimnames(mtx)
    mtx <- as.data.table(mtx)
    mtx[, names(mtx) := lapply(.SD, trans)]
    mtx <- as(mtx, "matrix")
    dimnames(mtx) <- dims
  }
  
  assays(scm)[[new_assay]] <- mtx
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Bins the ranges of an \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object. Typically, most assays will use either mean (for measurements) or sum (for counts). The transform is applied column-wise to optimize how HDF5 files access sample data. If HDF5 objects are used, transform functions should be  from \pkg{DelayedMatrixStats}.
#' 
#' In the output object, the number of CpGs in each region is saved in mcol(scm)$n_cpgs.
#' 
#' Reduced dimensionality data will be discarded.
#' @inheritParams generic_scMethrix_function
#' @param regions Granges; The regions from which to make the bins.
#' @param bin_size integer; The size of each bin. First bin will begin at the start position of the first genomic
#' region on the chromosome. If NULL, there will be one bin per region. Default 100000.
#' @param bin_by character; can create bins by # of base pairs "bp" or by # of CpG sites "cpg". Default "bp"
#' @param trans named vector of closures; The transforms for each assay in a named vector. Default NULL, meaning that 
#' operations for "counts" assay is sum(x, na.rm=TRUE), and for all other assays is mean(x, na.rm=TRUE)
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' regions <- GRanges(seqnames = c("chr1"), ranges = IRanges(1,200000000)) 
#' regions <- unlist(tile(regions,10))
#' bin_scMethrix(scMethrix_data, regions = regions)
#' @export
bin_scMethrix <- function(scm = NULL, regions = NULL, bin_size = 100000, bin_by = c("bp","cpg"), trans = NULL, 
                          overlap_type = c("within", "start", "end", "any", "equal"), h5_dir = NULL, verbose = TRUE, 
                          batch_size = 20, n_threads = 1, replace = FALSE) {
  #- Input Validation --------------------------------------------------------------------------
  yid <- NULL

  .validateExp(scm)
  .validateType(regions,c("granges","null"))
  .validateType(bin_size,c("integer","null"))
  bin_by <- .validateArg(bin_by,bin_scMethrix)
  sapply(trans, function (t) .validateType(t, c("function","null")))
  overlap_type <- .validateArg(overlap_type,subset_scMethrix)
  if (is_h5(scm)) .validateType(h5_dir,"string")
  .validateType(verbose,"boolean")
  .validateType(batch_size,"integer")
  n_threads <- .validateThreads(n_threads)
  .validateType(replace,"boolean")

  if (is.null(trans[["counts"]])) {
    trans <- c(trans, counts = function(x) sum(x,na.rm=T))}
 
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Binning experiment...")
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
    scm <- subset_scMethrix(scm, regions = regions) 
  } else { # If no region is specifed, use entire chromosomes
    if (verbose) message ("No regions specified; using whole chromosomes")
    regions = range(SummarizedExperiment::rowRanges(scm))
  }
  
  if (verbose) message("Checking ",length(regions)," input regions...",start_time())

  regions$rid <- paste0("rid_", 1:length(regions))
  
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) 
  
  if(nrow(overlap_indices) == 0){
    stop("No overlaps detected for input regions. No binning can be done.")
  }
  
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]

  if (!is.null(bin_size)) {
    
    #if (verbose) message("Binning by ",bin_by," with size of ",bin_size,"...")
    
    if (bin_by == "cpg") {
      
      rrng = GenomicRanges::GRanges()
      
      for(rid in 1:length(regions$rid)) {
        
        if (verbose) message("   Processing region ",rid,"/",length(regions$rid))
        
        idx <- which(overlap_indices[,yid == regions$rid[rid]])
        idx <- split_vector(idx,size=bin_size)

        for(i in idx) {
          gr <- range(rowRanges(scm[c(i[1],i[length(i)]),]))
          gr$n_cpgs <- length(i)
          rrng <- c(rrng,gr)
        }
      }
      
    } else if (bin_by == "bp") {

      rrng <- GenomicRanges::slidingWindows(regions,width=bin_size,step=bin_size)
      rrng <- unlist(as(rrng, "GRangesList"))
      
      # rrng <- bin_granges(rrng,bin_size = bin_size) #sort(unlist(tile(regions, width = bin_size))) #TODO: Should switch this to using RLE lookup
      
      idx <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
      idx <- idx[order(idx$subjectHits),]
      
      rrng <- rrng[unique(idx$subjectHits)]
      rrng$n_cpgs <- rle(idx$subjectHits)$lengths
    }
    
  } else { # If no bin_size is specified, use the entire region
    rrng <- regions
  }
  
  n_regions <- length(regions)
  rm(regions)
  
  assays <- list()
  
  rrng$rid <- paste0("rid_", 1:length(rrng))
  overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
  colnames(overlap_indices) <- c("xid", "yid")
  overlap_indices[,yid := paste0("rid_", yid)]
  
  if (verbose) message("Generated ",length(rrng)," bins in ",split_time())
  
  gc()

  # Function to process each bin
  
  for (name in SummarizedExperiment::assayNames(scm)) {

    if (verbose) message("Filling bins for the ",name," assay...")
    
    if (is.null(trans[[name]])) { # If no named vector is specified, default to mean
      op <- function(x) mean(x,na.rm=TRUE)#DelayedMatrixStats::colMeans2(x,na.rm=TRUE)
    } else {
      op <- trans[[name]]
    }
    
    if (is_h5(scm)) {

      DelayedArray::setAutoRealizationBackend("HDF5Array")
      
      cols <- split_vector(1:ncol(scm),size=batch_size)
      sink <- DelayedArray::AutoRealizationSink(c(length(unique(overlap_indices$yid)),ncol(scm)))
      grid <- DelayedArray::ArbitraryArrayGrid(list(length(unique(overlap_indices$yid)),cumsum(lengths(cols))))

      if (verbose) message("Generated ", length(cols), " chunks...")
      
      for (i in 1:length(cols)) {
        
        split_time()
        col <- cols[[i]]
        mtx <- as.data.table(get_matrix(scm,assay=name)[overlap_indices$xid,col])
        mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
        mtx <- mtx[,overlap_indices:=NULL]
        
        DelayedArray::write_block(block = as.matrix(mtx), viewport = grid[[as.integer(i)]], sink = sink)
        
        if (verbose) message("   Processed chunk ",i," (",split_time(),")")
      }
      
      assays[[name]] <- sink
      
      
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
      
    } else {

      ### Split by cols
      # worker <- function (mtx,overlap_indices,op) {
      #   mtx <- mtx[,lapply(.SD,op),by=(overlap_indices$yid)]
      #   mtx <- mtx[,overlap_indices:=NULL]
      #   return(mtx)
      # }
      # 
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

      mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] #TODO: Somehow missing rows if not subset, not sure why
      mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
      mtx <- mtx[,overlap_indices:=NULL]

      assays[[name]] <- as(mtx,class(get_matrix(scm,assay=name)))
    }
    
    if (verbose) message("Bins filled in ",split_time())
  }

  rrng <- rrng[which(rrng$rid %in% overlap_indices$yid)]
  
  rrng$rid <- NULL

  if (verbose) message("Rebuilding experiment...")
  
  if (is_h5(scm)) {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
                              replace = replace)  
  } else {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
  }
  
  if (verbose) message("Experiment binned for ", n_regions," regions containing ", length(m_obj)," total bins in ", stop_time())
  
  # parallel::stopCluster(cl)
  # rm(cl)
  
  return (m_obj)
}


























# 
# 
# 
# bin_scMethrix2 <- function(scm = NULL, regions = NULL, bin_size = 100000, bin_by = "bp", trans = NULL, 
#                           overlap_type = "within", h5_dir = NULL, verbose = TRUE, n_chunks = 1, n_threads = 1) {
#   
#   yid <- NULL
#   
#   if (is.null(trans)) {
#     trans <- c(counts = function(x) sum(x,na.rm=TRUE))
#   }
#   
#   if (!is(scm, "scMethrix")) {
#     stop("A valid scMethrix object needs to be supplied.", call. = FALSE)
#   }
#   
#   # if (is_h5(scm) && is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
#   
#   bin_by = match.arg(arg = bin_by, choices = c("bp","cpg"))
#   
#   if (!is.null(regions)) {
#     regions = cast_granges(regions)
#     scm <- subset_scMethrix(scm, regions = regions) 
#   } else { # If no region is specifed, use entire chromosomes
#     regions = range(rowRanges(scm))
#   }
#   
#   if (verbose) message("Subsetting for ",length(regions)," input regions...",start_time())
#   
#   regions$rid <- paste0("rid_", 1:length(regions))
#   
#   overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type)) 
#   
#   if(nrow(overlap_indices) == 0){
#     stop("No overlaps detected")
#   }
#   
#   colnames(overlap_indices) <- c("xid", "yid")
#   overlap_indices[,yid := paste0("rid_", yid)]
#   
#   if (!is.null(bin_size)) {
#     
#     if (verbose) message("Binning by ",bin_by," with size of ",bin_size,"...")
#     
#     if (bin_by == "cpg") {
#       
#       rrng = GRanges()
#       
#       for(rid in 1:length(regions$rid)) {
#         
#         if (verbose) message("   Processing region ",rid,"/",length(regions$rid))
#         
#         idx <- which(overlap_indices[,yid == regions$rid[rid]])
#         idx <- split_vector(idx,num = bin_size, by = "size")
#         
#         for(i in idx) {
#           gr <- range(rowRanges(scm[c(i[1],i[length(i)]),]))
#           gr$n_cpgs <- length(i)
#           rrng <- c(rrng,gr)
#         }
#       }
#       
#     } else if (bin_by == "bp") {
#       
#       rrng <- unlist(tile(regions, width = bin_size)) #TODO: Should switch this to using RLE lookup
#       
#       idx <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
#       
#       rrng <- rrng[unique(idx$subjectHits)]
#       rrng$n_cpgs <- rle(idx$subjectHits)$lengths
#     }
#     
#   } else { # If no bin_size is specified, use the entire region
#     rrng <- regions
#   }
#   
#   if (verbose) message("Generated ",length(rrng)," bins in ",split_time())
#   
#   assays <- list()
#   
#   rrng$rid <- paste0("rid_", 1:length(rrng))
#   overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, rrng, type = overlap_type))
#   colnames(overlap_indices) <- c("xid", "yid")
#   overlap_indices[,yid := paste0("rid_", yid)]
#   
#   # Function to process each bin
#     for (name in SummarizedExperiment::assayNames(scm)) {
#       
#       if (verbose) message("   Filling bins for the ",name," assay...")
#       
#       
#       if (is.null(trans[[name]])) { # If no named vector is specified, default to mean
#         op <- function(x) mean(x,na.rm=TRUE)
#       } else {
#         op <- trans[[name]]
#       }
#       
#       if (is_h5(scm)) {
#         
#         rid_list <- split_vector(rrng$rid,num = n_threads)
#         
#         
#         
#         
#         
#         
#         
#         cl <- parallel::makeCluster(n_threads)  
#         doParallel::registerDoParallel(cl) 
#         parallel::clusterEvalQ(cl, c(library(data.table)))
#         parallel::clusterExport(cl,list('overlap_indices','bin','mtx','yid'), envir = environment())
#         
#       } else {
#         
#         mtx <- data.table(get_matrix(scm,assay=name))[overlap_indices$xid,] #TODO: Somehow missing rows if not subset, not sure why
#         mtx <- mtx[,lapply(.SD,op),by=overlap_indices$yid]
#         mtx <- mtx[,overlap_indices:=NULL]
# 
#         assays[[name]] <- as(mtx,class(get_matrix(scm,assay=name)))
#         
#       }
#       
#       if (verbose) message("Bins filled in ",split_time())
#     }
#   
#   rrng$rid <- NULL
#   
#   if (verbose) message("Rebuilding experiment...")
#   
#   if (is_h5(scm)) {
#     m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = TRUE, 
#                               h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),
#                               replace = replace)  
#   } else {
#     m_obj <- create_scMethrix(assays = assays, rowRanges=rrng, is_hdf5 = FALSE, 
#                               genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData(scm),)
#   }
#   
#   if (verbose) message("Experiment binned for ",length(regions)," regions containing ",length(m_obj)," total bins in ",stop_time())
#   
#   # parallel::stopCluster(cl)
#   # rm(cl)
#   
#   return (m_obj)
# }
# 

#' Collapses multiple samples into a single sample by group
#' @details 
#' Multiple samples can be collapsed into a single meta-sample. Grouping for samples can be defined via colData. The collapse function can accept an arbitrary function for each assay on how to handle the collapsing (typically `mean` for scores, and `sum` for counts).
#' 
#' In the output object, `colData()` will contain a comma-delimited list of samples (`Samples`) that each group contains as well as the total number of CpGs in the group (`n_Samples`).
#' 
#' Reduced dimensionality data will be discarded.
#' 
#' @inheritParams generic_scMethrix_function
#' @param colname string; The colname from \code{colData(scm)} indicating which samples should be collapse together
#' @param trans named vector of closures; The transforms for each assay in a named vector. Default NULL, meaning that 
#' operations for "counts" assay is sum(x, na.rm=TRUE), and for all other assays is mean(x, na.rm=TRUE)
#' @param batch_size The number of CpGs to calculate at once.
#' \code{\link{IRanges}} package.
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' colData(scMethrix_data)["Cluster"] = c("X","X","Y","Y")
#' collapse_samples(scMethrix_data, colname = "Cluster")
#' @export
collapse_samples <- function(scm = NULL, colname = NULL, trans = NULL, h5_dir = NULL, batch_size = 100000, n_threads = 1, replace = FALSE, verbose = TRUE) {
  
  #- Input Validation --------------------------------------------------------------------------
  Group <- NULL
  
  .validateExp(scm)
  .validateType(colname,"string")
  sapply(trans, function (t) .validateType(t, c("function","null")))
  if (is_h5(scm)) .validateType(h5_dir,"string")
  .validateType(batch_size,"integer")
  n_threads <- .validateThreads(n_threads)
  .validateType(verbose,"boolean")
  .validateType(replace,"boolean")
  
  if (!(colname %in% colnames(colData(scm)))) 
    stop("Cannot find column `",colname,"` in colData (Avail: ",paste(names(colData(scm)),collapse=", "),")")
  
  #if (is.null(h5_dir) && is_h5(scm)) stop("Output directory must be specified")
  
  if (is.null(trans[["counts"]])) 
    trans <- c(trans, c(counts = function(x) rowSums(x,na.rm=TRUE)))#DelayedMatrixStats::rowSums2(x,na.rm=TRUE)))
  
 # if (any(sapply(trans, function (x) {!is(x, "function")}))) stop("Invalid operation in trans")
  
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Starting to collapse experiment...",start_time())
  
  assays <- list()
  overlaps_indicies <- data.table(Sample = sampleNames(scm), Group = factor(scm@colData[,colname]))

  for (name in SummarizedExperiment::assayNames(scm)) {
    
    if (verbose) message("   Collapsing samples for the ",name," assay...")
    
    if (is.null(trans[[name]])) { # If no named vector is specified, default to mean
      op <- function(x) rowMeans(x,na.rm=TRUE)#DelayedMatrixStats::rowMeans2(x,na.rm=TRUE)
    } else {
      op <- trans[[name]]
    }
    
    if (is_h5(scm)) {
      
      DelayedArray::setAutoRealizationBackend("HDF5Array")
      
      grps <- split(1:nrow(overlaps_indicies), overlaps_indicies$Group)
      cpgs <- split_vector(1:nrow(scm),size = batch_size)
      sink <- DelayedArray::AutoRealizationSink(c(length(rowRanges(scm)),uniqueN(overlaps_indicies$Group)))
      grid <- DelayedArray::RegularArrayGrid(dim(sink), spacings = c(length(rowRanges(scm)),1))

      if (verbose) message("Generated ", length(grps)*length(cpgs), " chunks")
      chunk = 0
      
      for (i in 1:length(grps)) {
        
        col <- NULL

        for (cpg in cpgs) {
          if (verbose) message("Processing chunk ", chunk <- chunk+1)
          mtx <- get_matrix(scm,assay=name)[cpg,grps[[i]]]
          mtx <- as.matrix(op(mtx))
          col <- rbind(col,mtx)
        }
        
        colnames(col) <- names(grps[i])
        DelayedArray::write_block(block = col, viewport = grid[[as.integer(i)]], sink = sink)

      }
      
      assays[[name]] <- sink
      
    } else {
      
      mtx <- data.table(get_matrix(scm,assay=name)) #TODO: Somehow missing rows if not subset, not sure why
      mtx <- data.table::setDT(lapply(split.default(mtx, overlaps_indicies$Group), op))[]

      assays[[name]] <- as(mtx,class(get_matrix(scm,assay=name)))
      
    }

  }
  
  colData <- data.frame(row.names = levels(overlaps_indicies$Group),
                        Samples= sapply(levels(overlaps_indicies$Group),function(grp) {
                          paste(overlaps_indicies[Group == grp]$Sample,collapse=",")}),
                        n_Samples = as.vector(table(overlaps_indicies$Group)))
  
  
  if (verbose) message("Rebuilding experiment...")
  
  if (is_h5(scm)) {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rowRanges(scm), is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData,
                              replace = replace)  
  } else {
    m_obj <- create_scMethrix(assays = assays, rowRanges=rowRanges(scm), is_hdf5 = FALSE, 
                              genome_name = scm@metadata$genome,desc = scm@metadata$desc,colData = colData,)
  }
  
  if (verbose) message("Experiment collapsed into ", nrow(colData)," sample groups in ",stop_time())
  
  return (m_obj)
  
}

#------------------------------------------------------------------------------------------------------------
#' Imputes the NA values of a \code{\link{scMethrix}} object.
#' @details Uses the inputted function to transform an assay in the \code{\link{scMethrix}} object
#' @param threshold The value for cutoff in the "score" assay to determine methylated or unmethylated status. 
#' Default = 50
#' @inheritParams generic_scMethrix_function
#' @return An \code{\link{scMethrix}} object
#' @examples
#' data('scMethrix_data')
#' @export
#' @import Melissa
#' @references Kapourani CA, Sanguinetti G (2019). “Melissa: Bayesian clustering and imputation of single cell methylomes.” Genome Biology, 20, 61. doi: 10.1186/s13059-019-1665-8.
impute_by_melissa <- function (scm, threshold = 50, assay = "score", new_assay = "impute") {
  
  #- Input Validation --------------------------------------------------------------------------
  . <- NULL
  
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    if (new_assay == "score") stop("Cannot overwrite the score assay")
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) warning("Imputation cannot be done on HDF5 data. Data will be cast as matrix for imputation.")
  
  scm <- transform_assay(scm, assay = assay, new_assay = "binary", trans = binarize)

  #- Function code -----------------------------------------------------------------------------
  # Convert Granges to genomic interval [-1,1]
  chrom_size <- scm@metadata$chrom_size
  chrom_start <- sapply(coverage(scm), function(x) {x@lengths[1]}+1)
  interval <- as.data.table(rowRanges(scm))[,.(seqnames,start)]
  
  for (i in 1: length(chrom_size)) {
    # The final square bracket at line below is added due to bug:
    # https://stackoverflow.com/questions/34667536/have-to-call-variable-twice-before-evaluated
    interval[seqnames==names(chrom_start)[i], interval := round((2*((start-chrom_start[i])/chrom_size[i]))-1, digits = 7)][]
  }
  
  mcols(scm) <- cbind(mcols(scm),interval = interval[,interval])
  
  # Put into Melissa met format
  loc <- as.matrix(rowRanges(scm)@seqnames)
  loc <- matrix(loc,dimnames = list(paste0(rowRanges(scm)@seqnames,":",interval$start)))
  cells <- list()
  
  for (n in 1:ncol(scm)) {
    
    met <- matrix(c(mcols(scm)$interval,as.matrix(get_matrix(scm[,n],assay="binary"))),
                  ncol = 2,dimnames = list(rownames(loc),NULL))
    mets <- lapply(rowRanges(scm)@seqnames@values, function(x) met[which(loc == x),,drop=FALSE])
    mets <- lapply(mets, function (m) m[-which(m[,2]==-1),,drop=FALSE])
    
    cells[[colnames(scm)[n]]] <- mets
  }
  
  # # Parameter options
  opts <- list()
  # # Load cell filenames
  # opts$met_files <- NULL
  # opts$cell_id <- colnames(scm)
  # opts$is_centre  <- is_centre   # Whether genomic region is already pre-centred
  # opts$is_window  <- is_window   # Use predefined window region
  # opts$upstream   <- upstream    # Upstream of centre
  # opts$downstream <- downstream  # Downstream of centre
  # opts$chrom_size <- chrom_size_file  # Chromosome size file
  # opts$chr_discarded <- chr_discarded # Chromosomes to discard
  # opts$cov        <- cov         # Regions with at least n CpGs
  # opts$sd_thresh  <- sd_thresh   # Variance of methylation within region
  # 
  melissa_obj <- structure(list(met = cells, anno_region = NULL, opts = opts),
                           class = "melissa_data_obj")
  
  # Do the imputation
  basis_obj <- BPRMeth::create_rbf_object(M = 4)
  test_obj <- Melissa::partition_dataset(melissa_obj)
  
  set.seed(123)
  melissa_obj <- Melissa::melissa(X = melissa_obj$met, basis = basis_obj,K = min(ncol(scm)-1,2))#,
                                  #vb_max_iter = 30, vb_init_nstart = 1, 
                                  #is_parallel = FALSE)
  # plot_melissa_profiles(melissa_obj = melissa_obj, region = 1, 
  #                       title = "Methylation profiles for region 25")
  
  
  
  imputation_obj <- Melissa::impute_test_met(obj = melissa_obj, test = test_obj)
                                    
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Generic imputation return function
#' @details Uses the specified imputation operation to evaluation an scMethrix object.
#' @param regions Granges; the regions to impute. Default is by chromosome.
#' @param type string/closure; the imputation to perform "kNN","iPCA",or "RF". Otherwise, a closure can be specified that returns the imputed matrix. Default = "kNN"
#' @param n_pc the range of principal components to check when using iPCA. Caution: this can be very time-intensive
#' @inheritParams generic_scMethrix_function
#' @inheritParams impute::impute.knn
#' @inheritParams missForest::missForest
#' @inheritParams missMDA::imputePCA
#' @return list; two \code{\link{scMethrix}} objects names 'training' and 'test'
#' @examples
#' data('scMethrix_data')
#' impute_regions(scMethrix_data)
#' @export
#' @references Hastie T, Tibshirani R, Narasimhan B, Chu G (2021). impute: impute: Imputation for microarray data. R package version 1.66.0.
#' @references Stekhoven, D. J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
#' @references Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. Analytical and Bioanalytical Chemistry, 5, 1241-1251.
#' @references Josse, J. and Husson, F. (2011). Selecting the number of components in PCA using cross-validation approximations. Computational Statistics and Data Analysis. 56 (6), pp. 1869-1879.
impute_regions <- function(scm = NULL, assay="score", new_assay = "impute", regions = NULL, n_chunks = 1, 
                               n_threads = 1, overlap_type=c("within", "start", "end", "any", "equal"), type=c("kNN","iPCA","RF"), verbose = TRUE, k=10, n_pc=2,...) {
 
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  assay <- .validateAssay(scm,assay)
  .validateType(new_assay,"string")
  .validateType(regions,c("granges","null"))
  .validateType(n_chunks,"integer")
  overlap_type <- .validateArg(overlap_type,impute_regions)
  if (!.validateType(type,"function",throws=F)){ 
    .validateType(type,"String")
    type = .validateArg(type,impute_regions)
  }
  .validateType(verbose,"boolean")
  .validateType(k,"integer")
  .validateType(n_pc,"integer")
  
  if (new_assay %in% SummarizedExperiment::assayNames(scm)) {
    if (new_assay == "score") stop("Cannot overwrite the score assay")
    warning("Name already exists in assay. It will be overwritten.", call. = FALSE)
  }
  
  if (is_h5(scm)) warning("Imputation cannot be done on HDF5 data. Data will be cast as matrix for imputation.")

  yid <- NULL
  
  #- Function code -----------------------------------------------------------------------------
  if (verbose) message("Starting imputation...",start_time())
  
  if (.validateType(type,"function",throws=F)) {
    op = type
  } else if (type == "kNN") {
    op <- function(mtx) impute::impute.knn(mtx, k = min(k,ncol(mtx)), 
                                           rowmax = 1.0, colmax = 1.0, maxp = 1500, ...)$data
  } else if (type == "iPCA") {
    if (length(n_pc) > 1) {
      warning("Caution: n_pc is given as range. This can be very time-intensive.")
      n_pc <- missMDA::estim_ncpPCA(as.matrix(get_matrix(scm,assay = assay)),ncp.min = n_pc[1], ncp.max = n_pc[2], 
                                    method.cv = "Kfold", verbose = TRUE)
      n_pc <- n_pc$ncp
    }
    
    op <- function(mtx) missMDA::imputePCA(mtx, ncp = n_pc, ...)$completeObs
  } else if (type == "RF") {
    op <- function(mtx) missForest::missForest(mtx, ...)$ximp
  } else {
    stop("Error in imputation. No valid algorithm specified. This should never be reached.")
  }
  
  if (!is.null(regions)) {
    regions = cast_granges(regions)
    scm <- subset_scMethrix(scm, regions = regions) 
    
    assays(scm)[[new_assay]] <- assays(scm)[[assay]]

    regions$rid <- paste0("rid_", 1:length(regions))
    overlap_indices <- as.data.table(GenomicRanges::findOverlaps(scm, regions, type = overlap_type))
    colnames(overlap_indices) <- c("xid", "yid")
    overlap_indices[,yid := paste0("rid_", yid)]
    rid_list <- regions$rid
    # rid_list <- split_vector(regions$rid,num = n_chunks)
    # rid_list <- lapply(rid_list, function(rids) {split_vector(rids,num = n_threads)})
  
    for (i in 1:length(rid_list)) {
  
      if (verbose) message("Parsing chunk ", i, " of ",length(rid_list))
      
      idx <- lapply(rid_list[i],FUN = function(rid) overlap_indices[overlap_indices$yid %in% rid]$xid)
      impute <- lapply(idx,function(i) {
        imputed <- op(get_matrix(scm[i,],assay)) #TODO: make this parallel
        if (!setequal(dim(imputed),c(length(idx),ncol(scm)))) stop("Error with imputation algorithm. Imputed matrix does not match the dimensions of the input matrix", call. = FALSE)
        return(imputed)
      })
      
      for (i in 1:length(idx)) assays(scm)[[new_assay]][idx[[i]],] <- impute[[i]]
    }
  } else {

    imputed <- op(as.matrix(get_matrix(scm,assay)))
    if (!setequal(dim(imputed),c(nrow(scm),ncol(scm)))) stop("Error with imputation algorithm. Imputed matrix does not match the dimensions of the input matrix", call. = FALSE)
    assays(scm)[[new_assay]] <- imputed
    
  }
  
  if (any(is.na(assay(scm,new_assay)))) warning("NAs still present in the new_assay. This should not happen.")
   
  if (verbose) message("Imputed in ",stop_time())
  
  return(scm)
}

#------------------------------------------------------------------------------------------------------------
#' Splits an scMethrix object into two for use as a training and test set
#' @details Typically used for teaching classification algorithms. The seed can be set for consistency.
#' @param training_prop numeric; The size of the training set as a proportion of the experiment (0 to 1)
#' For a range, the optimal value will be estimated; this is time-intensive.
#' @param seed string; value to use for sampling
#' @inheritParams generic_scMethrix_function
#' @return list; two \code{\link{scMethrix}} objects names 'training' and 'test'
#' @examples
#' data('scMethrix_data')
#' generate_training_set(scMethrix_data, training_prop = 0.2)
#' @export
generate_training_set <- function(scm = NULL, training_prop = 0.2, seed = "123") {
 
  #- Input Validation --------------------------------------------------------------------------
   .validateExp(scm)
  .validateType(training_prop,"numeric")
  .validateType(seed,"string")
  
  if (training_prop > 1 || training_prop < 0) stop("training_prop must in the range of [0,1]", call. = FALSE)
  
  #- Function code -----------------------------------------------------------------------------
  set.seed(seed)
  idx <- sort(sample(1:nrow(scm),floor(nrow(scm)*training_prop)))
  
  training <- scm[idx,]
  test <- scm[setdiff(1:nrow(scm),idx),]
  
  return(list(training = training,test = test))
}

#------------------------------------------------------------------------------------------------------------
#' Generates a random subset of CpG sites
#' @details From an \code{\link{scMethrix}} object, this will randomly select \code{n_cpgs} and create a new
#' object containing only those CpGs. This is typically used for approximation or visualization. The seed
#' can be specified for consistency. 
#' @param n_cpgs numeric; The number of CpGs to include
#' @param seed string; value to use for sampling
#' @inheritParams generic_scMethrix_function
#' @return scMethrix; an experiment with n_cpgs
#' @examples
#' data('scMethrix_data')
#' generate_random_subset(scMethrix_data,n_cpgs = round(nrow(scMethrix_data)/2))
#' @export
generate_random_subset <- function(scm = NULL, n_cpgs = 10000, seed = "123") {
  
  #- Input Validation --------------------------------------------------------------------------
  .validateExp(scm)
  .validateType(n_cpgs,"integer")
  .validateType(seed,"string")
  
  if (n_cpgs > nrow(scm) || n_cpgs < 1) {
    n_cpgs = max(1,n_cpgs)
    n_cpgs = min(nrow(scm),n_cpgs)
    warning("Invalid n_cpgs. Must be between 1 and ",nrow(scm),". Defaulted to ",n_cpgs)
  }
  
  #- Function code -----------------------------------------------------------------------------
  set.seed(seed)
  idx <- sort(sample(1:nrow(scm),n_cpgs))
  
  return(scm[idx,])
}

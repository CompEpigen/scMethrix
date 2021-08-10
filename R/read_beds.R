#' Versatile BedGraph reader.
#' @details Reads BED files and generates methylation matrices.
#' Optionally arrays can be serialized as on-disk HDFS5 arrays.
#' @param files BED files containing methylation d
#' @param ref_cpgs BED files containing list of CpG sites. Must be zero-based genome.
#' @param colData vector of input sample names
#' @param stranded Whether in input data is stranded. Default FALSE
#' @param genome_name Name of genome. Default hg19
#' @param batch_size integer; Max number of files to hold in memory at once. Default 20
#' @param n_threads number of threads to use. Default 1.
#' Be-careful - there is a linear increase in memory usage with number of threads. This option is does not work with Windows OS.
#' @param h5 Should the coverage and methylation matrices be stored as \code{\link{HDF5Array}}
#' @param h5_dir directory to store H5 based object
#' @param h5_temp temporary directory to store hdf5
#' @param desc Description of the experiment
#' @param verbose flag to output messages or not.
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param reads Manual input of reads. Typically used for testing.
#' @param replace Boolean flag for whether to delete the contents of h5_dir before saving
#' @param pipeline Default NULL. Currently supports "Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"
#' If not known use idx arguments for manual column assignments.
#' @param chr_idx integer; column index for chromosome in bedgraph files
#' @param start_idx integer; column index for start position in bedgraph files
#' @param end_idx integer; column index for end position in bedgraph files
#' @param beta_idx integer; column index for beta values in bedgraph files
#' @param M_idx integer; column index for read counts supporting Methylation in bedgraph files
#' @param U_idx integer; column index for read counts supporting Un-methylation in bedgraph files
#' @param strand_idx integer; column index for strand information in bedgraph files
#' @param cov_idx integer; column index for total-coverage in bedgraph files
#' @export
#' @return An object of class \code{\link{scMethrix}}
#' @rawNamespace import(data.table, except = c(shift, first, second))
#' @import SingleCellExperiment GenomicRanges tools
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
# Must generate an index CpG file first:
#   sort-bed [input files] | bedops --chop 1 --ec - > CpG_index

read_beds <- function(files = NULL, ref_cpgs = NULL, colData = NULL, genome_name = "hg19", batch_size = 20, n_threads = 0, 
                      h5 = FALSE, h5_dir = NULL, h5_temp = NULL, desc = NULL, verbose = TRUE,
                      zero_based = FALSE, reads = NULL, replace = FALSE, pipeline = NULL, stranded = FALSE,
                      chr_idx = NULL, start_idx = NULL, end_idx = NULL, beta_idx = NULL,
                      M_idx = NULL, U_idx = NULL, strand_idx = NULL, cov_idx = NULL) {

  if (verbose) message("File import starting...")
  
  if (is.null(files)) {
    stop("Missing input files [files]", call. = FALSE)
  }
  
  if (!all(grepl("\\.(bed|bedgraph)", files))) {
    stop("Input files must be of type .bed or .bedgraph", call. = FALSE)
  }
  
  # if (is.null(ref_cpgs) && h5) {
  #   stop("Reference CpGs must be provided for HDF5 format", call. = FALSE)
  # }
  
  if (n_threads > length(files)/2){
    n_threads <- min(n_threads,length(files)/2) # cannot have multiple threads with a single file being input
    warning("Too many threads specified. Each thread must have at least 2 files to process. 
            Defaulting to n_thread = ", n_threads)
  } else if (n_threads < 0) {
    n_threads <- 0
    warning("n_threads < 0. Defaulting to 0")
  }
  
  # Get the correct indexes of the input beds
  if (is.null(pipeline)) {
    if (verbose) message(paste0("BED column format:  Custom"))
    col_list <- parse_source_idx(chr_idx = chr_idx, start_idx = start_idx, end_idx = end_idx,
                                 beta_idx = beta_idx, cov_idx = cov_idx, strand_idx = strand_idx, 
                                 M_idx = M_idx, U_idx = U_idx, verbose = verbose)
  } else {
    pipeline <- match.arg(arg = pipeline, choices = c("Bismark_cov", "MethylDackel", "MethylcTools", "BisSNP", "BSseeker2_CGmap"))
    if (verbose) message(paste0("BED column format:  ", pipeline))
    
    col_list <- get_source_idx(protocol = pipeline)
    
    if (any(pipeline %in% c("Bismark_cov"))) {
      if (zero_based) {
        warning("*BismarkCov files are one based. You may want to re-run with zero_based=FALSE")
      }
    }
  }

  #Check the reference cpgs
  if (is.null(ref_cpgs)) {
    #Generate references from input files
    ref_cpgs <- read_index(files = files, col_list = col_list, n_threads = n_threads,
                           batch_size = batch_size, zero_based = zero_based)
  } else {
    # Check stranding of reference genome and add - strand if nec.
    if (stranded) {
      ref_plus <- data.table::copy(ref_cpgs)
      ref_plus[, `:=`(strand, "+")]
      ref_cpgs[, `:=`(start, start + 1)]
      ref_cpgs[, `:=`(strand, "-")]
      ref_cpgs <- data.table::rbindlist(list(ref_cpgs, ref_plus), use.names = TRUE)
      data.table::setkeyv(ref_cpgs, cols = c("chr", "start"))
      message(paste0("-CpGs stranded: ", format(nrow(ref_cpgs), big.mark = ",")), "(reference CpGs from both strands)")
      rm(ref_plus) 
      gc()
    }
  }
  
  col_list$max_value = max(read_bed_by_index(files = files[1], ref_cpgs,col_list=col_list)$beta,na.rm = TRUE)
  
  gc()
  
  if (h5) {
    
    if (is.null(h5_dir)) stop("Output directory must be specified", call. = FALSE)
    
    if(dir.exists(h5_dir) && !replace) stop("h5_dir already exists! Use 'replace=TRUE' to replace it. All 
                                            existing data in that directory will be deleted.") 
    
    #if (zero_based) {ref_cpgs[,2:3] <- ref_cpgs[,2:3]+1}
    if (is.null(reads)) reads <- read_hdf5_data(files = files, ref_cpgs = ref_cpgs, col_list = col_list, 
                                                n_threads = n_threads, h5_temp = h5_temp, batch_size = batch_size,
                                                zero_based = zero_based, verbose = verbose)
    message("Building scMethrix object")

    if(is.null(colData)) {
      colData <- data.frame()[1:(length(files)), ]
      row.names(colData) <- colnames(reads$score)
    } else {
      cd <- data.frame()[1:(length(files)), ]
      cd$Sample <- colnames(reads$score)
      cd <- merge(cd,colData,by="Sample")
      row.names(cd) <- cd$Sample
      cd$Sample <- NULL
      colData <- cd
    }
    
    # if (is.null(colData)) colData <- data.frame()[1:(length(files)), ]
    # row.names(colData) <- unlist(lapply(files,get_sample_name))
    
    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    chrom_size = sapply(coverage(ref_cpgs), function(x) {length(x)-x@lengths[1]})
    
    m_obj <- create_scMethrix(assays = reads, rowRanges=ref_cpgs, is_hdf5 = TRUE, 
                              h5_dir = h5_dir, genome_name = genome_name,desc = desc,colData = colData,
                              replace = replace,chrom_size = chrom_size)
    
    message("Object built!\n")

  } else {

    message("Reading in BED files") 
    
    reads <- read_mem_data(files, ref_cpgs, col_list = col_list, batch_size = batch_size, n_threads = n_threads,
                           zero_based = zero_based,verbose = verbose)
    
    message("Building scMethrix object")

    if(is.null(colData)) {
      colData <- data.frame()[1:(length(files)), ]
      row.names(colData) <- colnames(reads$score)
    } else {
      cd <- data.frame()[1:(length(files)), ]
      cd$Sample <- colnames(reads$score)
      cd <- merge(cd,colData,by="Sample")
      row.names(cd) <- cd$Sample
      cd$Sample <- NULL
      colData <- cd
    }

    ref_cpgs <- GenomicRanges::makeGRangesFromDataFrame(ref_cpgs)
    chrom_size = sapply(coverage(ref_cpgs), function(x) {length(x)-x@lengths[1]})
    
    m_obj <- create_scMethrix(assays = reads, 
                              rowRanges=ref_cpgs, is_hdf5 = FALSE, genome_name = genome_name, 
                              desc = desc, colData = colData, chrom_size = chrom_size )
    
  }
  
  gc()
  message("Object built!\n")
  return(m_obj)
}

#' Parse BED files for unique genomic coordinates
#' @details Create list of unique genomic regions from input BED files. Populates a list of batch_size+1 with 
#' the genomic coordinates from BED files, then runs \code{\link{unique}} when the list is full and keeps the running
#' results in the batch_size+1 position. Also indexes based on 'chr' and 'start' for later searching.
#' @param files List of BED files
#' @param n_threads integer for number of threads to use
#' @param batch_size Number of files to process before running unique. Default of 30.
#' @param zero_based Whether the input data is 0 or 1 based
#' @param col_list The column index object for the input BED files
#' @param verbose flag to output messages or not.
#' @return data.table containing all unique genomic coordinates
#' @import parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
#' @export
read_index <- function(files, col_list, n_threads = 0, zero_based = FALSE, batch_size = 200, verbose = TRUE) {

  # Parallel functionality
  if (n_threads != 0) {
    
    if (n_threads > detectCores(logical = TRUE)) {
      n_threads <- detectCores(logical = TRUE)-1
      warning("Too many threads. Defaulting to n_threads =",n_threads)
    }
    
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_index','start_time','split_time','stop_time','get_sample_name'))
    
    if (verbose) message("Generating CpG index")
    
    chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/n_threads)))
    
    rrng <- c(parallel::parLapply(cl,chunk_files,fun=read_index, 
                                  batch_size=round(batch_size/n_threads), col_list = col_list,
                                  n_threads = 0, zero_based = zero_based, verbose = FALSE))
    
    parallel::stopCluster(cl)
    
    rrng <- data.table::rbindlist(rrng)
    data.table::setkeyv(rrng, c("chr","start"))
    rrng <- unique(rrng)
    
    return(rrng)
  }
  
  # Single thread functionality  
  # Determine all unique CpG sites from input files
  if (verbose) message("Generating index",start_time())
  
  rrng <- vector(mode = "list", length = batch_size+1)

  for (i in 1:length(files)) {
    if (verbose) message("   Parsing: ",get_sample_name(files[i]),appendLF=FALSE)
    
    data <- data.table::fread(files[i], select = unname(col_list$col_idx[c("chr","start")]))
    
    # Concat the batch list if last element, otherwise save and iterate
    if (i%%batch_size == 0) {
      
      rrng[[batch_size]] <- data
      data <- data.table::rbindlist(rrng)
      data <- unique(data) #data <- distinct(data)# #data <- data[!duplicated(data),]
      rrng <- vector(mode = "list", length = batch_size+1)
      rrng[[batch_size+1]] <- data
    } else {
      rrng[[i%%batch_size]] <- data
    }
    if (verbose) message(" (",split_time(),")")
  }
  
  #Concat the final list
  rrng <- data.table::rbindlist(rrng)
  colnames(rrng) <- c("chr","start")
  data.table::setkeyv(rrng, c("chr","start"))
  rrng <- unique(rrng)
  
  # Remove the consecutive subsequent sites
  #rrng <- rrng[data.table::rowid(collapse::seqid(rrng$start)) %% 2 == 1]
  
  #Add the end position and offset if zero based
  if (zero_based) rrng$start <- rrng$start+1
  rrng$end <- rrng$start+1
  rrng <- data.table::setkeyv(rrng, c("chr","start"))
  if (verbose) message("Index generated! (",stop_time(),")")#, Avg.", (stop_time()/length(files)),")")
  
  return(rrng)
}

#' Parses BED files for methylation values using previously generated index genomic coordinates
#' @details Creates an NA-based vector populated with methlylation values from the input BED file in the
#' respective indexed genomic coordinates
#' @param files The BED file to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param zero_based Whether the input data is 0 or 1 based
#' @param col_list The column index object for the input BED files
#' @param strand_collapse Default FALSe
#' @param fill Fill the output matrix to match the ref_cpgs. This must be used for HDF5 input formats to ensure consistent
#' spacing for the grid. It is optional for in-memory formats.
#' @return data.table containing vector of all indexed methylation values for the input BED
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
#' @export
# read_bed_by_index <- function(file, ref_cpgs, col_list = NULL, zero_based=FALSE) {
# 
#   . <- chr <- start <- beta <- meth1 <- meth2 <- cov <- cov1 <- cov2 <- NULL
# 
#   # Format will be: chr | start | meth | cov
#   data <- suppressWarnings(data.table::fread(file, select = unname(col_list$col_idx), 
#                                              col.names = names(col_list$col_idx), key = c("chr", "start")))
# 
#   data[, `:=`(chr, as.character(chr))]
#   data[, `:=`(start, as.integer(start))]
#   
#   if (!is.null(col_list$col_idx["beta"]) && !is.null(col_list$max_value) && col_list$max_value != 1) {
#     data[,beta := (beta/col_list$max_value)]
#   }
#   
#   if (!is.null(col_list$fix_missing)) {
#     for (cmd in col_list$fix_missing) {
#       data[, eval(parse(text = cmd))]
#     }
#   }
# 
#   #if (zero_based) data[,start] <- data[,start]+1L
# 
#   if (col_list$has_cov) {
#     
#     #Do the search
#     sample <- cbind(data[.(ref_cpgs$chr, ref_cpgs$start)][,.(beta,cov)],
#                     data[.(ref_cpgs$chr, ref_cpgs$end)][,.(beta,cov)])
#     colnames(sample) <- c("beta1", "cov1", "beta2","cov2")
#     
#     #Get the meth values from start and end CpG by weighted mean for both reads (meth*cov)
#     sample[,meth1 := meth1 * cov1]
#     sample[,meth2 := meth2 * cov2]
#     sample[,beta := rowSums(.SD, na.rm = TRUE), .SDcols = c("beta1", "beta2")]
#     sample[,cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("cov1", "cov2")]
#     sample[cov == 0, cov := NA] #since above line evals NA+NA as 0
#     sample[,beta := beta / cov]
#     sample <- sample[,.(beta,cov)]
#     
#     s <- nrow(ref_cpgs)-sum(is.na(sample[,(beta)]))
#   } else {
#     #Get the methylation value as the mean of start and end site
#     sample <- rowMeans(cbind(
#       data[.(ref_cpgs$chr, ref_cpgs$start)][,(beta)], 
#       data[.(ref_cpgs$chr, ref_cpgs$end)][,(beta)]), na.rm=TRUE)
#     s <- nrow(ref_cpgs)-sum(is.na(sample))
#   }
#   
#   d <- nrow(data)
#   
#   if (s < 0.9*d) 
#   {warning(paste0("Only ",round(s/d*100,1) ,"% of CpG sites in '",basename(file),"' are present in ref_cpgs"))}
#   
#   # Save the sample name data
#   if (col_list$has_cov) {
#     sample <- data.table(sapply(sample, as.integer))
#     sample <- list(beta = sample[,.(beta)], cov = sample[,.(cov)])
#     names(sample$beta) <- get_sample_name(file)
#     names(sample$cov) <- get_sample_name(file)
#   } else {
#     sample <- data.table(as.integer(sample))
#     names(sample) <- get_sample_name(file)
#     sample <- list(beta = sample)
#   }
# 
#   return(sample)
# }

# read_bed_by_index2 <- function(files, ref_cpgs, col_list = NULL, zero_based=FALSE) {
#   
#   . <- chr <- start <- beta <- meth1 <- meth2 <- cov <- cov1 <- cov2 <- NULL
#   
#   reads <- lapply(files, function (file) {
#   
#     # Format will be: chr | start | meth | cov
#     data <- suppressWarnings(data.table::fread(file, select = unname(col_list$col_idx), 
#                                                col.names = names(col_list$col_idx), key = c("chr", "start")))
#     
#     if (!is.null(col_list$col_idx["beta"]) && !is.null(col_list$max_value) && col_list$max_value != 1) {
#       data[,beta := (beta/col_list$max_value)]
#     }
#     
#     if (!is.null(col_list$fix_missing)) {
#       for (cmd in col_list$fix_missing) {
#         data[, eval(parse(text = cmd))]
#       }
#     }
#     
#     #if (zero_based) data[,start] <- data[,start]+1L
#     
#     if (col_list$has_cov) {
#       
#       #Do the search
#       sample <- cbind(data[.(ref_cpgs$chr, ref_cpgs$start)][,.(beta,cov)],
#                       data[.(ref_cpgs$chr, ref_cpgs$end)][,.(beta,cov)])
#       colnames(sample) <- c("beta1", "cov1", "beta2","cov2")
#       
#       #Get the meth values from start and end CpG by weighted mean for both reads (meth*cov)
#       sample[,beta1 := beta1 * cov1]
#       sample[,beta2 := beta2 * cov2]
#       sample[,beta := rowSums(.SD, na.rm = TRUE), .SDcols = c("beta1", "beta2")]
#       sample[,cov := rowSums(.SD, na.rm = TRUE), .SDcols = c("cov1", "cov2")]
#       sample[cov == 0, cov := NA] #since above line evals NA+NA as 0
#       sample[,beta := beta / cov]
#       sample <- sample[,.(beta,cov)]
#       
#       s <- nrow(ref_cpgs)-sum(is.na(sample[,(beta)]))
#     } else {
#       #Get the methylation value as the mean of start and end site
#       sample <- rowMeans(cbind(
#         data[.(ref_cpgs$chr, ref_cpgs$start)][,(beta)], 
#         data[.(ref_cpgs$chr, ref_cpgs$end)][,(beta)]), na.rm=TRUE)
#       s <- nrow(ref_cpgs)-sum(is.na(sample))
#     }
#     
#     d <- nrow(data)
#     
#     if (s < 0.9*d) 
#     {warning(paste0("Only ",round(s/d*100,1) ,"% of CpG sites in '",basename(file),"' are present in ref_cpgs"))}
#     
#     # Save the sample name data
#     if (col_list$has_cov) {
#       sample <- data.table(sapply(sample, as.integer))
#       sample <- list(beta = sample[,.(beta)], cov = sample[,.(cov)])
#       names(sample$beta) <- get_sample_name(file)
#       names(sample$cov) <- get_sample_name(file)
#     } else {
#       sample <- data.table(as.integer(sample))
#       names(sample) <- get_sample_name(file)
#       sample <- list(beta = sample)
#     }
#     
#     return(sample)
#   })
#   
#   browser()
#   
#   if (length(reads) == 1) {
#     return (reads[[1]])
#   } else {
#     
#     return(list(beta = as.matrix(do.call(cbind, lapply(reads, `[[`, "beta"))),
#                 cov = as.matrix(do.call(cbind, lapply(reads, `[[`, "cov")))))
#   }
# }
# 
read_bed_by_index <- function(files, ref_cpgs, col_list = NULL, zero_based=FALSE, strand_collapse = FALSE, fill = TRUE) {

  . <- meths <- covs <- M <- U <- chr <- NULL

  for (i in 1:length(files)) {

    bed <- data.table::fread(files[[i]], select = unname(col_list$col_idx),
                             col.names = names(col_list$col_idx), key = c("chr", "start"))
    n_bed <- nrow(bed)
    
    #Scale beta to [0,1]
    if (!is.null(col_list$col_idx["beta"])) {
      if(!is.null(col_list$max_value) && col_list$max_value != 1) {
        bed[,beta := (beta/col_list$max_value)]
      } else {
        sample_row_idx = sample(x = seq_len(nrow(bed)), size = min(1000,nrow(bed)), replace = FALSE)
        max_beta = max(bed[sample_row_idx, beta], na.rm = TRUE)
        rm(sample_row_idx)
        bed[,beta := (beta/max_beta)]
      }
    }

    #Fill in other columns
    if (!is.null(col_list$fix_missing)) {
      for (cmd in col_list$fix_missing) {
        bed[, eval(parse(text = cmd))]
      }
    }
    
    bed[, c("end","strand"):=NULL] # TODO: this probably breaks stuff
    
    # Bring bedgraphs to 1-based cordinate
    if (zero_based) {
      bed[, `:=`(start, start + 1)]
    }
    
    # Merge with the ref_cpgs
    bed <- merge(ref_cpgs,bed)
    
    # If strand information needs to collapsed, bring start position of
    # crick strand to previous base (on watson base) and estimate new M, U
    # and beta values
    if (strand_collapse) {
      if (!all(c("M", "U") %in% names(bed))) stop("strand_collapse works only when M and U are available!")
      bed[, `:=`(start, ifelse(strand == "-", yes = start - 1, no = start))]
      bed = bed[, .(M = sum(M, na.rm = TRUE), U = sum(U, na.rm = TRUE)), .(chr, start)]
      bed[, `:=`(cov, M + U)]
      #bed[cov == 0, cov := NA]
      bed[, `:=`(beta, M/cov)]
      bed[, `:=`(strand, "+")]
    }

    meth <- bed[,c("chr","start","beta")]
    cov <- bed[,c("chr","start","cov")]

    setnames(meth, "beta", get_sample_name(files[[i]]))
    setnames(cov, "cov", get_sample_name(files[[i]]))

    if (is.null(meths) && is.null(covs)) {
      meths <- meth
      covs <- cov
    } else {
      meths <- merge(meths,meth,all=TRUE)
      covs <- merge(covs,cov,all=TRUE)
    }
  }
  
  if (fill) {
    suppressWarnings(meths <- merge(ref_cpgs,meths,all.x=TRUE)[,-c("end","strand")])
    if (col_list$has_cov) suppressWarnings(covs <- merge(ref_cpgs,covs,all.x=TRUE)[,-c("end","strand")])
  }
  
  if (col_list$has_cov) {
    return (list(beta = meths[,-c("chr","start")], cov = covs[,-c("chr","start")]))
  } else {
    return (list(beta = meths[,-c("chr","start")]))
  }
}

#' Writes values from input BED files into an in-disk \code{\link{HDF5Array}}
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param h5_temp The file location to store the \code{\link{RealizationSink}} object
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param col_list The column index object for the input BED files
#' @param verbose flag to output messages or not.
#' @param batch_size The number of file to hold in memory at once
#' @return List of \code{\link{HDF5Array}}. 1 is methylation, 2 is coverage. If no cov_idx is specified, 2 will be NULL
#' @import DelayedArray HDF5Array parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
read_hdf5_data <- function(files, ref_cpgs, col_list, batch_size = 20, n_threads = 0, h5_temp = NULL, zero_based = FALSE, 
                           verbose = TRUE) {
  
  if (verbose) message("Starting HDF5 object",start_time()) 
  
  if (is.null(h5_temp)) {h5_temp <- tempdir()}
  
  # Generate the realization sinks
  dimension <- as.integer(nrow(ref_cpgs))
  colData <- as.vector(unlist(lapply(files,get_sample_name)))
  
  M_sink <- HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                           dimnames = list(NULL,colData), type = "double",
                                           filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp),
                                           name = "M", level = 6)
  
  cov_sink <- if(!col_list$has_cov) NULL else 
    HDF5Array::HDF5RealizationSink(dim = c(dimension, length(files)),
                                   dimnames = list(NULL, colData), type = "integer",
                                   filepath = tempfile(pattern = "cov_sink_", tmpdir = h5_temp),
                                   name = "C", level = 6)
 
  # Determine the grids for the sinks
  if (n_threads == 0) {
    files <- split_vector(files,batch_size,by="size")
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(unlist(files))),
                                           spacings = c(dimension, length(files[[1]]))) 
  } else {
    files <- split_vector(files,ceiling(batch_size/n_threads),by="size")
    grid <- DelayedArray::RegularArrayGrid(refdim = c(dimension, length(unlist(files))),
                                           spacings = c(dimension, length(files[[1]]))) 
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_bed_by_index','get_sample_name'))
  }
  
  # Read data to the sinks
  for (i in 1:length(files)) {
    
    if (verbose) message("   Parsing: Chunk ", i,appendLF=FALSE)
    
    if (n_threads == 0) {
      
      bed <- read_bed_by_index(files = files[[i]], ref_cpgs = ref_cpgs, col_list = col_list, zero_based = zero_based)
      
      DelayedArray::write_block(block = as.matrix(bed[["beta"]]), viewport = grid[[i]], sink = M_sink)
      
      if (col_list$has_cov) DelayedArray::write_block(block = as.matrix(bed[["cov"]]),
                                                      viewport = grid[[i]], sink = cov_sink)
    } else {

      bed <- parallel::parLapply(cl,split_vector(files[[i]],n_threads,by="chunks"),
                                 fun=read_bed_by_index, ref_cpgs = ref_cpgs,
                                 col_list = col_list, zero_based = zero_based)
      
      DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "beta"))), 
                                viewport = grid[[i]], sink = M_sink)      
      
      if (!is.null(col_list$cov)) {
        DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "cov"))), 
                                  viewport = grid[[i]], sink = cov_sink)
      }
    }
    
    rm(bed)
    if (i%%5==0) gc()
    if (verbose) message(" (",split_time(),")")
  }
  
  if (n_threads != 0) parallel::stopCluster(cl)
  if (verbose) message("Object created in ",stop_time()) 
  
  if (col_list$has_cov) {
    reads = list(score = M_sink, counts = cov_sink)
  } else {
    reads = list(score = M_sink)
  }
  
  reads <- lapply(reads,function(x) as(x, "HDF5Array"))
  
  return(reads)
}

# read_hdf5_data2 <- function(files, ref_cpgs, col_list, n_threads = 0, n_chunks = 1, h5_temp = NULL, zero_based = FALSE, 
#                            verbose = TRUE) {
#   
#   if (verbose) message("Starting HDF5 object",start_time()) 
#   
#   if (is.null(h5_temp)) {h5_temp <- tempdir()}
#   
#   # Generate the realization sinks
#   n_cpg <- as.integer(nrow(ref_cpgs))
#   colData <- as.vector(unlist(lapply(files,get_sample_name)))
#   
#   M_sink <- HDF5Array::HDF5RealizationSink(dim = c(n_cpg, length(files)),
#                                            dimnames = list(NULL,colData), type = "integer",
#                                            filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp),
#                                            name = "M", level = 6)
#   
#   cov_sink <- if(!col_list$has_cov) NULL else 
#               HDF5Array::HDF5RealizationSink(dim = c(n_cpg, length(files)),
#                                              dimnames = list(NULL,colData), type = "integer",
#                                              filepath = tempfile(pattern = "cov_sink_", tmpdir = h5_temp),
#                                              name = "C", level = 6)
# 
#   files <- split_vector(files,n_chunks,by="chunks")
#   
#   grid <- DelayedArray::RegularArrayGrid(refdim = c(n_cpg, length(unlist(files))),
#                                            spacings = c(n_cpg, length(files[[1]]))) 
#   cl <- parallel::makeCluster(n_threads)  
#   doParallel::registerDoParallel(cl)  
#   parallel::clusterEvalQ(cl, c(library(data.table)))
#   parallel::clusterExport(cl,list('read_bed_by_index2','get_sample_name'))
#   on.exit(parallel::stopCluster(cl))
#   
#   # Read data to the sinks
#   for (i in 1:length(files)) {
#     
#     i <- as.integer(i) #TODO: figure out why this is necessary or grid[[i]] throws
#                        #each subscript must be a single integer when subsetting an RegularArrayGrid object with [['
#     
#     if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE)
# 
#     bed <- parallel::parLapply(cl,split_vector(files[[i]],n_threads,by="chunks"),fun=read_bed_by_index2, 
#                                  ref_cpgs = ref_cpgs, col_list = col_list, zero_based = zero_based)
# 
#     DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "beta"))), 
#                                 viewport = grid[[i]], sink = M_sink)      
#       
#     if (!is.null(col_list$cov)) {
#       DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "cov"))), 
#                                   viewport = grid[[i]], sink = cov_sink)
#     }
#     
#     
#     rm(bed)
#     gc()
#     if (verbose) message(" (",split_time(),")")
#   }
#   
#   if (verbose) message("Object created in ",stop_time()) 
#   
#   if (col_list$has_cov) {
#     reads = list(score = M_sink, counts = cov_sink)
#   } else {
#     reads = list(score = M_sink)
#   }
#   
#   reads <- lapply(reads,function(x) as(x, "HDF5Array"))
#   
#   return(reads)
# }
# 
# read_hdf5_data3 <- function(files, ref_cpgs, col_list, n_threads = 0, n_chunks = 1, h5_temp = NULL, zero_based = FALSE, 
#                             verbose = TRUE) {
#   
#   if (verbose) message("Starting HDF5 object",start_time()) 
#   
#   if (is.null(h5_temp)) {h5_temp <- tempdir()}
#   
#   # Generate the realization sinks
#   n_cpg <- as.integer(nrow(ref_cpgs))
#   colData <- as.vector(unlist(lapply(files,get_sample_name)))
#   
#   M_sink <- HDF5Array::HDF5RealizationSink(dim = c(n_cpg, length(files)),
#                                            dimnames = list(NULL,colData), type = "integer",
#                                            filepath = tempfile(pattern="M_sink_",tmpdir=h5_temp),
#                                            name = "M", level = 6)
#   
#   cov_sink <- if(!col_list$has_cov) NULL else 
#     HDF5Array::HDF5RealizationSink(dim = c(n_cpg, length(files)),
#                                    dimnames = list(NULL,colData), type = "integer",
#                                    filepath = tempfile(pattern = "cov_sink_", tmpdir = h5_temp),
#                                    name = "C", level = 6)
#   
#   files <- split_vector(files,n_chunks,by="chunks")
#   
#   grid <- DelayedArray::RegularArrayGrid(refdim = c(n_cpg, length(unlist(files))),
#                                          spacings = c(n_cpg, length(files[[1]]))) 
#   cl <- parallel::makeCluster(n_threads)  
#   doParallel::registerDoParallel(cl)  
#   parallel::clusterEvalQ(cl, c(library(data.table)))
#   parallel::clusterExport(cl,list('read_bed_by_index2','get_sample_name'))
#   on.exit(parallel::stopCluster(cl))
#   
#   # Read data to the sinks
#   for (i in 1:length(files)) {
#     
#     i <- as.integer(i) #TODO: figure out why this is necessary or grid[[i]] throws
#     #each subscript must be a single integer when subsetting an RegularArrayGrid object with [['
#     
#     if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE)
#     
#     bed <- parallel::parLapply(cl,split_vector(files[[i]],n_threads,by="chunks"),fun=read_bed_by_index2, 
#                                ref_cpgs = ref_cpgs, col_list = col_list, zero_based = zero_based)
#     
#     DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "beta"))), 
#                               viewport = grid[[i]], sink = M_sink)      
#     
#     if (!is.null(col_list$cov)) {
#       DelayedArray::write_block(block = as.matrix(do.call(cbind, lapply(bed, `[[`, "cov"))), 
#                                 viewport = grid[[i]], sink = cov_sink)
#     }
#     
#     
#     rm(bed)
#     gc()
#     if (verbose) message(" (",split_time(),")")
#   }
#   
#   if (verbose) message("Object created in ",stop_time()) 
#   
#   if (col_list$has_cov) {
#     reads = list(score = M_sink, counts = cov_sink)
#   } else {
#     reads = list(score = M_sink)
#   }
#   
#   reads <- lapply(reads,function(x) as(x, "HDF5Array"))
#   
#   return(reads)
# }

#' Writes values from input BED files into an in-memory \code{\link{matrix}}
#' @details Using the generated index for genomic coordinates, creates a NA-based dense matrtix of methylation
#' values for each BED file/sample. Each column contains the meth. values for a single sample.
#' @param files The BED files to parse
#' @param ref_cpgs The index of all unique coordinates from the input BED files
#' @param batch_size The number of files to hold in memory at once
#' @param n_threads The number of threads to use. 0 is the default thread with no cluster built.
#' @param zero_based Boolean flag for whether the input data is zero-based or not
#' @param col_list The column index object for the input BED files
#' @param verbose flag to output messages or not.
#' @return matrix of the methylation values for input BED files
#' @import parallel doParallel
#' @examples
#' \dontrun{
#' #Do Nothing
#' }
read_mem_data <- function(files, ref_cpgs, col_list, batch_size = 20, n_threads = 0, zero_based = FALSE, 
                          verbose = TRUE) {

  if (verbose) message("Reading BED data...",start_time()) 

  if (n_threads != 0) {
    # Parallel functionality
    if (verbose) message("Starting cluster with ",n_threads," threads.")
    
    cl <- parallel::makeCluster(n_threads)  
    doParallel::registerDoParallel(cl)  
    
    parallel::clusterEvalQ(cl, c(library(data.table)))
    parallel::clusterExport(cl,list('read_bed_by_index','start_time','split_time','stop_time','get_sample_name'))
    
    chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/n_threads)))
    
    reads <- c(parallel::parLapply(cl,files,fun=read_bed_by_index, ref_cpgs = ref_cpgs, 
                                   zero_based = zero_based, col_list = col_list))
    
    parallel::stopCluster(cl)
  
  } else {
    # Single thread functionality
    # if (verbose) message("   Parsing: Chunk ",i,appendLF=FALSE) #TODO: Get this workings
    
    files <- split_vector(files,batch_size,by="size")
    reads <- lapply(files,read_bed_by_index,ref_cpgs = ref_cpgs,zero_based = zero_based, col_list = col_list)
  }

  if (col_list$has_cov) {
    reads <- list(score = cbind(lapply(reads, `[[`, "beta")),
                  counts = cbind(lapply(reads, `[[`, "cov")))
  } else {
    reads <- list(score = cbind(lapply(reads, `[[`, "beta")))
  }
    
  reads <- lapply(reads,function(x) as.matrix(do.call(cbind, x)))
    
    # if (verbose) message(" (",split_time(),")")
  
  if (verbose) message("Data read in ",stop_time()) 
  
  return (reads)

}
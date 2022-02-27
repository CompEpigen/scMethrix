combineSE <- function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
  
  args <- unname(list(x, y, ...))
  
  # Check that each object has unique colnames
  colnames <- unlist(lapply(args, colnames))
  if (anyDuplicated(colnames)) {
    stop(paste0("Cannot combine ", class(args[[1]]), " objects with ", 
                "duplicate 'colnames'"))
  }
  
  # Check that each object has the same assays
  an <- lapply(args, assayNames)
  if (any(sapply(an, function(x, y) any(is.na(match(x, y))), 
                 y = an[[1]]))) {
    stop(paste0("All ", class(args[[1]]), " objects must have ", 
                "identical 'assayNames'"))
  }
  
  # Combine rowRanges or NAMES slot
  if (is(args[[1L]], "RangedSummarizedExperiment")) {
    # NOTE: rowRanges(x) includes elementMetadata(x) since 
    #       elementMetadata(x) = elementMetadata(rowRanges(x)).
    #       This is the reason for the use.mcols if-else block.
    if (use.mcols) {
      all_rowRanges <- do.call(c, lapply(args, rowRanges))
    } else {
      all_rowRanges <- do.call(c, lapply(args, function(x) {
        rr <- rowRanges(x)
        elementMetadata(rr) <- NULL
        rr
      }))
    }
    rowRanges <- unique(all_rowRanges)
    nr <- length(rowRanges)
  } else {
    if (any(vapply(args, function(x) is.null(x@NAMES), logical(1)))) {
      stop("Cannot combine ", class(args[[1]]), " objects with ", 
           "'NAMES' set to NULL")
    }
    all_NAMES <- do.call(c, lapply(args, function(x) x@NAMES))
    NAMES <- unique(all_NAMES)
    nr <- length(NAMES)
  }
  
  # Combine colData
  colData <- do.call(rbind, lapply(args, colData))
  
  # Combine elementMetadata
  if (use.mcols) {
    if (is(args[[1L]], "RangedSummarizedExperiment")) {
      elementMetadata <- mcols(rowRanges)
    } else {
      # IDEA: Create DataFrame with all_NAMES/all_rowRanges in one 
      #       column and elementMetadata in others, unique-ify, and 
      #       check that the number of unique rows equals nr.
      # WARNING: This will be slow for large SummarizedExperiment0 
      #          objects
      elementMetadata <- unique(
        cbind(DataFrame(all_NAMES), 
              do.call(rbind, lapply(args, elementMetadata))))
      # Now drop the all_NAMES column (assumes it is the first column)
      elementMetadata <- elementMetadata[, -c(1L), drop = FALSE]
      # Sanity check
      if (nrow(elementMetadata) != nr) {
        stop(paste0("'elementMetadata' must match across ", 
                    class(args[[1]]), " objects"))
      }
    }
  } else {
    elementMetadata <- DataFrame()
    elementMetadata@nrows <- nr
  }
  
  # Create assays of the correct dimension (fill with 'nomatch')
  # First, create the empty combined assay using the appropriate 
  # storage.mode (guessed from the storage.mode of the assay in the 
  # first sample).
  nomatch <- lapply(seq_along(an[[1]]), function(i) {
    storage.mode(nomatch) <- storage.mode(assay(args[[1]], 
                                                i, 
                                                withDimnames = FALSE))
    nomatch
  })
  assays <- lapply(nomatch, function(nm) {
    matrix(nm, nrow = nr, ncol = length(colnames))
  })
  names(assays) <- an[[1]]
  
  # NOTE: I suspect that there are faster and more efficient ways to 
  # combine the assays, perhaps at the C-level.
  if (is(args[[1L]], "RangedSummarizedExperiment")) {
    for (j in seq_along(args)) {
      ol <- findOverlaps(args[[j]], rowRanges, type = "equal")
      for (i in seq_along(assays)) {
        assays[[i]][subjectHits(ol),
                    match(colnames(args[[j]]), colnames)] <- 
          assay(args[[j]], i, withDimnames = FALSE)
      }
    }
  } else {
    for (j in seq_along(args)) {
      ol <- match(args[[j]]@NAMES, NAMES)
      for (i in seq_along(assays)) {
        assays[[i]][ol, match(colnames(args[[j]]), colnames)] <- 
          assay(args[[j]], i, withDimnames = FALSE)
      }
    }
  }
  assays <- Assays(assays)
  
  # Combine metadata
  metadata <- do.call(c, lapply(args, metadata))
  
  if (is(args[[1L]], "RangedSummarizedExperiment")) {
    # No need to replace elementMetadata slot since it is part of 
    # rowRanges.
    BiocGenerics:::replaceSlots(args[[1L]], 
                                rowRanges = rowRanges,
                                colData = colData, 
                                assays = assays,
                                metadata = metadata)
  } else {
    BiocGenerics:::replaceSlots(args[[1L]],
                                NAMES = NAMES,
                                colData = colData,
                                assays = assays,
                                metadata = metadata,
                                elementMetadata = elementMetadata)
  }
}
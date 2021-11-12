message("Test setup starting...")

# files <- c(system.file("extdata", "C1.bedgraph", package="scMethrix"),
#            system.file("extdata", "C2.bedgraph", package="scMethrix"),
#            system.file("extdata", "C3.bedgraph", package="scMethrix"),
#            system.file("extdata", "C4.bedgraph", package="scMethrix"))

files <- c("F:/scMethrix/inst/extdata/C1.bedgraph","F:/scMethrix/inst/extdata/C2.bedgraph",
           "F:/scMethrix/inst/extdata/C3.bedgraph","F:/scMethrix/inst/extdata/C4.bedgraph")

h5_dir <- paste0(tempdir(),"/sse")

colData <- data.frame(row.names = get_sample_name(files), Group = rep(1:2,2))

col_list <- parse_source_idx(chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)
scm.h5 <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
scm.mem <- read_beds(files,h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)

#scm.mem <- read_beds(files,h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)

mcols(scm.h5)$CpG <- mcols(scm.mem)$CpG <- 1:nrow(scm.mem)

n_cpg <- nrow(scm.mem)
n_samples <- ncol(scm.mem)

# These should match the error messages given in the respective functions in accessory_funcs.R
msg.validateExp <- "Invalid scMethrix"
msg.validateAssay <- "Invalid assay"
msg.validateArg <- "Invalid arg"
msg.validateType <- "Invalid type"
msg.validateValue <- "Invalid value"

# scMethrix_data <- scm.mem
# usethis::use_data(scMethrix_data,overwrite=TRUE)

message("Test setup completed")

imputation_test_helper <- function(func) {
  expect_error(imputation_test_helper("not scMethrix"),msg.validateExp)
  
  #invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  invisible(lapply(list(scm.mem), function(scm) {
    expect_error(func(scm,assay = "not an assay"))
    expect_error(func(scm,new_assay = "score"))
    expect_warning(func(scm,new_assay = "counts"))
    
    impute = func(scm,new_assay="impute")
    expect_true("impute" %in% SummarizedExperiment::assayNames(impute))
    
    sco <- get_matrix(impute,assay="score")
    imp <- get_matrix(impute,assay="impute")
    NAs <- which(is.na(sco))
    nonNAs <- which(!is.na(sco))
    
    expect_true(anyNA(sco) && !anyNA(imp))
    expect_equivalent(sco[nonNAs],imp[nonNAs])
    expect_false(all(sco[NAs] %in% imp[NAs]))
  }))
}

graph_test_helper <- function(scm, func, indiv_samples = TRUE, indiv_chr = FALSE, pheno = NULL, ...) {

  frmat <- function(x) sort(as.character(unique(x)))
  
  expect_error(func("not scMethrix"),msg.validateExp)
  
  if (is.null(pheno)) {
    plot <- func(scm = scm, verbose = FALSE, ...)
  } else {
    plot <- func(scm = scm, pheno = pheno, verbose = FALSE, ...)
  }

  expect_true("ggplot" %in% class(plot))
    
  if (indiv_samples) {
    expect_equivalent(frmat(plot$data$Sample),frmat(colnames(scm)))
  } else {
    expect_false("Sample" %in% colnames(plot$data))
  }
  
  if (indiv_chr) {
    expect_equal(frmat(plot$data$Chromosome), frmat(seqnames(rowRanges(scm))))
  } else {
    expect_false("Chromosome" %in% colnames(plot$data))
  }

  if (!is.null(pheno)) {
    expect_equal(frmat(plot$data$Pheno), frmat(colData(scm)[,pheno]))
  } else {
    if ("Pheno" %in% colnames(plot$data)) expect_equal(frmat(plot$data$Pheno),frmat(plot$data$Sample))
  }
  
  expect_error(print(plot),NA) # Checks if the plot is printable
  
  return(invisible(plot))
}

dim_red_graph_test_helper <- function(scm, func, color_anno = NULL, shape_anno=NULL, ...) {

  frmat <- function(x) sort(as.character(unique(x)))

  plot <- graph_test_helper(scm = scm, func = func, color_anno = color_anno, shape_anno = shape_anno, ...)

  if (!is.null(color_anno)) {
    expect_equal(frmat(plot$data$Color), frmat(colData(scm)[,color_anno]))
  } else {
    expect_false("Color" %in% colnames(plot$data))
  }
  
  if (!is.null(shape_anno)) {
    expect_equal(frmat(plot$data$Shape), frmat(colData(scm)[,shape_anno]))
    
  } else {
    expect_false("Shape" %in% colnames(plot$data))
  }
  
  return(invisible(plot))
}

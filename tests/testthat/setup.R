message("Test setup starting...")

# files <- c(system.file("extdata", "C1.bedgraph", package="scMethrix"),
#            system.file("extdata", "C2.bedgraph", package="scMethrix"),
#            system.file("extdata", "C3.bedgraph", package="scMethrix"),
#            system.file("extdata", "C4.bedgraph", package="scMethrix"))

files <- c("D:/Git/scMethrix/inst/extdata/C1.bedgraph","D:/Git/scMethrix/inst/extdata/C2.bedgraph",
           "D:/Git/scMethrix/inst/extdata/C3.bedgraph","D:/Git/scMethrix/inst/extdata/C4.bedgraph")

h5_dir <- tempfile("scm_h5_", tmpdir = tempdir())
colData <- DataFrame(row.names = get_sample_name(files), Group = rep(1:2,2))
refGenome <- "hg19"

col_list <- parse_source_idx(chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4, cov_idx = 5)
scm.h5  <- read_beds(files, is_h5 = TRUE,h5_dir = h5_dir, replace = TRUE, genome = refGenome,
                     chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
scm.mem <- read_beds(files, is_h5=FALSE, genome = refGenome,                     
                     chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)

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

# Checks to see if two experiments are the same, excluding H5 status
identical.scm <- function(scm1,scm2,exclude_is_h5 = F) {

  if (exclude_is_h5) {
    metadata(scm1) <- within(metadata(scm1), rm(is_h5))
    metadata(scm2) <- within(metadata(scm2), rm(is_h5))
  } 
  
  match = 
    identical(as.matrix(score(scm1)),as.matrix(score(scm2))) &&
    identical(colData(scm1),colData(scm2)) && 
    identical(rowData(scm1),rowData(scm2)) &&
    identical(metadata(scm1),metadata(scm2))
  
  return(match)
}

graph_test_helper <- function(scm, func, plot, indiv_samples = TRUE, indiv_chr = FALSE, pheno = NULL, 
                              expected_samples = NULL, samples = NULL, groups = NULL, ...) {

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

graph_test_helper2 <- function(plot, expected_x = NULL, expected_y = NULL, expected_label = NULL) {
  
  frmat <- function(x) sort(as.character(unique(x)))
  
  expect_true("ggplot" %in% class(plot))
  
  x_labs <- frmat(ggplot_build(plot)$layout$panel_params[[1]]$x$get_labels())
  y_labs <- frmat(ggplot_build(plot)$layout$panel_params[[1]]$y$get_labels())
  grp_labs <- frmat(ggplot_build(plot)$data[[1]]$label)
  
  if (!is.null(expected_x)) expect_equivalent(frmat(x_labs),frmat(expected_x))
  if (!is.null(expected_y)) expect_equivalent(frmat(y_labs),frmat(expected_y))
  if (!is.null(expected_label)) expect_equivalent(frmat(grp_labs),frmat(expected_label))
  
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

message("Test setup starting...")

files <- c(system.file("extdata", "C1.bedgraph", package="scMethrix"),
           system.file("extdata", "C2.bedgraph", package="scMethrix"),
           system.file("extdata", "C3.bedgraph", package="scMethrix"),
           system.file("extdata", "C4.bedgraph", package="scMethrix"))

files <- c("D:/Git/scMethrix/inst/extdata/C1.bedgraph","D:/Git/scMethrix/inst/extdata/C2.bedgraph",
           "D:/Git/scMethrix/inst/extdata/C3.bedgraph","D:/Git/scMethrix/inst/extdata/C4.bedgraph")

h5_dir <- paste0(tempdir(),"/sse")

col_list <- parse_source_idx(chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)
scm.h5 <- read_beds(files,h5=TRUE,h5_dir=h5_dir,replace=TRUE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)
scm.mem <- read_beds(files,h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5)

n_cpg <- nrow(scm.mem)
n_samples <- ncol(scm.mem)

# scMethrix_data <- scm.mem
# usethis::use_data(scMethrix_data,overwrite=TRUE)

message("Test setup completed")

imputation_test_helper <- function(func) {
  expect_error(func("not scMethrix"),"A valid scMethrix object needs to be supplied")
  
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

graph_test_helper <- function(func, indiv_samples = TRUE) {
  expect_error(func("not scMethrix"),"A valid scMethrix object needs to be supplied")
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    plot <- func(scm)
    expect_true("ggplot" %in% class(plot))
    if (indiv_samples) expect_equivalent(levels(plot$data$variable),colnames(scm))
  }))
}
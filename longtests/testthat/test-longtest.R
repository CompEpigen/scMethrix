test_that("longtest", {

  path <- paste0(h5_dir,"/HDF5mem")
  unlink(path, recursive = TRUE)
  suppressWarnings(dir.create(path,recursive=TRUE))
  
  scm.h5 <- read_beds(files,is_h5=TRUE,h5_dir=path,replace=TRUE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
  expect_true(validObject(scm.h5))
  
  scm.mem <- read_beds(files,is_h5=FALSE,chr_idx=1, start_idx=2, end_idx=3, beta_idx=4, cov_idx=5, colData = colData)
  expect_true(validObject(scm.mem))
  
  expect_true(identical.scm(h5,mem,exclude_is_h5=TRUE))
  
  exps <- invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  
    scm <- convert_scMethrix(scm)
    
    # Get all the plots
    plot_spar_start <- plot_sparsity(scm)
    plot_dens_start <- plot_density(scm)
    plot_cov_start <- plot_coverage(scm)
    plot_stat_start <- plot_stats(scm)
    metadata(scm) <- append(metadata(scm), list(plot_spar_start = plot_spar_start, plot_dens_start = plot_dens_start, 
                                                plot_cov_start = plot_cov_start, plot_stat_start = plot_stat_start))
    
    # Get stats
    scm <- get_rowdata_stats(scm)
    scm <- get_coldata_stats(scm)
    
    # Do coord opertions
    scm <- transform_assay(scm, new_assay="binarize1", trans=binarize)
    scm <- impute_regions(scm, new_assay = "impute1")
    scm <- bin_scMethrix(scm, bin_size = 10000000)
    scm <- subset_scMethrix(scm, samples=head(sampleNames(scm),length(sampleNames(scm))-1))
    scm <- subset_scMethrix(scm, regions=rowRanges(scm)[1:5], by="exclude")
    scm <- impute_regions(scm, new_assay = "impute2")
    scm <- transform_assay(scm, assay="impute2", new_assay="binarize2", trans=binarize)

    # Do dim reds
    #scm <- dim_red_scMethrix(scm,assay="impute2",type="tSNE", perplexity=10)
    scm <- dim_red_scMethrix(scm,assay="impute2",type="UMAP")
    scm <- dim_red_scMethrix(scm,assay="impute2",type="PCA")
    
    # Do clustering
    scm <- cluster_scMethrix(scm, assay= "impute2", type="hier",colname="Heir",n_clusters=2)
    scm <- cluster_scMethrix(scm, assay= "impute2", type="part",colname="Part",n_clusters=2)
    scm <- cluster_scMethrix(scm, assay= "impute2", type="model",colname="Model")
    
    # Get stats
    scm <- get_rowdata_stats(scm, assay = "impute2", suffix = "_impute2")
    scm <- get_coldata_stats(scm, assay = "impute2", suffix = "_impute2")
    
    # do more plots
    plot_spar_end <- plot_sparsity(scm, assay="impute2", phenotype="Heir")
    plot_dens_end <- plot_density(scm, assay="impute2", phenotype="Heir")
    plot_cov_end <- plot_coverage(scm, assay="impute2", phenotype="Heir")
    plot_stat_end <- plot_stats(scm, assay="impute2", phenotype="Heir")
    metadata(scm) <- append(metadata(scm), list(plot_spar_end = plot_spar_end, plot_dens_end = plot_dens_end, 
                                                plot_cov_end = plot_cov_end, plot_stat_end = plot_stat_end))
    
    scm <- convert_scMethrix(scm)
    
    return (scm)
  }))
  
  expect_true(identical.scm(exp[[1]],exp[[2]],exclude_is_h5=TRUE))
  
})

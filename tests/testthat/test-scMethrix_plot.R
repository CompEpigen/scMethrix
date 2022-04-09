#---- .prepare_plot_data ------------------------------------------------------------------------------------------------
test_that(".prepare_plot_data", {
  expect_error(get_region_summary("not scMethrix"),msg.validateExp)
  
  invisible(lapply(list(scm.mem,scm.h5), function(scm) {
    expect_error(.prepare_plot_data(scm, n_cpgs = "not an int"),msg.validateType)
    
    d <- .prepare_plot_data(scm,na.rm=F)
    expect_equal(dim(d),c(nrow(scm)*ncol(scm),3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_NAs <- nrow(d[is.na(Value)]) 
    d <- .prepare_plot_data(scm,)
    expect_equal(dim(d),c(nrow(scm)*ncol(scm)-n_NAs,3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_cpgs = 50
    d <- .prepare_plot_data(scm,n_cpgs = n_cpgs,na.rm=F)
    expect_equal(dim(d),c(n_cpgs*ncol(scm),3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
    n_NAs <- nrow(d[is.na(Value)]) 
    d <- .prepare_plot_data(scm,n_cpgs = n_cpgs)
    expect_equal(dim(d),c(n_cpgs*ncol(scm)-n_NAs,3))
    expect_equal(colnames(d),c("Sample","Value","Pheno"))
    
  }))
})

invisible(lapply(list(scm.mem,scm.h5), function(scm) {
  
#---- plot_violin ------------------------------------------------------------------------------------------------------
  test_that("plot_violin", {graph_test_helper(scm,plot_violin)})

#---- plot_density -----------------------------------------------------------------------------------------------------
  test_that("plot_density", {graph_test_helper(scm,plot_density)})

#---- plot_coverage ----------------------------------------------------------------------------------------------------
  test_that("plot_coverage", {
    graph_test_helper(scm, plot_coverage, type="histogram")
    graph_test_helper(scm, plot_coverage, type="density")
    graph_test_helper(scm, plot_coverage, type="histogram", pheno="Group")
    graph_test_helper(scm, plot_coverage, type="density",   pheno="Group")
  })
  
#---- plotDensity ------------------------------------------------------------------------------------------------------
  test_that("plotDensity", {
    
    sampNames <- sampleNames(scm)
    
    graph_test_helper2(plotDensity(scm, by = "Sample", type = "Density"), 
                       expected_label = sampNames)
    
    graph_test_helper2(plotDensity(scm, by = "Sample", type = "Violin"), 
                       expected_x = sampNames)                  
    
    graph_test_helper2(plotDensity(scm, by = "Sample", type = "Histogram"), 
                       expected_label = sampNames) 
                         
  })
  
  
  #---- plotStats --------------------------------------------------------------------------------------------------------
  test_that("plotStats", {
    
    chrNames <- levels(seqnames(scm))
    sampNames <- sampleNames(scm)
    
    graph_test_helper2(plotStats(scm, by = "Sample", collapse = FALSE), 
                       expected_x = sampNames, 
                       expected_label = chrNames)
                       
    graph_test_helper2(plotStats(scm, by = "Sample", collapse = TRUE), 
                       expected_x = sampNames, 
                       expected_label = "All")                  
                       
    graph_test_helper2(plotStats(scm, by = "Chromosome", collapse = FALSE), 
                       expected_x = chrNames, 
                       expected_label = sampNames)
    
    graph_test_helper2(plotStats(scm, by = "Chromosome", collapse = TRUE), 
                       expected_x = chrNames, 
                       expected_label = "All")                            
  })
  
  #---- plot_dim_red -----------------------------------------------------------------------------------------------------
  test_that("plot_dim_red", {
    invisible(lapply(list("PCA","tSNE","UMAP"), function(type) {
    
      scm.dimred <- suppressWarnings(impute_regions(scm))
      scm.dimred <- dim_red_scMethrix(scm.dimred, assay="impute", type = type, verbose = FALSE) 
      
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type)
      dim_red_graph_test_helper(scm.dimred, plot_dim_red, dim_red=type, color_anno="Group", shape_anno="Group")
    }))
  })
}))

#---- .getPalette ------------------------------------------------------------------------------------------------------
test_that(".getPalette", {
  
  expect_error(.getPalette(),NA)
  expect_error(.getPalette(nColors = 5, palette = "not a palette"),"Invalid palette")
  expect_error(.getPalette(nColors = "not an number"),msg.validateType)
  expect_error(.getPalette(nColors = -1),msg.validateValue)
  
  # Simple case
  nColors <- 5
  colors <- .getPalette(nColors = nColors)
  expect_equal(length(colors),nColors)
  expect_error(col2rgb(colors, NA))
  
  # Large case
  nColors <- 500
  colors <- .getPalette(nColors = nColors)
  expect_equal(length(colors),nColors)
  expect_error(col2rgb(colors, NA))
  
  # Qualitative
  nColors <- 5
  palette <- "Pastel 1"
  colors <- .getPalette(nColors = nColors, palette = palette)
  expect_equal(length(colors),nColors)
  expect_error(col2rgb(colors, NA))
  expect_equal(colors, colorspace::qualitative_hcl (nColors, palette = palette))
  
  # Sequential
  palette <- "Grays"
  colors <- .getPalette(nColors = nColors, palette = palette)
  expect_equal(length(colors),nColors)
  expect_error(col2rgb(colors, NA))
  expect_equal(colors, colorspace::sequential_hcl (nColors, palette = palette))
  
  # Diverging
  palette <- "Blue-Red"
  colors <- .getPalette(nColors = nColors, palette = palette)
  expect_equal(length(colors),nColors)
  expect_error(col2rgb(colors, NA))
  expect_equal(colors, colorspace::diverging_hcl (nColors, palette = palette))
})

#---- .getShapes -------------------------------------------------------------------------------------------------------
test_that(".getShapes", {
  
  expect_error(.getShapes(),NA)
  
  nShapes = 10
  shapes = .getShapes(nShapes)
  expect_true(.validateType(shapes,"integer"))
  expect_equal(length(shapes), nShapes)
  
  shapes <- .getShapes()
  d <- data.frame(p=shapes,i = 1:length(shapes)-1)
  plot <- ggplot() +
    scale_y_continuous(name="") +
    scale_x_continuous(name="") +
    scale_y_reverse() +
    scale_shape_identity() +
    geom_point(data=d, mapping=aes(x=i%%16, y=i%/%16, shape=p), size=5, fill="red") +
    geom_text(data=d, mapping=aes(x=i%%16, y=i%/%16+0.25, label=p), size=3) +
    theme_void()
  
  expect_error(print(plot),NA)
})

#---- .calcJitter ------------------------------------------------------------------------------------------------------
test_that(".calcJitter", {
  expect_equal(.calcJitter(-5),0)
  expect_equal(.calcJitter(100, max = 0.8),0.8)
  expect_equal(.calcJitter(5),0.5,tolerance = 0.5)
  expect_equal(.calcJitter(5,-10),0)
  expect_equal(.calcJitter(100,10),1)
})



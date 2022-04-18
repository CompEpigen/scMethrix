#---- scMethrix_report -------------------------------------------------------------------------------------------------
#' Creates a detailed interactive HTML summary report from an [`scMethrix-class`] object
#' @description Creates a detailed interactive html summary report from [`scMethrix-class`]  object.
#' If the directory contains required files (from previous run), it directly proceeds to generate html report.
#' @inheritParams generic_scMethrix_function
#' @param outputDir Output directory name where the files should be saved. Default = `tempdir()`
#' @param nCpG Number of CpGs to use for estimating beta value distribution. Default = `10000`
#' @param prefix If provided, the name of the report and the intermediate files will start with the prefix. Default = `""`
#' @return an interactive html report
#' @examples
#' \dontrun{
#' data('methrix_data')
#' methrix::methrix_report(meth = methrix_data)
#' }
#' @export
scMethrix_report <- function(scm, outputDir = tempdir(), prefix = NULL, nCpG = 10000, verbose = FALSE) {
  
  #---- Input validation ---------------------------------------------------
  .validateExp(scm)
  .validateType(outputDir,"directory")
  .validateType(prefix,c("string","null"))
  .validateType(nCpG,"integer")
  
  nCpG <- min(nrow(scm), nCpG)

  plotDensityCoverage <- plotCoverageAll <- plotCoverage <- NULL
  
  #---- Function code ------------------------------------------------------
  
  if (verbose) message("Generating report...", start_time())
  
  if (!dir.exists(outputDir)) 
      dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)
  
  #---- Fig 1 --------------------------------------------------------------
  if (verbose) message(paste0("   Step 1 of 5: Plotting density"))
  
  plotDensityScore <- plotDensity(scm, assay = "score")
  if (has_cov(scm)) plotDensityCoverage <- plotDensity(scm, assay = "counts")
  
  #---- Fig 2 --------------------------------------------------------------
  if (verbose) message(paste0("   Step 2 of 5: Reference CpGs Covered"))
  
  plotCpGsCoveredCountBySample <- plotStats(scm, assay = "score", collapse = FALSE, stat="Count", by = "Sample")
  plotCpGsCoveredFracBySample <-  plotStats(scm, assay = "score", collapse = FALSE, stat="Proportion", by = "Sample")
  
  plotCpGsCoveredCountBySampleAll <- plotStats(scm, assay = "score", collapse = TRUE, stat="Count", by = "Sample")
  plotCpGsCoveredFracBySampleAll <-  plotStats(scm, assay = "score", collapse = TRUE, stat="Proportion", by = "Sample")
  
  plotCpGsCoveredFracByChr <-     plotStats(scm, assay = "score", collapse = FALSE, stat="Proportion", by = "Chr")
  plotCpGsCoveredCountByChr <-    plotStats(scm, assay = "score", collapse = FALSE, stat="Count", by = "Chr")
  
  plotCpGsCoveredFracByChrAll <-     plotStats(scm, assay = "score", collapse = TRUE, stat="Proportion", by = "Chr")
  plotCpGsCoveredCountByChrAll <-    plotStats(scm, assay = "score", collapse = TRUE, stat="Count", by = "Chr")
  
  #---- Fig 3 --------------------------------------------------------------
  if (verbose) message(paste0("   Step 3 of 5: Methylation"))
  
  plotMethylationAll <- plotStats(scm, assay = "score", stat = "Mean", by = "Sample", collapse = TRUE)
  plotMethylation <- plotStats(scm, assay = "score", stat = "Mean", by = "Sample", collapse = FALSE)
  
  #---- Fig 4 --------------------------------------------------------------
  if (verbose) message(paste0("   Step 4 of 5: Coverage"))
  
  if (has_cov(scm)) {
    plotCoverageAll <- plotStats(scm, assay = "counts", stat = "Mean", by = "Sample", collapse = TRUE)
    plotCoverage <- plotStats(scm, assay = "counts", stat = "Mean", by = "Sample", collapse = FALSE)
  }
  
  #---- knit ---------------------------------------------------------------
  
  if (verbose) message(paste0("Knitting report"))
  
  md <- system.file("report", "scMethrix_summarize.Rmd", package = "scMethrix")
  outputFile <- "scMethrix_report.html"
  if (!is.null(prefix)) outputFile <- paste(prefix, outputFile, sep="_")

  rmarkdown::render(input = md, output_file = outputFile,
                    output_dir = outputDir, clean = TRUE, 
                    params = list(hasCov = has_cov(scm), 
                                  plotDensityScore = plotDensityScore,
                                  plotDensityCoverage = plotDensityCoverage,
                                  plotCpGsCoveredCountBySample = plotCpGsCoveredCountBySample,
                                  plotCpGsCoveredCountByChr = plotCpGsCoveredCountByChr,
                                  plotCpGsCoveredFracBySample = plotCpGsCoveredFracBySample,
                                  plotCpGsCoveredFracByChr = plotCpGsCoveredFracByChr,
                                  plotCpGsCoveredCountBySampleAll = plotCpGsCoveredCountBySampleAll,
                                  plotCpGsCoveredCountByChrAll = plotCpGsCoveredCountByChrAll,
                                  plotCpGsCoveredFracBySampleAll = plotCpGsCoveredFracBySampleAll,
                                  plotCpGsCoveredFracByChrAll = plotCpGsCoveredFracByChrAll,
                                  plotMethylationAll = plotMethylationAll,
                                  plotMethylation = plotMethylation,
                                  plotCoverageAll = plotCoverageAll,
                                  plotCoverage = plotCoverage))
  
  stop_time()
}
  
  
  
  
#   fig1 <- suppressWarnings(normalizePath(file.path(outputDir, paste0(prefix, "_MC_per_chr.tsv"))))
#   
#   if (file.exists(fig1)) {
#     message("File already present. Skipping step 1..")
#   } else {
#     if (recal_stats) {
#       per_chr_stat <- getStats(m = meth, per_chr = TRUE)
#       gc()
#     } else {
#       per_chr_stat <- metadata(meth)$descriptive_stats$chr_stat
#       colnames(per_chr_stat)[which(colnames(per_chr_stat) == "chr")] <- "Chromosome"
#     }
#     data.table::fwrite(x = per_chr_stat, file = fig1, sep = "\t")
#   }
#   
#   # Global methylation/Coverage
#   message(paste0("Step 2 of 5: Global methylation/Coverage per sample\n"))
#   
#   fig2 <- suppressWarnings(normalizePath(file.path(outputDir, paste0(prefix, "_global_MC_per_samp.tsv"))))
# 
#   if (file.exists(fig2)) {
#     message("File already present. Skipping step 2..")
#   } else {
#     if (recal_stats) {
#       genome_stat <- getStats(m = meth, per_chr = FALSE)
#       gc()
#     } else {
#       genome_stat <- metadata(meth)$descriptive_stats$genome_stat
#     }
#     data.table::fwrite(x = genome_stat, file = fig2, sep = "\t")
#   }
#   
#   # n CpGs covered per chromomse
#   message(paste0("Step 3 of 5: Reference CpGs covered per chromosome"))
#   
#   fig3 <- suppressWarnings(normalizePath(file.path(outputDir, paste0(prefix, "_n_covered_per_chr.tsv"))))
#   
#   contig_nCpGs <- metadata(meth)$ref_CpG
#   colnames(contig_nCpGs) <- c("chr", "total_CpGs")
#   if (file.exists(fig3)) {
#     message("File already present. Skipping step 3..")
#   } else {
#     if (recal_stats) {
#       cov_tbl <- lapply(seq_len(ncol(meth)), function(i) {
#         data.table::as.data.table(as.data.frame(table(rowData(meth)[which(!is.na(get_matrix(m = meth,
#                                                                                             "M")[, i])), "chr"])))
#       })
#       names(cov_tbl) <- colnames(meth)
#       gc()
#       cov_tbl <- data.table::rbindlist(l = cov_tbl, use.names = TRUE,
#                                        fill = TRUE, idcol = "Sample_Name")
#       colnames(cov_tbl) <- c("Sample_Name", "chr", "n_covered")
#       non_cov_tbl <- merge(cov_tbl, contig_nCpGs, by = "chr",
#                            all.x = TRUE)
#       #non_cov_tbl[, `:=`(n_covered, total_CpGs - n_non_covered)]
#       #non_cov_tbl[, `:=`(n_non_covered, NULL)]
#     } else {
#       non_cov_tbl <- data.table::melt(metadata(meth)$descriptive_stats$n_cpgs_covered,
#                                       id.vars = "chr")
#       colnames(non_cov_tbl) <- c("chr", "Sample_Name", "n_covered")
#       non_cov_tbl <- merge(non_cov_tbl, contig_nCpGs, by = "chr",
#                            all.x = TRUE)
#     }
#     data.table::fwrite(x = non_cov_tbl, file = fig3, sep = "\t")
#   }
#   
#   
#   # Common CpGs covered by all samples
#   message(paste0("Step 4 of 5: Common reference CpGs covered across all samples"))
#   
#   fig4 <- suppressWarnings(normalizePath(file.path(outputDir, paste0(prefix, "_n_covered_by_all_samples.tsv"))))
#   
#   if (file.exists(fig4)) {
#     message("File already present. Skipping step 4..")
#   } else {
#     # na_vec = apply(get_matrix(meth, type = 'M'), MARGIN = 1, anyNA)
#     if (is_h5(meth)) {
#       na_vec <- DelayedMatrixStats::rowAnyNAs(get_matrix(meth, type = "M"))
#     } else {
#       na_vec <- matrixStats::rowAnyNAs(get_matrix(meth, type = "M"))
#     }
#     
#     mf_chr_summary <- as.data.frame(table(rowData(x = meth)[, "chr"],
#                                           na_vec))
#     
#     if (!plot_beta_dist) {
#       rm(na_vec)
#       gc()
#     }
#     
#     data.table::setDT(x = mf_chr_summary)
#     mf_chr_summary <- mf_chr_summary[na_vec == FALSE]
#     mf_chr_summary[, `:=`(na_vec, NULL)]
#     colnames(mf_chr_summary) <- c("chr", "n_CpG")
#     mf_chr_summary <- merge(metadata(meth)$ref_CpG, mf_chr_summary)
#     colnames(mf_chr_summary)[2] <- c("total_CpGs")
#     mf_chr_summary[, `:=`(fract_CpG, n_CpG/total_CpGs)]
#     data.table::fwrite(x = mf_chr_summary, file = fig4, sep = "\t")
#   }
#   
#   # Density plot data
#   message(paste0("Step 5 of 5: beta value distribution"))
#   if (plot_beta_dist) {
#     dens_files <- list.files(path = outputDir, pattern = "*_density\\.tsv\\.gz$")
#     if (length(dens_files) == nrow(colData(meth))) {
#       message("Files already present. Skipping step 5..")
#     } else {
#       if (!exists(x = "na_vec")) {
#         na_vec <- apply(get_matrix(meth, type = "M"), MARGIN = 1,
#                         anyNA)
#       }
#       
#       row_idx <- sample(which(na_vec == FALSE), size = min(beta_nCpG,
#                                                            length(which(na_vec == FALSE))), replace = FALSE)
#       
#       lapply(X = seq_len(nrow(colData(meth))), FUN = function(i) {
#         
#         i_dens <- density(get_matrix(meth, type = "M")[row_idx,
#                                                        i], na.rm = TRUE)
#         if (!is.null(prefix)){
#           if (!dir.exists(paste0(outputDir, "/", prefix, "/"))) {
#             dir.create(path = paste0(outputDir, "/", prefix, "/"), showWarnings = FALSE, recursive = TRUE)
#           }
#           data.table::fwrite(x = data.table::data.table(x = i_dens$x,
#                                                         y = i_dens$y), file = paste0(outputDir, "/", prefix, "/", rownames(colData(x = meth))[i],
#                                                                                      "_density.tsv.gz"), sep = "\t")
#         } else {
#           data.table::fwrite(x = data.table::data.table(x = i_dens$x,
#                                                         y = i_dens$y), file = paste0(outputDir, "/", rownames(colData(x = meth))[i],
#                                                                                      "_density.tsv.gz"), sep = "\t")
#         }
#       })
#       rm(na_vec)
#     }
#     gc()
#   }
#   
#   fig5 <- suppressWarnings(normalizePath(file.path(outputDir, paste0(prefix, "_contig_lens.tsv"))))
# 
#   data.table::fwrite(x = metadata(meth)$chrom_sizes, file = fig5, sep = "\t")
#   
#   message(paste0("Knitting report"))
#   md <- system.file("report", "summarize_methrix.Rmd", package = "methrix")
#   
#   output_file <- paste0(prefix, "_methrix_reports.html")
#   
#   rmarkdown::render(input = md, output_file = output_file,
#                     outputDir = outputDir, clean = TRUE, params = list(prefix = prefix, n_covered_tsv = fig3,
#                                                                          n_covered_by_all_samples_tsv = fig4, mc_per_chr_stat = fig1,
#                                                                          mc_per_sample_stat = fig2, chr_lens = fig5))
#   
#   browseURL(url = paste0(outputDir, "/", output_file))
#   
#   message(data.table::timetaken(started.at = start_proc_time))
# }

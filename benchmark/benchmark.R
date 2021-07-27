
save_dir <- "D:/Git/sampleData/benchmark/"

### Overall function benchmarks #######################################

scm <- scm.big.mem

bench.old <- microbenchmark(
  "old"= old <<- bin_scMethrix_st(scm),times = 1,unit = "s")

bench.new <- microbenchmark(
  "new"= bin_scMethrix_test(scm),
  times = 1,unit = "s"
)

old <- bin_scMethrix_st(scm)
new <- bin_scMethrix_test(scm)

b <- microbenchmark(
  "old" = read_index(files,col_list),
  "new" = read_index2(files,col_list),
  times = 1,unit = "ms")

bench <- microbenchmark(
  ### Operations
  "get_metadata_stats"={get_metadata_stats(scm)},
  "remove_assay"={remove_assay(scm,assay="counts")},
  "merge_scMethrix"={merge_scMethrix(scm[1:(nrow(scm)/2)-1,],scm[((nrow(scm)/2)+1):nrow(scm),],by="row")},
  "export_bed"={export_bed(scm,path=paste0(tempdir(),"/exp"))},
  # "get_region_summary"={get_region_summary(scm,n_chunks = 8,n_threads=8)},
  "get_matrix"={get_matrix(scm)},
  "convert_scMethrix"={convert_scMethrix(scm,h5_dir=paste0(tempdir(),"/out"))},
  "subset_scMethrix"={subset_scMethrix(scm,regions=rowRanges(scm)[1:floor(nrow(scm)/2),])},
  "get_stats"={get_stats(scm)},
  "remove_uncovered"={remove_uncovered(scm)},
  "mask_by_coverage"={mask_by_coverage(scm,low_threshold=2,avg_threshold = 1)},
  "mask_by_sample"={mask_by_sample(scm,low_threshold=2)},
  
  ### Transformations
  "transform_assay"={transform_assay(scm,assay="score",new_assay="binarize",trans=binarize)},
  "bin_scMethrix"={bin_scMethrix(scm,n_threads = 8)},
  #"impute_by_melissa"={impute_by_melissa(scm)},
  #"impute_by_iPCA"={impute_by_iPCA(scm)},
  #"impute_by_RF"={impute_by_RF(scm)},
  #"impute_by_kNN"={impute_by_kNN(scm)},
  "generate_training_set"={generate_training_set(scm)},
  "generate_random_subset"={generate_random_subset(scm)},
  # ""={},
  # ""={},
  # ""={},
  # ""={},
  # ""={},
  times = 1,unit = "s"
)

### Benchmark the indexing #######################################
read.index <- microbenchmark(
  # "base" = {
  #   for (i in 1:length(files)) {
  #     read.delim(files[i],header=FALSE,sep="\t")
  #   }},
  # 
  "fread" = {for (i in 1:length(files)) {data.table::fread(files[i], header = FALSE, select = c(1:2))}},
  "b = 10" = {read_index(files,batch_size=10)},
  "b = 25" = {read_index(files,batch_size=25)},
  "b = 50" = {read_index(files,batch_size=50)},
  "b = 100" = {read_index(files,batch_size=100)},
  "b = 200" = {index <<- read_index(files,batch_size=200)}, 
  times = 1,unit = "s")
read.index$name <- "read_index (batch test)"
saveRDS(read.index, file = paste0(save_dir,"read.index.rds"))
saveRDS(index, file = paste0(save_dir,"index.rds"))

# graph_benchmark(read.index,xlabel="",unit="m")

read.parallel.index <- microbenchmark(
  "read.index" = {read_index(files,batch_size=200)},
  "n.c = 1" = {read_index(files,batch_size=200,n_threads = 1)},
  "n.c = 2" = {read_index(files,batch_size=100,n_threads = 2)},
  "n.c = 4" = {read_index(files,batch_size=50,n_threads = 4)},
  "n.c = 8" = {read_index(files,batch_size=25,n_threads = 8)},
  times = 1,unit = "s")
read.parallel.index$name <- "read_parallel_index (core test)\nbatch = cores*files = 200"
saveRDS(read.parallel.index, file = paste0(save_dir,"read.parallel.index.rds"))

read.subset <- microbenchmark(
  "subset" = {ref_cpgs <<- subset_ref_cpgs(mm10_cpgs,index)},
  times = 1,unit = "m")
read.subset$name <- "subset to mm10 genome"
saveRDS(read.subset, file = paste0(save_dir,"read.subset.rds"))
saveRDS(ref_cpgs, file = paste0(save_dir,"ref_cpgs.rds"))

# graph_benchmark(read.parallel.index,xlabel="",unit="s")

bench.index <- graph_benchmark(rbind(read.subset),xlabel="",unit="m")
bench.index

### Benchmark the genome vs subset ###############################

read.compare <- microbenchmark(
  "mm10" =      {read_hdf5_data(files,mm10_cpgs)},
  "gen_cpgs" =  {read_hdf5_data(files,ref_cpgs)},
  times = 1,unit = "s",
  setup = {
    #mm10_cpgs <- methrix::extract_CPGs(ref_genome = "BSgenome.Mmusculus.UCSC.mm10")
    #mm10_cpgs <- mm10_cpgs$cpgs[,1:3]
    #index <- read_index(files,batch_size=200)
    #ref_cpgs <<- subset_ref_cpgs(mm10_cpgs,index)
  })
read.compare$name <- "mm10 vs subset"
saveRDS(read.compare, file = paste0(save_dir,"read.compare.rds"))

bench.compare <- graph_benchmark(read.compare,xlabel="",unit="m")
bench.compare

### Benchmark the reading ########################################

read.data <- microbenchmark(
  "read.h5" = {read_hdf5_data(files,ref_cpgs)},
  "n.c = 1" = {read_hdf5_data(files,ref_cpgs,n_threads=1)},
  "n.c = 2" = {read_hdf5_data(files,ref_cpgs,n_threads=2)},
  "n.c = 4" = {read_hdf5_data(files,ref_cpgs,n_threads=4)},
  "n.c = 8" = {reads <<- read_hdf5_data(files,ref_cpgs,n_threads=8)},
  times = 1,unit = "s")
read.data$name <- "read_parallel_bed_by_index (core test)"
saveRDS(read.data, file = paste0(save_dir,"read.data.rds"))
saveRDS(reads, file = paste0(save_dir,"reads.rds"))

bench.read <- graph_benchmark(read.data,xlabel="",unit="m")
bench.read

#data <- rbind(read.index,read.parallel.index,read.data)
#saveRDS(data, file = "benchmark.rds")

### Benchmark the reading coverage ########################################

read.coverage <- microbenchmark(
  "w/o coverage" = {read_hdf5_data(files,ref_cpgs,n_threads=8)},
  "with coverage" = {read_hdf5_data(files,ref_cpgs,n_threads=8,cov_idx = c(5,6))},
  times = 1,unit = "s")
read.coverage$name <- "read_parallel_bed_by_index (core test)"
saveRDS(read.coverage, file = paste0(save_dir,"read.coverage"))

bench.coverage <- graph_benchmark(read.coverage,xlabel="",unit="m")
bench.coverage

#data <- rbind(read.index,read.parallel.index,read.data)
#saveRDS(data, file = "benchmark.rds")


### Benchmark object creation ########################################

write.object <- microbenchmark(
  "read.bed" = {read_beds(files,ref_cpgs,reads = reads,h5 = TRUE,h5_dir=dir,verbose=TRUE,replace=TRUE)},
  times = 1,unit = "s",
   setup = {
     dir <- paste0(tempdir(),"\\bench")
 })

write.object$name <- "write_HDF5"
saveRDS(write.object, file = paste0(save_dir,"benchmark - write.rds"))

bench.write <- graph_benchmark(write.object,xlabel="",unit="m")
bench.write

data <- rbind(read.index,read.parallel.index,read.subset,read.data,write.object)
saveRDS(data, file = paste0(save_dir,"benchmark.rds"))

### The graphing functions #######################################

CairoWin()

graph_benchmark <- function (bench, xlabel=NULL,unit="ns") {
  
  if(is.null(xlabel)) xlabel = "expr"
  
  bench <- as.data.frame(bench)
  
  if (unit != "ns") {
    
    bench$time <- switch(unit,
      s = {bench$time/1e9},
      m = {bench$time/60e9},
      h = {bench$time/3600e9},
      stop("Unrecognized unit")
    )
  }

  stats <- aggregate(. ~ expr + name, bench, function(x) c(mean = mean(x), max = max(x)))
  stats <- do.call(data.frame, stats)
  stats$time.mean <- sprintf(stats$time.mean, fmt = '%#.2f')
  
  ggplot(bench, aes(x=expr, y=time)) + 
    geom_boxplot() + 
    labs(x = xlabel, y = paste0("CPU time [",unit,"]")) +
    geom_text(data = stats, aes(label = time.mean, y = time.max),vjust = -0.5) + 
    theme_bw() + expand_limits(y = max(bench$time)*1.05) + facet_grid(~name,scales="free")
  
}
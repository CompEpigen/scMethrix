
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

# graph_benchmark(read.index,xlabel="",unit="s")

read.parallel.index <- microbenchmark(
  "read.index" = {read_index(files,batch_size=200)},
  "n.c = 1" = {read_index(files,batch_size=200,n_threads = 1)},
  "n.c = 2" = {read_index(files,batch_size=100,n_threads = 2)},
  "n.c = 4" = {read_index(files,batch_size=50,n_threads = 4)},
  "n.c = 8" = {read_index(files,batch_size=25,n_threads = 8)},
  times = 1,unit = "s")
read.parallel.index$name <- "read_parallel_index (core test)\nbatch = cores*files = 200"

read.subset <- microbenchmark(
  "subset" = {ref_cpgs <<- subset_ref_cpgs(mm10_cpgs,index)},
  times = 1,unit = "s")
read.subset$name <- "subset to mm10 genome"

# graph_benchmark(read.parallel.index,xlabel="",unit="s")

bench.index <- graph_benchmark(rbind(read.index,read.parallel.index,read.subset),xlabel="",unit="s")
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

### Benchmark the reading ########################################

read.data <- microbenchmark(
  "read.h5" = {read_hdf5_data(files,ref_cpgs)},
  "n.c = 1" = {read_hdf5_data(files,ref_cpgs,n_threads=1)},
  "n.c = 2" = {read_hdf5_data(files,ref_cpgs,n_threads=2)},
  "n.c = 4" = {read_hdf5_data(files,ref_cpgs,n_threads=4)},
  "n.c = 8" = {reads <<- read_hdf5_data(files,ref_cpgs,n_threads=8)},
  times = 1,unit = "s")
read.data$name <- "read_parallel_bed_by_index (core test)"

bench.read <- graph_benchmark(rbind(read.data),xlabel="",unit="m")
bench.read

#data <- rbind(read.index,read.parallel.index,read.data)
#saveRDS(data, file = "benchmark.rds")

### Benchmark object creation ########################################

write.object <- microbenchmark(
  "read.bed" = {read_beds(files,ref_cpgs,reads = reads,h5 = TRUE,h5_dir=dir,verbose=TRUE,replace=TRUE)},
  times = 1,unit = "s")
# ,
#   setup = {
#     #dir <- paste0(tempdir(),"\\bench")
#     })

write.object$name <- "read_bed (core test)"

bench.write <- graph_benchmark(write.object,xlabel="",unit="m")
bench.write

data <- rbind(read.index,read.parallel.index,read.subset,read.data,write.object)
saveRDS(data, file = "benchmark.rds")

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
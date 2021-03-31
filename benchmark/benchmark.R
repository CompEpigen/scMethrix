library(microbenchmark)
library(measurements)
library(magrittr)
library(doParallel)
library(parallel)
library(Cairo)

### Benchmark the indexing #######################################

read.index <- microbenchmark(
  # "base" = {
  #   for (i in 1:length(files)) {
  #     read.delim(files[i],header=FALSE,sep="\t")
  #   }},
  # 
  "fread" = {for (i in 1:length(files)) {data.table::fread(files[i], header = FALSE, select = c(1:2))}},
  "b = 10" = {read_index(files,batch_size=10,verbose = FALSE)},
  "b = 25" = {read_index(files,batch_size=25,verbose = FALSE)},
  "b = 50" = {read_index(files,batch_size=50,verbose = FALSE)},
  "b = 100" = {read_index(files,batch_size=100,verbose = FALSE)}, 
  "b = 200" = {read_index(files,batch_size=250,verbose = FALSE)}, 
  times = 1,unit = "s")
read.index$name <- "read_index (batch test)"


read.parallel.index <- microbenchmark(
  "b = 200" = {read_index(files,batch_size=200,verbose = FALSE)},
  "n.c = 1" = {read_parallel_index(files,batch_size=200,verbose = FALSE,no_cores = 1)},
  "n.c = 2" = {read_parallel_index(files,batch_size=100,verbose = FALSE,no_cores = 2)},
  "n.c = 4" = {read_parallel_index(files,batch_size=50,verbose = FALSE,no_cores = 4)},
  "n.c = 8" = {read_parallel_index(files,batch_size=25,verbose = FALSE,no_cores = 8)},
  times = 1,unit = "s")
read.parallel.index$name <- "read_parallel_index (core test)\nbatch = cores*files = 200"

g1 <- graph_benchmark(rbind(read.index,read.parallel.index),xlabel="Batch size/Number of cores",unit="s")
g1

### Benchmark the reading ########################################


index <- read_parallel_index(files,batch_size=50,verbose = FALSE,no_cores = 4)
read_parallel_bed_by_index <- function(files, index, no_cores = 1) {
  
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl)  
  
  clusterEvalQ(cl, c(library(data.table)))
  
  foo <- function(files,index) {
    for (i in 1:length(files)) {  
      read_bed_by_index(files[i],index)
    }
  }
  
  clusterExport(cl,list('foo',"read_bed_by_index","get_sample_name"))
  
  chunk_files <- split(files, ceiling(seq_along(files)/(length(files)/(no_cores))))
  c(parLapply(cl,chunk_files,fun=foo,index=index))
  
  stopCluster(cl)
}

read.data <- microbenchmark(
  
  "raw" = {for (i in 1:length(files)) {read_bed_by_index(files[i],index)}},
  "n.c = 1" = {read_parallel_bed_by_index(files,index,no_cores=1)},
  "n.c = 2" = {read_parallel_bed_by_index(files,index,no_cores=2)},
  "n.c = 4" = {read_parallel_bed_by_index(files,index,no_cores=4)},
  "n.c = 8" = {read_parallel_bed_by_index(files,index,no_cores=8)},
  times = 1,unit = "s")
read.data$name <- "read_parallel_bed_by_index (core test)"

graph_benchmark(read.data,xlabel="Number of cores",unit="s")

beep()

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
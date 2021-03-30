library(microbenchmark)
library(measurements)
library(magrittr)
library(doParallel)
library(parallel)
library(Cairo)

read.index <- microbenchmark(
  # "base" = {
  #   for (i in 1:length(files)) {
  #     read.delim(files[i],header=FALSE,sep="\t")
  #   }},
  # 
   "fread" = {
     for (i in 1:length(files)) {
       data.table::fread(files[i], header = FALSE, select = c(1:2))
     }},
  
  "batch = 10" = {
    read_index(files,batch_size=10,verbose = FALSE)
  },
  
  "batch = 25" = {
    read_index(files,batch_size=25,verbose = FALSE)
  },
  
  "batch = 50" = {
    read_index(files,batch_size=50,verbose = FALSE)
  },
  
  "batch = 100" = {
    read_index(files,batch_size=100,verbose = FALSE)
  },  
  
  "batch = 200" = {
    read_index(files,batch_size=250,verbose = FALSE)
  },  

times = 1,unit = "s")

read.parallel.index.50 <- microbenchmark(
  # "base" = {
  #   for (i in 1:length(files)) {
  #     read.delim(files[i],header=FALSE,sep="\t")
  #   }},
  # 
  "no.cores = 1" = {
    read_index(files,batch_size=200,verbose = FALSE)
  },
  
  "no.cores = 2" = {
    read_parallel_index(files,batch_size=100,verbose = FALSE,no_cores = 2)
  },
  
  "no.cores = 4" = {
    read_parallel_index(files,batch_size=50,verbose = FALSE,no_cores = 4)
  },
  
  "no.cores = 8" = {
    read_parallel_index(files,batch_size=25,verbose = FALSE,no_cores = 8)
  },
 
  times = 1,unit = "s")

##################################################################

CairoWin()

graph_benchmark(xx,xlabel="read.index (n = 500)",unit="s")

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

  stats <- aggregate(. ~ expr, bench, function(x) c(mean = mean(x), max = max(x)))
  stats <- do.call(data.frame, stats)
  stats$time.mean <- sprintf(stats$time.mean, fmt = '%#.2f')

  ggplot(bench, aes(x=expr, y=time)) + 
    geom_boxplot() + 
    labs(x = xlabel, y = paste0("CPU time [",unit,"]")) +
    geom_text(data = stats, aes(label = time.mean, y = time.max),vjust = -0.5) + 
    theme_bw() + expand_limits(y = max(bench$time)*1.05)
  
}

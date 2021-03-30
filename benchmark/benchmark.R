library(microbenchmark)
library(measurements)
library(magrittr)



#### Benchmark the index generation
read.index <- microbenchmark(
  "base" = {
    for (i in 1:length(files)) {
      read.delim(files[i],header=FALSE,sep="\t")
    }},
  
  
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
  
  "batch = 100" = {
    read_index(files,batch_size=100,verbose = FALSE)
  },  

times = 5,unit = "s")

read.index.plot <- graph_benchmark(read.index, "read.index","s")

graph_benchmark <- function (bench, xlabel,unit=NULL) {

  bench <- as.data.frame(bench)
  
  if (!is.null(unit)) {
    
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
    theme_classic()
  
}


beep()


index <- microbenchmark(
  index.1 = (index.1 <<- read_index(files = files[1:1])),
  index.10 = (index.10 <<- read_index(files = files[1:10])),
  index.100 = (index.100 <<- read_index(files = files[1:100])),
  index.1000 = (index.1000 <<- read_index(files = files[1:1000])),
  times = 1,unit="s"
)

read <- microbenchmark(
  read.1 = (read.1 <<- read_hdf5_data(files = files[1:1], index = index.1)),
  read.10 = (read.10 <<- read_hdf5_data(files = files[1:10], index = index.10)),
  read.100 = (read.100 <<- read_hdf5_data(files = files[1:100], index = index.100)),
  read.1000 = (read.1000 <<- read_hdf5_data(files = files[1:1000], index = index.1000)),
  times = 1,unit="s"
)

write.h5 <- microbenchmark(
  write.h5.1 = read_beds(files = files[1:1],h5=TRUE, index = index.1, reads = read.1),
  write.h5.10 = read_beds(files = files[1:10],h5=TRUE, index = index.10, reads = read.10),
  # write.h5.100 = read_beds(files = files[1:100],h5=TRUE, index = index.100, reads = read.100),
  # write.h5.1000 = read_beds(files = files[1:1000],h5=TRUE, index = index.1000, reads = read.1000),
  times = 1,unit="s"
)

write.mem <- microbenchmark(
  write.mem.1 = read_beds(files = files[1:1]),
  write.mem.10 = read_beds(files = files[1:10]),
  write.mem.100 = read_beds(files = files[1:100]),
  write.mem.1000 = read_beds(files = files[1:1000]),
  times = 1,unit="s"
)




convert_to_unit <- function(x,
                            unit=c("ns", "us", "ms", "s", "t",
                                   "hz", "khz", "mhz", "eps", "f")) {
  unit <- match.arg(unit)
  
  switch (unit,
          t=unit <- sprintf ("%ss", find_prefix(x * 1e-9,
                                                minexp = -9, maxexp = 0, mu = FALSE)),
          f=unit <- sprintf ("%shz", find_prefix(1e9 / x,
                                                 minexp =  0, maxexp = 6, mu = FALSE))
  )
  unit <- tolower(unit)
  switch (unit,
          ns  ={attr(x, "unit") <- "nanoseconds"           ; unclass(x      )},
          us  ={attr(x, "unit") <- "microseconds"          ; unclass(x / 1e3)},
          ms  ={attr(x, "unit") <- "milliseconds"          ; unclass(x / 1e6)}, 
          s   ={attr(x, "unit") <- "seconds"               ; unclass(x / 1e9)},
          eps ={attr(x, "unit") <- "evaluations per second"; unclass(1e9 / x)},
          hz  ={attr(x, "unit") <- "hertz"                 ; unclass(1e9 / x)},
          khz ={attr(x, "unit") <- "kilohertz"             ; unclass(1e6 / x)},
          mhz ={attr(x, "unit") <- "megahertz"             ; unclass(1e3 / x)},
          stop("Unknown unit '", unit, "'.")
  )
}







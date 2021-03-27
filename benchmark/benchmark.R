index <- microbenchmark(
  index.1 = (index.1 <<- read_index(files = files[1:1])),
  index.10 = (index.10 <<- read_index(files = files[1:10])),
  index.100 = (index.100 <<- read_index(files = files[1:100])),
  index.1000 = (index.1000 <<- read_index(files = files[1:1000])),
  times = 1,unit="s"
)

read <- microbenchmark(
  index.1 = read_hdf5_data(files = files[1:1], index = index.1),
  index.10 = read_hdf5_data(files = files[1:10], index = index.10),
  index.100 = read_hdf5_data(files = files[1:100], index = index.100),
  index.1000 = read_hdf5_data(files = files[1:1000], index = index.1000),
  times = 1,unit="s"
)

read.h5 <- microbenchmark(
  index.1 = read_beds(files = files[1:1],h5=TRUE),
  index.10 = read_beds(files = files[1:10],h5=TRUE),
  index.100 = read_beds(files = files[1:100],h5=TRUE),
  index.1000 = read_beds(files = files[1:1000],h5=TRUE),
  times = 1,unit="s"
)

read.mem <- microbenchmark(
  index.1 = read_beds(files = files[1:1]),
  index.10 = read_beds(files = files[1:10]),
  index.100 = read_beds(files = files[1:100]),
  index.1000 = read_beds(files = files[1:1000]),
  times = 1,unit="s"
)












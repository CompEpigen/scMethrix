### Setup
setwd("D:/Documents/School/Thesis/methrix/methrix_data_generation")
files <- list.files (getwd(),full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]
hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
chk_op <- function (vals,op) {
  suffs <- sub(".*\\.", "", names(vals))
  evals <- lapply(unique(suffs),function(suff) {
    browser()
    equiv <- vals[which(suffs == suff)]
    return(op(equiv[[1]],equiv[[2]]))
  })
  
  return (all(unlist(evals)))
}
  
chk_equal <- function (vals) chk_op(vals, function(x,y) all.equal(x,y))
chk_length <- function (vals) chk_op(vals, function(x,y) length(x)==length(y))

### Reading in
scm = read_beds(files=files,h5=FALSE, ref_cpgs = hg19_cpgs$cpgs[,1:3], 
                  chr_idx=1, start_idx=2, M_idx=3, U_idx=4, stranded = FALSE, zero_based = FALSE)
m = read_bedgraphs(files, h5= FALSE, ref_cpgs = hg19_cpgs, 
                   chr_idx=1, start_idx=2, M_idx=3, U_idx=4,stranded = FALSE, zero_based = FALSE)
m <- methrix::remove_uncovered(m)

expect_equal(length(scm),length(m))
expect_equal(rownames(colData(m)),rownames(colData(scm)))

# microbenchmark(
#   scm = ,
#   m = ,
#   times = 3,
#   check = chk_len
# )

### Get matrix
microbenchmark(
  scm = get_matrix(scm),
  m = methrix::get_matrix(m),
  scm.loci = get_matrix(scm,add_loci = TRUE),
  m.loci = methrix::get_matrix(m,add_loci = TRUE),
  scm.gr = get_matrix(scm,add_loci = TRUE,in_granges=TRUE),
  m.gr = methrix::get_matrix(m,add_loci = TRUE,in_granges=TRUE),
  times = 1,
  check = chk_equal
)

### Subset
microbenchmark(
  setup = {regions = data.table(chr = 'chr21', start = 1, end =  30000000)},
  scm.con = subset_scMethrix(scm, contigs = 'chr21'),
    m.con = subset_methrix(m, contigs = 'chr21'),
  scm.sam = subset_scMethrix(scm, samples = 'C1'),
  m.sam = subset_methrix(m, samples = 'C1'),
  scm.reg = subset_scMethrix(scm, regions = regions),
  m.reg = subset_methrix(m, regions = regions),
  times = 1,
  check = chk_length
)

### Coverage Filter
microbenchmark(
  scm.cov = {
    scm.cov <- mask_by_coverage(scm,low_threshold = 1, avg_threshold = 5000)
    scm.cov <- mask_by_sample(scm.cov,low_threshold = 3)
    scm.cov <- remove_uncovered(scm.cov)
  },
  m.cov = methrix::coverage_filter(m, cov_thr = 1, min_samples = 3),
  times = 1,
  check = chk_length
)

### Region filter
microbenchmark(
  setup = {regions = data.table(chr = 'chr21', start = 27867971, end =  27868103)},
  scm.filt = subset_scMethrix(scm,regions = regions,by="exclude"),
  m.filt = methrix::region_filter(m,regions = regions),
  times = 1,
  check = chk_length
)

### Mask methrix
microbenchmark(
  scm.cov = {
    scm.cov <- mask_by_coverage(scm,low_threshold = 1, avg_threshold = 5000)
  },
  m.cov = methrix::mask_methrix(m, low_count = 1, high_quantile  = NULL),
  times = 1,
  check = chk_length
)

### Get region summary
microbenchmark(
  setup = {regions = data.table(chr = c('chr21','chr21'), start = c(27867971,27868110), end =  c(27868103,27868900))},
  scm.sum = get_region_summary(scm,regions = regions),
  m.sum = methrix::get_region_summary(m,regions = regions),
  times = 1,
  check = chk_equal
)

### Get stats
microbenchmark(
  scm.perchr <- getStats(scm),
  m.perchr <- methrix::get_stats(m),
  scm.all <- getStats(scm),
  m.all <- methrix::get_stats(m),
  times = 1,
  check = chk_equal
)

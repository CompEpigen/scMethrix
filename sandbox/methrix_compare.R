### Setup
setwd("D:/Documents/School/Thesis/methrix/methrix_data_generation")
files <- list.files (getwd(),full.names = TRUE)
files <- files[grepl(".*bedgraph$", files,ignore.case = TRUE)]
hg19_cpgs = methrix::extract_CPGs(ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')
chk_len <- function (vals) length(vals[[1]]) == length(vals[[2]]) #TODO make this check for pairs of objects with same suffix
chk_eq <- function (vals) all.equal(vals[[1]],vals[[2]])

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
  times = 3,
  check = chk_eq
)

### Subset
microbenchmark(
  scm.sub = subset_scMethrix(scm, contigs = 'chr21'),
    m.sub = subset_methrix(m, contigs = 'chr21'),
  times = 3,
  check = chk_len
)

### Coverage Filter
microbenchmark(
  scm.cov = {scm <- mask_by_coverage(scm,low_threshold = 1,avg_threshold = 5)
  scm.cov <- mask_by_sample(scm.cov,low_threshold = 3)
  scm.cov <- remove_uncovered(scm.cov)
    },
  m.cov = methrix::coverage_filter(m, cov_thr = 1, min_samples = 3),
  times = 3,
  check = chk_len
)

### Region filter
microbenchmark(
  setup = {regions = data.table(chr = 'chr21', start = 27867971, end =  27868103)},
  scm.filt = subset_scMethrix(scm,regions = regions,by="exclude"),
  m.filt = methrix::region_filter(m,regions = regions),
  times = 3,
  check = chk_len
)

### Get region summary
microbenchmark(
  setup = {regions = data.table(chr = c('chr21','chr21'), start = c(27867971,27868110), end =  c(27868103,27868900))},
  scm.sum = get_region_summary(scm,regions = regions),
  m.sum = methrix::get_region_summary(m,regions = regions),
  times = 3,
  check = chk_eq
)

### Get stats
microbenchmark(
  scm.perchr <- get_stats(scm),
  m.perchr <- methrix::get_stats(m),
  scm.all <- get_stats(scm),
  m.all <- methrix::get_stats(m),
  times = 3,
  check = chk_eq
)



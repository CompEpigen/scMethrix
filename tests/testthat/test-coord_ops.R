#---- cast_granges -----------------------------------------------------------------------------------------------------
test_that("cast_granges",{
  expect_error(cast_granges("not a Granges"),"Invalid input class")
  
  gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equal(gr,cast_granges(gr))
  
  df <- data.frame(chr=as.character(gr@seqnames),start=gr@ranges@start,end=gr@ranges@width)
  expect_equal(gr,cast_granges(df))
})

#---- bin_granges ------------------------------------------------------------------------------------------------------
test_that("bin_granges",{
  expect_error(bin_granges(gr="not granges"),msg.validateType)
  regions <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(1,100))
  expect_equal(length(bin_granges(regions,bin_size=10)),10) 
  expect_equal(reduce(bin_granges(regions,bin_size=10)),regions)
})

#---- subset_ref_cpgs --------------------------------------------------------------------------------------------------
test_that("subset_ref_cpgs",{
  
  ref_cpgs = data.frame(chr="chr1",start=(1:5*2-1), end=(1:5*2))
  gen_cpgs = data.frame(chr="chr1",start=(6:10*2-1), end=(6:10*2))
  gen_cpgs = rbind(ref_cpgs[1:3,],gen_cpgs)
  rows <- sample(nrow(gen_cpgs))
  gen_cpgs <- gen_cpgs[rows, ] #randomize the order
  
  sub_cpgs = subset_ref_cpgs(ref_cpgs, gen_cpgs)
  
  expect_equal(ref_cpgs[1:3,],sub_cpgs)
  
})

#---- .getGRchrStats ---------------------------------------------------------------------------------------------------
test_that(".getGRchrStats",{
  scm <- scm.mem
  gr <- rowRanges(scm)
  idx <- .getGRchrStats(gr)
  setDT(idx)

  expect_equal(dim(idx),c(length(GenomeInfoDb::seqlevels(gr)),7))
  
  chr1 <- rowRanges(subset_scMethrix(scm,contigs = "chr1"))
  chr2 <- rowRanges(subset_scMethrix(scm,contigs = "chr2"))
  
  expect_equal(idx[Chromosome == "chr1",Sites],length(chr1))
  expect_equal(idx[Chromosome == "chr2",Sites],length(chr2))
  
  expect_equal(idx[Chromosome == "chr1",Start.loci],start(range(chr1)))
  expect_equal(idx[Chromosome == "chr2",Start.loci],start(range(chr2)))
  
  expect_equal(idx[Chromosome == "chr1",End.loci],end(range(chr1)))
  expect_equal(idx[Chromosome == "chr2",End.loci],end(range(chr2)))
  
  expect_equal(idx[Chromosome == "chr1",Start.idx],1)
  expect_equal(idx[Chromosome == "chr2",Start.idx],length(chr1)+1)
  
  expect_equal(idx[Chromosome == "chr1",End.idx],length(chr1))
  expect_equal(idx[Chromosome == "chr2",End.idx],length(chr1)+length(chr2))
  
  expect_equal(idx[Chromosome == "chr1",Width],end(range(chr1))-start(range(chr1)))
  expect_equal(idx[Chromosome == "chr2",Width],end(range(chr2))-start(range(chr2)))
})


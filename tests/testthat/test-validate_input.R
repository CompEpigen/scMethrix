
test_that(".validateExp",{
  expect_error(.validateExp("not scMethrix"),msg.validateExp)
  expect_true(.validateExp(scm.mem))
  expect_true(.validateExp(scm.h5))
})

test_that(".validateAssay",{
  expect_error(.validateAssay("not scMethrix"),msg.validateExp)
  expect_equivalent(.validateAssay(scm.mem,assay="score"),"score")
  expect_equivalent(.validateAssay(scm.mem,assay="sco"),"score")
  expect_error(.validateAssay(scm.mem,assay="not an assay"),msg.validateAssay)
})


test_that(".validateArg",{
  func <- function(var = c("banana","banjo")) {}
  
  var = "banana"
  expect_equivalent(.validateArg(var,func),"banana")    
  var = "banjo"
  expect_equivalent(.validateArg(var,func),"banjo")  
  var = "bAnA"
  expect_equivalent(.validateArg(var,func),"banana")  
  expect_error(.validateArg(var,func,ignore.case = F), msg.validateArg)  
  var = c("banana","banjo")
  expect_equivalent(.validateArg(var,func),"banana")    
  var = "ban"
  expect_error(.validateArg(var,func), msg.validateArg)  
  var = "bad input"
  expect_error(.validateArg(var,func), msg.validateArg) 
  
  #TODO: test input for argument list
})

test_that(".validateType",{
  
  expect_error(.validateType(input = "an input"),           "No valid type specified")
  expect_error(.validateType(input = "an input",            type = "not a type"),msg.validateArg)
  
  expect_true (.validateType(input = 10,                     type = "integer"))
  expect_true (.validateType(input = c(10,20,30),            type = "integer"))
  expect_true (.validateType(input = 10,                     type = "INT"))
  expect_error(.validateType(input = "not an int",           type = "integer"),msg.validateType)
  expect_error(.validateType(input = list(10,"not an int"),  type = "integer"),msg.validateType)
  
  expect_true (.validateType(input = 10.5,                   type = "numeric"))
  expect_true (.validateType(input = c(10.5,20.5,30.5),      type = "numeric"))
  expect_error(.validateType(input = "not an num",           type = "numeric"),msg.validateType)
  expect_error(.validateType(input = list(10.5,"not an num"),type = "numeric"),msg.validateType)
  
  expect_true (.validateType(input = "A",                    type = "character"))
  expect_true (.validateType(input = c("A","B"),             type = "character"))
  expect_error(.validateType(input = "not a char",           type = "character"),msg.validateType)
  expect_error(.validateType(input = c("A","not a char"),    type = "character"),msg.validateType)
  expect_error(.validateType(input = list("A",1),            type = "character"),msg.validateType)
  
  expect_true (.validateType(input = "str1",                 type = "string"))
  expect_true (.validateType(input = c("str1","str2"),       type = "string"))
  expect_error(.validateType(input = 0,                      type = "string"),msg.validateType)
  expect_error(.validateType(input = list("str1",0),         type = "string"),msg.validateType)
  
  # expect_true (.validateType(input = c(1,2,3),               type = "vector"))
  # expect_error(.validateType(input = 1,                      type = "vector"),msg.validateType)
  
  expect_true (.validateType(input = list(1,2),              type = "list"))
  expect_true (.validateType(input = list(list(1,2),list(1,2)), type = "list"))
  expect_error(.validateType(input = "not lst",              type = "list"),msg.validateType)
  expect_true (.validateType(input = list(list(1,2),"not lst"), type = "list"),msg.validateType)
  
  expect_true (.validateType(input = TRUE,                   type = "boolean"))
  expect_true (.validateType(input = list(TRUE,TRUE),        type = "boolean"))
  expect_error(.validateType(input = "not a bool",           type = "boolean"),msg.validateType)  
  expect_error(.validateType(input = list(T,"not a bool"),   type = "boolean"),msg.validateType) 
  expect_true (.validateType(input = TRUE,                   type = "logical"))
  expect_error(.validateType(input = "not a bool",           type = "logical"),msg.validateType)
  
  tmpfile <- tempfile()
  file.create(tmpfile)
  expect_true (.validateType(input = tmpfile,                type = "file"))
  expect_true (.validateType(input = list(tmpfile,tmpfile),  type = "file"))
  expect_error(.validateType(input = "not a file",           type = "file"),msg.validateType)
  
  expect_true (.validateType(input = tempdir(),              type = "directory"))
  expect_error(.validateType(input = "not an directory",     type = "directory"),msg.validateType)
  
  gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), IRanges(1:10, width=10:1))
  expect_true (.validateType(input = GRanges(),              type = "GRanges"))
  expect_true (.validateType(input = gr,                     type = "GRanges"))
  expect_true (.validateType(input = list(gr,gr),            type = "GRanges"))
  expect_true(.validateType(input = list(gr,list(gr,gr)), type = "GRanges"))
  expect_error(.validateType(input = "not an Granges",       type = "GRanges"),msg.validateType)
  expect_error(.validateType(input = list(gr,"not an Granges"), type = "GRanges"),msg.validateType)
  expect_error(.validateType(input = list(gr,list(gr,"not an Granges")), type = "GRanges"),msg.validateType)
  
  
  expect_true (.validateType(input = sum,                    type = "function"))
  expect_true (.validateType(input = list(sum,sum),          type = "function"))
  expect_true (.validateType(input = function(x) x+1,        type = "function"))
  expect_error(.validateType(input = "not a function",       type = "function"),msg.validateType)
  
  expect_true (.validateType(input = NULL,                   type = "null"))
  expect_true (.validateType(input = list(NULL,NULL),        type = "null"))
  expect_error(.validateType(input = "not null",             type = "null"),msg.validateType)
  
})

# scMethrix 0.1.0.0003 (19-11-21)

* Combined the masking function into a single function with additional capabilities
* Added more save/load functionality
* Changed reduce_cpgs to now affect entire experiment (now reduce_scMethrix)
* Changed ambiguity in colData that caused previous bugs
* Added stats function for CpGs and samples
* Improved plot outputs
* Added conversion from GRset
* Updated the vignette with new dataset and new functions

# scMethrix 0.1.0.0002 (02-11-21)

* Fixed fill option for import
* Fixed writeBlock HDF5 bug for non-multiple batch sizes
* Added rtracklayer::liftOver support
* Updated vignettes

# scMethrix 0.1.0.0001 (31-10-21)

* Fixed bugs reported by maurerv (#11): 
+ Fixed colData inconsistency in read_beds
+ Added extract_CPGs
+ Fixed scope of .validateValues
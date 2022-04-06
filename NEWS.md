# scMethrix 0.4.0.0000 (3-15-22)
-   Expanded export functions
-   Downgraded most imports to suggests
-   Updated the vignette with new dataset and new functions
-   Added R-CMD check
-   Added HTML report

# scMethrix 0.3.0.0000 (1-12-22)

-   Simplified constructor
-   Corrected imputation for loci with no reads
-   Simplified conversion between HDF5 and in-memory storage
-   Expanded plotting functionality
-   Added consistency checks (validity, HDF5, etc)
-   Improved import performance
-   Added long tests
-   Added support for .idat files

# scMethrix 0.2.0.0000 (19-11-21)

-   Combined the masking function into a single function with additional capabilities
-   Added more save/load functionality
-   Changed reduce_cpgs to now affect entire experiment (now reduce_scMethrix)
-   Changed ambiguity in colData that caused previous bugs
-   Added stats function for CpGs and samples
-   Improved plot outputs
-   Added conversion from GRset

# scMethrix 0.1.0.0002 (02-11-21)

-   Fixed fill option for import
-   Fixed writeBlock HDF5 bug for non-multiple batch sizes
-   Added rtracklayer::liftOver support
-   Updated vignettes

# scMethrix 0.1.0.0001 (31-10-21)

-   Fixed bugs reported by maurerv (\#11):
-   Fixed colData inconsistency in read_beds
-   Added extract_CPGs
-   Fixed scope of .validateValues

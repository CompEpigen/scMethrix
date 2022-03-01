# scMethrix <a href="https://compepigen.github.io/scMethrix/"><img src="vignettes/package_logo.png" align="right" height="150"/></a>

## Fast and efficient summarization of generic BedGraph and .IDAT files from single cell bisulfite sequencing

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/CompEpigen/scMethrix/issues)  
<!--[![R-CMD-check](https://github.com/CompEpigen/scMethrix/workflows/R-CMD-check/badge.svg)](https://github.com/CompEpigen/scMethrix/actions) --> [![License: GPL3](https://img.shields.io/badge/license-GPL3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![minimal R version: 3.6](https://img.shields.io/badge/R%3E%253D-3.6-6666ff.svg)](https://www.r-project.org/) [![Package Version = 0.3.0](https://img.shields.io/badge/Package%20version-0.3.0-orange.svg?style=flat-square)](https://github.com/CompEpigen/scMethrix/blob/master/NEWS) [![Last-changedate](https://img.shields.io/badge/last%20change-2022--01--29-yellowgreen.svg)](https://github.com/CompEpigen/scMethrix/blob/master/NEWS)

<!-- [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FCompEpigen%2FscMethrix&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com) -->

## Introduction

`scMethrix` is an extension of [`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html), and focusing on single-cell methylation data. It provides set of functions which allows easy importing of [BedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)-like files or raw binary intensity data ([.IDAT](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_array_analysis_workflows.pdf)) generated from methylation microarrays. These files are first aggregated into an experiment object, then processed using numerous internal functions or other downstream analysis tools, and can exported into numerous common file formats for further downstream processing.

A detailed example data analysis is provided in our [vignettes](https://compepigen.github.io/scMethrix/articles/x01_load_data.html). Further examples using external tools are also provided.

## Package overview

<img src="https://github.com/CompEpigen/scMethrix/blob/97b06165946e6febd08435b2564cedf8d1dfe43b/vignettes/package_summary.PNG?raw=true" height="65%" width="65%"/>

## Installation

<!--
``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
#Installing stable version from BioConductor
BiocManager::install("methrix")
```
-->

``` r
#Installing developmental version from GitHub
BiocManager::install("CompEpigen/scMethrix")
```

<!--
***NOTE***

Installation from BioConductor requires the BioC and R versions to be up-to-date. This arises from the restrictions imposed by BioConductor community which might cause package incompatibilities with the earlier versions of R (e.g, R \< 4.0). In that case, installing from GitHub might be easier since it is much more merciful with regards to versions.
-->

## Updates

See [here](https://github.com/CompEpigen/scMethrix/blob/master/NEWS)

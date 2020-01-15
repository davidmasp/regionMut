
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regionMut <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

<!-- badges: end -->

The goal of regionMut is to perform an statistical analysis of mutation
enrichment in regions of interest such as genes. Regions will be
represented in channels of bed files.

Strand regions can be used only once.

## Installation

The installation process is divided in the R package which basically
hosts all the functionallity and a small launcher or executable that can
be used to obtain CLI functionality.

Currently it only tests in UNIX systems with bash. It could also be
installed in the WSL.

### R package

<!--
You can install the released version of regionMut from 
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("regionMut")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# Note if private you need the PAT in your .Renviron file
devtools::install_github("davidmasp/regionMut@develop")
```

### Executables

Although you can in theory run regionMut as R package, regionMut ships
with some executables that work as standalone programs. They can be
found in the [`exec/`](exec/) folder and they can be installed in the
system using the [`regionmut`](regionmut) file which is written in bash.

To install the executables (aka launchers) you can either clone the
repository and link the executable to a file in your path or download
the executable as stand alone and place it somewhere in your path. This
last option is exemplified here:

    cd somewhere/in/your/PATH
    URL=""
    wget ${URL}
    
    # test
    regionmut region -h

The launcher won’t drastically change so to update the package and its
functionalities you will need to only update the R package itself.

## Usage

The normal use of regionmut is by the CLI so the 2 steps below should be
performed before this.

### Input files

to fill

An example of how the input table should look is found
[here](inst/testdata/test_bins.tsv)

| master\_name | bins     | file\_name                                      | strand |
| :----------- | :------- | :---------------------------------------------- | :----- |
| bed1         | a        | inst/testdata/test\_bed1\_bina.bed              | \*     |
| bed1         | b        | inst/testdata/test\_bed1\_binb.bed              | \*     |
| bed1         | c        | inst/testdata/test\_bed1\_binc.bed              | \*     |
| bed2         | K        | inst/testdata/test\_bed2\_binK.bed              | \*     |
| bed2         | L        | inst/testdata/test\_bed2\_binL.bed              | \*     |
| bed2         | M        | inst/testdata/test\_bed2\_binM.bed              | \*     |
| bed3         | S\_plus  | inst/testdata/test\_bed3\_binS\_featOnPlus.bed  | \+     |
| bed3         | S\_minus | inst/testdata/test\_bed3\_binS\_featOnMinus.bed | \-     |

### Arguments

to fill

### Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(regionMut)
## basic example code
```
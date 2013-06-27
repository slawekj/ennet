ENNET algorithm
=====

This is an implementation of ENNET algorithm for Gene Regulatory Network inference from mRNA expression data, in form of an R package.

Requirements:
- R environment, the package was implemented and tested using R version 2.13.1.
- foreach, plyr R packages installed.
- g++ compiler and standard libraries.


Installation:
- Download and unpack the source code to the download folder.
- Install the package by invoking the following command in the installation folder, e.g. a folder containing README.md file:

$R CMD build ennet && R CMD INSTALL ennet

Usage:
- Once the R package is installed, please refer to the reference manual, e.g. inside R console type:

\>library(ennet)

\>?ennet

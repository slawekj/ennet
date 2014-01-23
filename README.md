ENNET algorithm
=====

This is an implementation of ENNET algorithm for Gene Regulatory Network inference from mRNA expression data, in form of a R package. See http://www.biomedcentral.com/1752-0509/7/106 for more information.

Requirements:
- R environment, the package was implemented and tested using R version 2.13.1.
- "foreach", "plyr" R packages installed. In order to install the packages go to R console, and type:

\>install.packages("foreach")

\>install.packages("plyr")

- g++ compiler and standard libraries.


Installation:
- Download and unpack the source code to the download folder.
- Install the package by invoking the following command in the installation folder, e.g. a folder containing README.md file, from command line:

$R CMD build ennet && R CMD INSTALL ennet

Usage:
- Once the R package is installed, please refer to the reference manual, e.g. inside R console type:

\>library(ennet)

\>?ennet

LICENSE
=====

Copyright (C) 2013  Janusz Slawek

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program, see LICENSE.

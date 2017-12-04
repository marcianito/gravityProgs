Repository which allows a visualization and user interface for other R packages
======================================================================

## Description

This is an R-shiny App, allowing a web-browser based visualization and usage of other of my R packages.
Supported so far is: [gravityInf](http://github.com/marcianito/gravityInf).
Soon to follow: [UmbrellaEffect](http://github.com/marcianito/UmbrellaEffect).

This repository includes R-shiny Apps, which are run locally on your computer.
It is therefore necessary to have a functional R-base installation and the below mentioned R-package dependencies installed.
Once everything is installed, it can be used with only one command: 
`./runShiny <NameOfAppToLaunch>` 

(Windows launch will follow).

Help files are included in the GUI as hover over.

For bug fixes, comments or further development please contact: mreich@posteo.de.

## Usage

From within the repository folder, execute:

`./runShiny <NameOfAppToLaunch>` 

Options so far for "NameOfAppToLaunch" are:

1. Infiltration
2. UmbrellaEffect

## Dependencies

### Computationally
* r-base version 3.3.1
* [gravityInf](http://github.com/marcianito/gravityInf)
* [UmbrellaEffect](http://github.com/marcianito/UmbrellaEffect)
* shiny R-packages: shiny, shinyjs, shinyBS. shinythemes
* other R-packages necessary for and listed at my other R-packages

**Warning**: depending on your model discretization, in both space and time, it might be
necessary to run this analysis on a cluster (or have at least a high performance machine).

### Data-wise

Depending on your package use of interest, have a look at their needs.


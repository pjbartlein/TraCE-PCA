# TraCE-PCA

## Large-scale PCA/eigenvector, and rolling correlation and PCA of TraCE-21ka data

This repository contains R-code for performing large-scale PCA eigenvector analysis using the singular-value decomposition as implemented by the `RSpectra` package, using TraCE-21ka data as an example.  In addition, rolling correlation and PCA analysis using the `roll` package and `rollapply()` function in the `zoo` package is illustratred using a set of area-average time series for the mid-continent of North America.  THe correlations and compoent loadings are plotted using the 'corrplot' and 'qgraph' packages.

The contents of the various folders are as follows

`/animations`  -- `ImageMagick` scripts for animating qgraphs of the rolling correlations  
`/data` -- Data sets  
`/figs`-- Figures  
`/ncl` -- `NCL` code for plotting component loading maps and scores time series  
`/R` -- `R` code


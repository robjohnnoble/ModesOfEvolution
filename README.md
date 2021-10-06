# Data and code for the paper "Spatial structure governs the mode of tumour evolution"

The `ImageAnalysis` folder contains code and data related to image analysis of tumour histology slides.

Within `ModelOutput`, each folder contains the input and output files from an instance of the [demon model](https://github.com/robjohnnoble/demon_model), which can be analysed using functions in the [demonanalysis R package](https://github.com/robjohnnoble/demonanalysis). For example, the following R code will draw a Muller plot, tumour clone grid, and clonal origin times plot:
``` r
install.packages("devtools")
library(devtools)
  
install_github("robjohnnoble/demonanalysis")
library(demonanalysis)

plot_figure2("ModelOutput/FissionModelK2046")
```

`ModesFigures.R` contains code for plotting other figures from the paper.

`RealTumourTreesData` contains data for constructing the evolutionary trees of real tumours, sourced from published studies.
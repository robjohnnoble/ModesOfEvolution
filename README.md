# Example data from paper "Spatial structure governs the mode of tumour evolution"

Each folder contains the input and output files from an instance of the [demon model](https://github.com/robjohnnoble/demon_model).

Output can be analysed using functions in the [demonanalysis R package](https://github.com/robjohnnoble/demonanalysis).

For example, the following R code will draw a Muller plot, tumour clone grid, and clonal origin times plot:
``` r
install.packages("devtools")
library(devtools)
  
install_github("robjohnnoble/demonanalysis")
library(demonanalysis)

plot_figure2("FissionModel")
```
